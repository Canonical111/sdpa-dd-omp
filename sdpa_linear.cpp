/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-03: fatal eigensolver failure exits non-zero; rdpotf2_ returns on all paths. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-05: the sparse Schur Cholesky reports a non-positive pivot as FAILURE instead of silently zeroing it and returning success. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23: threads the sparse Schur-complement
   Cholesky's k1 update loop, behind measured work gates and an explicit mode switch. Ported from
   sdpa-gmp-omp, where the same threading was derived and reviewed; see
   review/DD-PORT-PLAN-2026-08-23.md in the recipe repo for why it transfers and what does not.

   WHY IT IS RACE-FREE WITHOUT ATOMICS, and bit-identical at any thread count. For a fixed pivot
   i, each k1 owns one target row j = column_index[k1]; column indices ascend and are unique
   within a segment, so distinct k1 write disjoint index ranges. Every read comes from pivot row
   i, which this loop never writes. Each destination is therefore touched by exactly one k1 per
   pivot, the k2 order inside it is unchanged, and a barrier per pivot preserves the i order --
   so every destination sees the identical sequence of subtractions in the identical order. That
   is a statement about the loop's dependence structure, not about the arithmetic, which is why
   it holds at double-double exactly as it does at 256 bits.

   WHAT THE gmp VERSION NEEDS AND THIS ONE DOES NOT. There, `sp_ele[k3] -= tmp * tmp2` builds a
   __gmp_expr temporary -- an mpf_init2 + malloc + free around one multiply and one subtract --
   so it carries a per-thread mpf_class scratch and asserts a uniform-precision invariant over
   the factor. `dd_real` is a 16-byte value type: the same expression allocates nothing and has
   no runtime precision, so the scratch, its plumbing and that invariant do not exist here.

   THE GATES ARE UNCALIBRATED FOR dd. A dd multiply-add is roughly an order of magnitude cheaper
   than an mpf one at 256 bits, so break-even sits at MORE updates, not fewer, and gmp's values
   must not be copied. The defaults below are deliberately conservative floors against pathology;
   calibrate on the target machine before claiming any of them is a tuned optimum. */

#include <sdpa_linear.h>
#include <sdpa_dataset.h>

#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <climits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace sdpa {

namespace {

// Conservative floors, NOT tuned optima -- see the notice above.
#ifndef SDPA_DD_MIN_SPCHOL_WORK
#define SDPA_DD_MIN_SPCHOL_WORK 200000ULL /* matched updates in one pivot */
#endif
#ifndef SDPA_DD_MIN_SPCHOL_WIDTH
#define SDPA_DD_MIN_SPCHOL_WIDTH 32ULL /* target rows in one pivot */
#endif

// Pivot bookkeeping every thread computes identically, so the per-pivot branch cannot diverge
// across the team and no shared array is needed.
inline int spchol_width(int a1, int a2) { return a2 - a1 - 1; }
inline unsigned long long spchol_work(int a1, int a2) {
    const unsigned long long w = (unsigned long long)spchol_width(a1, a2);
    return w * (w + 1ULL) / 2ULL; // matched updates under suffix containment
}

// One k1 column's updates. Identical arithmetic and identical order to the original loop body.
inline void spchol_k1(SparseMatrix &aMat, int *diagonalIndex, int k1, int indexA2) {
    const dd_real a = aMat.sp_ele[k1]; // row i is read-only here, so a copy is safe and cheap
    int k3 = diagonalIndex[aMat.column_index[k1]];
    const int indexB2 = diagonalIndex[aMat.column_index[k1] + 1];
    for (int k2 = k1; k2 < indexA2; ++k2) {
        const dd_real b = aMat.sp_ele[k2];
        const int tmp3 = aMat.column_index[k2];
        // k3 is a MONOTONE cursor across k2 -- never reset. That keeps the search amortised and
        // is correct only because column indices ascend within a segment.
        for (; k3 < indexB2; ++k3) {
            if (aMat.column_index[k3] == tmp3) {
                aMat.sp_ele[k3] -= a * b;
                k3++;
                break;
            }
        }
    }
}

// Bounded env override of a compile-time default. Strict: a malformed value is fatal rather than
// silently ignored, and it is parsed ONCE at the entry point so no code path can skip validation.
unsigned long long spchol_gate(const char *name, unsigned long long dflt) {
    const char *e = getenv(name);
    if (e == NULL || e[0] == '\0') {
        return dflt;
    }
    errno = 0;
    char *end = NULL;
    const unsigned long long v = strtoull(e, &end, 10);
    if (end == e || *end != '\0' || errno == ERANGE) {
        rError(name << " must be a non-negative integer (got \"" << e << "\")");
    }
    return v;
}

// unset/auto -> gated; serial -> never thread; force -> thread every pivot (test/benchmark).
enum SpcholMode { SPCHOL_AUTO, SPCHOL_SERIAL, SPCHOL_FORCE };
SpcholMode spchol_mode() {
    const char *e = getenv("SDPA_SPCHOL_MODE");
    if (e == NULL || e[0] == '\0' || strcmp(e, "auto") == 0) {
        return SPCHOL_AUTO;
    }
    if (strcmp(e, "serial") == 0) {
        return SPCHOL_SERIAL;
    }
    if (strcmp(e, "force") == 0) {
        return SPCHOL_FORCE;
    }
    // Strict: an unrecognised mode is refused, never treated as the default. A silent fallback
    // would let a typo run the very path the caller was trying to select against.
    rError("SDPA_SPCHOL_MODE must be auto, serial or force (got \"" << e << "\")");
    return SPCHOL_AUTO;
}

// Parse-and-validate every spchol tunable exactly once per process, on the first Cholesky of
// EITHER route. Called from both getCholesky overloads, because a problem whose bMat routes dense
// never reaches the sparse one -- and a validator only one path can reach is not a validator.
// (Caught by its own negative test: with this inside the sparse overload alone,
// SDPA_SPCHOL_MODE=bogus was silently accepted on a dense-route problem.)
struct SpcholCfg {
    SpcholMode mode;
    unsigned long long gate_work;
    unsigned long long gate_width;
    bool log;
};

const SpcholCfg &spchol_cfg() {
    static bool done = false;
    static SpcholCfg c;
    if (!done) {
        c.mode = spchol_mode();
        c.gate_work = spchol_gate("SDPA_DD_MIN_SPCHOL_WORK", SDPA_DD_MIN_SPCHOL_WORK);
        c.gate_width = spchol_gate("SDPA_DD_MIN_SPCHOL_WIDTH", SDPA_DD_MIN_SPCHOL_WIDTH);
        c.log = (getenv("SDPA_SPCHOL_LOG") != NULL);
        done = true;
    }
    return c;
}

} // namespace

dd_real Lal::getMinEigen(DenseMatrix &lMat, DenseMatrix &xMat, DenseMatrix &Q, Vector &out, Vector &b, Vector &r, Vector &q, Vector &qold, Vector &w, Vector &tmp, Vector &diagVec, Vector &diagVec2, Vector &workVec) {
    dd_real alpha, beta, value;
    dd_real min = 1.0e+51, min_old = 1.0e+52, min_min = 1.0e+50;
    dd_real error = 1.0e+10;

    int nDim = xMat.nRow;
    int k = 0, kk = 0;

    diagVec.initialize(min_min);
    diagVec2.setZero();
    q.setZero();
    r.initialize(MONE);
    beta = sqrt((dd_real)nDim); // norm of "r"

    // nakata 2004/12/12
    while (k < nDim && k < sqrt((dd_real)nDim) + 10 && beta > 1.0e-16 &&
           (abs(min - min_old) > (1.0e-5) * abs(min) + (1.0e-8)
            // && (fabs(min-min_old) > (1.0e-3)*fabs(min)+(1.0e-6)
            || abs(error * beta) > (1.0e-2) * abs(min) + (1.0e-4))) {
        // rMessage("k = " << k);
        qold.copyFrom(q);
        value = MONE / beta;
        Lal::let(q, '=', r, '*', &value);

        // w = (lMat^T)*q
        w.copyFrom(q);
        Rtrmv("Lower", "Transpose", "NotUnit", nDim, lMat.de_ele, nDim, w.ele, 1);

        Lal::let(tmp, '=', xMat, '*', w);
        w.copyFrom(tmp);
        Rtrmv("Lower", "NoTranspose", "NotUnit", nDim, lMat.de_ele, nDim, w.ele, 1);
        // w = lMat*xMat*(lMat^T)*q
        // rMessage("w = ");
        // w.display();
        Lal::let(alpha, '=', q, '.', w);
        diagVec.ele[k] = alpha;
        Lal::let(r, '=', w, '-', q, &alpha);
        Lal::let(r, '=', r, '-', qold, &beta);
        // rMessage("r = ");
        // r.display();

        if (kk >= sqrt((dd_real)k) || k == nDim - 1 || k > sqrt((dd_real)nDim + 9)) {
            kk = 0;
            out.copyFrom(diagVec);
            b.copyFrom(diagVec2);
            out.ele[nDim - 1] = diagVec.ele[k];
            b.ele[nDim - 1] = 0.0;

            // rMessage("out = ");
            // out.display();
            // rMessage("b = ");
            // b.display();

            mplapackint info;
            int kp1 = k + 1;
            Rsteqr("I_withEigenvalues", kp1, out.ele, b.ele, Q.de_ele, Q.nRow, workVec.ele, info);
            if (info < 0) {
                rError(" rLanczos :: bad argument " << -info << " Q.nRow = " << Q.nRow << ": nDim = " << nDim << ": kp1 = " << kp1);
            } else if (info > 0) {
                rMessage(" rLanczos :: cannot converge " << info);
                break;
            }

            // rMessage("out = ");
            // out.display();
            // rMessage("Q = ");
            // Q.display();

            min_old = min;
#if 0
      min = 1.0e+50;
      error = 1.0e+10;
      for (int i=0; i<k+1; ++i) {
	if (min>out.ele[i]){
	  min = out.ele[i];
	  error = Q.de_ele[k+Q.nCol*i];
	}
      }
#else
            // out have eigen values with ascending order.
            min = out.ele[0];
            error = Q.de_ele[k];
#endif

        } // end of 'if ( kk>=sqrt(k) ...)'
        // printf("\n");

        Lal::let(value, '=', r, '.', r);
        beta = sqrt(value);
        diagVec2.ele[k] = beta;
        ++k;
        ++kk;
    } // end of while
    // rMessage("k = " << k);
    return min - abs(error * beta);
}

dd_real Lal::getMinEigenValue(DenseMatrix &aMat, Vector &eigenVec, Vector &workVec) {
    // aMat is rewritten.
    // aMat must be symmetric.
    // eigenVec is the space of eigen values
    // and needs memory of length aMat.nRow
    // workVec is temporary space and needs
    // 3*aMat.nRow-1 length memory.
    mplapackint N = aMat.nRow;
    mplapackint LWORK, info;
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        LWORK = 3 * N - 1;
        // "N" means that we need not eigen vectors
        // "L" means that we refer only lower triangular.
        Rsyev("NonVectors", "Lower", N, aMat.de_ele, N, eigenVec.ele, workVec.ele, LWORK, info);
        if (info != 0) {
            if (info < 0) {
                rMessage("getMinEigenValue:: info is mistaken " << info);
            } else {
                rMessage("getMinEigenValue:: cannot decomposition");
            }
            exit(EXIT_FAILURE);
            return 0.0;
        }
        return eigenVec.ele[0];
        // Eigen values are sorted by ascending order.
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return 0.0;
}

bool Lal::getInnerProduct(dd_real &ret, Vector &aVec, Vector &bVec) {
    int N = aVec.nDim;
    if (N != bVec.nDim) {
        rError("getInnerProduct:: different memory size");
    }
    ret = Rdot(N, aVec.ele, 1, bVec.ele, 1);

    return _SUCCESS;
}

bool Lal::getInnerProduct(dd_real &ret, BlockVector &aVec, BlockVector &bVec) {
    if (aVec.nBlock != bVec.nBlock) {
        rError("getInnerProduct:: different memory size");
    }
    bool total_judge = _SUCCESS;
    ret = 0.0;
    dd_real tmp_ret;
    for (int l = 0; l < aVec.nBlock; ++l) {
        bool judge = getInnerProduct(tmp_ret, aVec.ele[l], bVec.ele[l]);
        ret += tmp_ret;
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

bool Lal::getInnerProduct(dd_real &ret, DenseMatrix &aMat, DenseMatrix &bMat) {
    if (aMat.nRow != bMat.nRow || aMat.nCol != bMat.nCol) {
        rError("getInnerProduct:: different memory size");
    }
    int length;
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        length = aMat.nRow * aMat.nCol;
        ret = Rdot(length, aMat.de_ele, 1, bMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::getInnerProduct(dd_real &ret, SparseMatrix &aMat, DenseMatrix &bMat) {
    if (aMat.nRow != bMat.nRow || aMat.nCol != bMat.nCol) {
        rError("getInnerProduct:: different memory size");
    }
    int length;
    int amari, shou;

    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        // Attension: in SPARSE case, only half elements
        // are stored. And bMat must be DENSE case.
        ret = 0.0;
// rMessage("aMat.NonZeroCount == " << aMat.NonZeroCount);
#if 0
    for (int index=0; index<aMat.NonZeroCount; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      dd_real value = aMat.sp_ele      [index];
      // rMessage("i=" << i << "  j=" << j);
      if (i==j) {
	ret+= value*bMat.de_ele[i+bMat.nRow*j];
      } else {
	ret+= value*(bMat.de_ele[i+bMat.nRow*j]
		     + bMat.de_ele[j+bMat.nRow*i]);

      }
    }
#else
        amari = aMat.NonZeroCount % 4;
        shou = aMat.NonZeroCount / 4;
        for (int index = 0; index < amari; ++index) {
            int i = aMat.row_index[index];
            int j = aMat.column_index[index];
            dd_real value = aMat.sp_ele[index];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                ret += value * bMat.de_ele[i + bMat.nRow * j];
            } else {
                ret += value * (bMat.de_ele[i + bMat.nRow * j] + bMat.de_ele[j + bMat.nRow * i]);
            }
        }
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = aMat.row_index[index];
            int j1 = aMat.column_index[index];
            dd_real value1 = aMat.sp_ele[index];
            dd_real ret1;
            // rMessage("i=" << i << "  j=" << j);
            if (i1 == j1) {
                ret1 = value1 * bMat.de_ele[i1 + bMat.nRow * j1];
            } else {
                ret1 = value1 * (bMat.de_ele[i1 + bMat.nRow * j1] + bMat.de_ele[j1 + bMat.nRow * i1]);
            }
            int i2 = aMat.row_index[index + 1];
            int j2 = aMat.column_index[index + 1];
            dd_real value2 = aMat.sp_ele[index + 1];
            dd_real ret2;
            // rMessage("i=" << i << "  j=" << j);
            if (i2 == j2) {
                ret2 = value2 * bMat.de_ele[i2 + bMat.nRow * j2];
            } else {
                ret2 = value2 * (bMat.de_ele[i2 + bMat.nRow * j2] + bMat.de_ele[j2 + bMat.nRow * i2]);
            }
            int i3 = aMat.row_index[index + 2];
            int j3 = aMat.column_index[index + 2];
            dd_real value3 = aMat.sp_ele[index + 2];
            dd_real ret3;
            // rMessage("i=" << i << "  j=" << j);
            if (i3 == j3) {
                ret3 = value3 * bMat.de_ele[i3 + bMat.nRow * j3];
            } else {
                ret3 = value3 * (bMat.de_ele[i3 + bMat.nRow * j3] + bMat.de_ele[j3 + bMat.nRow * i3]);
            }
            int i4 = aMat.row_index[index + 3];
            int j4 = aMat.column_index[index + 3];
            dd_real value4 = aMat.sp_ele[index + 3];
            dd_real ret4;
            // rMessage("i=" << i << "  j=" << j);
            if (i4 == j4) {
                ret4 = value4 * bMat.de_ele[i4 + bMat.nRow * j4];
            } else {
                ret4 = value4 * (bMat.de_ele[i4 + bMat.nRow * j4] + bMat.de_ele[j4 + bMat.nRow * i4]);
            }
            // ret += ret1;
            // ret += ret2;
            // ret += ret3;
            // ret += ret4;
            ret += (ret1 + ret2 + ret3 + ret4);
        }
#endif
        break;
    case SparseMatrix::DENSE:
        length = aMat.nRow * aMat.nCol;
        ret = Rdot(length, aMat.de_ele, 1, bMat.de_ele, 1);
        break;
    }
    return _SUCCESS;
}

bool Lal::getCholesky(DenseMatrix &retMat, DenseMatrix &aMat) {
    (void)spchol_cfg(); // validate the spchol tunables even when this problem routes dense
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("getCholesky:: different memory size");
    }
    int length, shou, amari;
    mplapackint info;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
#if 1
        Rpotrf("Lower", retMat.nRow, retMat.de_ele, retMat.nRow, info);
#else
        info = choleskyFactorWithAdjust(retMat);
#endif
        if (info != 0) {
            rMessage("cannot cholesky decomposition");
            rMessage("Could you try with smaller gammaStar?");
            return FAILURE;
        }
// Make matrix as lower triangular matrix
#if 0
    for (int j=0; j<retMat.nCol; ++j) {
      for (int i=0; i<j; ++i) {
	retMat.de_ele[i+retMat.nCol*j] = 0.0;
      }
    }
#else
        for (int j = 0; j < retMat.nCol; ++j) {
            shou = j / 4;
            amari = j % 4;
            for (int i = 0; i < amari; ++i) {
                retMat.de_ele[i + retMat.nCol * j] = 0.0;
            }
            for (int i = amari, count = 0; count < shou; ++count, i += 4) {
                retMat.de_ele[i + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 1 + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 2 + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 3 + retMat.nCol * j] = 0.0;
            }
        }
#endif
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

// nakata 2004/12/01
// modified 2008/05/20    "aMat.sp_ele[indexA1] = 0.0;"
// aMat = L L^T
//
// 2026-08-05 ("C7"): report failure instead of silently patching a non-positive pivot.
//
// Upstream set a negative pivot to 0.0, kept going, and returned `true`
// unconditionally, so a Schur complement that is not positive definite was reported
// to the caller as a successful factorisation. The zeroed pivot then propagates: the
// scaling loop below multiplies that whole column by 0.0, so the returned "L" is not
// a factor of anything and the search direction computed from it is silently wrong.
// A pivot of exactly 0.0 was worse still -- it fell into the `else` arm and computed
// 1.0/sqrt(0.0) = +inf, which poisons the rest of the factorisation with inf/NaN.
// Both cases are now a reported FAILURE, which is what the DENSE twin
// choleskyFactorWithAdjust already does (`info > 0` -> rMessage + FAILURE, below),
// and what the CALLER already expects: Newton::compute_DyVec does
//     bool ret = Lal::getCholesky(sparse_bMat, diagonalIndex);
//     if (ret == FAILURE) { return FAILURE; }
// The plumbing was there; only this callee was incapable of ever using it.
//
// Numerically inert on any input whose Schur complement stays positive definite:
// neither branch is reachable unless a pivot is <= 0, and this function's arithmetic
// is otherwise untouched. Proven bit-identical with patches/regress.sh (10 problems,
// including truss6, the only problem in a 93-problem census that reaches this
// function at all).
namespace {

// The untouched serial factorisation. Kept as its own routine so that SDPA_SPCHOL_MODE=serial,
// a one-member team, and a build without OpenMP all run THE SAME code rather than an emulation
// of it inside a parallel region.
bool spchol_serial(SparseMatrix &aMat, int *diagonalIndex) {
    const int nDIM = aMat.nRow;
    for (int i = 0; i < nDIM; ++i) {
        const int a1 = diagonalIndex[i], a2 = diagonalIndex[i + 1];
        if (!(aMat.sp_ele[a1] > 0.0)) {
            // Non-positive (or NaN) pivot: the Schur complement is not positive
            // definite, so there is no Cholesky factor to return. Say so.
            rMessage("sparse cholesky miss condition :: not positive definite"
                     << " :: pivot " << i << " = " << aMat.sp_ele[a1]);
            return FAILURE;
        }
        aMat.sp_ele[a1] = 1.0 / sqrt(aMat.sp_ele[a1]); // inverse diagonal
        for (int k1 = a1 + 1; k1 < a2; ++k1) {
            aMat.sp_ele[k1] *= aMat.sp_ele[a1];
        }
        for (int k1 = a1 + 1; k1 < a2; ++k1) {
            spchol_k1(aMat, diagonalIndex, k1, a2);
        }
    }
    return _SUCCESS;
}

} // namespace

bool Lal::getCholesky(SparseMatrix &aMat, int *diagonalIndex) {
    if (aMat.type != SparseMatrix::SPARSE) {
        rError("Lal::getCholesky aMat is not sparse format");
    }

    // Validated on the first Cholesky of either route (see spchol_cfg), so a malformed value is
    // refused whether or not this problem's bMat routes sparse.
    const SpcholCfg &cfg = spchol_cfg();
    const SpcholMode mode = cfg.mode;
    const unsigned long long gate_work = cfg.gate_work;
    const unsigned long long gate_width = cfg.gate_width;
    const bool want_log = cfg.log;

#ifdef _OPENMP
    const int team = (mode == SPCHOL_SERIAL) ? 1 : omp_get_max_threads();
#else
    const int team = 1;
#endif
    if (team < 2) {
        if (want_log) {
            rMessage("spchol: serial (team=" << team << ", mode=" << (int)mode << ")");
        }
        return spchol_serial(aMat, diagonalIndex);
    }

#ifdef _OPENMP
    const int nDIM = aMat.nRow;
    const bool force = (mode == SPCHOL_FORCE);
    bool fail = false;
    int failed_pivot = -1;
    bool one_thread = false, one_thread_r = false;
    unsigned long long pivots_par = 0, pivots_seq = 0;
    int actual = 1, max_width = 0;

#pragma omp parallel num_threads(team) default(none)                                        \
    shared(aMat, diagonalIndex, fail, failed_pivot, one_thread, one_thread_r,               \
           pivots_par, pivots_seq, actual, max_width)                                        \
    firstprivate(nDIM, force, gate_work, gate_width)
    {
        // THE team, read inside the region that owns it. num_threads() is a REQUEST: OMP_DYNAMIC
        // or a resource limit may return fewer threads, including one, and OpenMP does not
        // promise two consecutive regions get the same team. Probing with a throwaway region and
        // dispatching on its answer would be an unwarranted inference about this one.
        //
        // The branch is COLLECTIVE: every thread evaluates the same team size, so none takes a
        // different path and the worksharing constructs below are encountered by all of them in
        // the same order. With one thread, that thread runs the untouched serial routine.
        const int here = omp_get_num_threads();
        if (here < 2) {
            one_thread = true;
            one_thread_r = spchol_serial(aMat, diagonalIndex);
        } else {
#pragma omp single
            { actual = here; }

            for (int i = 0; i < nDIM; ++i) {
                const int a1 = diagonalIndex[i], a2 = diagonalIndex[i + 1];
                // Only the pivot test and its inverse square root are inherently one-thread work.
#pragma omp single
                {
                    if (!(aMat.sp_ele[a1] > 0.0)) {
                        fail = true;
                        failed_pivot = i;
                    } else {
                        aMat.sp_ele[a1] = 1.0 / sqrt(aMat.sp_ele[a1]);
                    }
                } // implicit barrier: publishes fail and the inverted diagonal
                if (fail)
                    break; // identical in every thread, so the worksharing sequence stays identical

                // Row scaling, shared out rather than left on one thread: it is O(L_i) per pivot
                // and O(nnz) over the factorisation, so leaving it serial puts an Amdahl floor
                // under the routine that only bites once the update loop is fast. Elementwise and
                // independent -- each k1 writes its own element and reads only the already
                // published diagonal -- so this is bit-identical.
#pragma omp for
                for (int k1 = a1 + 1; k1 < a2; ++k1)
                    aMat.sp_ele[k1] *= aMat.sp_ele[a1];
                // implicit barrier: the scaled row must be complete before any update reads it

                // Same predicate in every thread -- derived only from values all of them have.
                const bool par = force || (spchol_work(a1, a2) >= gate_work &&
                                           (unsigned long long)spchol_width(a1, a2) >= gate_width);
                if (par) {
#pragma omp single
                    {
                        pivots_par++;
                        if (spchol_width(a1, a2) > max_width)
                            max_width = spchol_width(a1, a2);
                    }
#pragma omp for schedule(dynamic)
                    for (int k1 = a1 + 1; k1 < a2; ++k1)
                        spchol_k1(aMat, diagonalIndex, k1, a2);
                    // implicit barrier: pivot i completes before i+1 begins. NEVER nowait.
                } else {
                    // Too small to share out. One thread does it; the end barrier is retained so
                    // the team stays in lockstep across pivots.
#pragma omp single
                    {
                        pivots_seq++;
                        for (int k1 = a1 + 1; k1 < a2; ++k1)
                            spchol_k1(aMat, diagonalIndex, k1, a2);
                    }
                }
            }
        }
    }

    if (one_thread) {
        if (want_log) {
            rMessage("spchol: requested team=" << team << " but received 1; ran serial");
        }
        return one_thread_r;
    }
    if (want_log) {
        // Non-vacuity evidence: a caller can see whether any pivot actually ran threaded, rather
        // than inferring it from a wall-clock difference.
        rMessage("spchol: team=" << actual << " pivots_threaded=" << pivots_par
                                << " pivots_serial=" << pivots_seq << " max_width=" << max_width
                                << " gate_work=" << gate_work << " gate_width=" << gate_width);
    }
    if (fail) {
        rMessage("sparse cholesky miss condition :: not positive definite"
                 << " :: pivot " << failed_pivot);
        return FAILURE;
    }
    return _SUCCESS;
#endif
}

bool Lal::getInvLowTriangularMatrix(DenseMatrix &retMat, DenseMatrix &aMat) {
    // Make inverse with refference only to lower triangular.
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("getCholesky:: different memory size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        retMat.setIdentity();
        /* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: Rtrsm -> Rtrsm_omp. This solve
           against an identity RHS is the larger half of the Cholesky-inverse phase (42.0% of
           wall on gpp124-1, 39.7% on arch0 at 8 threads) and was entirely serial. Rtrsm_omp
           splits it over the columns of retMat, which leaves each output element's arithmetic
           and its order untouched, so the result is bit-identical to Rtrsm; it falls back to
           Rtrsm below its work gate and for any case but Left/Lower/NoTranspose. Deliberately
           NOT done by threading Rtrsm itself, which would also thread Rpotrf's inner solves.
           See git log. */
        Rtrsm_omp("Left", "Lower", "NoTraspose", "NonUnitDiagonal", aMat.nRow, aMat.nCol, MONE, aMat.de_ele, aMat.nRow, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::getSymmetrize(DenseMatrix &aMat) {
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        if (aMat.nRow != aMat.nCol) {
            rError("getSymmetrize:: different memory size");
        }
        for (int index = 0; index < aMat.nRow - 1; ++index) {
            int index1 = index + index * aMat.nRow + 1;
            int index2 = index + (index + 1) * aMat.nRow;
            int length = aMat.nRow - 1 - index;
            // aMat.de_ele[index1] += aMat.de_ele[index2]
            Raxpy(length, MONE, &aMat.de_ele[index2], aMat.nRow, &aMat.de_ele[index1], 1);
            // aMat.de_ele[index1] /= 2.0
            dd_real half = 0.5;
            Rscal(length, half, &aMat.de_ele[index1], 1);
            // aMat.de_ele[index2] = aMat.de_ele[index1]
            Rcopy(length, &aMat.de_ele[index1], 1, &aMat.de_ele[index2], aMat.nRow);
        }
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::getTranspose(DenseMatrix &retMat, DenseMatrix &aMat) {
    if (aMat.nRow != aMat.nCol) {
        rError("getTranspose:: different memory size");
        // Of course, a non-symmetric matrix has
        // its transposed matrix,
        // but in this algorithm we have to make
        // transposed matrix only when symmetric matrix.
    }
    retMat.copyFrom(aMat);
    switch (aMat.type) {
    case DenseMatrix::DENSE:
#if 0
    for (int i=0; i<aMat.nRow; ++i) {
      for (int j=0; j<=i; ++j) {
	int index1 = i+aMat.nCol*j;
	int index2 = j+aMat.nCol*i;
	retMat.de_ele[index1] = aMat.de_ele[index2];
	retMat.de_ele[index2] = aMat.de_ele[index1];
      }
    }
#else
        for (int i = 0; i < aMat.nRow; ++i) {
            int shou = (i + 1) / 4;
            int amari = (i + 1) / 4;
            for (int j = 0; j < amari; ++j) {
                int index1 = i + aMat.nCol * j;
                int index2 = j + aMat.nCol * i;
                retMat.de_ele[index1] = aMat.de_ele[index2];
                retMat.de_ele[index2] = aMat.de_ele[index1];
            }
            for (int j = amari, counter = 0; counter < shou; ++counter, j += 4) {
                int index1 = i + aMat.nCol * j;
                int index_1 = j + aMat.nCol * i;
                retMat.de_ele[index1] = aMat.de_ele[index_1];
                retMat.de_ele[index_1] = aMat.de_ele[index1];
                int index2 = i + aMat.nCol * (j + 1);
                int index_2 = (j + 1) + aMat.nCol * i;
                retMat.de_ele[index2] = aMat.de_ele[index_2];
                retMat.de_ele[index_2] = aMat.de_ele[index2];
                int index3 = i + aMat.nCol * (j + 2);
                int index_3 = (j + 2) + aMat.nCol * i;
                retMat.de_ele[index3] = aMat.de_ele[index_3];
                retMat.de_ele[index_3] = aMat.de_ele[index3];
                int index4 = i + aMat.nCol * (j + 3);
                int index_4 = (j + 3) + aMat.nCol * i;
                retMat.de_ele[index4] = aMat.de_ele[index_4];
                retMat.de_ele[index_4] = aMat.de_ele[index4];
            }
        }
#endif
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

int Lal::rdpotf2_(char *uplo, int *n, double *a, int *lda, int *info) {
#if USE_DOUBLE
    int nRow = *lda;
    for (int j = 0; j < *n; ++j) {
        double ajj = a[j + nRow * j] - ddot_f77(&j, &a[j], lda, &a[j], lda);

        // Here is point.(start)
        if (ajj <= (float)-1.0e-6) {
            a[j + j * nRow] = ajj;
            *info = j + 1;
            return 0;
        }
        if (ajj <= (float)1.0e-14) {
            ajj = 1e100;
            a[j + j * nRow] = ajj;
        } else {
            ajj = sqrt(ajj);
            a[j + j * nRow] = ajj;
        }
        // Here is point.(end)

        if (j < *n - 1) {
            int i = *n - 1 - j;
            dgemv_f77("No transpose", &i, &j, &DMONE, &a[j + 1], lda, &a[j], lda, &DONE, &a[(j + 1) + nRow * j], &IONE, strlen("No transpose"));
            double d1 = 1.0 / ajj;
            dscal_f77(&i, &d1, &a[(j + 1) + nRow * j], &IONE);
        }
    }
#endif
    return 0;
}

int Lal::rdpotrf_(char *uplo, int *n, double *a, int *lda, int *info) {
#if USE_DOUBLE
    // This funciton makes Cholesky factorization
    // in only case Lower Triangular.
    // That is, A will be L*L**T, not U**T*U.
    int nRow = *lda;
    *info = 0;

    int nb = ilaenv_f77(&IONE, "DPOTRF", "L", n, &IMONE, &IONE, &IMONE, strlen("DPOTRF"), strlen("L"));
    if (nb <= 1 || nb >= *n) {
        // Here is point.
        rdpotf2_(uplo, n, a, lda, info);
    } else {

        for (int j = 0; j < *n; j += nb) {
            int jb = min(nb, *n - j);
            dsyrk_f77("Lower", "No transpose", &jb, &j, &DMONE, &a[j], lda, &DONE, &a[j + nRow * j], lda, strlen("Lower"), strlen("No transpose"));
            // Here is point.
            rdpotf2_("Lower", &jb, &a[j + nRow * j], lda, info);
            if (*info != 0) {
                *info = *info + j - 1;
                return 0;
            }
            if (j + jb <= *n - 1) {
                int i = *n - j - jb;
                dgemm_f77("No transpose", "Transpose", &i, &jb, &j, &DMONE, &a[j + jb], lda, &a[j], lda, &DONE, &a[(j + jb) + nRow * j], lda, strlen("No transpose"), strlen("Transpose"));
                dtrsm_f77("Right", "Lower", "Transpose", "Non-unit", &i, &jb, &DONE, &a[j + nRow * j], lda, &a[(j + jb) + nRow * j], lda, strlen("Right"), strlen("Lower"), strlen("Transpose"), strlen("Non-unit"));
            }
        }
    }
#endif
    return 0;
}

bool Lal::choleskyFactorWithAdjust(DenseMatrix &aMat) {
    mplapackint info = 0;
#if 1
    // aMat.display();
    TimeStart(START1);
    Rpotrf("Lower", aMat.nRow, aMat.de_ele, aMat.nRow, info);
    TimeEnd(END1);
    // rMessage("Schur colesky  ::"  << TimeCal(START1,END1));
    // aMat.display();
#elif 1
    dpotrf_f77("Lower", &aMat.nRow, aMat.de_ele, &aMat.nRow, &info, strlen("Lower"));
#else
    rdpotrf_("Lower", &aMat.nRow, aMat.de_ele, &aMat.nRow, &info);
#endif
    if (info < 0) {
        rMessage("cholesky argument is wrong " << -info);
    } else if (info > 0) {
        rMessage("cholesky miss condition :: not positive definite"
                 << " :: info = " << info);

        return FAILURE;
    }
    return _SUCCESS;
#if 0
  dd_real ZERO_DETECT = 1.0e-3;
  dd_real NONZERO = 1.0e-7;
  // no idea version
  // if Cholesky factorization failed, then exit soon.
  int info = 1; // info == 0 means success
  int start = 0;
  while (start<aMat.nRow) {
    int N = aMat.nRow - start;
    dpotf2_("Lower",&N,&aMat.de_ele[start+start*aMat.nRow],
	    &aMat.nRow,&info);
    if (info <=0) {
      // rMessage("Cholesky is very nice");
      break;
    }
    start += (info-1); // next target
    dd_real wrong = aMat.de_ele[start+start*aMat.nRow];
    if (wrong < -ZERO_DETECT) {
      rMessage("cholesky adjust position " << start);
      rMessage("cannot cholesky decomposition"
	       " with adjust " << wrong);
      return FAILURE;
    }
    aMat.de_ele[start+start*aMat.nRow] = NONZERO;
    if (start<aMat.nRow-1) {
      // improve the right down element of 0
      for (int j=1; j<=aMat.nRow-1-start; ++j) {
	dd_real& migi  = aMat.de_ele[start+(start+j)*aMat.nRow];
	dd_real& shita = aMat.de_ele[(start+j)+start*aMat.nRow];
	dd_real& mishi = aMat.de_ele[(start+j)+(start+j)*aMat.nRow];
	// rMessage(" mishi = " << mishi);
	if (mishi < NONZERO) {
	  // rMessage(" mishi < NONZERO ");
	  mishi = NONZERO;
	  migi  = NONZERO * 0.1;
	  shita = NONZERO * 0.1;
	} else if (migi*shita > NONZERO*mishi) {
	  // rMessage(" migi*migi > NONZERO*mishi ");
	  migi  = sqrt(NONZERO*mishi) * 0.99;
	  shita = sqrt(NONZERO*mishi) * 0.99;
	}
      }
    }
    rMessage("cholesky adjust position " << start);
  }
  if (info < 0) {
    rError("argument is something wrong " << info);
  }
  return _SUCCESS;
#endif
}

bool Lal::solveSystems(Vector &xVec, DenseMatrix &aMat, Vector &bVec) {
    // aMat must have done Cholesky factorized.
    if (aMat.nCol != xVec.nDim || aMat.nRow != bVec.nDim || aMat.nRow != aMat.nCol) {
        rError("solveSystems:: different memory size");
    }
    if (aMat.type != DenseMatrix::DENSE) {
        rError("solveSystems:: matrix type must be DENSE");
    }
    xVec.copyFrom(bVec);
    Rtrsv("Lower", "NoTranspose", "NonUnit", aMat.nRow, aMat.de_ele, aMat.nCol, xVec.ele, 1);
    Rtrsv("Lower", "Transpose", "NonUnit", aMat.nRow, aMat.de_ele, aMat.nCol, xVec.ele, 1);
    return _SUCCESS;
}

// nakata 2004/12/01
bool Lal::solveSystems(Vector &xVec, SparseMatrix &aMat, Vector &bVec) {
#define TUNEUP 0
#if TUNEUP
    if (aMat.nCol != xVec.nDim || aMat.nRow != bVec.nDim || aMat.nRow != aMat.nCol) {
        printf("A.row:%d A.col:%d x.row:%d b.row:%d\n", aMat.nCol, aMat.nRow, xVec.nDim, bVec.nDim);
        rError("solveSystems(sparse):: different memory size");
    }
    int length;
    int amari, shou, counter;

    switch (aMat.type) {
    case SparseMatrix::SPARSE:
#endif
        // Attension: in SPARSE case, only half elements
        // are stored. And bMat must be DENSE case.
        // rMessage("aMat.NonZeroCount == " << aMat.NonZeroCount);
        xVec.copyFrom(bVec);
#if TUNEUP

        shou = aMat.NonZeroCount / 4;
        amari = aMat.NonZeroCount % 4;
        int i, j;
        dd_real value;

        for (int index = 0; index < amari; ++index) {
            int i = aMat.row_index[index];
            int j = aMat.column_index[index];
            dd_real value = aMat.sp_ele[index];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
        }

        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            i = aMat.row_index[index];
            j = aMat.column_index[index];
            value = aMat.sp_ele[index];
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
            i = aMat.row_index[index + 1];
            j = aMat.column_index[index + 1];
            value = aMat.sp_ele[index + 1];
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
            i = aMat.row_index[index + 2];
            j = aMat.column_index[index + 2];
            value = aMat.sp_ele[index + 2];
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
            i = aMat.row_index[index + 3];
            j = aMat.column_index[index + 3];
            value = aMat.sp_ele[index + 3];
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
        }

        for (int index = aMat.NonZeroCount - 1; index >= aMat.NonZeroCount - amari; --index) {
            i = aMat.row_index[index];
            j = aMat.column_index[index];
            value = aMat.sp_ele[index];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[j] -= value * xVec.ele[i];
            }
        }

        for (int index = aMat.NonZeroCount - amari - 1, counter = 0; counter < shou; ++counter, index -= 4) {
            i = aMat.row_index[index];
            j = aMat.column_index[index];
            value = aMat.sp_ele[index];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[i] -= value * xVec.ele[j];
            }
            i = aMat.row_index[index - 1];
            j = aMat.column_index[index - 1];
            value = aMat.sp_ele[index - 1];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[i] -= value * xVec.ele[j];
            }
            i = aMat.row_index[index - 2];
            j = aMat.column_index[index - 2];
            value = aMat.sp_ele[index - 2];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[i] -= value * xVec.ele[j];
            }
            i = aMat.row_index[index - 3];
            j = aMat.column_index[index - 3];
            value = aMat.sp_ele[index - 3];
            // rMessage("i=" << i << "  j=" << j);
            if (i == j) {
                xVec.ele[i] *= value;
            } else {
                xVec.ele[i] -= value * xVec.ele[j];
            }
        }
#else
    for (int index = 0; index < aMat.NonZeroCount; ++index) {
        int i = aMat.row_index[index];
        int j = aMat.column_index[index];
        dd_real value = aMat.sp_ele[index];
        // rMessage("i=" << i << "  j=" << j);
        if (i == j) {
            xVec.ele[i] *= value;
        } else {
            xVec.ele[j] -= value * xVec.ele[i];
        }
    }
    for (int index = aMat.NonZeroCount - 1; index >= 0; --index) {
        int i = aMat.row_index[index];
        int j = aMat.column_index[index];
        dd_real value = aMat.sp_ele[index];
        value = aMat.sp_ele[index];
        // rMessage("i=" << i << "  j=" << j);
        if (i == j) {
            xVec.ele[i] *= value;
        } else {
            xVec.ele[i] -= value * xVec.ele[j];
        }
    }
#endif
#if TUNEUP
        break;
    case SparseMatrix::DENSE:
        xVec.copyFrom(bVec);
        F77_FUNC(dtrsv, DTRSV)("Lower", "NoTranspose", "NonUnit", &aMat.nRow, aMat.de_ele, &aMat.nCol, xVec.ele, &IONE);
        F77_FUNC(dtrsv, DTRSV)("Lower", "Transpose", "NonUnit", &aMat.nRow, aMat.de_ele, &aMat.nCol, xVec.ele, &IONE);
        return _SUCCESS;
    }
#endif
    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix &retMat, DenseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nRow || bMat.nCol != retMat.nCol || retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
            // attension::scalar is loval variable.
        }
        Rgemm("NoTranspose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);

        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix &retMat, SparseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nRow || bMat.nCol != retMat.nCol) {
        rError("multiply :: different matrix size");
    }
    retMat.setZero();
    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            for (int index = 0; index < aMat.NonZeroCount; ++index) {
                int i = aMat.row_index[index];
                int j = aMat.column_index[index];
                dd_real value = aMat.sp_ele[index];
                if (i != j) {
#define MULTIPLY_NON_ATLAS 0
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[i + retMat.nRow * t] += value * bMat.de_ele[j + bMat.nRow * t];
                        retMat.de_ele[j + retMat.nRow * t] += value * bMat.de_ele[i + bMat.nRow * t];
                    }
#else
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + i, retMat.nRow);
                    Raxpy(bMat.nCol, value, bMat.de_ele + i, bMat.nRow, retMat.de_ele + j, retMat.nRow);
#endif
                } else {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[j + retMat.nRow * t] += value * bMat.de_ele[j + bMat.nRow * t];
                    }
#else
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + j, retMat.nRow);
#endif
                }
            }    // end of 'for index'
        } else { // scalar!=NULL
            for (int index = 0; index < aMat.NonZeroCount; ++index) {
                int i = aMat.row_index[index];
                int j = aMat.column_index[index];
                dd_real value = aMat.sp_ele[index] * (*scalar);
                if (i != j) {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[i + retMat.nRow * t] += value * bMat.de_ele[j + bMat.nRow * t];
                        retMat.de_ele[j + retMat.nRow * t] += value * bMat.de_ele[i + bMat.nRow * t];
                    }
#else
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + i, retMat.nRow);
                    Raxpy(bMat.nCol, value, bMat.de_ele + i, bMat.nRow, retMat.de_ele + j, retMat.nRow);
#endif
                } else {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[j + retMat.nRow * t] += value * bMat.de_ele[j + bMat.nRow * t];
                    }
#else
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + j, retMat.nRow);
#endif
                }
            } // end of 'for index'
        }     // end of 'if (scalar==NULL)
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            scalar = &MONE;
            // attension:: scalar is local variable.
        }
        Rgemm("NoTranspose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);
        break;

    } // end of switch

    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix &retMat, DenseMatrix &aMat, SparseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nRow || bMat.nCol != retMat.nCol) {
        rError("multiply :: different matrix size");
    }
    retMat.setZero();
    switch (bMat.type) {
    case SparseMatrix::SPARSE:
        // rMessage("Here will be faster by atlas");
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            for (int index = 0; index < bMat.NonZeroCount; ++index) {
                int i = bMat.row_index[index];
                int j = bMat.column_index[index];
                dd_real value = bMat.sp_ele[index];
                if (i != j) {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[t + retMat.nRow * j] += value * aMat.de_ele[t + aMat.nRow * i];
                        retMat.de_ele[t + retMat.nRow * i] += value * aMat.de_ele[t + aMat.nRow * j];
                    }
#else
                    Raxpy(bMat.nCol, value, &aMat.de_ele[aMat.nRow * j], 1, &retMat.de_ele[retMat.nRow * i], 1);
                    Raxpy(bMat.nCol, value, &aMat.de_ele[aMat.nRow * i], 1, &retMat.de_ele[retMat.nRow * j], 1);
#endif
                } else {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[t + retMat.nRow * j] += value * aMat.de_ele[t + aMat.nRow * j];
                    }
#else
                    Raxpy(bMat.nCol, value, &aMat.de_ele[aMat.nRow * j], 1, &retMat.de_ele[retMat.nRow * j], 1);
#endif
                }
            }    // end of 'for index'
        } else { // scalar!=NULL
            for (int index = 0; index < bMat.NonZeroCount; ++index) {
                int i = bMat.row_index[index];
                int j = bMat.column_index[index];
                dd_real value = bMat.sp_ele[index] * (*scalar);
                if (i != j) {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[t + retMat.nCol * j] += value * aMat.de_ele[t + bMat.nCol * i];
                        retMat.de_ele[t + retMat.nCol * i] += value * aMat.de_ele[t + bMat.nCol * j];
                    }
#else
                    Raxpy(bMat.nCol, value, aMat.de_ele + (aMat.nRow * j), 1, retMat.de_ele + (retMat.nRow * i), 1);
                    Raxpy(bMat.nCol, value, aMat.de_ele + (aMat.nRow * i), 1, retMat.de_ele + (retMat.nRow * j), 1);
#endif
                } else {
#if MULTIPLY_NON_ATLAS
                    for (int t = 0; t < bMat.nCol; ++t) {
                        retMat.de_ele[t + retMat.nCol * j] += value * aMat.de_ele[t + aMat.nCol * j];
                    }
#else
                    Raxpy(bMat.nCol, value, aMat.de_ele + (aMat.nRow * j), 1, retMat.de_ele + (retMat.nRow * j), 1);
#endif
                }
            } // end of 'for index'
        }     // end of 'if (scalar==NULL)
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            scalar = &MONE;
            // attension: scalar is local variable.
        }
        Rgemm("NoTranspose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);
        break;
    } // end of switch

    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix &retMat, DenseMatrix &aMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("multiply :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
        Rscal(length, *scalar, retMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::multiply(Vector &retVec, Vector &aVec, dd_real *scalar) {
    if (retVec.nDim != aVec.nDim) {
        rError("multiply :: different vector size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    Rcopy(retVec.nDim, aVec.ele, 1, retVec.ele, 1);
    Rscal(retVec.nDim, *scalar, retVec.ele, 1);
    return _SUCCESS;
}

bool Lal::multiply(BlockVector &retVec, BlockVector &aVec, dd_real *scalar) {
    if (retVec.nBlock != aVec.nBlock) {
        rError("multiply:: different memory size");
    }
    bool total_judge = _SUCCESS;
    for (int l = 0; l < aVec.nBlock; ++l) {
        bool judge = multiply(retVec.ele[l], aVec.ele[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

bool Lal::multiply(Vector &retVec, DenseMatrix &aMat, Vector &bVec, dd_real *scalar) {
    if (retVec.nDim != aMat.nRow || aMat.nCol != bVec.nDim || bVec.nDim != retVec.nDim) {
        rError("multiply :: different matrix size");
    }
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
        }
        Rgemv("NoTranspose", aMat.nRow, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bVec.ele, 1, 0.0, retVec.ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::tran_multiply(DenseMatrix &retMat, DenseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nCol || aMat.nRow != bMat.nRow || bMat.nCol != retMat.nCol || retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
            // scalar is local variable
        }
        // The Point is the first argument is "Transpose".
        Rgemm("Transpose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele, aMat.nCol, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }

    return _SUCCESS;
}

bool Lal::multiply_tran(DenseMatrix &retMat, DenseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nCol || bMat.nRow != retMat.nRow || retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
        }
        // The Point is the first argument is "NoTranspose".
        Rgemm("NoTranspose", "Transpose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nCol, 0.0, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::plus(Vector &retVec, Vector &aVec, Vector &bVec, dd_real *scalar) {
    if (retVec.nDim != aVec.nDim || aVec.nDim != bVec.nDim) {
        rError("plus :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    if (retVec.ele != aVec.ele) {
        Rcopy(retVec.nDim, aVec.ele, 1, retVec.ele, 1);
    }
    Raxpy(retVec.nDim, *scalar, bVec.ele, 1, retVec.ele, 1);
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix &retMat, DenseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow || retMat.nCol != bMat.nCol || retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("plus :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        if (retMat.de_ele != aMat.de_ele) {
            Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
        }
        Raxpy(length, *scalar, bMat.de_ele, 1, retMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix &retMat, SparseMatrix &aMat, DenseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow || retMat.nCol != bMat.nCol) {
        rError("plus :: different matrix size");
    }
    // ret = (*scalar) * b
    if (multiply(retMat, bMat, scalar) == FAILURE) {
        return FAILURE;
    }
    int length;
    // ret += a
    int shou, amari;
    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
#if 0
    for (int index=0; index<aMat.NonZeroCount; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      dd_real value = aMat.sp_ele      [index];
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
#else
        shou = aMat.NonZeroCount / 4;
        amari = aMat.NonZeroCount % 4;
        for (int index = 0; index < amari; ++index) {
            int i = aMat.row_index[index];
            int j = aMat.column_index[index];
            dd_real value = aMat.sp_ele[index];
            if (i != j) {
                retMat.de_ele[i + retMat.nCol * j] += value;
                retMat.de_ele[j + retMat.nCol * i] += value;
            } else {
                retMat.de_ele[i + retMat.nCol * i] += value;
            }
        } // end of 'for index'
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = aMat.row_index[index];
            int j1 = aMat.column_index[index];
            dd_real value1 = aMat.sp_ele[index];
            if (i1 != j1) {
                retMat.de_ele[i1 + retMat.nCol * j1] += value1;
                retMat.de_ele[j1 + retMat.nCol * i1] += value1;
            } else {
                retMat.de_ele[i1 + retMat.nCol * i1] += value1;
            }
            int i2 = aMat.row_index[index + 1];
            int j2 = aMat.column_index[index + 1];
            dd_real value2 = aMat.sp_ele[index + 1];
            if (i2 != j2) {
                retMat.de_ele[i2 + retMat.nCol * j2] += value2;
                retMat.de_ele[j2 + retMat.nCol * i2] += value2;
            } else {
                retMat.de_ele[i2 + retMat.nCol * i2] += value2;
            }
            int i3 = aMat.row_index[index + 2];
            int j3 = aMat.column_index[index + 2];
            dd_real value3 = aMat.sp_ele[index + 2];
            if (i3 != j3) {
                retMat.de_ele[i3 + retMat.nCol * j3] += value3;
                retMat.de_ele[j3 + retMat.nCol * i3] += value3;
            } else {
                retMat.de_ele[i3 + retMat.nCol * i3] += value3;
            }
            int i4 = aMat.row_index[index + 3];
            int j4 = aMat.column_index[index + 3];
            dd_real value4 = aMat.sp_ele[index + 3];
            if (i4 != j4) {
                retMat.de_ele[i4 + retMat.nCol * j4] += value4;
                retMat.de_ele[j4 + retMat.nCol * i4] += value4;
            } else {
                retMat.de_ele[i4 + retMat.nCol * i4] += value4;
            }
        } // end of 'for index'
#endif
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        length = retMat.nRow * retMat.nCol;
        Raxpy(length, 1.0, aMat.de_ele, 1, retMat.de_ele, 1);
        break;
    } // end of switch
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix &retMat, DenseMatrix &aMat, SparseMatrix &bMat, dd_real *scalar) {
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow || retMat.nCol != bMat.nCol) {
        rError("plus :: different matrix size");
    }
    // ret = a
    if (retMat.copyFrom(aMat) == FAILURE) {
        return FAILURE;
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length, shou, amari;
    // ret += (*scalar) * b
    switch (bMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
#if 0
    for (int index=0; index<bMat.NonZeroCount; ++index) {
      int        i = bMat.row_index   [index];
      int        j = bMat.column_index[index];
      dd_real value = bMat.sp_ele      [index] * (*scalar);
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
#else
        shou = bMat.NonZeroCount / 4;
        amari = bMat.NonZeroCount % 4;
        for (int index = 0; index < amari; ++index) {
            int i = bMat.row_index[index];
            int j = bMat.column_index[index];
            dd_real value = bMat.sp_ele[index] * (*scalar);
            if (i != j) {
                retMat.de_ele[i + retMat.nCol * j] += value;
                retMat.de_ele[j + retMat.nCol * i] += value;
            } else {
                retMat.de_ele[i + retMat.nCol * i] += value;
            }
        } // end of 'for index'
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = bMat.row_index[index];
            int j1 = bMat.column_index[index];
            dd_real value1 = bMat.sp_ele[index] * (*scalar);
            if (i1 != j1) {
                retMat.de_ele[i1 + retMat.nCol * j1] += value1;
                retMat.de_ele[j1 + retMat.nCol * i1] += value1;
            } else {
                retMat.de_ele[i1 + retMat.nCol * i1] += value1;
            }
            int i2 = bMat.row_index[index + 1];
            int j2 = bMat.column_index[index + 1];
            dd_real value2 = bMat.sp_ele[index + 1] * (*scalar);
            if (i2 != j2) {
                retMat.de_ele[i2 + retMat.nCol * j2] += value2;
                retMat.de_ele[j2 + retMat.nCol * i2] += value2;
            } else {
                retMat.de_ele[i2 + retMat.nCol * i2] += value2;
            }
            int i3 = bMat.row_index[index + 2];
            int j3 = bMat.column_index[index + 2];
            dd_real value3 = bMat.sp_ele[index + 2] * (*scalar);
            if (i3 != j3) {
                retMat.de_ele[i3 + retMat.nCol * j3] += value3;
                retMat.de_ele[j3 + retMat.nCol * i3] += value3;
            } else {
                retMat.de_ele[i3 + retMat.nCol * i3] += value3;
            }
            int i4 = bMat.row_index[index + 3];
            int j4 = bMat.column_index[index + 3];
            dd_real value4 = bMat.sp_ele[index + 3] * (*scalar);
            if (i4 != j4) {
                retMat.de_ele[i4 + retMat.nCol * j4] += value4;
                retMat.de_ele[j4 + retMat.nCol * i4] += value4;
            } else {
                retMat.de_ele[i4 + retMat.nCol * i4] += value4;
            }
        } // end of 'for index'
#endif
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        length = retMat.nRow * retMat.nCol;
        Raxpy(length, *scalar, bMat.de_ele, 1, retMat.de_ele, 1);
        break;
    } // end of switch
    return _SUCCESS;
}

bool Lal::plus(BlockVector &retVec, BlockVector &aVec, BlockVector &bVec, dd_real *scalar) {
    if (retVec.nBlock != aVec.nBlock || retVec.nBlock != bVec.nBlock) {
        rError("plus:: different nBlock size");
    }
    bool total_judge = _SUCCESS;
    for (int l = 0; l < retVec.nBlock; ++l) {
        bool judge = plus(retVec.ele[l], aVec.ele[l], bVec.ele[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

// ret = a '*' (*scalar)
bool Lal::let(Vector &retVec, const char eq, Vector &aVec, const char op, dd_real *scalar) {
    switch (op) {
    case '*':
        return multiply(retVec, aVec, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '*' (*scalar)
bool Lal::let(BlockVector &retVec, const char eq, BlockVector &aVec, const char op, dd_real *scalar) {
    switch (op) {
    case '*':
        return multiply(retVec, aVec, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(Vector &retVec, const char eq, Vector &aVec, const char op, Vector &bVec, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retVec, aVec, bVec, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retVec, aVec, bVec, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' '*' 't' 'T' b*(*scalar)
bool Lal::let(DenseMatrix &retMat, const char eq, DenseMatrix &aMat, const char op, DenseMatrix &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
        return multiply(retMat, aMat, bMat, scalar);
        break;
    case 't':
        // ret = aMat**T * bMat
        return tran_multiply(retMat, aMat, bMat, scalar);
        break;
    case 'T':
        // ret = aMat * bMat**T
        return multiply_tran(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool Lal::let(DenseMatrix &retMat, const char eq, SparseMatrix &aMat, const char op, DenseMatrix &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
        return multiply(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool Lal::let(DenseMatrix &retMat, const char eq, DenseMatrix &aMat, const char op, SparseMatrix &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
        return multiply(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = aMat '*' '/' bVec
bool Lal::let(Vector &rVec, const char eq, DenseMatrix &aMat, const char op, Vector &bVec) {
    switch (op) {
    case '*':
        return multiply(rVec, aMat, bVec, NULL);
        break;
    case '/':
        // ret = aMat^{-1} * bVec;
        // aMat is positive definite
        // and already colesky factorized.
        return solveSystems(rVec, aMat, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// nakata 2004/12/01
// ret = aMat '*' '/' bVec
bool Lal::let(Vector &rVec, const char eq, SparseMatrix &aMat, const char op, Vector &bVec) {
    switch (op) {
    case '/':
        // ret = aMat^{-1} * bVec;
        // aMat is positive definite
        // and already colesky factorized.
        return solveSystems(rVec, aMat, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, Vector &aVec, const char op, Vector &bVec) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aVec, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, DenseMatrix &aMat, const char op, DenseMatrix &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, DenseMatrix &aMat, const char op, SparseMatrix &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, bMat, aMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, SparseMatrix &aMat, const char op, DenseMatrix &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, BlockVector &aVec, const char op, BlockVector &bVec) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aVec, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

/////////////////////////////////////////////////////////////////////////

bool Lal::getInnerProduct(dd_real &ret, DenseLinearSpace &aMat, DenseLinearSpace &bMat) {
    bool total_judge = _SUCCESS;
    ret = 0.0;
    dd_real tmp_ret;

    // for SDP
    if (aMat.SDP_nBlock != bMat.SDP_nBlock) {
        rError("getInnerProduct:: different memory size");
    }
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::getInnerProduct(tmp_ret, aMat.SDP_block[l], bMat.SDP_block[l]);
        ret += tmp_ret;
        if (judge == FAILURE) {
            rMessage(" something failed");
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  if (aMat.SOCP_nBlock != bMat.SOCP_nBlock) {
    rError("getInnerProduct:: different memory size");
  }
  for (int l=0; l<aMat.SOCP_nBlock; ++l) {
    bool judge = Lal::getInnerProduct(tmp_ret,aMat.SOCP_block[l],bMat.SOCP_block[l]);
    ret += tmp_ret;
    if (judge == FAILURE) {
      rMessage(" something failed");
      total_judge = FAILURE;
    }
  }
#endif

    // for LP
    if (aMat.LP_nBlock != bMat.LP_nBlock) {
        rError("getInnerProduct:: different memory size");
    }
    for (int l = 0; l < aMat.LP_nBlock; ++l) {
        tmp_ret = aMat.LP_block[l] * bMat.LP_block[l];
        ret += tmp_ret;
    }

    return total_judge;
}

bool Lal::getInnerProduct(dd_real &ret, SparseLinearSpace &aMat, DenseLinearSpace &bMat) {
    bool total_judge = _SUCCESS;
    ret = 0.0;
    dd_real tmp_ret;

    // for SDP
    for (int l = 0; l < aMat.SDP_sp_nBlock; ++l) {
        int index = aMat.SDP_sp_index[l];
        bool judge = Lal::getInnerProduct(tmp_ret, aMat.SDP_sp_block[l], bMat.SDP_block[index]);
        ret += tmp_ret;
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  for (int l=0; l<aMat.SOCP_sp_nBlock; ++l) {
    int index = aMat.SOCP_sp_index[l];
    bool judge = Lal::getInnerProduct(tmp_ret,aMat.SOCP_sp_block[l],bMat.SOCP_block[index]);
    ret += tmp_ret;
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
#endif

    for (int l = 0; l < aMat.LP_sp_nBlock; ++l) {
        int index = aMat.LP_sp_index[l];
        tmp_ret = aMat.LP_sp_block[l] * bMat.LP_block[index];
        ret += tmp_ret;
    }

    return total_judge;
}

bool Lal::multiply(DenseLinearSpace &retMat, DenseLinearSpace &aMat, dd_real *scalar) {
    bool total_judge = _SUCCESS;

    // for SDP
    if (retMat.SDP_nBlock != aMat.SDP_nBlock) {
        rError("multiply:: different memory size");
    }
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::multiply(retMat.SDP_block[l], aMat.SDP_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  if (retMat.SOCP_nBlock!=aMat.SOCP_nBlock) {
    rError("multiply:: different memory size");
  }
  for (int l=0; l<aMat.SOCP_nBlock; ++l) {
    bool judge = Lal::multiply(retMat.SOCP_block[l],aMat.SOCP_block[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
#endif

    // fo LP
    if (retMat.LP_nBlock != aMat.LP_nBlock) {
        rError("multiply:: different memory size");
    }
    for (int l = 0; l < aMat.LP_nBlock; ++l) {
        retMat.LP_block[l] = aMat.LP_block[l] * (*scalar);
    }

    return total_judge;
}

bool Lal::plus(DenseLinearSpace &retMat, DenseLinearSpace &aMat, DenseLinearSpace &bMat, dd_real *scalar) {
    bool total_judge = _SUCCESS;

    // for SDP
    if (retMat.SDP_nBlock != aMat.SDP_nBlock || retMat.SDP_nBlock != bMat.SDP_nBlock) {
        rError("plus:: different nBlock size");
    }
    for (int l = 0; l < retMat.SDP_nBlock; ++l) {
        bool judge = Lal::plus(retMat.SDP_block[l], aMat.SDP_block[l], bMat.SDP_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  if (retMat.SOCP_nBlock!=aMat.SOCP_nBlock 
      || retMat.SOCP_nBlock!=bMat.SOCP_nBlock) {
    rError("plus:: different nBlock size");
  }
  for (int l=0; l<retMat.SOCP_nBlock; ++l) {
    bool judge = Lal::plus(retMat.SOCP_block[l],aMat.SOCP_block[l],
			   bMat.SOCP_block[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
#endif

    // for LP
    if (retMat.LP_nBlock != aMat.LP_nBlock || retMat.LP_nBlock != bMat.LP_nBlock) {
        rError("plus:: different nBlock size");
    }
    for (int l = 0; l < retMat.LP_nBlock; ++l) {
        if (scalar == NULL) {
            retMat.LP_block[l] = aMat.LP_block[l] + bMat.LP_block[l];
        } else {
            retMat.LP_block[l] = aMat.LP_block[l] + bMat.LP_block[l] * (*scalar);
        }
    }

    return total_judge;
}

// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
bool Lal::plus(DenseLinearSpace &retMat, SparseLinearSpace &aMat, DenseLinearSpace &bMat, dd_real *scalar) {
    bool total_judge = _SUCCESS;

    // for SDP
    for (int l = 0; l < aMat.SDP_sp_nBlock; ++l) {
        int index = aMat.SDP_sp_index[l];
        bool judge = Lal::plus(retMat.SDP_block[index], aMat.SDP_sp_block[l], bMat.SDP_block[index], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  for (int l=0; l<aMat.SOCP_sp_nBlock; ++l) {
    int index = aMat.SOCP_sp_index[l];
    bool judge = Lal::plus(retMat.SOCP_block[index],aMat.SOCP_sp_block[l],
			   bMat.SOCP_block[index],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
#endif

    // for LP
    for (int l = 0; l < aMat.LP_sp_nBlock; ++l) {
        int index = aMat.LP_sp_index[l];
        if (scalar == NULL) {
            retMat.LP_block[index] = aMat.LP_sp_block[l] + bMat.LP_block[index];
        } else {
            retMat.LP_block[index] = aMat.LP_sp_block[l] + bMat.LP_block[index] * (*scalar);
        }
    }

    return total_judge;
}

// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
bool Lal::plus(DenseLinearSpace &retMat, DenseLinearSpace &aMat, SparseLinearSpace &bMat, dd_real *scalar) {
    bool total_judge = _SUCCESS;

    // for SDP
    for (int l = 0; l < bMat.SDP_sp_nBlock; ++l) {
        int index = bMat.SDP_sp_index[l];
        bool judge = Lal::plus(retMat.SDP_block[index], aMat.SDP_block[index], bMat.SDP_sp_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for SOCP
#if 0
  for (int l=0; l<bMat.SOCP_sp_nBlock; ++l) {
    int index = bMat.SOCP_sp_index[l];
    bool judge = Lal::plus(retMat.SOCP_block[index],aMat.SOCP_block[index],
			   bMat.SOCP_sp_block[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
#endif

    // for LP
    for (int l = 0; l < bMat.LP_sp_nBlock; ++l) {
        int index = bMat.LP_sp_index[l];
        if (scalar == NULL) {
            retMat.LP_block[index] = aMat.LP_block[index] + bMat.LP_sp_block[l];
        } else {
            retMat.LP_block[index] = aMat.LP_block[index] + bMat.LP_sp_block[l] * (*scalar);
        }
    }

    return total_judge;
}

bool Lal::getSymmetrize(DenseLinearSpace &aMat) {
    bool total_judge = _SUCCESS;
    // for SDP
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::getSymmetrize(aMat.SDP_block[l]);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

bool Lal::getTranspose(DenseLinearSpace &retMat, DenseLinearSpace &aMat) {
    // for SDP
    if (retMat.SDP_nBlock != aMat.SDP_nBlock) {
        rError("getTranspose:: different memory size");
    }
    bool total_judge = _SUCCESS;
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::getTranspose(retMat.SDP_block[l], aMat.SDP_block[l]);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

// ret = a '*' (*scalar)
bool Lal::let(DenseLinearSpace &retMat, const char eq, DenseLinearSpace &aMat, const char op, dd_real *scalar) {
    switch (op) {
    case '*':
        return multiply(retMat, aMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(DenseLinearSpace &retMat, const char eq, DenseLinearSpace &aMat, const char op, DenseLinearSpace &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(DenseLinearSpace &retMat, const char eq, SparseLinearSpace &aMat, const char op, DenseLinearSpace &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(DenseLinearSpace &retMat, const char eq, DenseLinearSpace &aMat, const char op, SparseLinearSpace &bMat, dd_real *scalar) {
    dd_real minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, DenseLinearSpace &aMat, const char op, DenseLinearSpace &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, SparseLinearSpace &aMat, const char op, DenseLinearSpace &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(dd_real &ret, const char eq, DenseLinearSpace &aMat, const char op, SparseLinearSpace &bMat) {
    switch (op) {
    case '.':
        return getInnerProduct(ret, bMat, aMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

} // namespace sdpa
