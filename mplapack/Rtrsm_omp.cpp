/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: new file. Column-parallel Rtrsm for the
   single case the solver uses (Left / Lower / NoTranspose); every other case is handed to the
   serial Rtrsm. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: work gate calibrated on thanos and pi;
   MPLAPACK_OMP_MIN_TRSM_WORK_8T renamed to MPLAPACK_OMP_MIN_TRSM_WORK because the measured
   break-even does not track the team size. Measured crossover is n = 10-13 (thanos) and
   n = 12-18 (pi), not the n >= 64 previously assumed. See mplapack_omp_tuning.h and git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: second handled case added,
   Right / Lower / Transpose, which is Rpotrf's block-column panel solve. It is parallelised
   over the m ROWS of B, not the columns -- see "Why rows for Right" below. The Left case and
   its gate constant are untouched: the same operands take the same branch and produce the same
   bits as before this change. See git log. */

/*  Why a separate entry point rather than a pragma inside Rtrsm.cpp:
    Rtrsm is also called from Rpotrf2's recursion and from Rpotrf's UPPER path, at shapes for
    which nothing has been measured. Threading the generic kernel would thread those too.
    Exactly three call sites reach this file's symbol: Lal::getInvLowTriangularMatrix (Left) and
    Rpotrf's lower blocked panel (Right).

    Why columns for Left:
    for side == "Left" the outermost loop of every branch of Rtrsm runs over the n columns of B,
    and every read and write inside it is b[.. + (j-1)*ldb] for that one j; a is read-only. So
    the columns of B are independent and the arithmetic within a column is untouched by the
    split -- which is what makes the result bit-identical to the serial kernel at any thread
    count.

    Why ROWS for Right -- the axis flips, and getting it wrong is silent:
    the Right/Lower/Transpose branch of Rtrsm (Rtrsm.cpp:263-285) is

        for (k = 1; k <= n; k++) {
            if (nounit) { temp = one / a[(k-1)+(k-1)*lda];
                          for (i=1..m) b[(i-1)+(k-1)*ldb] = temp * b[(i-1)+(k-1)*ldb]; }
            for (j = k+1; j <= n; j++)
                if (a[(j-1)+(k-1)*lda] != zero) {
                    temp = a[(j-1)+(k-1)*lda];
                    for (i=1..m) b[(i-1)+(j-1)*ldb] = b[(i-1)+(j-1)*ldb] - temp*b[(i-1)+(k-1)*ldb];
                }                                             ^^^^^^^^^^^^^^^^^^^ column k, not j
            if (alpha != one) for (i=1..m) b[(i-1)+(k-1)*ldb] = alpha * b[(i-1)+(k-1)*ldb];
        }

    Column j is updated FROM column k for every k < j, so the columns are a dependent chain --
    splitting Right over columns compiles, converges, and returns wrong answers. The index i, by
    contrast, appears only as "for all i" applied elementwise: no operation in the branch mixes
    two different values of i, and a is read-only. So for side == "Right" it is the m rows of B
    that are independent, and a contiguous row block carries its whole k/j chain in order. Hence
    bit-identity at any thread count, same as the Left case.

    Schedule for Right: every row costs exactly the same, n*(n+1)/2 multiply-adds, so an even
    contiguous partition is already perfectly balanced and there is nothing for dynamic
    scheduling to fix. Blocked rather than one-row-per-iteration so that the innermost loop stays
    long and unit-stride (B is column-major) and the a[] zero-tests are not repeated m times. */

#include <mpblas_dd.h>
#include "mplapack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Rtrsm_omp(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb) {
    const dd_real zero = 0.0;
    //
    //     The two cases this file implements. Each test is positive (Mlsame against the
    //     wanted letter), never "not the other one", so a malformed argument string falls
    //     through to the serial kernel instead of being silently treated as one of them.
    //
    //       left  : B := alpha*inv(A)*B,      A lower triangular, not transposed  (column-split)
    //       right : B := alpha*B*inv(A**T),   A lower triangular                  (ROW-split)
    //
    //     Rtrsm treats "T" and "C" identically for a real type, so both are accepted.
    //
    bool left_lower_notrans = Mlsame_dd(side, "L") && Mlsame_dd(uplo, "L") && Mlsame_dd(transa, "N") && (Mlsame_dd(diag, "N") || Mlsame_dd(diag, "U"));
    bool right_lower_trans = Mlsame_dd(side, "R") && Mlsame_dd(uplo, "L") && (Mlsame_dd(transa, "T") || Mlsame_dd(transa, "C")) && (Mlsame_dd(diag, "N") || Mlsame_dd(diag, "U"));
    //
    //     Argument errors and the quick-return/alpha==0 paths are left to Rtrsm so that
    //     Mxerbla reporting lives in exactly one place. Rtrsm's nrowa is m for side=="L" and
    //     n for side=="R", so the lda test differs between the two cases; do not merge them.
    //
    bool valid = (m > 0) && (n > 0) && (ldb >= std::max((mplapackint)1, m)) && (alpha != zero);
    bool parallel_left = left_lower_notrans && valid && (lda >= std::max((mplapackint)1, m));
    bool parallel_right = right_lower_trans && valid && (lda >= std::max((mplapackint)1, n));
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    bool have_team = (nthreads > 1) && !omp_in_parallel();
    //
    //     Left: unchanged from the calibrated B3 form -- same work expression, same constant,
    //     same width axis. The same operands take the same branch as before the Right case
    //     existed.
    //
    parallel_left = parallel_left && have_team && ((double)m * (double)m * (double)n >= MPLAPACK_OMP_TRI_WORK(MPLAPACK_OMP_MIN_TRSM_WORK, nthreads)) && (n >= MPLAPACK_OMP_MIN_TRI_WIDTH);
    //
    //     Right: the work is m rows times n*(n+1)/2 multiply-adds each, i.e. m*n*n to within a
    //     factor of two -- the same convention the Left gate uses. The width axis is m, because
    //     m is what gets split. Both differ from Left, and so does the call site (Rpotrf's
    //     panel, not the Cholesky inverse), so it gets its own constant rather than borrowing
    //     MIN_TRSM_WORK.
    //
    parallel_right = parallel_right && have_team && ((double)m * (double)n * (double)n >= MPLAPACK_OMP_POTRF_WORK(MPLAPACK_OMP_MIN_POTRF_TRSM_WORK, nthreads)) && (m >= MPLAPACK_OMP_MIN_POTRF_WIDTH);
#else
    /* Without OpenMP there is nothing to gate; everything goes to the serial kernel. */
    parallel_left = false;
    parallel_right = false;
#endif
    if (!parallel_left && !parallel_right) {
        Rtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
        return;
    }
#ifdef _OPENMP
    const dd_real one = 1.0;
    bool nounit = Mlsame_dd(diag, "N");
    if (parallel_right) {
        //
        //           Form  B := alpha*B*inv( A**T ).   (right, lower, transpose)
        //
        //     Rows are split into one contiguous block per thread; every row costs the same, so
        //     the even partition below is exactly balanced and no dynamic schedule is needed.
        //     (t*m)/nblk is used rather than a ceil-divided chunk width so that no block comes
        //     out empty and the last one is not left oversized. nblk <= m by construction.
        //
        //     The body is Rtrsm.cpp:264-284 verbatim with `1 <= i <= m` narrowed to the block's
        //     row range; keep it a transcription, including `b = b - temp*b` rather than `-=`.
        //
        mplapackint nblk = (mplapackint)nthreads;
        if (nblk > m) {
            nblk = m;
        }
#pragma omp parallel for schedule(static, 1) num_threads((int)nblk)
        for (mplapackint t = 0; t < nblk; t = t + 1) {
            mplapackint i0 = 1 + (t * m) / nblk;
            mplapackint i1 = ((t + 1) * m) / nblk;
            mplapackint i = 0;
            mplapackint j = 0;
            mplapackint k = 0;
            dd_real temp = 0.0;
            for (k = 1; k <= n; k = k + 1) {
                if (nounit) {
                    temp = one / a[(k - 1) + (k - 1) * lda];
                    for (i = i0; i <= i1; i = i + 1) {
                        b[(i - 1) + (k - 1) * ldb] = temp * b[(i - 1) + (k - 1) * ldb];
                    }
                }
                for (j = k + 1; j <= n; j = j + 1) {
                    if (a[(j - 1) + (k - 1) * lda] != zero) {
                        temp = a[(j - 1) + (k - 1) * lda];
                        for (i = i0; i <= i1; i = i + 1) {
                            b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - temp * b[(i - 1) + (k - 1) * ldb];
                        }
                    }
                }
                if (alpha != one) {
                    for (i = i0; i <= i1; i = i + 1) {
                        b[(i - 1) + (k - 1) * ldb] = alpha * b[(i - 1) + (k - 1) * ldb];
                    }
                }
            }
        }
        //
        //     End of Rtrsm_omp (right).
        //
        return;
    }
    //
    //           Form  B := alpha*inv( A )*B.   (lower, no transpose)
    //
    //     schedule(dynamic,1), not the default static split: against the identity RHS this
    //     call is given, column j costs O((m-j)^2) because the leading zeros of the identity
    //     are skipped, so a contiguous block split hands the first thread roughly a third of
    //     the total work. Measured on thanos (8 threads): contiguous 2.57-3.10x vs
    //     dynamic,1 3.9-6.5x.
    //
    //     num_threads() is clamped to n so a narrow B never creates idle team members.
    //
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads < (int)n ? nthreads : (int)n)
    for (mplapackint j = 1; j <= n; j = j + 1) {
        mplapackint i = 0;
        mplapackint k = 0;
        if (alpha != one) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = alpha * b[(i - 1) + (j - 1) * ldb];
            }
        }
        for (k = 1; k <= m; k = k + 1) {
            if (b[(k - 1) + (j - 1) * ldb] != zero) {
                if (nounit) {
                    b[(k - 1) + (j - 1) * ldb] = b[(k - 1) + (j - 1) * ldb] / a[(k - 1) + (k - 1) * lda];
                }
                for (i = k + 1; i <= m; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[(k - 1) + (j - 1) * ldb] * a[(i - 1) + (k - 1) * lda];
                }
            }
        }
    }
#endif
    //
    //     End of Rtrsm_omp .
    //
}
