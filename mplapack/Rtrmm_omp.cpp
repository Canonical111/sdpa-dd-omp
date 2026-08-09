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

/* MODIFIED from upstream (BSD 2-clause; this file is BSD-licensed MPLAPACK, not GPL -- the original copyright notice above is retained and the change is recorded here), 2026-08-04: new file. Column-parallel Rtrmm for the
   single case the solver uses (Left / Lower / Transpose); every other case is handed to the
   serial Rtrmm. The parallel body is a verbatim copy of that branch of Rtrmm.cpp, including the
   zero-skip this fork added there on 2026-08-03 -- the two must be kept in step. See git log. */
/* MODIFIED from upstream (BSD 2-clause; this file is BSD-licensed MPLAPACK, not GPL -- the original copyright notice above is retained and the change is recorded here), 2026-08-04: work gate calibrated on thanos and pi;
   MPLAPACK_OMP_MIN_TRMM_WORK_8T renamed to MPLAPACK_OMP_MIN_TRMM_WORK because the measured
   break-even does not track the team size. Measured crossover is n = 10-12 (thanos) and
   n = 12-16 (pi), not the n >= 45 previously assumed. See mplapack_omp_tuning.h and git log. */

/*  Why a separate entry point rather than a pragma inside Rtrmm.cpp:
    Rtrmm is also called ~30 times by Rlarfb, always with side == "Right" and on the small k-wide
    workspace of the QR/tridiagonal reductions. Threading the generic kernel would thread those
    too, with no measurement behind them. Only Jal::getInvCholAndInv calls this file's symbol.

    Why columns:
    for side == "Left" the outermost loop of every branch of Rtrmm runs over the n columns of B,
    and every read and write inside it is b[.. + (j-1)*ldb] for that one j; a is read-only. So
    the columns of B are independent and the arithmetic within a column is untouched by the
    split -- which is what makes the result bit-identical to the serial kernel at any thread
    count. This is NOT true for side == "Right", where the "B*A" loops read b[..+(k-1)*ldb] from
    other columns and it is the m ROWS that are independent instead. This file implements Left
    only and refuses anything else, so that distinction cannot be got wrong by accident. */

#include <mpblas_dd.h>
#include "mplapack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Rtrmm_omp(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb) {
    const dd_real zero = 0.0;
    //
    //     The only case this file implements: B := alpha*A**T*B with A lower triangular,
    //     applied from the left. Each test is positive (Mlsame against the wanted letter),
    //     never "not the other one", so a malformed argument string falls through to the
    //     serial kernel instead of being silently treated as this case. "C" is accepted
    //     alongside "T" because Rtrmm itself treats them identically for a real matrix.
    //
    bool handled = Mlsame_dd(side, "L") && Mlsame_dd(uplo, "L") && (Mlsame_dd(transa, "T") || Mlsame_dd(transa, "C")) && (Mlsame_dd(diag, "N") || Mlsame_dd(diag, "U"));
    //
    //     Argument errors and the quick-return/alpha==0 paths are left to Rtrmm so that
    //     Mxerbla reporting lives in exactly one place.
    //
    bool valid = (m > 0) && (n > 0) && (lda >= std::max((mplapackint)1, m)) && (ldb >= std::max((mplapackint)1, m)) && (alpha != zero);
    //
    bool parallel = handled && valid;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    parallel = parallel && (nthreads > 1) && !omp_in_parallel() && ((double)m * (double)m * (double)n >= MPLAPACK_OMP_TRI_WORK(MPLAPACK_OMP_MIN_TRMM_WORK, nthreads)) && (n >= MPLAPACK_OMP_MIN_TRI_WIDTH);
#else
    /* Without OpenMP there is nothing to gate; everything goes to the serial kernel. */
    parallel = false;
#endif
    if (!parallel) {
        Rtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
        return;
    }
#ifdef _OPENMP
    bool nounit = Mlsame_dd(diag, "N");
    //
    //           Form  B := alpha*A**T*B.   (lower)
    //
    //     schedule(dynamic,1) for the same reason as Rtrsm_omp: the B this call is given is
    //     inv(L), lower triangular, so the zero-skip makes column j cost O((m-j)^2) and a
    //     contiguous block split is badly unbalanced. Measured on thanos (8 threads):
    //     contiguous 2.57-3.10x vs dynamic,1 3.9-6.5x.
    //
    //     num_threads() is clamped to n so a narrow B never creates idle team members.
    //
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads < (int)n ? nthreads : (int)n)
    for (mplapackint j = 1; j <= n; j = j + 1) {
        mplapackint i = 0;
        mplapackint k = 0;
        dd_real temp = 0.0;
        for (i = 1; i <= m; i = i + 1) {
            temp = b[(i - 1) + (j - 1) * ldb];
            if (nounit) {
                temp = temp * a[(i - 1) + (i - 1) * lda];
            }
            for (k = i + 1; k <= m; k = k + 1) {
                if (b[(k - 1) + (j - 1) * ldb] != zero) {
                    temp += a[(k - 1) + (i - 1) * lda] * b[(k - 1) + (j - 1) * ldb];
                }
            }
            b[(i - 1) + (j - 1) * ldb] = alpha * temp;
        }
    }
#endif
    //
    //     End of Rtrmm_omp .
    //
}
