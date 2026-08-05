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

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: new file. Column-parallel Rsyrk for the
   single case Rpotrf's lower path uses (Lower / NoTranspose); every other case is handed to the
   serial Rsyrk. See git log. */

/*  WHY A SEPARATE SYMBOL rather than a pragma inside Rsyrk.cpp:
    Rsyrk is also called from Rpotrf2 (mplapack/Rpotrf2.cpp:123,140), which is the RECURSIVE
    unblocked Cholesky -- it calls itself, so a pragma in the generic kernel would fork at every
    level of that recursion and for every diagonal block of every SDP block, shapes for which
    nothing has been measured. Rsyr2k/Rlarfb reach the generic kernel too. Only Rpotrf's lower
    blocked path calls this file's symbol.

    WHY COLUMNS:
    in the Lower/NoTranspose branch of Rsyrk (mplapack/Rsyrk.cpp:157-177) the outer loop runs
    over the n columns of C; iteration j reads only `a` (read-only for the whole call) and
    writes only c[j..n, j]. Two different j never touch the same element of C, and the order of
    the accumulation WITHIN a column -- the `l` loop, which is the reduction axis -- is
    untouched by the split. That is what makes the result bit-identical to serial Rsyrk at any
    thread count and any schedule.

    The `l` (=k) axis is the reduction. Splitting it would reorder a dd_real sum and change the
    last bits; this fork publishes a determinism guarantee (1-thread and 8-thread output are
    compared in CI), so k is never split.

    SCHEDULE: column j costs k*(n-j+1) multiply-adds, so the cost profile is a descending ramp.
    A contiguous static split hands thread 0 of 8 a fraction
        sum_{j<=n/8}(n-j+1) / (n(n+1)/2) -> (15/128)/(1/2) = 0.234
    of the total, capping the speedup at 4.27x however many threads are used. schedule(dynamic,1)
    removes that cap. It costs nothing in determinism: the columns are independent, so which
    thread runs which column cannot change a single bit.

    LATENT HAZARD -- READ BEFORE EDITING Rgemm_NT_omp.cpp.
    This file deliberately does NOT implement Rsyrk by delegating to Rgemm_NT_omp (which is
    already threaded and gated), for two reasons: Rgemm would write the strictly-upper triangle
    of C, which the caller's storage does not own, and it would do twice the arithmetic. But a
    future editor may be tempted, and the temptation is a trap worth recording, because the
    delegation is bit-identical ONLY while Rgemm_NT_omp keeps its present inner-loop form

        C[i + j * ldc] += temp * A[i + l * lda];        temp = alpha * B[j + l * ldb];

    which is elementwise identical to this file's

        c[(i-1) + (j-1)*ldc] += temp * a[(i-1) + (l-1)*lda];   temp = alpha * a[(j-1)+(l-1)*lda];

    with B == A. Netlib's OTHER gemm form -- accumulate TEMP from zero over l, then write
    C = alpha*TEMP + beta*C -- computes the same value in a different order and would NOT be
    bit-identical in dd_real. Rgemm_NT_omp's zero-skip and this file's zero-skip also test the
    same element, and that too is load-bearing. That kernel was already edited once by this fork
    (the A2 "restored netlib dgemm zero-skip" change), so this is a live risk, not a theoretical
    one -- and it is live TODAY regardless of this file, because Rpotrf's lower path calls
    Rgemm("No transpose","Transpose",...) at Rpotrf.cpp:139 and its output therefore depends on
    Rgemm_NT_omp's accumulation order. The invariant is pinned by the CI cell
    "kernel accumulation-form invariant" (.github/workflows/build.yml) driving
    tests/kernel_form_invariant.cpp, which fails if either kernel's form drifts. */

#include <mpblas_dd.h>
#include "mplapack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Rsyrk_omp(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real const beta, dd_real *c, mplapackint const ldc) {
    const dd_real zero = 0.0;
    //
    //     The only case this file implements: C := alpha*A*A**T + beta*C, lower triangle.
    //     Each test is positive (Mlsame against the wanted letter), never "not the other
    //     one", so a malformed argument string falls through to the serial kernel instead
    //     of being silently treated as this case.
    //
    bool handled = Mlsame_dd(uplo, "L") && Mlsame_dd(trans, "N");
    //
    //     Argument errors, the quick-return paths and the alpha==0 path are left to Rsyrk so
    //     that Mxerbla reporting lives in exactly one place. nrowa == n for trans=="N".
    //
    bool valid = (n > 0) && (k > 0) && (alpha != zero) && (lda >= std::max((mplapackint)1, n)) && (ldc >= std::max((mplapackint)1, n));
    //
    bool parallel = handled && valid;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    //
    //     Gate on WORK, not on one dimension: the panel Rsyrk of a Cholesky is short and fat
    //     (n = jb <= 64) but k grows to the whole factored prefix, so a gate on n alone would
    //     judge the first panel of a 1600x1600 factorisation the same as the last. The work
    //     expression k*n*n is the multiply-add count k*n*(n+1)/2 to within a factor of two,
    //     the same convention MPLAPACK_OMP_MIN_TRSM_WORK uses.
    //
    parallel = parallel && (nthreads > 1) && !omp_in_parallel() && ((double)k * (double)n * (double)n >= MPLAPACK_OMP_POTRF_WORK(MPLAPACK_OMP_MIN_POTRF_SYRK_WORK, nthreads)) && (n >= MPLAPACK_OMP_MIN_POTRF_WIDTH);
#else
    /* Without OpenMP there is nothing to gate; everything goes to the serial kernel. */
    parallel = false;
#endif
    if (!parallel) {
        Rsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
        return;
    }
#ifdef _OPENMP
    /* Declared here rather than beside `zero` above: without OpenMP the whole body below is
       compiled out and an unused `one` trips -Werror=unused-but-set-variable, which this
       project's CI enables. */
    const dd_real one = 1.0;
    //
    //           Form  C := alpha*A*A**T + beta*C.   (lower, no transpose)
    //
    //     num_threads() is clamped to n so a narrow panel never creates idle team members.
    //     The body below is Rsyrk.cpp:158-176 verbatim, including the `+=` and the zero-skip;
    //     any deviation would break bit-identity, so keep it a transcription.
    //
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads < (int)n ? nthreads : (int)n)
    for (mplapackint j = 1; j <= n; j = j + 1) {
        mplapackint i = 0;
        mplapackint l = 0;
        dd_real temp = 0.0;
        if (beta == zero) {
            for (i = j; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = zero;
            }
        } else if (beta != one) {
            for (i = j; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
            }
        }
        for (l = 1; l <= k; l = l + 1) {
            if (a[(j - 1) + (l - 1) * lda] != zero) {
                temp = alpha * a[(j - 1) + (l - 1) * lda];
                for (i = j; i <= n; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                }
            }
        }
    }
#endif
    //
    //     End of Rsyrk_omp .
    //
}
