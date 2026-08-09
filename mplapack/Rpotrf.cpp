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

#include <mpblas_dd.h>
#include <mplapack_dd.h>

void Rpotrf(const char *uplo, mplapackint const n, dd_real *a, mplapackint const lda, mplapackint &info) {
    bool upper = false;
    mplapackint nb = 0;
    mplapackint j = 0;
    mplapackint jb = 0;
    const dd_real one = 1.0;
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    upper = Mlsame_dd(uplo, "U");
    if (!upper && !Mlsame_dd(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla_dd("Rpotrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Determine the block size for this environment.
    //
    nb = iMlaenv_dd(1, "Rpotrf", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code.
        //
        Rpotrf2(uplo, n, a, lda, info);
    } else {
        //
        //        Use blocked code.
        //
        if (upper) {
            //
            //           Compute the Cholesky factorization A = U**T*U.
            //
            for (j = 1; j <= n; j = j + nb) {
                //
                //              Update and factorize the current diagonal block and test
                //              for non-positive-definiteness.
                //
                jb = std::min(nb, n - j + 1);
                Rsyrk("Upper", "Transpose", jb, j - 1, -one, &a[(j - 1) * lda], lda, one, &a[(j - 1) + (j - 1) * lda], lda);
                Rpotrf2("Upper", jb, &a[(j - 1) + (j - 1) * lda], lda, info);
                if (info != 0) {
                    goto statement_30;
                }
                if (j + jb <= n) {
                    //
                    //                 Compute the current block row.
                    //
                    Rgemm("Transpose", "No transpose", jb, n - j - jb + 1, j - 1, -one, &a[(j - 1) * lda], lda, &a[((j + jb) - 1) * lda], lda, one, &a[(j - 1) + ((j + jb) - 1) * lda], lda);
                    Rtrsm("Left", "Upper", "Transpose", "Non-unit", jb, n - j - jb + 1, one, &a[(j - 1) + (j - 1) * lda], lda, &a[(j - 1) + ((j + jb) - 1) * lda], lda);
                }
            }
            //
        } else {
            //
            //           Compute the Cholesky factorization A = L*L**T.
            //
            for (j = 1; j <= n; j = j + nb) {
                //
                //              Update and factorize the current diagonal block and test
                //              for non-positive-definiteness.
                //
                jb = std::min(nb, n - j + 1);
                /* MODIFIED from upstream (BSD 2-clause; this file is BSD-licensed MPLAPACK, not GPL -- the original copyright notice above is retained and the change is recorded here), 2026-08-04: Rsyrk -> Rsyrk_omp and
                   Rtrsm -> Rtrsm_omp in the LOWER path only. These are the two panel kernels
                   that were still serial here; the Rgemm below already threads via
                   Rgemm_NT_omp. Rsyrk_omp splits the jb output columns, Rtrsm_omp splits the
                   n-j-jb+1 ROWS of B (Right side -- the columns are a dependent chain there);
                   both leave the arithmetic of each independent unit and its order untouched,
                   so the factorisation is bit-identical to the serial one at any thread count.
                   Each falls back to its serial kernel below its work gate and for every case
                   it does not implement. The UPPER path above is deliberately left alone: the
                   solver never asks for it, so it has no measurement behind it. See git log. */
                Rsyrk_omp("Lower", "No transpose", jb, j - 1, -one, &a[(j - 1)], lda, one, &a[(j - 1) + (j - 1) * lda], lda);
                Rpotrf2("Lower", jb, &a[(j - 1) + (j - 1) * lda], lda, info);
                if (info != 0) {
                    goto statement_30;
                }
                if (j + jb <= n) {
                    //
                    //                 Compute the current block column.
                    //
                    Rgemm("No transpose", "Transpose", n - j - jb + 1, jb, j - 1, -one, &a[((j + jb) - 1)], lda, &a[(j - 1)], lda, one, &a[((j + jb) - 1) + (j - 1) * lda], lda);
                    Rtrsm_omp("Right", "Lower", "Transpose", "Non-unit", n - j - jb + 1, jb, one, &a[(j - 1) + (j - 1) * lda], lda, &a[((j + jb) - 1) + (j - 1) * lda], lda);
                }
            }
        }
    }
    goto statement_40;
//
statement_30:
    info += j - 1;
//
statement_40:;
    //
    //     End of Rpotrf
    //
}
