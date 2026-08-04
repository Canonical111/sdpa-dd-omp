/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_dd.h,v 1.10 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_DD_H_
#define _MPBLAS_DD_H_

#include <qd/dd_real.h>
#include <mplapack_config.h>
#include <mplapack_utils_dd.h>

dd_real Rdot(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Rcopy(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Raxpy(mplapackint const n, dd_real const da, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Rscal(mplapackint const n, dd_real const da, dd_real *dx, mplapackint const incx);
bool Mlsame_dd(const char *a, const char *b);
void Mxerbla_dd(const char *srname, int info);
void Rswap(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
dd_real Rnrm2(mplapackint const n, dd_real *x, mplapackint const incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rsymv(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rsyr2(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *y, mplapackint const incy, dd_real *a, mplapackint const lda);
void Rger(mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *y, mplapackint const incy, dd_real *a, mplapackint const lda);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: declared the two column-parallel
   Left-side triangular kernels. See git log.
   These are NOT drop-in replacements for Rtrmm/Rtrsm at every call site by policy: they exist
   so that the Cholesky-inverse phase can be threaded WITHOUT also threading Rtrsm inside
   Rpotrf or Rtrmm inside Rlarfb. Call them only from Lal::getInvLowTriangularMatrix and
   Jal::getInvCholAndInv. Each parallelises the Left-side case it implements over the columns
   of B and delegates every other case, and every sub-threshold call, to the serial kernel;
   results are bit-identical to it at any thread count. */
void Rtrmm_omp(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
void Rtrsm_omp(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real const beta, dd_real *c, mplapackint const ldc);

#endif
