/*
 * Copyright (c) 2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rdot.cpp,v 1.5 2010/08/07 05:50:09 nakatamaho Exp $
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *
 * $Id: Rdot.cpp,v 1.5 2010/08/07 05:50:09 nakatamaho Exp $

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
Based on http://www.netlib.org/blas/ddot.f
Rdot forms the dot product of two vectors.
*/

/* MODIFIED from upstream (BSD 2-clause; this file is BSD-licensed MPLAPACK, not GPL -- the original copyright notice above is retained and the change is recorded here), 2026-07-31: Rdot runs serially in the original summation order. See git log. */
#include <mpblas_dd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

dd_real Rdot_omp(mplapackint n, dd_real *dx, mplapackint incx, dd_real *dy, mplapackint incy) {
    mplapackint ix = 0;
    mplapackint iy = 0;
    mplapackint i;
    dd_real temp;

    temp = 0.0;

    if (incx < 0)
        ix = (-n + 1) * incx;
    if (incy < 0)
        iy = (-n + 1) * incy;

    temp = 0.0;
    if (incx == 1 && incy == 1) {
// no reduction for multiple precision
        /* SERIAL, in the original sequential order, bit for bit.
           Upstream combined per-thread partial sums under a critical section, so the
           summation order followed whichever thread won the lock. In double-double that
           perturbs low bits and steers the trajectory: theta2 at 24 threads gave
           51/65/80/56/50 iterations across 5 runs.

           A fixed-chunk parallel reduction was tried as the fix and REJECTED on measurement:
           over 8 problems it was never faster, and it pushed arch4 onto a 40%-longer
           trajectory (101 vs 72 iterations). Running dot serially is both faster and a
           stronger guarantee -- the result matches a serial run exactly, not merely across
           thread counts. The rejected experiment is kept, with its data, in
           patches/experiments/rdot_chunked.py; it is not part of the shipped series. */
        for (i = 0; i < n; i++) {
            temp += dx[i] * dy[i];
        }
    } else {
        for (i = 0; i < n; i++) {
            temp += dx[ix] * dy[iy];
            ix = ix + incx;
            iy = iy + incy;
        }
    }
    return temp;
}
