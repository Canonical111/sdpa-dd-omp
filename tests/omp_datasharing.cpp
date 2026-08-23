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

/* NEW FILE (GPLv2, matching the solver sources it guards), 2026-08-15: fixture for the
 * OpenMP data-sharing portability leg. Not compiled into the solver; built only by CI.
 * See review/GCC8-PORTABILITY-FIX-PLAN.md and git log. Ported verbatim from
 * sdpa-gmp-omp 2026-08-23: the fixture has no includes and names no numeric type, so it is
 * identical across the dd/gmp/qd forks by construction rather than by coincidence.
 *
 * Fixture for the OpenMP data-sharing portability leg in .github/workflows/build.yml.
 *
 * Deliberately has NO #include of any kind, so it compiles with nothing installed, in any
 * container, on any compiler that understands -fopenmp. That is the whole point: the leg that
 * needs an end-of-life compiler must not also need a working package manager.
 *
 * WHY A FIXTURE AND NOT JUST THE REAL SOURCES. The real sources are also checked, by whichever
 * leg can obtain GMP. But a leg that only ever compiles code that WORKS cannot tell you it is
 * still able to detect code that does not. This file compiles three ways, and the leg requires
 * the shipped spelling to succeed AND each broken spelling to be REJECTED by the compiler that
 * is supposed to reject it. If a future toolchain quietly stops enforcing one of these rules,
 * that shows up as a failed negative control rather than as silence.
 *
 * THE RULE BEING PINNED. Through OpenMP 4.5 a const-qualified variable with no mutable member
 * is PREDETERMINED shared, and naming it in a data-sharing clause is an error. OpenMP 5.0
 * deleted that rule, so under default(none) it MUST be named. GCC <= 8 implements the first,
 * GCC >= 9 and clang the second. The released solver shipped the spelling GCC 8 rejects, and
 * nothing caught it: see review/GCC8-PORTABILITY-FIX-PLAN.md.
 *
 *   SPELLING_BAD_SHARED_CONST   const named in shared()  -> must FAIL on GCC <= 8
 *   SPELLING_BAD_UNLISTED       const not named at all   -> must FAIL on clang
 *   (default)                   const in firstprivate()  -> must COMPILE everywhere
 *
 * `n` is initialised from a function call on purpose. A const initialised from a compile-time
 * constant is folded and never acquires a data-sharing attribute at all, so a fixture written
 * that way compiles under every spelling on every compiler and tests nothing -- which is
 * exactly what the first version of this reproducer did.
 */

int sink[4];
int runtime_value();

int omp_datasharing_fixture() {
    const int n = runtime_value();
#if defined(SPELLING_BAD_SHARED_CONST)
#pragma omp parallel default(none) shared(sink, n)
#elif defined(SPELLING_BAD_UNLISTED)
#pragma omp parallel default(none) shared(sink)
#else
#pragma omp parallel default(none) shared(sink) firstprivate(n)
#endif
    {
        sink[0] = n;
    }
    return sink[0];
}
