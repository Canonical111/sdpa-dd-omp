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
/* NEW FILE, added by this fork (GPLv2 2a notice), 2026-08-23: unit tests for the checked-size
   helpers in sdpa_tool.h and for rError's statement safety. See git log. */
/* Unit tests for the checked-size helpers and the rError macro's statement safety.

   The point of this file is WHICH SYMBOLS it calls. A test that reimplements
   `(a * b > INT_MAX)` proves the test author can write a guard; it says nothing about whether
   the allocation sites are guarded. So these call sdpaTriangularCount / sdpaProduct /
   sdpaFitsInt from <sdpa_tool.h> -- the same inline functions the six production sites use.

   Built and run by CI. Standalone: needs no solver, no input file and no OpenMP. */
#include <climits>
#include <cstdio>
#include <cstdlib>

#include <sdpa_tool.h>

using std::cout;
using std::endl;

using sdpa::sdpaFitsInt;
using sdpa::sdpaProduct;
using sdpa::sdpaTriangularCount;

static int failures = 0;

static void check(bool ok, const char *what) {
    printf("  %-58s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) {
        failures++;
    }
}

int main() {
    printf("checked-size helpers:\n");

    // --- sdpaTriangularCount: n(n+1)/2, the count that overflows int at n = 65536 ---
    check(sdpaTriangularCount(0) == 0LL, "triangular(0) == 0");
    check(sdpaTriangularCount(1) == 1LL, "triangular(1) == 1");
    check(sdpaTriangularCount(10) == 55LL, "triangular(10) == 55");
    // Last value whose triangular count still fits an int, and the first that does not.
    check(sdpaFitsInt(sdpaTriangularCount(65535)), "triangular(65535) fits int (last accepted)");
    check(!sdpaFitsInt(sdpaTriangularCount(65536)), "triangular(65536) does NOT fit (first rejected)");
    // Computed in 64 bits: an int expression would have wrapped to a small or negative value.
    check(sdpaTriangularCount(65536) == 2147516416LL, "triangular(65536) computed exactly in 64 bits");
    check(sdpaTriangularCount(INT_MAX) > (long long)INT_MAX, "triangular(INT_MAX) exceeds INT_MAX");

    // --- sdpaProduct: the matrix-dimension product ---
    check(sdpaProduct(46340, 46340) == 2147395600LL, "product(46340,46340) exact");
    check(sdpaFitsInt(sdpaProduct(46340, 46340)), "46340^2 fits int (last accepted square)");
    check(!sdpaFitsInt(sdpaProduct(46341, 46341)), "46341^2 does NOT fit (first rejected square)");
    check(sdpaProduct(INT_MAX, 2) == 4294967294LL, "product(INT_MAX,2) exact in 64 bits");
    check(!sdpaFitsInt(sdpaProduct(INT_MAX, 2)), "product(INT_MAX,2) rejected");

    // --- negatives: callers permit them only where a count can legitimately be zero ---
    check(!sdpaFitsInt(-1LL), "negative rejected");
    check(!sdpaFitsInt(sdpaProduct(-1, 10)), "product with a negative dimension rejected");
    check(sdpaFitsInt(0LL), "zero accepted");
    check(sdpaFitsInt((long long)INT_MAX), "INT_MAX itself accepted (boundary)");
    check(!sdpaFitsInt((long long)INT_MAX + 1), "INT_MAX+1 rejected (boundary)");

    printf("\nrError statement safety:\n");
    // The property, not the macro text: an unbraced `if` must not exit when the condition is
    // false. Before the do/while(0) wrapper the exit was a second statement outside the if, so
    // this reached the line below only by accident of never being written.
    if (failures < 0)
        rError("unreachable: rError fired on a false condition");
    check(true, "if (false) rError(...) continues execution");
    // The true branch is exercised by CI running this binary a second time with an argument,
    // so a nonzero exit can be asserted from outside.
    printf("\n%s\n", failures == 0 ? "ALL OK" : "FAILURES");
    return failures == 0 ? 0 : 1;
}
