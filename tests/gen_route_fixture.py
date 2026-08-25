#!/usr/bin/env python3
"""Emit small problems that sit ON the bMat route-chooser's boundaries.

Finding L. The chooser has four gates and two policies. Since 2026-08-24 the DEFAULT is the
0.40 fill-derived policy -- selected by unset, `auto` or `fill`, which are synonyms -- and
`legacy` is the released 0.25 aggregate policy. (When this file was written the polarity was the
other way round: `auto` was 0.25 and `fill` was the opt-in one.)

Until this fixture existed, the only route assertion in CI compared example1's PRINTED OUTPUT
between unset and `auto`. That cannot detect a policy being switched underneath the default,
because on a problem where both routes reach the same answer the output is identical by design.
The route itself has to be the thing asserted, and it is observable: SDPA_BMAT_LOG=1 prints
`bmat: gateN ... -> DENSE|SPARSE`.

Cases, chosen so each one straddles a gate rather than sitting comfortably inside it. The gate
comparisons are strict `>`, so the equality cases are the interesting half:

  m100    m = 100 exactly    gate 1 wants m > 100, so 100 must NOT pass
  m101    m = 101            the first m that does
  nb5     nBlock = 5         gate 1 wants nBlock > 5, so 5 must NOT pass
  nb6     nBlock = 6         the first nBlock that does
  wide    one block holds EXACTLY half the constraints -- gate 2 wants > 0.5*m, so this must NOT
          trip it. It is also THE CASE THAT MATTERS: `legacy` then stops at gate 3 (aggregate
          11,430 over its 10,000 cutoff) while the current default carries on to gate 4 and
          routes SPARSE. Without a fixture where the two policies disagree, "unset == auto" is
          satisfied by every possible route table and proves nothing.
  wider   one more constraint in that block, which DOES trip gate 2
  dense   a problem dense enough that both policies decline it at gate 2
  band    a banded problem whose ordered fill is small while its aggregate density is not. It
          does NOT diverge -- both policies route it DENSE -- and that is the point: they reach
          the same verdict through DIFFERENT cutoffs on the same aggregate (49,500 against the
          default's 36,000 and `legacy`'s 22,500), which pins the cutoff wiring where the route
          alone would not distinguish them.

Usage: gen_route_fixture.py <case>
"""
import sys

def emit(m, blocks, touch):
    """touch(k) -> list of (block_index0, row1, col1, value) for constraint k (1-based)."""
    print(m)
    print(len(blocks))
    print(" ".join(str(b) for b in blocks))
    print(" ".join("1.0" for _ in range(m)))
    for b, n in enumerate(blocks):
        for i in range(1, n + 1):
            print("0 %d %d %d 1.0" % (b + 1, i, i))
    for k in range(1, m + 1):
        for (b, r, c, v) in touch(k):
            print("%d %d %d %d %s" % (k, b + 1, r, c, v))

def block_diagonal(m, nblock, bsize):
    per = m // nblock
    def touch(k):
        b = min((k - 1) // per, nblock - 1)
        d = ((k - 1) % bsize) + 1
        return [(b, d, d, "1.0")]
    return touch

def banded(m, nblock, bsize, span):
    """Constraint k touches `span` consecutive blocks: fill stays low, aggregate does not."""
    per = max(m // nblock, 1)
    def touch(k):
        b0 = min((k - 1) // per, nblock - 1)
        out = []
        for t in range(span):
            b = (b0 + t) % nblock
            d = (((k - 1) + t) % bsize) + 1
            out.append((b, d, d, "1.0"))
        return out
    return touch

CASES = {
    # gate 1 boundaries
    "m100":  lambda: emit(100, [8] * 12, block_diagonal(100, 12, 8)),
    "m101":  lambda: emit(101, [8] * 12, block_diagonal(101, 12, 8)),
    "nb5":   lambda: emit(200, [8] * 5,  block_diagonal(200, 5, 8)),
    "nb6":   lambda: emit(200, [8] * 6,  block_diagonal(200, 6, 8)),
    # gate 2 boundary: one block carries exactly half the constraints, then half plus one
    "wide":  lambda: emit(200, [12] * 8, _half(200, 8, 12, 0)),
    "wider": lambda: emit(200, [12] * 8, _half(200, 8, 12, 1)),
    # policy divergence -- `wide` is the case where legacy and the current default disagree
    "dense": lambda: emit(160, [10] * 8, banded(160, 8, 10, 8)),
    "band":  lambda: emit(300, [10] * 20, banded(300, 20, 10, 6)),
}

def _half(m, nblock, bsize, extra):
    big = m // 2 + extra
    def touch(k):
        d = ((k - 1) % bsize) + 1
        if k <= big:
            return [(0, d, d, "1.0")]
        b = 1 + ((k - big - 1) % (nblock - 1))
        return [(b, d, d, "1.0")]
    return touch

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in CASES:
        sys.stderr.write("usage: %s {%s}\n" % (sys.argv[0], ",".join(sorted(CASES))))
        sys.exit(2)
    CASES[sys.argv[1]]()
