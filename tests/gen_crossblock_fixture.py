#!/usr/bin/env python3
"""Emit a sparse SDPA problem whose Schur entries are shared BETWEEN blocks.

This fixture exists because tests/gen_spchol_fixture.py structurally cannot fail the test that
matters here. That one gives every constraint exactly one block, so each Schur entry is written
from exactly one block, and the threaded assembly's cross-block accumulation is never exercised
-- which is why a real data race in the released assembly (a `nowait` that let one worker enter
block l+1 while another was still in block l, both accumulating into the same sparse_bMat slot)
passed every stream-identity check the fork had.

So here each constraint touches TWO blocks. A pair of constraints that co-occur in two different
blocks then accumulates into the same Schur entry from both, which is the dependence the barrier
protects.

The gates, in the order they apply:

  gate 1  m > 100 and nBlock > 5                      -> m = 480, nBlock = 24
  gate 2  no block's constraint count exceeds 0.5*m   -> 40 per block, against 240
  gate 4  ordered fill below the cutoff               -> routes SPARSE under the default `auto`
  assembly gate  pairs >= SDPA_BMAT_ASM_MIN_PAIRS     -> 24 * 820 = 19,680, against 8,000

It routes sparse and threads its assembly with NO environment variables set, so the regression it
guards is on the default path rather than behind a forced mode.

THE STRIDE MUST BE 1, and the first version of this fixture got that wrong. Sharing a Schur
entry between two blocks is necessary but not sufficient: the two blocks must also be IN FLIGHT
AT THE SAME TIME for the race to land. Blocks are walked in order, so with a stride of 7 the
sharing blocks are seven apart and the workers -- who stay within about one block of each other
-- never overlap them. That fixture ran the released, racing binary six times and got six
identical answers. With adjacent blocks a single fast worker rolling into block b+1 while a slow
one finishes block b is enough, and the same binary gives six different answers.

Generated rather than committed so its structure is auditable arithmetic rather than a blob.
"""
NBLOCK = 24
BSIZE = 9             # 9*10/2 = 45 distinct symmetric entries, and a block needs 40 of them
STRIDE = 1            # ADJACENT blocks: see the note below on why this must be 1
PER_BLOCK_FIRST = 20  # -> 40 constraints per block (20 as first block, 20 as second)

# Distinct entries, diagonal first. Each constraint gets its OWN entry in each block it touches,
# so no two A_i are equal and the Schur complement stays nonsingular. The first version of this
# fixture reused entries and produced a singular Schur complement that failed Cholesky at
# iteration 0 -- a fixture that cannot complete an iteration cannot test the assembly.
ENTRIES = [(d, d) for d in range(1, BSIZE + 1)]
ENTRIES += [(r, c) for c in range(1, BSIZE + 1) for r in range(1, c)]
assert len(ENTRIES) >= 2 * PER_BLOCK_FIRST, "block too small to give every constraint its own entry"

m = NBLOCK * PER_BLOCK_FIRST   # 480
print(m)
print(NBLOCK)
print(" ".join(str(BSIZE) for _ in range(NBLOCK)))
print(" ".join("1.0" for _ in range(m)))

# F0 = identity on every block: strictly feasible at X = I.
for b in range(1, NBLOCK + 1):
    for i in range(1, BSIZE + 1):
        print("0 %d %d %d 1.0" % (b, i, i))

# Constraint k lives in b0 = k // PER_BLOCK_FIRST and in b1 = (b0 + STRIDE) % NBLOCK. In its
# FIRST block it occupies entry j; in its SECOND it occupies entry PER_BLOCK_FIRST + j. A block
# therefore sees entries 0..19 from its "first" constraints and 20..39 from its "second" ones --
# 40 distinct patterns, no duplicates.
#
# Two constraints that co-occur in two different blocks contribute to the SAME Schur entry from
# both blocks. That cross-block accumulation is the whole point of this fixture.
k = 0
for b0 in range(NBLOCK):
    for j in range(PER_BLOCK_FIRST):
        k += 1
        b1 = (b0 + STRIDE) % NBLOCK
        r0, c0 = ENTRIES[j]
        r1, c1 = ENTRIES[PER_BLOCK_FIRST + j]
        print("%d %d %d %d 1.0" % (k, b0 + 1, r0, c0))
        print("%d %d %d %d 1.0" % (k, b1 + 1, r1, c1))
