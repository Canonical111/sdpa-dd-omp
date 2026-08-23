#!/usr/bin/env python3
"""Emit a block-diagonal SDPA sparse problem that REACHES the sparse Schur Cholesky.

The bundled example1 is m=10, so it exits at gate 1 (m <= 100) and the sparse path is never
touched -- which is why CI needs its own fixture rather than reusing it. Requirements, in the
order the gates apply:

  gate 1  m > 100 and nBlock > 5
  gate 2  no block's constraint count exceeds 0.5*m
  gate 3  aggregate density below the cutoff
  gate 4  ordered fill below 0.40

A block-diagonal problem where each constraint touches exactly one small block satisfies all
four: the Schur complement inherits the block structure, so both the aggregate pattern and its
ordered fill stay far below any of the cutoffs, and no single block carries many constraints.

Deliberately NOT a physically meaningful SDP -- it exists to exercise a code path, and it is
generated rather than committed as data so its structure is auditable arithmetic rather than an
opaque blob.
"""
NBLOCK = 24      # > 5, and enough blocks that the Schur complement is genuinely sparse
BSIZE = 8        # block order: sets how wide a Schur pivot can get (see WIDTH below)
PER_BLOCK = 8    # constraints per block -> m = 192 > 100 (gate 1), and 8 <= 0.5*m (gate 2)

# WIDTH MATTERS, and the first version of this fixture got it wrong. With 40 blocks of 3
# constraints the Schur complement's widest pivot was 2, so the admission guard capped the team
# at 2 and -- as the per-worker counters then showed -- exactly ONE worker ever updated the
# factor, while `pivots_threaded` happily reported 120. A fixture that cannot occupy more than
# one worker cannot support a `workers_used >= 2` assertion, so the constraints per block are
# raised until the pivots are wide enough for several workers to share one.

m = NBLOCK * PER_BLOCK
print(m)
print(NBLOCK)
print(" ".join(str(BSIZE) for _ in range(NBLOCK)))
# c vector: every constraint has the same right-hand side
print(" ".join("1.0" for _ in range(m)))

# F0 = identity on every block, so the problem is strictly feasible at X = I.
for b in range(1, NBLOCK + 1):
    for i in range(1, BSIZE + 1):
        print("0 %d %d %d 1.0" % (b, i, i))

# Constraint k touches ONLY block (k // PER_BLOCK): one diagonal entry, plus one off-diagonal
# to give the block interior structure. Each constraint therefore contributes to exactly one
# block, which is what keeps the Schur complement block-diagonal and sparse.
k = 0
for b in range(1, NBLOCK + 1):
    for j in range(PER_BLOCK):
        k += 1
        d = (j % BSIZE) + 1
        print("%d %d %d %d 1.0" % (k, b, d, d))
        if d < BSIZE:
            print("%d %d %d %d 0.5" % (k, b, d, d + 1))
