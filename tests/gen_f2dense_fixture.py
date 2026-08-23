#!/usr/bin/env python3
"""Generate an SDPA fixture that REACHES the hasF2Gcal consumption path.

The path is taken only when, for a pair (i, j) in one i-group:
  * useFormula[i] == F2, and
  * A[j]'s block is stored DENSE, and
  * the pair is not the first such pair in the group.

So the fixture needs constraint matrices dense enough that SparseMatrix::changeToDense
converts them (>= 20% of n*n entries), several constraints in one block so a group has more
than one pair, and a Kappa small enough that the cost heuristic prefers F2 over F1 --
F2 wins when Kappa*up <= n^2, which is why this is exercised with -k 0.3.

Every SDPLIB and bootstrap problem in dd's benchkit reports f2dense_total = 0, so without a
purpose-built fixture the branch cannot be reached at all. Verify with SDPA_BMAT_ASM_CENSUS=1.
"""
import sys

n = 6      # block order: n*n = 36 entries, small so n^2 stays under Kappa*up
m = 5      # constraints: group i then holds (m - i + 1) pairs

# MIXED sparsity is the point. Pairs within group i are (i,i), (i,i+1), ..., so the FIRST pair
# of the group has A_j = A_i. If A_i is SPARSE, calF2 takes its sparse branch, which never
# touches hasF2Gcal -- so G = X*F is NEVER computed for that group. A LATER pair with a dense
# A_j then reads the indeterminate flag; if it reads true it skips computing G and contracts
# against whatever work2 still holds from a previous group. That is the harmful shape, and an
# all-dense fixture cannot produce it: there the first pair always computes G, so both readings
# of the indeterminate flag happen to agree.
SPARSE_CONSTRAINTS = {1}   # A_1 sparse, A_2..A_m dense

print(m)
print(1)
print(n)
print(" ".join("1.0" for _ in range(m)))

def dense_block(scale, shift):
    """A fully dense symmetric upper triangle -- forces changeToDense."""
    out = []
    for r in range(1, n + 1):
        for c in range(r, n + 1):
            v = scale * (1.0 + 0.1 * ((r * 7 + c * 3 + shift) % 5))
            if r == c:
                v += 2.0 * scale        # diagonally dominant: stays positive definite
            out.append((r, c, v))
    return out

# F0 (the constant term) must make the initial point feasible; a scaled identity-plus-dense
# keeps every A_k in range while remaining positive definite.
for r, c, v in dense_block(1.0, 0):
    print(f"0 1 {r} {c} {v:.15e}")
def sparse_block(scale, shift):
    """Few enough entries that changeToDense leaves it SPARSE (< 20% of n*n)."""
    out = [(r, r, 2.0 * scale) for r in range(1, n + 1)]     # diagonal only: n of n*n entries
    out.append((1, 2, 0.5 * scale * (1 + 0.1 * shift)))
    return out

for k in range(1, m + 1):
    blk = sparse_block(0.1, k) if k in SPARSE_CONSTRAINTS else dense_block(0.1, k)
    for r, c, v in blk:
        print(f"{k} 1 {r} {c} {v:.15e}")
