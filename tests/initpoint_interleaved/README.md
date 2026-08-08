# Interleaved LP/SDP initial-point regression (review2 finding 2)

Block structure `{-1, 2}`: an LP block FOLLOWED by an SDP block. The sparse
(`-is init_sparse.ini-s`) and dense (`-id init_dense.ini`) files describe the
IDENTICAL initial point, so the two runs must produce identical iteration
counts and objectives.

The pre-fix dense reader consumed every compacted SDP block first and the
flattened LP part afterwards -- the file, however, is written in the original
bLOCKsTRUCT order. With this structure the LP scalar was consumed as the first
SDP entry, the point read was not the point written, and the run either aborted
("initial point is not positive definite") or silently returned nothing
(iteration 0, noINFO, exit 0 -- observed on the pre-fix binary with exactly
these files).

CI runs both spellings and diffs iterations + objectives; any difference fails.
Note the dense flag is `-id`, not `-ii` -- unknown flags are silently ignored,
which is how the first version of this very test managed to test nothing.
