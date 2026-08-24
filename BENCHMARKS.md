# Benchmarks — sdpa-dd-omp

Fork base: upstream `6eaad8d9`. Comparison points: **sdpa-dd 7.1.2** (the 2009 release,
single-threaded by construction) and **upstream master** built with OpenMP.

## Methodology

External wall-clock seconds, **median of 3 repeats**, spread reported where it exceeds
rounding. Runs are pinned to physical cores on Linux (`taskset` + `OMP_PROC_BIND=true
OMP_PLACES=cores`, CPU set recorded per row); macOS exposes no per-process affinity, so those
runs use OpenMP binding plus one unrecorded warmup. The parameter file is passed explicitly
and its SHA-256 recorded. Every repeat is a row in the raw TSVs published alongside this
document; a run is counted only if it exited cleanly and every field parsed. Iteration counts
and objectives are checked across repeats and across configurations -- a wall-time ratio
between builds with different iteration counts measures path length, not speed, and is
flagged.

Solver-internal timers are elapsed time in sdpa-dd/sdpa-gmp; upstream sdpa-qd's timer
reports process CPU time (summed over threads), which this fork fixes -- all qd numbers here
use external wall time and the corrected clock.

## pi — i9-13900K (8P+16E, 24 physical cores), Ubuntu 24.04, gcc 13

## pi-i9-13900k — dd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

## Thread scaling against upstream — 2026-08-23

Upstream `nakatamaho/sdpa-dd` at **`6eaad8d`**, the exact commit this fork branched from and still
its `master`, rebuilt from its own recipe with `--enable-openmp=yes`. Same pinning, same
parameter file, idle host (i9-13900K, 24 physical cores). This isolates the fork's changes and
nothing else.

### The large sparse problems — where upstream gains almost nothing

`main loop time` over a fixed 4-iteration budget, medians of 2:

| `dE4` (m=7401, routes **sparse**) | 1 | 2 | 4 | 8 | 16 | 24 | 1→24 |
|---|---:|---:|---:|---:|---:|---:|---:|
| upstream `6eaad8d` | 49.19 | 49.12 | 47.76 | 47.13 | 46.91 | 46.93 | **1.05×** |
| **this fork** | 47.08 | 47.29 | 24.38 | 14.11 | 8.31 | **6.44** | **7.31×** |
| fork advantage | 1.04× | 1.04× | 1.96× | 3.34× | 5.64× | **7.29×** | |

**Upstream gets 4.8% out of 24 cores on `dE4`.** Not a regression — nothing at all, because `dE4`
routes sparse and upstream threads *no part* of the sparse path: neither the Schur-complement
Cholesky nor its assembly. Both are threaded here, which is the whole of the difference.

| `dE3` (m=6067, routes **dense** at default) | 1 | 2 | 4 | 8 | 16 | 24 | 1→24 |
|---|---:|---:|---:|---:|---:|---:|---:|
| upstream `6eaad8d` | 434.62 | 408.85 | 211.26 | 121.36 | 70.81 | 60.84 | 7.14× |
| **this fork** | 176.82 | 168.20 | 88.26 | 48.92 | 27.69 | **22.93** | 7.71× |
| fork advantage | **2.46×** | 2.43× | 2.39× | 2.48× | 2.56× | 2.65× | |

`dE3` is the honest counter-example: upstream **does** scale here, 7.14×, because it routes dense
and the dense path *is* threaded upstream. The fork's lead at 24 threads is 2.65× — and **13.5×**
with `SDPA_BMAT_MODE=fill`, which takes the sparse route upstream cannot thread (4.49 s).

The 1-thread column is worth its own note: **2.46× before any threading enters**, from the
`Rgemm` zero-skip and the FP-contraction pin.

### Two problems from a user, and what they show

Both from `input/sdpaquestions/`, run with that folder's own `paramdd.sdpa`. These are the
sharpest demonstrations of two of this fork's changes, because they are real inputs that upstream
handles badly rather than synthetic fixtures.

#### `max_custom_scaled.dat-s` (m=16) — upstream is 24.5× slower on 24 cores than on 1

`main loop time`, full solve, best of 3:

| threads | 1 | 2 | 4 | 8 | 16 | 24 | 1→24 |
|---|---:|---:|---:|---:|---:|---:|---:|
| upstream `6eaad8d` | 0.0144 | 0.0896 | 0.1194 | 0.1540 | 0.2386 | 0.3545 | **0.04× — 24.5× slower** |
| **this fork** | 0.0119 | 0.0120 | 0.0119 | 0.0120 | 0.0120 | 0.0120 | **0.99× — flat** |

Upstream degrades **monotonically at every step**: each doubling of the thread count makes it
worse. The problem is tiny (7 blocks of order 4, 3, 10, 5, 6, 2, 2), so the fork/join costs far
more than the work it distributes, and upstream enters the parallel region regardless.

This fork is **flat to four decimal places** — 0.0119 s at 1 thread, 0.0120 s at 24 — because the
work/width gate declines to thread what cannot repay itself. At 24 threads the fork is **29.6×
faster**, entirely by *not* threading.

Upstream's answer also **drifts** here: **16 distinct objective values across 18 runs**, spanning
`-2.2168869933897821e-02` to `-2.2168869946397239e-02` — a relative spread of 5.6e-10, i.e. the
last ~10 digits of a double-double result are not reproducible. This fork returns **one** value in
all 18 runs, at every thread count.

*(The user's own recorded runs show the same effect end to end: 0.600 s total multi-core against
0.019 s single-core, same 43 iterations, same answer.)*

#### `_max.dat-s` — upstream cannot read it

Not a performance result. This file carries **601-digit mantissas**, and QD's decimal converter
overflows internally at ≥309 digits while *returning success*, so the values enter the SDP
matrices as NaN:

| | iterations | phase | objValPrimal |
|---|---:|---|---|
| upstream, as given | **0** | `dFEAS` | `-0.0` |
| upstream, same problem truncated to 30 digits | 43 | `pdOPT` | `-2.2168869950759854e-02` |
| **this fork**, as given | **43** | `pdOPT` | `-2.2168869945991022e-02` |

Upstream reports `cholesky miss condition :: not positive definite` at iteration 0 — an
infeasibility diagnosis for a problem that is perfectly feasible. Truncating the mantissas to 30
significant digits, still far beyond double-double's ~32-digit resolution so the *values* are
unchanged, makes upstream solve it correctly: **the cause is token length alone.**

This is the defect fixed by the reader work already in this fork, and the measurement above is an
independent confirmation of it against a fresh upstream build at the fork's own base commit.

---

### Small problems — where upstream's threading is actively harmful

Full solves, best of 2:

| problem | upstream 1 thr | upstream 24 thr | upstream 1→24 | fork 1→24 |
|---|---:|---:|---:|---:|
| `control1` | 0.018 s | 0.203 s | **0.09× — 11× slower** | **1.00×** |
| `truss5` | 0.800 s | 2.674 s | **0.30× — 3.3× slower** | 2.29× |
| `theta1` | 0.229 s | 0.238 s | 0.96× — slower | 2.00× |
| `arch0` | 13.044 s | 7.162 s | 1.82× | 6.33× |
| `gpp100` | 0.830 s | 0.305 s | 2.72× | 4.57× |

Adding threads makes upstream **11× slower** on `control1` and **3.3× slower** on `truss5`: the
work is far too small to repay a fork/join, and upstream enters the parallel region anyway.

The fork's `control1` row is the fix, visible directly — **1.00×, flat to four decimal places
across all 24 threads.** The work/width gating declines to thread what cannot pay for itself, so
adding cores costs exactly nothing rather than an order of magnitude.

**So "upstream scales negatively" is true, but only for small problems.** On `gpp100`, `arch0` and
`dE3` it scales perfectly respectably. Its distinctive failure on the large sparse problems is not
negative scaling at all — it is *flatness*, and that is the larger gap.

Raw rows: `dd_upstream_vs_fork.tsv`, `dd_small_scaling.tsv` in the review artifacts.

---

| problem | m | v712 | upstream1 | upstream24 | upstream8P | optimized1 | optimized8P | optimized24 | optimized24arena1 |
|---|---|---|---|---|---|---|---|---|---|
| control1 ‡ | 21 | 0.010 | 0.010 | 0.170 ±24% | 0.070 ±29% | 0.010 | 0.010 | 0.010 | 0.010 |
| gpp100 | 101 | 0.720 | 0.830 | 0.320 ±3% | 0.290 ±3% | 0.750 | 0.260 | 0.240 | 0.240 |
| control2 | 66 | 0.180 | 0.190 | 0.310 ±6% | 0.170 ±6% | 0.190 | 0.070 | 0.070 | 0.070 |
| gpp124-1 | 125 | 1.380 ±1% | 1.590 | 0.490 | 0.480 | 1.430 ±1% | 0.470 | 0.430 ±2% | 0.430 |
| theta1 ‡ | 104 | 0.140 | 0.230 | 0.240 ±4% | 0.170 ±6% | 0.180 | 0.120 | 0.120 | 0.120 |
| mcp100 ‡ | 100 | 0.750 | 1.280 ±1% | 0.600 ±12% | 0.530 ±15% | 0.810 ±2% | 0.410 | 0.400 | 0.400 |
| mcp124-1 ‡ | 124 | 1.060 | 2.130 | 1.090 ±34% | 0.950 ±46% | 1.150 ±1% | 0.590 ±3% | 0.590 ±2% | 0.590 |
| truss5 ‡ | 208 | 0.850 | 0.800 | 2.680 ±13% | 1.510 ±28% | 0.790 | 0.470 | 0.490 | 0.490 |
| arch0 ‡ | 174 | 10.220 | 13.050 | 7.370 ±4% | 7.280 ±5% | 9.950 | 2.690 ±1% | 2.470 ±1% | 2.450 |
| arch2 | 174 | 11.300 | 14.450 | 8.780 | 8.690 | 11.460 | 2.880 ±6% | 2.550 | 2.550 |
| arch4 ‡ | 174 | 11.290 | 14.650 | 8.860 ±3% | 8.790 ±6% | 11.620 | 2.900 | 2.570 | 2.570 |
| arch8 | 174 | 10.490 | 13.440 | 8.140 ±2% | 8.070 | 10.640 | 2.660 ±1% | 2.350 | 2.350 |
| control5 ‡ | 351 | 11.060 | 11.540 | 5.130 ±1% | 4.920 | 11.210 | 2.410 ±1% | 1.900 | 1.910 |
| control6 | 496 | 26.330 | 27.800 | 11.190 | 10.770 | 27.070 | 5.330 ±1% | 3.960 ±2% | 3.980 ±1% |
| theta2 ‡ | 498 | 4.310 | 5.050 | 2.240 ±26% | 2.250 ±5% | 4.600 | 2.100 ±1% | 2.130 | 2.130 |
| qap8 ‡ | 529 | 1.880 ±1% | 1.850 | 1.090 ±6% | 1.020 ±4% | 1.850 | 0.770 ±4% | 0.790 | 0.800 |
| control7 | 666 | 57.310 | 60.310 | 21.720 | 22.200 | 58.780 | 10.900 ±1% | 7.770 | 7.830 ±3% |
| qap9 | 748 | 4.600 | 4.730 | 2.190 ±5% | 2.270 ±11% | 4.720 | 1.660 | 1.640 ±1% | 1.650 ±1% |
| theta3 ‡ | 1106 | 44.380 | 51.150 | 15.210 ±36% | 13.920 ±19% | 48.900 | 14.510 | 13.650 | 13.660 |
| control8 ‡ | 861 | 111.310 | 118.680 | 38.910 ±2% | 41.550 ±2% | 116.010 | 20.700 | 14.660 ±2% | 14.460 ±2% |
| **total** | | **309.6** | **343.8** | **136.7** | **135.9** | **322.1** | **71.9** | **58.8** | **58.7** |

‡ = `v712` and `optimized24` take **different numbers of iterations** on this problem (control1 77 vs 71, theta1 43 vs 50, mcp100 53 vs 54, mcp124-1 49 vs 50, truss5 63 vs 62, arch0 75 vs 72, arch4 71 vs 72, control5 60 vs 59, theta2 54 vs 56, qap8 20 vs 19, theta3 63 vs 68, control8 60 vs 61). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**optimized24 vs v712, per iteration: 5.13x** (5.319 s -> 1.036 s per iteration, summed over the 20 of 20 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **3.10x**.

End-to-end wall time totals 309.6 s -> 58.8 s (5.27x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 12 of 20 problems (‡ above), so it mixes speed with path length.
Restricted to the 8 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 112.3 s -> 19.0 s = **5.91x** (gpp100, control2, gpp124-1, arch2, arch8, control6, control7, qap9).

### Integrity

- arch0/upstream24: iteration count VARIES across repeats (72/73/75)
- arch0/upstream8P: iteration count VARIES across repeats (72/73/75)
- arch4/upstream24: iteration count VARIES across repeats (71/72/73)
- arch4/upstream8P: iteration count VARIES across repeats (71/72/75)
- arch8/upstream24: iteration count VARIES across repeats (65/66)
- control1/upstream24: iteration count VARIES across repeats (69/74/88)
- control1/upstream8P: iteration count VARIES across repeats (70/81/90)
- control2/upstream24: iteration count VARIES across repeats (60/61/63)
- control2/upstream8P: iteration count VARIES across repeats (60/61/63)
- control5/upstream24: iteration count VARIES across repeats (59/60)
- control8/upstream24: iteration count VARIES across repeats (60/61)
- control8/upstream8P: iteration count VARIES across repeats (60/61)
- gpp100/upstream24: iteration count VARIES across repeats (23/24)
- gpp100/upstream24: objective VARIES across repeats
- gpp100/upstream8P: iteration count VARIES across repeats (23/24)
- gpp100/upstream8P: objective VARIES across repeats
- gpp124-1/upstream24: objective VARIES across repeats
- gpp124-1/upstream8P: objective VARIES across repeats
- mcp100/upstream24: iteration count VARIES across repeats (51/53/58)
- mcp100/upstream8P: iteration count VARIES across repeats (51/56/58)
- mcp124-1/upstream24: iteration count VARIES across repeats (55/70/79)
- mcp124-1/upstream8P: iteration count VARIES across repeats (63/66/93)
- qap8/upstream24: iteration count VARIES across repeats (20/21)
- qap8/upstream24: objective VARIES across repeats
- qap8/upstream8P: iteration count VARIES across repeats (20/21)
- qap8/upstream8P: objective VARIES across repeats
- qap9/upstream24: iteration count VARIES across repeats (19/20)
- qap9/upstream24: objective VARIES across repeats
- qap9/upstream8P: iteration count VARIES across repeats (19/21)
- qap9/upstream8P: objective VARIES across repeats
- theta1/upstream24: iteration count VARIES across repeats (44/45/46)
- theta1/upstream8P: iteration count VARIES across repeats (44/46/47)
- theta2/upstream24: iteration count VARIES across repeats (43/48/56)
- theta2/upstream8P: iteration count VARIES across repeats (51/52/53)
- theta3/upstream24: iteration count VARIES across repeats (51/65/75)
- theta3/upstream8P: iteration count VARIES across repeats (57/68)
- truss5/upstream24: iteration count VARIES across repeats (58/59/66)
- truss5/upstream8P: iteration count VARIES across repeats (60/68/79)
- control1: ITERATION COUNT differs BETWEEN configs (71, 77) -- compare per-iteration cost, not wall time
- gpp100: objective differs BETWEEN configs: -4.4943074465425575e+01, -4.4943074465425582e+01
- gpp124-1: objective differs BETWEEN configs: -7.3429256740235607e+00, -7.3429256740235651e+00
- theta1: ITERATION COUNT differs BETWEEN configs (43, 50) -- compare per-iteration cost, not wall time
- mcp100: ITERATION COUNT differs BETWEEN configs (53, 54) -- compare per-iteration cost, not wall time
- mcp124-1: ITERATION COUNT differs BETWEEN configs (49, 50) -- compare per-iteration cost, not wall time
- truss5: ITERATION COUNT differs BETWEEN configs (62, 63) -- compare per-iteration cost, not wall time
- arch0: ITERATION COUNT differs BETWEEN configs (72, 75) -- compare per-iteration cost, not wall time
- arch4: ITERATION COUNT differs BETWEEN configs (71, 72) -- compare per-iteration cost, not wall time
- control5: ITERATION COUNT differs BETWEEN configs (59, 60) -- compare per-iteration cost, not wall time
- theta2: ITERATION COUNT differs BETWEEN configs (54, 56) -- compare per-iteration cost, not wall time
- qap8: objective differs BETWEEN configs: -7.5694312988433944e+02, -7.5695155348092976e+02
- qap8: ITERATION COUNT differs BETWEEN configs (19, 20) -- compare per-iteration cost, not wall time
- qap9: objective differs BETWEEN configs: -1.4099194834660793e+03, -1.4099198052029226e+03
- theta3: ITERATION COUNT differs BETWEEN configs (63, 68) -- compare per-iteration cost, not wall time
- control8: ITERATION COUNT differs BETWEEN configs (60, 61) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | v712 | upstream1 | upstream24 | upstream8P | optimized1 | optimized8P | optimized24 | optimized24arena1 |
|---|---|---|---|---|---|---|---|---|
| control1 | 2.5 | 5.0 | 4.6 | 4.7 | 5.0 | 4.6 | 4.5 | 4.7 |
| gpp100 | 5.0 | 7.5 | 7.0 | 7.0 | 7.5 | 7.0 | 7.1 | 7.0 |
| control2 | 2.7 | 5.1 | 5.0 | 5.0 | 5.0 | 5.0 | 5.0 | 4.8 |
| gpp124-1 | 6.2 | 8.5 | 8.3 | 8.3 | 8.7 | 8.2 | 8.2 | 8.3 |
| theta1 | 3.2 | 5.5 | 5.5 | 5.3 | 5.5 | 5.6 | 5.7 | 5.3 |
| mcp100 | 4.7 | 6.8 | 6.8 | 6.7 | 6.8 | 6.8 | 6.8 | 6.8 |
| mcp124-1 | 5.7 | 8.0 | 8.0 | 7.7 | 8.0 | 7.7 | 8.0 | 8.0 |
| truss5 | 4.2 | 6.8 | 6.7 | 6.5 | 7.0 | 6.5 | 6.8 | 6.5 |
| arch0 | 8.2 | 10.7 | 10.2 | 10.3 | 10.7 | 10.2 | 10.3 | 10.2 |
| arch2 | 8.2 | 10.5 | 10.5 | 10.2 | 10.8 | 10.3 | 10.5 | 10.5 |
| arch4 | 8.2 | 10.7 | 10.6 | 10.2 | 10.7 | 10.2 | 10.5 | 10.5 |
| arch8 | 8.2 | 10.5 | 10.5 | 10.3 | 10.7 | 10.3 | 10.5 | 10.6 |
| control5 | 6.0 | 8.5 | 8.0 | 8.0 | 8.4 | 8.3 | 8.3 | 9.8 |
| control6 | 9.0 | 11.2 | 11.0 | 11.0 | 11.3 | 11.2 | 11.3 | 12.5 |
| theta2 | 8.7 | 11.2 | 10.7 | 10.7 | 11.2 | 10.8 | 11.8 | 11.3 |
| qap8 | 8.0 | 10.2 | 10.2 | 10.1 | 10.5 | 10.3 | 10.5 | 11.3 |
| control7 | 13.2 | 15.5 | 15.3 | 15.3 | 15.5 | 15.5 | 15.5 | 17.8 |
| qap9 | 13.2 | 15.7 | 15.3 | 15.3 | 15.7 | 15.5 | 15.5 | 17.5 |
| theta3 | 26.5 | 29.0 | 28.5 | 28.5 | 28.8 | 28.5 | 30.8 | 29.1 |
| control8 | 19.5 | 21.5 | 21.5 | 21.5 | 21.7 | 21.8 | 21.8 | 25.2 |


## Mac — Apple M1 Max (8P+2E), macOS

Reference columns (`v712`, `upstream1`, `upstream8`): Homebrew gcc-15, carried over unchanged.
Fork columns (`optimized1`, `optimized8`): re-measured with Homebrew gcc 16.1.0 on the
FP-contraction-pinned build. The compiler differs between the two halves of this table, so
fork-vs-reference ratios here carry a compiler change as well as the fork's own changes.

## mac-m1max — dd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | v712 | upstream1 | upstream8 | optimized1 | optimized8 |
|---|---|---|---|---|---|---|
| control1 | 21 | 0.040 ±25% | 0.050 ±40% | 0.820 ±49% | 0.030 | 0.030 ±67% |
| gpp100 | 101 | 0.920 ±2% | 1.000 ±5% | 0.800 ±10% | 0.580 ±3% | 0.220 ±18% |
| control2 ‡ | 66 | 0.220 ±5% | 0.220 | 0.810 ±11% | 0.160 | 0.090 |
| gpp124-1 | 125 | 1.680 ±1% | 1.860 ±3% | 1.200 | 1.070 ±5% | 0.360 ±19% |
| theta1 ‡ | 104 | 0.200 ±5% | 0.270 | 0.680 ±7% | 0.160 ±6% | 0.140 ±14% |
| mcp100 ‡ | 100 | 1.010 ±1% | 1.930 ±1% | 1.570 ±8% | 0.660 ±5% | 0.330 ±6% |
| mcp124-1 ‡ | 124 | 2.330 ±3% | 3.390 ±3% | 2.770 ±27% | 0.870 ±1% | 0.400 ±2% |
| truss5 ‡ | 208 | 1.130 | 0.950 ±1% | 8.350 ±2% | 0.650 ±5% | 0.570 ±9% |
| arch0 | 174 | 12.180 | 15.650 | 5.280 ±4% | 9.060 ±2% | 2.080 ±3% |
| arch2 | 174 | 14.140 | 17.720 | 5.630 ±4% | 10.870 ±1% | 2.310 ±4% |
| arch4 ‡ | 174 | 15.230 | 17.880 | 5.820 ±14% | 10.910 ±1% | 2.280 ±1% |
| arch8 ‡ | 174 | 12.990 ±1% | 16.480 | 5.310 ±5% | 9.910 | 2.060 ±1% |
| control5 ‡ | 351 | 12.070 ±1% | 12.670 ±1% | 4.060 ±1% | 7.840 ±1% | 1.600 ±3% |
| control6 | 496 | 29.500 ±1% | 30.870 ±2% | 7.580 ±3% | 19.260 ±1% | 3.760 ±11% |
| theta2 ‡ | 498 | 4.700 ±3% | 6.270 ±4% | 3.250 ±11% | 3.300 ±1% | 1.050 ±7% |
| qap8 | 529 | 2.140 ±1% | 2.260 ±3% | 1.320 ±17% | 1.460 ±1% | 0.410 ±12% |
| control7 | 666 | 62.820 ±1% | 65.790 ±2% | 13.740 ±1% | 40.910 ±1% | 6.800 ±4% |
| qap9 ‡ | 748 | 5.700 ±2% | 5.570 | 2.600 ±10% | 3.520 | 0.810 ±11% |
| theta3 ‡ | 1106 | 58.810 ±1% | 52.400 | 16.270 ±20% | 33.040 ±1% | 6.720 ±29% |
| control8 ‡ | 861 | 122.180 ±1% | 128.870 | 24.510 ±6% | 81.100 ±1% | 12.900 ±5% |
| **total** | | **360.0** | **382.1** | **112.4** | **235.4** | **44.9** |

‡ = `v712` and `optimized8` take **different numbers of iterations** on this problem (control2 61 vs 64, theta1 44 vs 50, mcp100 55 vs 54, mcp124-1 78 vs 50, truss5 59 vs 62, arch4 76 vs 72, arch8 65 vs 66, control5 60 vs 59, theta2 51 vs 56, qap9 20 vs 19, theta3 74 vs 68, control8 60 vs 61). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**optimized8 vs v712, per iteration: 7.78x** (6.024 s -> 0.774 s per iteration, summed over the 20 of 20 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **4.69x**.

End-to-end wall time totals 360.0 s -> 44.9 s (8.01x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 12 of 20 problems (‡ above), so it mixes speed with path length.
Restricted to the 8 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 123.4 s -> 16.0 s = **7.73x** (control1, gpp100, gpp124-1, arch0, arch2, control6, qap8, control7).

### Integrity

- arch0/upstream8: iteration count VARIES across repeats (72/75)
- arch4/upstream8: iteration count VARIES across repeats (72/81)
- arch8/upstream8: iteration count VARIES across repeats (65/66)
- control1/upstream8: iteration count VARIES across repeats (103/69/70)
- control2/upstream8: iteration count VARIES across repeats (61/63)
- control5/upstream8: iteration count VARIES across repeats (59/60)
- control8/upstream8: iteration count VARIES across repeats (60/61)
- gpp100/upstream8: objective VARIES across repeats
- gpp124-1/upstream8: objective VARIES across repeats
- mcp100/upstream8: iteration count VARIES across repeats (53/57/59)
- mcp124-1/upstream8: iteration count VARIES across repeats (57/76/77)
- qap8/upstream8: iteration count VARIES across repeats (19/21)
- qap8/upstream8: objective VARIES across repeats
- qap9/upstream8: iteration count VARIES across repeats (19/21)
- qap9/upstream8: objective VARIES across repeats
- theta1/upstream8: iteration count VARIES across repeats (46/47)
- theta2/upstream8: iteration count VARIES across repeats (51/54/56)
- theta3/upstream8: iteration count VARIES across repeats (51/58)
- truss5/upstream8: iteration count VARIES across repeats (61/62)
- control1: ITERATION COUNT differs BETWEEN configs (69, 71) -- compare per-iteration cost, not wall time
- gpp100: objective differs BETWEEN configs: -4.4943074465425575e+01, -4.4943074465425589e+01
- control2: ITERATION COUNT differs BETWEEN configs (59, 61, 64) -- compare per-iteration cost, not wall time
- gpp124-1: objective differs BETWEEN configs: -7.3429256740235607e+00, -7.3429256740235651e+00
- theta1: ITERATION COUNT differs BETWEEN configs (44, 45, 50) -- compare per-iteration cost, not wall time
- mcp100: ITERATION COUNT differs BETWEEN configs (54, 55, 66) -- compare per-iteration cost, not wall time
- mcp124-1: ITERATION COUNT differs BETWEEN configs (50, 62, 78) -- compare per-iteration cost, not wall time
- truss5: ITERATION COUNT differs BETWEEN configs (59, 62, 63) -- compare per-iteration cost, not wall time
- arch4: ITERATION COUNT differs BETWEEN configs (72, 76) -- compare per-iteration cost, not wall time
- arch8: ITERATION COUNT differs BETWEEN configs (65, 66) -- compare per-iteration cost, not wall time
- control5: ITERATION COUNT differs BETWEEN configs (59, 60) -- compare per-iteration cost, not wall time
- theta2: ITERATION COUNT differs BETWEEN configs (51, 56, 59) -- compare per-iteration cost, not wall time
- qap8: objective differs BETWEEN configs: -7.5694293990998108e+02, -7.5694312382795681e+02, -7.5694312988433944e+02
- qap9: objective differs BETWEEN configs: -1.4099198052029226e+03, -1.4099198799866915e+03, -1.4099345258838855e+03
- qap9: ITERATION COUNT differs BETWEEN configs (19, 20) -- compare per-iteration cost, not wall time
- theta3: ITERATION COUNT differs BETWEEN configs (62, 68, 74) -- compare per-iteration cost, not wall time
- control8: ITERATION COUNT differs BETWEEN configs (60, 61) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | v712 | upstream1 | upstream8 | optimized1 | optimized8 |
|---|---|---|---|---|---|
| control1 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| gpp100 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| control2 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| gpp124-1 | 8.5 | 8.5 | 8.6 | 8.4 | 8.4 |
| theta1 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| mcp100 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| mcp124-1 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| truss5 | 8.5 | 8.5 | 8.5 | 8.4 | 8.4 |
| arch0 | 9.2 | 9.5 | 9.6 | 9.3 | 9.5 |
| arch2 | 9.2 | 9.5 | 9.7 | 9.4 | 9.5 |
| arch4 | 9.3 | 9.5 | 9.7 | 9.4 | 9.5 |
| arch8 | 9.3 | 9.5 | 9.7 | 9.4 | 9.6 |
| control5 | 8.5 | 8.5 | 8.7 | 8.4 | 8.6 |
| control6 | 10.3 | 10.5 | 11.9 | 10.4 | 11.9 |
| theta2 | 9.6 | 9.9 | 11.0 | 9.8 | 10.8 |
| qap8 | 9.3 | 9.5 | 11.6 | 9.4 | 11.1 |
| control7 | 15.4 | 15.7 | 17.4 | 15.6 | 17.5 |
| qap9 | 14.2 | 14.5 | 16.8 | 14.3 | 16.7 |
| theta3 | 27.7 | 27.9 | 30.1 | 27.8 | 30.0 |
| control8 | 21.8 | 22.0 | 24.6 | 22.0 | 24.8 |


### The fork columns pin floating-point contraction, and now match x86-64 exactly

The two fork columns above (`optimized1`, `optimized8`) were re-measured; `v712`,
`upstream1` and `upstream8` are carried over unchanged from the previously published Mac run,
because re-building upstream with this fork's flags would change what "upstream" means in the
comparison.

This fork now compiles with `-ffp-contract=off` by default. On aarch64 the compiler was
otherwise free to fuse `a*b+c` into a single FMA, which is *more* accurate per operation but
breaks the two-sum/two-product identities double-double arithmetic is built on, so the solver
took a different sequence of iterates on this machine than on x86-64. With contraction pinned,
the M1 trajectory is **identical to x86-64 on 20 of 20 problems** —
20/20 on iteration count and 20/20 on the objective,
compared as decimal strings, not to a tolerance. Pinning moved the M1 iteration count on
9 problems (`control1`, `control2`, `control5`, `mcp100`, `mcp124-1`, `theta1`, `theta2`, `theta3`, `truss5`) and the objective
on 3 (`gpp124-1`, `qap8`, `qap9`); those moves are the
M1 converging onto the x86-64 path, not away from it. Per-problem record:
[`bench/b4_mac_rebaseline.tsv`](bench/b4_mac_rebaseline.tsv), whose `pinned_M1_iters` and
`pinned_M1_obj` columns are re-checked against every benchmark row when this document is
generated.

Consequence for the table: the fork columns and the reference columns no longer walk the same
path on the same problems, which is why the headline above is stated per iteration. The pin's
own runtime cost is **not** separable from the compiler change in this table — both fork
columns carry both — so no figure is quoted for it here.

Termination phase is now recorded per run, which an iterations-and-objective table cannot
show. Across both fork configurations and all three repeats it is constant per problem:
**noINFO** on 12 (`arch0`, `arch2`, `arch4`, `arch8`, `control1`, `control2`, `control5`, `control6`, `control7`, `control8`, `qap8`, `qap9`); **pFEAS** on 6 (`mcp100`, `mcp124-1`, `theta1`, `theta2`, `theta3`, `truss5`); **pdINF** on 2 (`gpp100`, `gpp124-1`). `gpp100` (23 iterations) and `gpp124-1` (24 iterations) terminate in `pdINF` — primal-dual infeasible — and stop early; their wall times are correspondingly small and are not fast *solves*.

Parameter file SHA-256 (16-hex prefix), identical on every re-measured row: `0abc064f632dd6ff`.

## thanos — AMD EPYC 7232P (8 physical cores), Ubuntu

The `fork*` binary here is the one produced by the README build instructions, from a fresh
clone of the published repository — this table validates the installation guide's output,
not a hand-configured tree.

## thanos-epyc7232p — dd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | v712 | upstream1 | upstream8 | fork1 | fork8 |
|---|---|---|---|---|---|---|
| control1 ‡ | 21 | 0.110 ±9% | 0.120 | 0.180 ±39% | 0.100 ±10% | 0.110 |
| gpp100 | 101 | 2.620 ±13% | 3.190 | 1.220 ±25% | 2.860 | 1.240 ±2% |
| control2 | 66 | 0.650 | 0.740 | 0.540 ±4% | 0.690 | 0.300 |
| gpp124-1 | 125 | 5.060 | 6.130 | 1.760 ±31% | 5.470 | 2.260 ±1% |
| theta1 ‡ | 104 | 0.510 | 0.840 | 0.590 ±8% | 0.640 | 0.450 |
| mcp100 ‡ | 100 | 2.810 | 4.880 | 1.880 ±27% | 3.010 | 1.700 |
| mcp124-1 ‡ | 124 | 4.410 | 8.170 | 3.850 ±33% | 4.260 | 2.420 |
| truss5 ‡ | 208 | 3.110 | 2.880 | 3.520 ±15% | 2.820 | 1.760 ±7% |
| arch0 | 174 | 36.560 | 50.410 | 30.770 ±1% | 38.070 | 10.690 ±3% |
| arch2 | 174 | 42.500 | 56.040 | 36.180 | 43.950 | 11.310 ±4% |
| arch4 ‡ | 174 | 42.450 | 56.780 | 36.720 ±2% | 44.510 | 11.370 ±5% |
| arch8 | 174 | 39.510 | 52.070 | 33.560 | 40.830 | 10.530 ±5% |
| control5 ‡ | 351 | 39.770 | 43.670 | 17.990 ±2% | 41.840 | 8.470 ±7% |
| control6 | 496 | 95.180 | 105.820 | 39.530 ±1% | 101.650 | 19.080 ±3% |
| theta2 ‡ | 498 | 14.990 | 19.440 | 11.000 ±49% | 17.480 | 8.390 ±2% |
| qap8 ‡ | 529 | 6.930 | 7.140 | 4.030 ±2% | 7.020 | 3.110 ±7% |
| control7 | 666 | 208.020 | 230.790 | 83.440 | 222.140 | 40.160 ±1% |
| qap9 ‡ | 748 | 19.760 | 18.540 | 9.290 ±14% | 18.240 | 7.560 ±2% |
| theta3 ‡ | 1106 | 131.040 | 201.840 | 66.660 ±11% | 190.910 | 73.730 |
| control8 | 861 | 411.300 | 455.220 | 155.070 ±2% | 439.560 | 77.330 |
| **total** | | **1107.3** | **1324.7** | **537.8** | **1226.0** | **292.0** |

‡ = `v712` and `fork8` take **different numbers of iterations** on this problem (control1 77 vs 71, theta1 43 vs 50, mcp100 55 vs 54, mcp124-1 56 vs 50, truss5 63 vs 62, arch4 71 vs 72, control5 60 vs 59, theta2 51 vs 56, qap8 20 vs 19, qap9 22 vs 19, theta3 50 vs 68). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**fork8 vs v712, per iteration: 3.82x** (19.475 s -> 5.096 s per iteration, summed over the 20 of 20 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **2.56x**.

End-to-end wall time totals 1107.3 s -> 292.0 s (3.79x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 11 of 20 problems (‡ above), so it mixes speed with path length.
Restricted to the 9 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 841.4 s -> 172.9 s = **4.87x** (gpp100, control2, gpp124-1, arch0, arch2, arch8, control6, control7, control8).

### Integrity

- arch0/upstream8: iteration count VARIES across repeats (72/73)
- arch4/upstream8: iteration count VARIES across repeats (72/73)
- control1/upstream8: iteration count VARIES across repeats (109/69/82)
- control2/upstream8: iteration count VARIES across repeats (60/62)
- control5/upstream8: iteration count VARIES across repeats (59/60)
- control8/upstream8: iteration count VARIES across repeats (60/61)
- gpp100/upstream8: objective VARIES across repeats
- gpp124-1/upstream8: objective VARIES across repeats
- mcp100/upstream8: iteration count VARIES across repeats (52/54/58)
- mcp124-1/upstream8: iteration count VARIES across repeats (53/66/78)
- qap8/upstream8: iteration count VARIES across repeats (19/20)
- qap8/upstream8: objective VARIES across repeats
- qap9/upstream8: iteration count VARIES across repeats (19/21)
- qap9/upstream8: objective VARIES across repeats
- theta1/upstream8: iteration count VARIES across repeats (45/47/48)
- theta2/upstream8: iteration count VARIES across repeats (50/62/82)
- theta3/upstream8: iteration count VARIES across repeats (52/54/58)
- truss5/upstream8: iteration count VARIES across repeats (59/62/67)
- control1: ITERATION COUNT differs BETWEEN configs (71, 77) -- compare per-iteration cost, not wall time
- gpp100: objective differs BETWEEN configs: -4.4943074465425575e+01, -4.4943074465425582e+01
- gpp124-1: objective differs BETWEEN configs: -7.3429256740235580e+00, -7.3429256740235607e+00
- theta1: ITERATION COUNT differs BETWEEN configs (43, 50) -- compare per-iteration cost, not wall time
- mcp100: ITERATION COUNT differs BETWEEN configs (54, 55) -- compare per-iteration cost, not wall time
- mcp124-1: ITERATION COUNT differs BETWEEN configs (50, 56) -- compare per-iteration cost, not wall time
- truss5: ITERATION COUNT differs BETWEEN configs (62, 63) -- compare per-iteration cost, not wall time
- arch4: ITERATION COUNT differs BETWEEN configs (71, 72) -- compare per-iteration cost, not wall time
- control5: ITERATION COUNT differs BETWEEN configs (59, 60) -- compare per-iteration cost, not wall time
- theta2: ITERATION COUNT differs BETWEEN configs (51, 56) -- compare per-iteration cost, not wall time
- qap8: objective differs BETWEEN configs: -7.5694312988433944e+02, -7.5694487097376361e+02
- qap8: ITERATION COUNT differs BETWEEN configs (19, 20) -- compare per-iteration cost, not wall time
- qap9: objective differs BETWEEN configs: -1.4099198052029226e+03, -1.4099353243156127e+03
- qap9: ITERATION COUNT differs BETWEEN configs (19, 22) -- compare per-iteration cost, not wall time
- theta3: ITERATION COUNT differs BETWEEN configs (50, 68) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | v712 | upstream1 | upstream8 | fork1 | fork8 |
|---|---|---|---|---|---|
| control1 | 11.7 | 11.6 | 11.7 | 11.8 | 11.8 |
| gpp100 | 11.8 | 11.8 | 11.7 | 11.7 | 11.7 |
| control2 | 11.8 | 11.7 | 11.7 | 11.7 | 11.7 |
| gpp124-1 | 11.7 | 11.7 | 11.8 | 11.7 | 11.8 |
| theta1 | 11.8 | 11.7 | 11.8 | 11.7 | 11.6 |
| mcp100 | 11.7 | 11.7 | 11.7 | 11.7 | 11.7 |
| mcp124-1 | 11.7 | 11.8 | 11.8 | 11.8 | 11.7 |
| truss5 | 11.8 | 11.8 | 11.7 | 11.8 | 11.7 |
| arch0 | 11.7 | 11.7 | 11.8 | 11.7 | 11.8 |
| arch2 | 11.8 | 11.8 | 11.7 | 11.8 | 11.7 |
| arch4 | 11.8 | 11.8 | 11.8 | 11.8 | 11.8 |
| arch8 | 11.8 | 11.8 | 11.9 | 11.9 | 11.8 |
| control5 | 11.8 | 11.7 | 11.8 | 11.7 | 11.8 |
| control6 | 11.7 | 11.8 | 11.7 | 11.8 | 11.7 |
| theta2 | 11.8 | 11.8 | 11.9 | 11.8 | 12.0 |
| qap8 | 11.8 | 11.8 | 11.8 | 11.7 | 11.8 |
| control7 | 13.2 | 15.4 | 15.1 | 15.6 | 16.4 |
| qap9 | 13.3 | 15.6 | 15.7 | 15.6 | 16.5 |
| theta3 | 26.2 | 28.6 | 28.7 | 28.7 | 29.9 |
| control8 | 19.5 | 21.9 | 21.2 | 21.7 | 22.5 |


Two figures against 7.1.2 appear above because mpack (7.1.2) and mplapack (master) round
differently and take different iteration counts on 12 of 20 problems. Against
upstream the right reading is: **serial** upstream and optimized trajectories match on
20/20 problems on pi and 20/20 on thanos, but only
11/20 on the Mac, where pinning FP contraction moved the fork onto the
x86-64 trajectory and therefore *off* the trajectory this machine's own unpinned upstream
takes (`control1`, `control2`, `control5`, `mcp100`, `mcp124-1`, `theta1`, `theta2`, `theta3`, `truss5`). And **threaded upstream is nondeterministic** (the integrity
sections above flag every instance), so the headline vs-upstream ratios are observed
end-to-end speedups rather than strictly same-trajectory comparisons. The optimized build's
own trajectory is identical at every thread count on every machine here.

> **`v7.1.3-omp.2` was not, and this is worth stating where the reproducibility claims are made.**
> That release carried a data race in the threaded sparse assembly: on a problem routed sparse
> whose blocks share Schur-complement entries — SDPLIB `truss6` at default settings — it returned
> a different answer on nearly every run. The tables above are unaffected (their problems either
> route dense, fall below the threading gate, or were re-measured with repetition), but the
> *claim* was false for that release and is only true again from `v7.1.3-omp.3`. Note also what
> the claim's shape cost: stated as identity **across thread counts**, it directed every check at
> comparisons between thread counts, and a race that moves both sides passes those. The
> reproducibility evidence below now includes repetition at ONE thread count.

## Reproducibility evidence

The optimized build is **thread-count independent**: an FNV-1a checksum over every scalar of
the final `xMat`/`yVec`/`zMat` (normalised representation; -0.0 folded to +0.0) is identical
at 1, 2, 4, 8 and 24 threads on all 12 problems tested -- see `bench/statehash_pi.tsv` for
the full record with source commit and parameter checksum. The record was measured on the
**published tip** (`6ebb025`) built from a clean clone, and all 60 cells are hash-identical
to the record taken on the pre-publication series. The instrument was validated
against a build known to differ, so agreement is informative rather than vacuous.

Unpatched upstream is **not** reproducible when threaded: across the same benchmark the
harness flagged 84 anomalies on pi and 20 on the Mac -- every one in an upstream threaded
configuration, including three problems whose *objective* differs run to run
(`gpp100`, `gpp124-1`, `qap8`).

Raw data: [`bench/pi_dd_v2.tsv`](bench/pi_dd_v2.tsv),
[`bench/mac_dd_v3.tsv`](bench/mac_dd_v3.tsv) (fork columns re-measured; reference columns
carried over from [`bench/mac_dd_v2.tsv`](bench/mac_dd_v2.tsv), which is retained),
[`bench/b4_mac_rebaseline.tsv`](bench/b4_mac_rebaseline.tsv),
[`bench/dd_v2_thanos.tsv`](bench/dd_v2_thanos.tsv),
[`bench/statehash_pi.tsv`](bench/statehash_pi.tsv).
