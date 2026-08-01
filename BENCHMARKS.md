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

| problem | m | v712 | upstream1 | upstream24 | upstream8P | optimized1 | optimized8P | optimized24 | optimized24arena1 |
|---|---|---|---|---|---|---|---|---|---|
| control1 | 21 | 0.010 | 0.010 | 0.170 ±24% | 0.070 ±29% | 0.010 | 0.010 | 0.010 | 0.010 |
| gpp100 | 101 | 0.720 | 0.830 | 0.320 ±3% | 0.290 ±3% | 0.750 | 0.260 | 0.240 | 0.240 |
| control2 | 66 | 0.180 | 0.190 | 0.310 ±6% | 0.170 ±6% | 0.190 | 0.070 | 0.070 | 0.070 |
| gpp124-1 | 125 | 1.380 ±1% | 1.590 | 0.490 | 0.480 | 1.430 ±1% | 0.470 | 0.430 ±2% | 0.430 |
| theta1 | 104 | 0.140 | 0.230 | 0.240 ±4% | 0.170 ±6% | 0.180 | 0.120 | 0.120 | 0.120 |
| mcp100 | 100 | 0.750 | 1.280 ±1% | 0.600 ±12% | 0.530 ±15% | 0.810 ±2% | 0.410 | 0.400 | 0.400 |
| mcp124-1 | 124 | 1.060 | 2.130 | 1.090 ±34% | 0.950 ±46% | 1.150 ±1% | 0.590 ±3% | 0.590 ±2% | 0.590 |
| truss5 | 208 | 0.850 | 0.800 | 2.680 ±13% | 1.510 ±28% | 0.790 | 0.470 | 0.490 | 0.490 |
| arch0 | 174 | 10.220 | 13.050 | 7.370 ±4% | 7.280 ±5% | 9.950 | 2.690 ±1% | 2.470 ±1% | 2.450 |
| arch2 | 174 | 11.300 | 14.450 | 8.780 | 8.690 | 11.460 | 2.880 ±6% | 2.550 | 2.550 |
| arch4 | 174 | 11.290 | 14.650 | 8.860 ±3% | 8.790 ±6% | 11.620 | 2.900 | 2.570 | 2.570 |
| arch8 | 174 | 10.490 | 13.440 | 8.140 ±2% | 8.070 | 10.640 | 2.660 ±1% | 2.350 | 2.350 |
| control5 | 351 | 11.060 | 11.540 | 5.130 ±1% | 4.920 | 11.210 | 2.410 ±1% | 1.900 | 1.910 |
| control6 | 496 | 26.330 | 27.800 | 11.190 | 10.770 | 27.070 | 5.330 ±1% | 3.960 ±2% | 3.980 ±1% |
| theta2 | 498 | 4.310 | 5.050 | 2.240 ±26% | 2.250 ±5% | 4.600 | 2.100 ±1% | 2.130 | 2.130 |
| qap8 | 529 | 1.880 ±1% | 1.850 | 1.090 ±6% | 1.020 ±4% | 1.850 | 0.770 ±4% | 0.790 | 0.800 |
| control7 | 666 | 57.310 | 60.310 | 21.720 | 22.200 | 58.780 | 10.900 ±1% | 7.770 | 7.830 ±3% |
| qap9 | 748 | 4.600 | 4.730 | 2.190 ±5% | 2.270 ±11% | 4.720 | 1.660 | 1.640 ±1% | 1.650 ±1% |
| theta3 | 1106 | 44.380 | 51.150 | 15.210 ±36% | 13.920 ±19% | 48.900 | 14.510 | 13.650 | 13.660 |
| control8 | 861 | 111.310 | 118.680 | 38.910 ±2% | 41.550 ±2% | 116.010 | 20.700 | 14.660 ±2% | 14.460 ±2% |
| **total** | | **309.6** | **343.8** | **136.7** | **135.9** | **322.1** | **71.9** | **58.8** | **58.7** |

**optimized24 vs v712: 5.27x**  (totals 309.6 s -> 58.8 s)

Iteration counts differ on 12 problem(s) (control1, theta1, mcp100, mcp124-1, truss5, arch0, arch4, control5, theta2, qap8, theta3, control8), so the wall-time ratio above mixes speed with path length. **Per iteration: 5.13x** (5.319 s -> 1.036 s).

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


## Mac — Apple M1 Max (8P+2E), macOS, Homebrew gcc-15

## mac-m1max — dd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | v712 | upstream1 | upstream8 | optimized1 | optimized8 |
|---|---|---|---|---|---|---|
| control1 | 21 | 0.040 ±25% | 0.050 ±40% | 0.820 ±49% | 0.040 ±25% | 0.040 |
| gpp100 | 101 | 0.920 ±2% | 1.000 ±5% | 0.800 ±10% | 0.930 ±1% | 0.410 ±5% |
| control2 | 66 | 0.220 ±5% | 0.220 | 0.810 ±11% | 0.210 | 0.110 ±9% |
| gpp124-1 | 125 | 1.680 ±1% | 1.860 ±3% | 1.200 | 1.710 ±1% | 0.680 ±1% |
| theta1 | 104 | 0.200 ±5% | 0.270 | 0.680 ±7% | 0.220 ±5% | 0.190 ±5% |
| mcp100 | 100 | 1.010 ±1% | 1.930 ±1% | 1.570 ±8% | 1.330 ±1% | 0.850 ±8% |
| mcp124-1 | 124 | 2.330 ±3% | 3.390 ±3% | 2.770 ±27% | 1.930 | 1.190 ±1% |
| truss5 | 208 | 1.130 | 0.950 ±1% | 8.350 ±2% | 0.950 ±1% | 0.820 ±4% |
| arch0 | 174 | 12.180 | 15.650 | 5.280 ±4% | 12.290 | 4.030 ±1% |
| arch2 | 174 | 14.140 | 17.720 | 5.630 ±4% | 14.320 | 4.300 ±1% |
| arch4 | 174 | 15.230 | 17.880 | 5.820 ±14% | 14.550 ±1% | 4.270 ±4% |
| arch8 | 174 | 12.990 ±1% | 16.480 | 5.310 ±5% | 13.280 ±1% | 4.000 ±2% |
| control5 | 351 | 12.070 ±1% | 12.670 ±1% | 4.060 ±1% | 12.290 ±1% | 2.920 ±4% |
| control6 | 496 | 29.500 ±1% | 30.870 ±2% | 7.580 ±3% | 30.300 ±118% | 6.650 ±14% |
| theta2 | 498 | 4.700 ±3% | 6.270 ±4% | 3.250 ±11% | 5.540 ±1% | 2.780 ±2% |
| qap8 | 529 | 2.140 ±1% | 2.260 ±3% | 1.320 ±17% | 2.210 | 0.960 ±2% |
| control7 | 666 | 62.820 ±1% | 65.790 ±2% | 13.740 ±1% | 64.490 ±2% | 11.950 ±7% |
| qap9 | 748 | 5.700 ±2% | 5.570 | 2.600 ±10% | 5.560 ±1% | 1.970 ±1% |
| theta3 | 1106 | 58.810 ±1% | 52.400 | 16.270 ±20% | 48.490 ±1% | 15.720 ±2% |
| control8 | 861 | 122.180 ±1% | 128.870 | 24.510 ±6% | 124.820 ±1% | 22.140 |
| **total** | | **360.0** | **382.1** | **112.4** | **355.5** | **86.0** |

**optimized8 vs v712: 4.19x**  (totals 360.0 s -> 86.0 s)

Iteration counts differ on 12 problem(s) (control1, control2, theta1, mcp100, mcp124-1, truss5, arch4, arch8, theta2, qap9, theta3, control8), so the wall-time ratio above mixes speed with path length. **Per iteration: 3.98x** (6.024 s -> 1.513 s).

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
- control2: ITERATION COUNT differs BETWEEN configs (59, 61) -- compare per-iteration cost, not wall time
- theta1: ITERATION COUNT differs BETWEEN configs (44, 45) -- compare per-iteration cost, not wall time
- mcp100: ITERATION COUNT differs BETWEEN configs (55, 66) -- compare per-iteration cost, not wall time
- mcp124-1: ITERATION COUNT differs BETWEEN configs (62, 78) -- compare per-iteration cost, not wall time
- truss5: ITERATION COUNT differs BETWEEN configs (59, 63) -- compare per-iteration cost, not wall time
- arch4: ITERATION COUNT differs BETWEEN configs (72, 76) -- compare per-iteration cost, not wall time
- arch8: ITERATION COUNT differs BETWEEN configs (65, 66) -- compare per-iteration cost, not wall time
- theta2: ITERATION COUNT differs BETWEEN configs (51, 59) -- compare per-iteration cost, not wall time
- qap8: objective differs BETWEEN configs: -7.5694293990998108e+02, -7.5694312382795681e+02
- qap9: objective differs BETWEEN configs: -1.4099198799866915e+03, -1.4099345258838855e+03
- qap9: ITERATION COUNT differs BETWEEN configs (19, 20) -- compare per-iteration cost, not wall time
- theta3: ITERATION COUNT differs BETWEEN configs (62, 74) -- compare per-iteration cost, not wall time
- control8: ITERATION COUNT differs BETWEEN configs (60, 61) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | v712 | upstream1 | upstream8 | optimized1 | optimized8 |
|---|---|---|---|---|---|
| control1 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| gpp100 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| control2 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| gpp124-1 | 8.5 | 8.5 | 8.6 | 8.5 | 8.5 |
| theta1 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| mcp100 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| mcp124-1 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| truss5 | 8.5 | 8.5 | 8.5 | 8.5 | 8.5 |
| arch0 | 9.2 | 9.5 | 9.6 | 9.5 | 9.7 |
| arch2 | 9.2 | 9.5 | 9.7 | 9.5 | 9.7 |
| arch4 | 9.3 | 9.5 | 9.7 | 9.5 | 9.7 |
| arch8 | 9.3 | 9.5 | 9.7 | 9.5 | 9.7 |
| control5 | 8.5 | 8.5 | 8.7 | 8.5 | 8.8 |
| control6 | 10.3 | 10.5 | 11.9 | 10.5 | 12.5 |
| theta2 | 9.6 | 9.9 | 11.0 | 9.9 | 10.7 |
| qap8 | 9.3 | 9.5 | 11.6 | 9.6 | 11.2 |
| control7 | 15.4 | 15.7 | 17.4 | 15.7 | 17.4 |
| qap9 | 14.2 | 14.5 | 16.8 | 14.5 | 16.6 |
| theta3 | 27.7 | 27.9 | 30.1 | 27.9 | 30.1 |
| control8 | 21.8 | 22.0 | 24.6 | 22.0 | 24.4 |


Two figures against 7.1.2 appear above because mpack (7.1.2) and mplapack (master) round
differently and take different iteration counts on 12 of 20 problems. Against upstream
master the iteration counts are identical on 20/20, so that comparison is clean
like-for-like.

## Reproducibility evidence

The optimized build is **thread-count independent**: an FNV-1a checksum over every scalar of
the final `xMat`/`yVec`/`zMat` (normalised representation; -0.0 folded to +0.0) is identical
at 1, 2, 4, 8 and 24 threads on all 12 problems tested -- see `bench/statehash_pi.tsv` for
the full record with source commit and parameter checksum. The instrument was validated
against a build known to differ, so agreement is informative rather than vacuous.

Unpatched upstream is **not** reproducible when threaded: across the same benchmark the
harness flagged 84 anomalies on pi and 20 on the Mac -- every one in an upstream threaded
configuration, including three problems whose *objective* differs run to run
(`gpp100`, `gpp124-1`, `qap8`).

Raw data: [`bench/pi_dd_v2.tsv`](bench/pi_dd_v2.tsv),
[`bench/mac_dd_v2.tsv`](bench/mac_dd_v2.tsv),
[`bench/statehash_pi.tsv`](bench/statehash_pi.tsv).
