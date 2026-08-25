# sdpa-dd-omp

An OpenMP-threaded, reproducible fork of **sdpa-dd**, a semidefinite-programming solver that works
in double-double arithmetic — roughly 32 significant digits, carried as an unevaluated sum of two
IEEE doubles.

Upstream builds with OpenMP, but threading it as shipped is not merely unprofitable. On small
problems more cores make it monotonically **slower** — 24.5× slower on 24 cores than on one, on a
user's m=16 input. On large sparse problems it threads *neither* of the two regions that dominate
the runtime, so 24 cores buy it **4.9%**. And its threaded runs are not reproducible: one SDPLIB
problem ran 69, 74 and 88 interior-point iterations across three identical invocations.

This fork threads the regions that actually dominate, makes threaded results bit-identical at any
thread count, fixes a decimal reader that silently turned high-precision input into NaN, and
re-derives the rule that decides how the Schur complement is factored.

Fork of [nakatamaho/sdpa-dd](https://github.com/nakatamaho/sdpa-dd) at `6eaad8d9` — see that
repository for upstream's own documentation. Reported upstream; not adopted there, which is why
this fork exists.

| | |
|---|---|
| **Build and run it** | [INSTALL.md](INSTALL.md) — build, verify, choose a thread count, troubleshoot |
| **Why it works this way** | [doc/technical.pdf](doc/technical.pdf) — mechanisms, the derivation behind the factorisation rule, every environment variable, the exit-status contract, and what is *not* established (source: [doc/technical.tex](doc/technical.tex)) |
| **Full benchmark tables** | [BENCHMARKS.md](BENCHMARKS.md), raw per-repeat data in [`bench/`](bench/) |

## What was improved

**Threading the regions that dominate.** The `bMat` constraint loop, the **sparse Schur-complement
Cholesky**, the **sparse Schur-complement assembly**, the Cholesky panel kernels and the X/Z
inverse-Cholesky triangulars are all threaded here. The last two of those five are the ones that
matter on large problems and are entirely serial upstream — which is the whole explanation for the
`dE4` row in the tables below.

Both the assembled Schur complement and the finished Cholesky factor are **bit-identical at any
thread count by construction**, not by luck: for a fixed pivot each worker owns one destination
row, every read comes from the pivot row, and no sum is ever reassociated.

The backward triangular solve is the one region deliberately left serial. It *is* a reassociated
sum, and no arrangement avoids that. `sdpa-gmp-omp` threads it behind an opt-in flag by fixing the
chunk count independently of the team size, which yields one reproducible answer at every thread
count — but a *different* answer from the serial reduction, which is exactly why it is opt-in
there. dd could do the same with a `qd_real` accumulator; that would likewise be a new answer, not
the old one computed faster. This fork keeps the pass serial rather than spend bit-identity with
the serial computation on a fraction of one phase. [doc/technical.pdf](doc/technical.pdf) §5.5 has
the argument.

**Threading no longer makes small problems slower.** Upstream enters its parallel regions
regardless of problem size — it has a threshold mechanism, but it is disabled behind `if (0)` and
cannot simply be switched on, because the `_ref` bodies it dispatches to are absent and would not
link. On a user's `max_custom_scaled.dat-s` (m=16) upstream degrades at *every* thread count,
ending **24.5× slower on 24 cores than on 1**; this fork is flat to four decimal places, because
the work/width gate declines to thread what cannot repay itself.

**Hardware FMA for the bundled QD on Apple Silicon — about 1.6×.** QD's own `configure` runs its
FMA probe only under `case $host in powerpc*-*-*)`, so the default `--enable-fma=auto` silently
resolves to *none* on aarch64: `QD_FMA`/`QD_FMS` stay undefined and `two_prod` — the most-executed
primitive in the whole solver — compiles the Dekker two-way split instead of `p=a*b; err=fma(a,b,-p)`.
On Apple Silicon that is **72 instructions against 8** (`dd_real::operator*`: 71 against 13), the
difference including two conditional branches for the split's overflow guard. Enabling it moves the
solver's main loop by **1.59× in aggregate** at one thread, 1.21–1.63× per problem, against a
repeat spread of 0.1–2.8%.

It is also **bit-identical and strictly more correct**, both measured rather than assumed: 16
SDPLIB problems across 96 runs give 0/16 hash mismatches in all four comparisons with identical
iteration counts, and near `DBL_MAX` the Dekker split overflows and returns `inf`/`nan` for a
finite product where FMA does not. The gate is deliberately **aarch64-only** — on x86-64 without
`-mfma`, GCC lowers `__builtin_fma` to a libm *call*, which is slower than the split, and this
project builds generic x86-64 on purpose so binaries stay bit-identical across AMD and Intel nodes.
x86-64 is provably untouched: 74/74 objects identical after stripping, `.text` byte-identical, and
all 60 published state-hash cells regenerated and matched.

**High-precision input is read correctly.** QD's decimal converter overflows internally at 309 or
more mantissa digits and *returns success while producing NaN*, so a 600-digit spelling of an
ordinary number entered the SDP matrices as NaN and surfaced as a bogus `cholesky miss condition`
at iteration 0 — the solver reporting an infeasible problem when the problem was fine. On a user's
600-digit input, upstream terminates after **0 iterations**; this fork reaches **pdOPT in 43
iterations**. Truncating that same input to 45 digits makes *upstream* solve it correctly, which
isolates the cause to token length rather than to arithmetic.

**The same answers on arm64 and x86-64.** GCC and clang default to `-ffp-contract=fast`, and
`dd_real` arithmetic is header-inline — so *the compiler* decides whether the compensation terms of
the double-double algorithms get folded into an FMA, and folding there deletes the very rounding
error `two_sum` exists to capture. Baseline x86-64 has no FMA to fold into and arm64 does, so the
two platforms ran different arithmetic. This fork pins `-ffp-contract=off` by default across every
SDPA, `mplapack` and bundled-QD translation unit; an M1 build then matches the x86-64 reference on
**20 of 20** problems. Hardware FMA is untouched (`__builtin_fma` still compiles to one `fmsub`),
x86-64 object code is byte-identical with and without the flag, and the cost is ~1% per iteration.

**Reproducibility, measured rather than asserted.** A checksum over every scalar of the final
`xMat`/`yVec`/`zMat` is identical at 1/2/4/8/24 threads, and CI compares the *assembled matrix and
the factor themselves* — not the printed solution, which carries 17 digits where a `dd_real`
carries 32, so a factor can change without any printed field moving. Repeated runs at **one**
thread count are checked too; that is the check a race passes the between-thread-count comparison
by moving both sides, and its absence is exactly what let `v7.1.3-omp.2` ship.

**A re-derived factorisation choice (2026-08-24).** Four gates decide dense vs sparse for the Schur
complement. Gate 3 is a cheap pre-screen on the *aggregate* sparsity pattern; gate 4 is the real
test, on the *ordered fill*. Gate 3 screened at 0.25·m² while gate 4 rejected at 0.40·m² — so the
pre-screen was the stricter test and overruled the real one. Demoting gate 3 to gate 4's own
constant is provable rather than optimistic: symbolic factorisation only *adds* entries, so
aggregate > F implies fill > F, which makes 0.40 the **largest sound cutoff** and 0.25 conservatism
with nothing behind it. **One tunable replaces two.** No SDPLIB problem changes route; on a
221-instance bootstrap census, 167 stop taking a route that cost 2.6–6.7× in time and 1.8–2.5× in
memory. Not "sparse always wins": `truss5` has an ordered fill of 1.0 and dense is genuinely faster
there — both rules send it dense, which is the point of testing fill rather than assuming.

**Two smaller inherited defects, fixed.** `mplapack` 2.0.1 dropped netlib `dgemm`'s zero-skip guard,
so multiplying by zero still paid full multiprecision cost — restored in `Rgemm_NN_omp.cpp` and
`Rgemm_NT_omp.cpp`, worth up to 2× serial. And the reduced per-formula timers reported
worker-seconds where elapsed time was expected, so a region looked *more* expensive the more it was
parallelised.

**Reporting that does not lie.** Upstream exits 0 on every path, including fatal errors, so a
crashed run is indistinguishable from a solved one in any harness that checks exit codes. This fork
distinguishes them, and reports a recoverable numerical failure as `solveStatus = PARTIAL` naming
the failing iteration. The full contract is in [doc/technical.pdf](doc/technical.pdf) §10.

## The benchmarks

Three views, and they measure different things:

- **20 SDPLIB problems on three machines** — the broad-coverage set, and the one to start from,
  though every problem in it takes a dense `bMat` and so misses the path below entirely;
- **the large sparse problems** — where the difference actually lives, because they are the only
  ones that reach the threaded sparse Cholesky and assembly;
- **small problems** — where upstream's threading is actively harmful and this fork's gating is the
  whole story.

[BENCHMARKS.md](BENCHMARKS.md) is the full dossier and says which table to use for what.

### 20 SDPLIB problems (m = 21…1106), three machines

External wall clock, median of 3 pinned repeats. The headline is quoted **per iteration** — because
the FP-contraction pin changed iteration counts on 12 of the 20 problems, and a wall-clock ratio
across differing iteration counts mixes speed with path length.

| | EPYC 7232P (8 cores) | i9-13900K (24 cores) | M1 Max (8P+2E) |
|---|---:|---:|---:|
| sdpa-dd 7.1.2 (2009), single-threaded | 1107.3 s | 309.6 s | 360.0 s |
| upstream master, threaded, at its best | 537.8 s | 136.7 s | 112.4 s |
| **this fork**, at its best | **292.0 s** @8 thr | **58.8 s** @24 thr | **44.9 s** @8 thr |
| **vs sdpa-dd 7.1.2, per iteration** | **3.82×** | **5.13×** | **7.78×** |

The per-iteration row is the one to quote; the wall-clock rows above it are given because they are
what a user actually waits for, but a ratio between them mixes speed with path length. Per-machine
detail, per-problem rows and the iteration counts are in [BENCHMARKS.md](BENCHMARKS.md).

**Threaded upstream is nondeterministic** — `control1` across three repeats: 69/74/88 iterations at
24 threads on the i9, 69/82/109 at 8 threads on the EPYC, where this fork runs 71/71/71 on both —
so the vs-upstream figures are observed end-to-end speedups rather than strictly same-trajectory
comparisons.

### The large sparse problems, against upstream

Measured 2026-08-25 on current `main` against upstream `6eaad8d` — the commit this fork branched
from and still its `master` — each rebuilt from its own recipe, i9-13900K (24 physical cores),
**medians of three**, the two builds interleaved cell by cell. Both problems take the fork's
**default** route:

| | upstream, 1 thread | upstream, 24 threads | **this fork, 24 threads** | fork vs upstream |
|---|---:|---:|---:|---:|
| `dE4` (m=7401, routes sparse) | 49.20 s | 46.92 s — **1.05×** | **6.42 s** | **≈7.3×** |
| `dE3` (m=6067, routes sparse since 2026-08-24) | 434.62 s | 60.80 s — 7.15× | **4.57 s** | **≈13×** |

`dE4` is the striking column: **24 cores buy upstream 4.9%**, because it routes sparse and upstream
threads neither the sparse Schur-complement Cholesky nor its assembly. Both are threaded here.

Against this fork's own historical baseline, measured the same way under a fixed four-iteration
protocol, `dE4` runs in about **6.4 s** against **46.329 s** for the unmodified fork —
about a **7.2×** improvement.

These are **fixed-budget comparisons, not time-to-convergence results** — neither problem converges
at double-double precision under either tolerance tested (the shipped `epsilonStar=1e-30` and a
relaxed `1e-20`). Raw rows and full provenance:
[bench/dd-port3-2026-08-24/](bench/dd-port3-2026-08-24/).

### Small problems

Upstream scales **negatively**: 11× slower on `control1` and 3.3× slower on `truss5` at 24 threads
than at 1. This fork's work gating keeps `control1` flat to four decimal places. Full curves in
[BENCHMARKS.md](BENCHMARKS.md).

## License

GPL v2, unchanged from upstream (`COPYING`). Copyright remains with the original SDPA authors; the
patches in this fork are contributed under the same license. Every modified or added source file
carries an in-file, dated change notice naming its own licence: GPLv2 §2a for the SDPA sources, and
the 2-clause BSD terms for `mplapack/`, which carries no GNU licence at all.
