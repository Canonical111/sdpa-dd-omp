# Installing sdpa-dd-omp
> **Runtime knobs:** every environment variable this fork adds, with whether it can change
> results, is documented in [RUNTIME.md](RUNTIME.md).


Double-double (~32 significant digits) SDP solver. Linux and macOS, x86-64 and aarch64.

## The short version

```bash
git clone https://github.com/Canonical111/sdpa-dd-omp.git
cd sdpa-dd-omp
bash .claude/skills/install-sdpa-omp/scripts/install.sh
```

That configures, builds, and **verifies** — it refuses to report success unless the bundled
example solves to the expected objective and the OpenMP runtime is actually linked. Add
`--serial` for a build with no OpenMP, or `--prefix=DIR` to install the binary.

The script is the tested path: CI executes it on every push, in both modes, so a regression in it
cannot ship quietly.

## Doing it by hand

```bash
autoreconf -fi                      # upstream ships no configure
./configure --enable-openmp=yes
make -j"$(nproc)"                   # nproc is fine HERE (build parallelism, not solver threads)
```

**Always run `autoreconf -fi`**, even on a tree you have built before: a `configure` that survived
a `git pull` can be stale relative to `configure.ac`, and reusing it silently configures
yesterday's build.

macOS needs a real GCC (Apple clang does not ship OpenMP):

```bash
brew install gcc
./configure --enable-openmp=yes CXX=g++-14 CC=gcc-14   # match your brew version
```

`nproc` is safe for `make -j` but is **wrong for the solver** — see thread counts below.

## Verifying the build

The bundled example is a real check, not a smoke test — but it is tiny, so it cannot demonstrate
threading:

```bash
./sdpa_dd -ds example1.dat-s -o out.result -p param.sdpa
grep -E 'phase.value|objValPrimal' out.result
```

Expect `phase.value = pdOPT` and `objValPrimal = -4.1899999999999999e+01`. A different objective
means something is wrong with the build; a different *iteration count* alone can be a compiler
difference and is covered by the FP-contraction pin described in the README.

### The property worth checking yourself

This fork guarantees **bit-identical results at any thread count**. That is a real, testable
claim — but `example1` cannot test it. It is m=10, so it routes dense, never reaches the threaded
Cholesky or the threaded assembly, and would print the same answer on a build with no threading
at all. Use a problem that actually threads:

```bash
python3 tests/gen_crossblock_fixture.py > xb.dat-s      # m=480, routes sparse, threads
for t in 1 4 8; do
  OMP_NUM_THREADS=$t ./sdpa_dd -ds xb.dat-s -o t$t.result -p param.sdpa
done
diff <(grep -E 'objValPrimal|objValDual|phase.value|Iteration =' t1.result) \
     <(grep -E 'objValPrimal|objValDual|phase.value|Iteration =' t8.result)
```

**And run it twice at the SAME thread count**, which is the check that matters most:

```bash
for r in 1 2 3 4 5; do
  ./sdpa_dd -ds xb.dat-s -o r$r.result -p param.sdpa >/dev/null
  grep -h '^objValPrimal' r$r.result
done | sort -u | wc -l          # must print 1
```

That second form is not belt-and-braces. `v7.1.3-omp.2` contained a data race in the threaded
assembly that produced a different answer on nearly every run; comparing thread counts would not
have caught it, because both sides of the comparison were wrong in different ways. It is fixed,
and both forms above are in CI.

No output from the `diff`, and `1` from the second, means the guarantee holds on your build. **Two settings deliberately change results and
must not be varied during this check:** `SDPA_BMAT_MODE` (a different factorisation route) and
`SDPA_SPCHOL_MODE=force` only insofar as it changes *which* pivots thread — the factor itself stays
bit-identical, but leave both unset for a clean comparison.

## Choosing a thread count

Match `OMP_NUM_THREADS` to **physical cores**, never SMT threads. On Linux, pin:

```bash
OMP_PROC_BIND=true OMP_PLACES=cores taskset -c 0-23 \
  OMP_NUM_THREADS=24 ./sdpa_dd -ds problem.dat-s -o out.result -p param.sdpa
```

**On a two-socket machine, measure before committing a long run.** Compare one socket's
physical-core count against the whole machine's, pinning each run to exactly the cores it should
use. Crossing the socket boundary helps on some machines and costs badly on others, and the answer
is problem-dependent as well as machine-dependent.

**On a hybrid CPU (P+E cores), benchmark rather than assume.** On this fork's own i9-13900K
measurements all 24 physical cores beat the 8 P-cores alone on the five-problem total (58.8 s
against 71.9 s) — but that is one machine and one problem set. See BENCHMARKS.md.

**How much threading buys depends entirely on which route your problem takes.** Both columns below
come from **one `v7.1.3-omp.3` binary in one session** on pi (i9-13900K, 24 physical cores),
medians of three interleaved repeats — dE3 (m=6067) and dE4 (m=7401):

| | 1 thread | 24 threads | 1→24 |
|---|---:|---:|---:|
| dE4, sparse route (its default) | 46.612 s | **6.399 s** | **7.28×** |
| dE3, dense route (its default) | 176.357 s | 22.645 s | 7.79× |
| dE3, sparse route (`SDPA_BMAT_MODE=fill`) | 31.370 s | **4.475 s** | 7.01× |

`fill` is also about **5× faster than dE3's own default route** at 24 threads, on 2.3× less memory
— see below.

> **One caveat, because the data does not support more precision than this.** The `dE3`-with-`fill`
> 24-thread cell is short (~4.5 s) and jittery: six runs across two campaigns span 4.34–4.90 s, a
> **12%** spread, against ≤2.5% for every other cell and ≤0.1% for the 1-thread cells. So the
> `fill`-versus-default advantage is *about 5×*, not 5.06× or 4.92×; both figures have been quoted
> and neither deserves three significant figures. The dE4 and dE3-dense rows are steady and do
> support theirs.

Raw rows and full provenance — build hash, compiler, input and parameter hashes, affinity, route
and overlap policy per cell: [bench/dd-port3-2026-08-24/dd_v3_scaling_1t_24t.tsv](bench/dd-port3-2026-08-24/dd_v3_scaling_1t_24t.tsv).

These are `main loop time` over a **fixed 4-iteration budget** — per-iteration solver cost, not
time to solution. The assembly race fix costs `dE4` +0.7% and `dE3`-via-`fill` +2.2% in time, and —
since the setup pass that decides whether two blocks may overlap was reduced from one `int` to
**two bits** per stored Schur entry — **nothing measurable in memory**: `dE4` is 382.0 MB against
382.2 MB before the fix, and `dE3`-via-`fill` 277.5 MB against 277.5 MB.

Neither problem converged at double-double precision under **either of the two tolerance settings
tested** — the shipped `epsilonStar=1e-30` and a relaxed `1e-20`. These runs therefore compare a
fixed four-iteration computation budget rather than time to convergence; a looser parameter choice
is a different experiment, and nothing here rules it out.

## Route selection, and the one knob worth knowing

`SDPA_BMAT_MODE` selects how the Schur complement is factored: `auto` (default, unchanged),
`fill`, `dense`, `sparse`. All are strictly parsed — a typo is refused rather than silently
treated as `auto`.

**`fill` is opt-in and can be a large win on large sparse problems.** On dE3 it takes the route
`auto` declines: **4.5 s against 22.6 s, and 277 MB against 671 MB** — about 5× faster on 2.4× less
memory. `SDPA_BMAT_LOG=1` prints which gate decided and why.

That is not one lucky problem. A route census identified **seven distinct structures** on which the
two policies disagree, covering all 167 affected instances, and `fill` has now been measured on
**every one of them**:

| m | `auto` (dense) | `fill` (sparse) | faster | less memory |
|---:|---:|---:|---:|---:|
| 2,439 | 1.92 s / 116 MB | 0.61 s / 50 MB | 3.17× | 2.30× |
| 4,489 | 9.22 s / 392 MB | 3.59 s / 213 MB | 2.57× | 1.84× |
| 5,278 | 14.50 s / 519 MB | 3.88 s / 234 MB | 3.74× | 2.22× |
| 6,067 | 22.97 s / 671 MB | 4.50 s / 277 MB | 5.11× | 2.42× |
| 8,359 | 66.20 s / 1,289 MB | 15.13 s / 587 MB | 4.38× | 2.20× |
| 10,614 | 150.47 s / 2,032 MB | 22.36 s / 831 MB | 6.73× | 2.45× |
| 11,227 | 175.81 s / 2,358 MB | 47.94 s / 1,183 MB | 3.67× | 1.99× |

**It remains opt-in anyway, and the reason is worth understanding before you set it.**
`SDPA_BMAT_MODE` is *result-changing*: dense and sparse are different factorisation routes and
follow different iterate trajectories. The figures above are per-iteration costs over a fixed
four-iteration budget, so they cannot tell you whether the two routes need the same *number* of
iterations to reach a tolerance — and these instances do not converge at double-double precision
under the tolerances tested, so no time-to-solution comparison exists to quote. What is established
is that on per-iteration cost and peak memory, `fill` dominates on every structure that
distinguishes the policies.

Raw rows: [bench/dd-port3-2026-08-24/dd_fill_seven_structures.tsv](bench/dd-port3-2026-08-24/dd_fill_seven_structures.tsv).

`SDPA_BMAT_MAX_GB` caps the dense allocation and is enforced on every route to dense, so a problem
that would silently need more memory than you have fails with a number instead.

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `configure` reports no OpenMP support | Apple clang, or missing `-fopenmp` | install GCC and pass `CXX=g++-14 CC=gcc-14` |
| built, but threading changes nothing | binary is serial, or threads unpinned | check `nm ./sdpa_dd \| grep GOMP`; pin with `taskset` + `OMP_PLACES=cores` |
| more threads is *slower* | crossed a socket boundary, or oversubscribed SMT | use one socket's physical cores; never count SMT threads |
| `SDPA_* must be ...` on startup | a typo in an env knob | the parser is strict by design; fix the value |
| stale `configure` after `git pull` | autotools inputs newer than `configure` | `autoreconf -fi` (the installer does this unconditionally) |
| iteration count differs across compilers | FP contraction | see the README's `-ffp-contract=off` note |

Every environment variable, the exit-status contract, and what each knob can and cannot change
are documented in [RUNTIME.md](RUNTIME.md), which is the maintained reference for all 25 of them.

## A successful build leaves the tree almost clean

`autoreconf -fi` refreshes the tracked `INSTALL` file, and on macOS the QD build replaces
`config.guess`/`config.sub`. Those are the only expected working-tree changes. If restoring them,
run `git diff` on each first and confirm the changes are only those refreshes — never
blanket-restore.
