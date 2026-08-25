# Building sdpa-dd-omp

Double-double (~32 significant digits) SDP solver. Linux and macOS, x86-64 and aarch64.

Verified by execution on Ubuntu x86-64 and macOS Apple Silicon, and by CI on every push. Every trap
in the troubleshooting table is a real failure that was hit, not a hypothetical.

## The short version

```bash
git clone https://github.com/Canonical111/sdpa-dd-omp.git
cd sdpa-dd-omp
bash .claude/skills/install-sdpa-omp/scripts/install.sh
```

That configures, builds, and **verifies** — it refuses to report success unless the bundled example
solves to the expected objective and the OpenMP runtime is actually linked. Add `--serial` for a
build with no OpenMP (legitimate and CI-tested, not a failure mode), or `--prefix=DIR` to install
the binary.

The script is the tested path: CI executes it on every push, in both modes, so a regression in it
cannot ship quietly. That same directory is packaged as a Claude Code skill and is discovered
automatically in a clone, if you would rather an agent do it.

## Doing it by hand

**Prerequisites.** A C++ compiler with OpenMP, and the Autotools — `autoconf`, `automake`,
`libtool`. Upstream ships no `configure`, so `autoreconf` is the *first* command below, not an
optional refresh; neither GCC nor the Xcode command-line tools provide it.

```bash
autoreconf -fi                      # upstream ships no configure
./configure --enable-openmp=yes
make -j"$(nproc)"                   # nproc is fine HERE (build parallelism, not solver threads)
```

**Always run `autoreconf -fi`**, even on a tree you have built before: a `configure` that survived a
`git pull` can be stale relative to `configure.ac`, and reusing it silently configures yesterday's
build.

`nproc` is safe for `make -j` but is **wrong for the solver** — see [Choosing a thread
count](#choosing-a-thread-count).

### macOS (Apple Silicon)

Apple clang does not ship OpenMP, so plain `./configure` produces a **silent serial binary**. A real
GCC is required — and so are the Autotools, which macOS does not provide:

```bash
brew install gcc autoconf automake libtool
GCC=$(ls "$(brew --prefix gcc)"/bin/gcc-[0-9]* | head -1)   # resolve the version you have --
GXX=$(ls "$(brew --prefix gcc)"/bin/g++-[0-9]* | head -1)   # brew may have just upgraded it
autoreconf -fi
./configure --enable-openmp=yes CC="$GCC" CXX="$GXX"
make -j8
otool -L sdpa_dd | grep omp                                # proof the build is threaded
```

The two `ls` lines are the same version discovery the installer script does. Hard-coding
`gcc-14`/`g++-14` works until the next `brew upgrade` renames the binary, so resolve it rather than
spell it.

## Verifying the build

### It solves, and gets the right answer

```bash
./sdpa_dd -ds example1.dat-s -o out.result -p param.sdpa
grep -E 'phase.value|objValPrimal' out.result
```

Expect `phase.value = pdOPT` and `objValPrimal = -4.1899999999999999e+01`. A different objective
means something is wrong with the build. A different *iteration count* alone can be a compiler
difference, and is what the FP-contraction pin exists to remove — see
[doc/technical.pdf](doc/technical.pdf) §4.

### The property worth checking yourself

This fork guarantees **bit-identical results at any thread count**. That is a real, testable claim —
but `example1` cannot test it. It is m=10, so it routes dense, never reaches the threaded Cholesky
or the threaded assembly, and would print the same answer on a build with no threading at all. Use a
problem that actually threads:

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
assembly that produced a different answer on nearly every run, and comparing thread counts would not
have caught it, because both sides of the comparison were wrong in different ways. It is fixed, and
both forms above are in CI. The full account is in [doc/technical.pdf](doc/technical.pdf) §6.

No output from the `diff`, and `1` from the second, means the guarantee holds on your build. **Two
settings deliberately change results and must not be varied during this check:** `SDPA_BMAT_MODE`
(a different factorisation route) and `SDPA_SPCHOL_MODE=force` — only insofar as it changes *which*
pivots thread; the factor itself stays bit-identical. Leave both unset for a clean comparison.

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
against 71.9 s) — but that is one machine and one problem set. See [BENCHMARKS.md](BENCHMARKS.md).

### How much threading buys depends on which route your problem takes

Both columns below come from **one binary in one session** on an i9-13900K (24 physical cores),
medians of three interleaved repeats — `dE3` (m=6067) and `dE4` (m=7401):

| | 1 thread | 24 threads | 1→24 | peak RSS @24 |
|---|---:|---:|---:|---:|
| dE4, sparse route (its default) | 46.555 s | **6.390 s** | **7.29×** | 382 MB |
| dE3, sparse route (**its default since 2026-08-24**) | 31.361 s | **4.417 s** | 7.10× | **277 MB** |
| dE3, dense route (`SDPA_BMAT_MODE=legacy`) | 176.287 s | 22.526 s | 7.83× | 671 MB |

Measured on current `main`. Across all 18 runs there are exactly **two distinct objectives, one per
problem**: identical at 1 and 24 threads, and identical between dE3's dense and sparse routes. Raw
rows and full provenance — build hash, compiler, input and parameter hashes, affinity, route and
overlap policy per cell:
[bench/dd-port3-2026-08-24/dd_current_scaling.tsv](bench/dd-port3-2026-08-24/dd_current_scaling.tsv).

`dE3`'s default route changed on 2026-08-24. It now takes the sparse route: **about 5× faster on
2.4× less memory** than the route it used to take. `SDPA_BMAT_MODE=legacy` restores the old one.

> **One caveat, because the data does not support more precision than this.** The `dE3` sparse
> 24-thread cell is short (~4.5 s) and jittery: six runs across two campaigns span 4.34–4.90 s, a
> **12%** spread, against ≤2.5% for every other cell and ≤0.1% for the 1-thread cells. So that
> advantage is *about 5×*, not 5.06× or 4.92×; both figures have been quoted and neither deserves
> three significant figures. The dE4 and dE3-dense rows are steady and do support theirs.

> **The two memory figures are `main`'s, not the latest release's.** The two-bit overlap map is not
> in `v7.1.3-omp.3`, which holds one `int` per Schur entry and uses 405.8 MB and 294.2 MB on the
> same two cells — about 6% more. Everything else here, including every timing, applies to both.

These are `main loop time` over a **fixed 4-iteration budget** — per-iteration solver cost, not time
to solution. Neither problem converged at double-double precision under
**either of the two tolerance settings tested**: the shipped `epsilonStar=1e-30` and a relaxed
`1e-20`. A looser parameter choice is a different experiment, and nothing here rules it out.

## Running it — where the reference lives

Everything about *running* the solver rather than building it is in
[doc/technical.pdf](doc/technical.pdf):

- **all 25 `SDPA_*` environment variables**, what each does, and which of them can change a computed
  value — §9, "Environment variables";
- **the exit-status contract**, every code from 0 to 3, including `solveStatus = PARTIAL` and
  `failureIteration` — §10;
- why the route chooser exists, what its four gates decide, and the proof that gate 3 is not a
  heuristic — §7, "The factorisation choice";
- the reproducibility oracles and their negative controls — §8.

**In normal use none of it is needed:** set `OMP_NUM_THREADS` as above and leave every `SDPA_*`
variable alone. Three are worth knowing about:

| variable | when you want it |
|---|---|
| `SDPA_BMAT_MODE=legacy` | reproduce a dd result from before 2026-08-24 bit-for-bit |
| `SDPA_BMAT_LOG=1` | see which factorisation route was chosen, and which gate decided |
| `SDPA_BMAT_MAX_GB=N` | cap the dense allocation, so a problem that would need more memory than you have fails with a number instead of thrashing |

`SDPA_BMAT_MODE` is the only one of the 25 that can change a computed result, and it does so by
design: dense and sparse are different factorisation routes and follow different iterate
trajectories. All values are strictly parsed — a typo is refused rather than silently treated as
the default.

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `configure` reports no OpenMP support | Apple clang, or missing `-fopenmp` | install GCC and pass `CXX=g++-14 CC=gcc-14` |
| built, but threading changes nothing | binary is serial, or threads unpinned | check `nm ./sdpa_dd \| grep GOMP`; pin with `taskset` + `OMP_PLACES=cores` |
| more threads is *slower* | crossed a socket boundary, or oversubscribed SMT | use one socket's physical cores; never count SMT threads |
| `SDPA_* must be ...` on startup | a typo in an env knob | the parser is strict by design; fix the value |
| stale `configure` after `git pull` | autotools inputs newer than `configure` | `autoreconf -fi` (the installer does this unconditionally) |
| iteration count differs across compilers | FP contraction | see [doc/technical.pdf](doc/technical.pdf) §4; `--enable-fp-contract=fast` recovers the old behaviour |
| `make` stops inside SPOOLES | the in-tree `Make.inc` patch has regressed | please report it — do **not** repair `Make.inc` by hand, because the next `make` re-extracts SPOOLES and discards the edit |

## A successful build leaves the tree almost clean

`autoreconf -fi` refreshes the tracked `INSTALL` file, and on macOS the QD build replaces
`config.guess`/`config.sub`. Those are the only expected working-tree changes. If restoring them,
run `git diff` on each first and confirm the changes are only those refreshes — never
blanket-restore.
