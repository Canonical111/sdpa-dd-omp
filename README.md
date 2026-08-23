# sdpa-dd-omp

OpenMP-threaded, reproducible fork of sdpa-dd. ("-omp" and not "parallel": SDPARA is the
SDPA group's own MPI solver, and "mp" already means multiprecision in this family.)

Fork of [nakatamaho/sdpa-dd](https://github.com/nakatamaho/sdpa-dd) at commit
`6eaad8d9` (the original upstream README is preserved as
[README-UPSTREAM.md](README-UPSTREAM.md)), carrying four patches that fix performance and
reproducibility regressions. **Build it:** [INSTALL.md](INSTALL.md). The issues were reported to the upstream maintainer and have not
been adopted there; this fork exists so the fixes are usable.

## What changed (4 commits on top of upstream)

| commit | change |
|---|---|
| 1 | restore netlib dgemm's zero-skip in `Rgemm_NN/NT_omp` — dropped in mplapack 2.0.1, up to 2× serial regression |
| 2 | gate OpenMP on work / width / nesting; run `Rdot` serially — up to 13× on small problems, removes a nondeterminism |
| 3 | thread the Schur-complement (`bMat`) constraint loop — the dominant gain |
| 4 | report the reduced per-formula timers as worker-seconds, not elapsed time |

The four patches above are the original port. The fork has since gained substantially more --
a ten-item correctness batch (exit statuses, container ownership, input bounds, an
uninitialised-read fix), record-bounded input parsing with line-numbered diagnostics, threaded
Cholesky panel kernels and X/Z inverse-Cholesky triangulars behind measured work gates,
hardware FMA for the bundled QD on aarch64, and the FP-contraction pin described below. Each
change carries its evidence in its commit message; [BENCHMARKS.md](BENCHMARKS.md) holds the
regenerated tables.

Every modified or added source file carries an in-file, dated change notice naming its own
licence: GPLv2 §2a for the SDPA sources, and the 2-clause BSD terms for `mplapack/`, which
carries no GNU licence at all. The
complete, always-current list is `git diff --stat <upstream-base>..HEAD` — an enumeration
here went stale twice and is deliberately not repeated. `git log` has the rationale per
change.

## Measured results

20 SDPLIB problems (m = 21…1106), external wall clock, median of 3 pinned repeats:

| | EPYC 7232P (8 cores) | i9-13900K (24 cores) | M1 Max (8P+2E) |
|---|---|---|---|
| sdpa-dd 7.1.2 (2009) | 1107.3 s | 309.6 s | 360.0 s |
| upstream master, threaded | 537.8 s | 136.7 s | 112.4 s |
| **this fork** (2026-08-02 campaign) | **292.0 s** | **58.8 s** | **86.0 s** |

> The table above is the ORIGINAL release campaign, kept as history. The current
> tables live in [BENCHMARKS.md](BENCHMARKS.md), regenerated 2026-08-07 after the
> threading and contraction work: the same optimized8 column now totals 44.9 s and
> the headline is quoted PER ITERATION (7.78x vs v712), because the contraction pin
> changed trajectories on 9 of 20 problems and a wall ratio across different
> iteration counts mixes speed with path length.
| vs upstream master (end-to-end) | **1.84×** | **2.33×** | **1.31×** |

Serial upstream and optimized trajectories matched on 20/20 problems, on all three machines, in that campaign (pre-contraction-pin; see BENCHMARKS.md for the current comparison).
**Threaded upstream is nondeterministic** (e.g. `control1` across three repeats: 69/74/88
iterations at 24 threads on the i9, 69/82/109 at 8 threads on the EPYC; this fork: 71/71/71
on both), so the vs-upstream figures are observed end-to-end speedups rather than strictly
same-trajectory comparisons. The EPYC fork binary is the one built by this README's own
instructions from a fresh clone. Full tables, methodology and raw per-repeat data:
[BENCHMARKS.md](BENCHMARKS.md) and [`bench/`](bench/).

Results are **thread-count independent**: a checksum over every scalar of the final
`xMat`/`yVec`/`zMat` is identical at 1/2/4/8/24 threads. Unpatched upstream is *not*
reproducible when threaded (iteration counts and even objectives vary run to run).

Results are also **platform independent between arm64 and x86-64**, which upstream is not.
GCC and clang default to `-ffp-contract=fast` for C++, and `dd_real` arithmetic is
header-inline, so without a pin *the compiler* decides whether the compensation terms of the
double-double algorithms get folded into an FMA — and folding there is not an identity, because
it deletes the very rounding error `two_sum` exists to capture. On baseline x86-64 there is no
FMA to fold into, so those binaries behave as if contraction were off; on arm64 `fmadd` is
baseline ISA, so they do not. The two platforms then run different arithmetic.

This fork therefore builds with `-ffp-contract=off` **by default**, applied to every SDPA and
MPLAPACK translation unit *and* to the bundled QD. With it, an M1 build matches the published
x86-64 reference on **20 of 20** problems — iteration counts and full-precision objectives —
where the unpinned build matches on none at the state-hash level. Hardware FMA is unaffected:
`-ffp-contract=off` does not disable an explicit `__builtin_fma`, so QD's `two_prod` keeps its
single `fmsub` and the aarch64 speedup is retained in full. Cost is ~1% per iteration at 1
thread and ~1.7% at 8. x86-64 output is unchanged — all 76 translation units compile to
byte-identical `.text` with and without the flag.

Pass `--enable-fp-contract=fast` to recover the old, platform-dependent behaviour. Note that
iteration counts published for arm64 *before* this default are counts of the unpinned solver;
see [`bench/b4_mac_rebaseline.tsv`](bench/b4_mac_rebaseline.tsv) for the exact deltas.

Evidence in this repository: [BENCHMARKS.md](BENCHMARKS.md) (generated tables, methodology),
[`bench/pi_dd_v2.tsv`](bench/pi_dd_v2.tsv) and [`bench/mac_dd_v2.tsv`](bench/mac_dd_v2.tsv)
(per-repeat raw data), [`bench/statehash_pi.tsv`](bench/statehash_pi.tsv) (thread-count
independence record). The investigation history, harness and patch generators live in the
companion repository (to be published).

## Building

```bash
autoreconf -fi                       # upstream ships configure.ac only, no configure
./configure --enable-openmp=yes
make -j$(nproc)
```

Verified on a fresh Ubuntu 24.04 clone (gcc 13) and by CI on every push. Plain `make` is
expected to work: upstream's `.POSIX:` Make.inc selects `c99`, which rejects the flags and used
to stop the build with a `struct timezone` error, so the tree now applies
`external/spooles/patches/patch-Make.inc` itself and passes the configured compiler on the
SPOOLES make line. CI asserts that the patched flags reached the real compile lines, so a
regression here fails the build rather than being papered over by hand. **If** `make` still
stops inside SPOOLES, that fix has regressed — please report it; do not repair `Make.inc`
manually, because the next `make` re-extracts SPOOLES and discards the edit.

Prefer it automated? This repository packages an agent skill —
[`.claude/skills/install-sdpa-omp/`](.claude/skills/install-sdpa-omp/) — that performs the
whole installation (compiler detection, the SPOOLES rescue, OpenMP verification, smoke test)
on Linux and macOS, verified on x86-64 and Apple Silicon. Claude Code discovers it
automatically in a clone; `bash .claude/skills/install-sdpa-omp/scripts/install.sh` also
works standalone.

### macOS (Apple Silicon) — verified on an M1 Max

**Do not run plain `./configure`**: `/usr/bin/gcc` is Apple clang, which has no OpenMP —
configure reports `option to support OpenMP... unsupported` and *silently builds a serial
binary*. Homebrew GCC is required, and SPOOLES needs it spelled out too:

```bash
brew install gcc autoconf automake libtool
GCC=$(ls $(brew --prefix gcc)/bin/gcc-[0-9]* | head -1)   # resolve the current version
GXX=$(ls $(brew --prefix gcc)/bin/g++-[0-9]* | head -1)   # (brew install may have just upgraded it)
autoreconf -fi
./configure CC="$GCC" CXX="$GXX" --enable-openmp=yes
make -j8          # the SPOOLES compiler/flags fix is applied in-tree; no manual Make.inc edit
otool -L sdpa_dd | grep gomp    # must print libgomp -- the proof the build is threaded
```

Run with `OMP_NUM_THREADS=<physical cores>`, pinned (`taskset`/`OMP_PLACES=cores`).
Serial builds (`--enable-openmp=no`) are supported and CI-checked.

## License

GPL v2, unchanged from upstream (`COPYING`). Copyright remains with the original
SDPA authors; the patches in this fork are contributed under the same license.

## Exit status

Scripts that loop over problems can rely on the exit code:

| outcome | exit |
|---|---|
| solver ran to any stopping condition (`pdOPT`, `pdFEAS`, `pFEAS`, `dFEAS`, `pdINF`, `pINF_dFEAS`, `pFEAS_dINF`, `pUNBD`, `dUNBD`, `noINFO`) | **0** |
| iteration limit reached | **0** |
| infeasibility / unboundedness detected | **0** |
| malformed input, unreadable file, invalid parameter | **1** with a diagnostic (line-numbered for data files) |
| numerical failure with nothing valid to print -- no iteration completed, or the updated `X`/`Z` left the cone **and the rollback could not be refactored** | **2**, `solveStatus = FAILURE` in the result file, no solution section |
| recoverable late failure after `k` good iterations -- a Schur factorisation failure, or an `X`/`Z` update that left the cone and was **rolled back**; the last valid iterate is printed and labelled | **3**, `solveStatus = PARTIAL`, `failureIteration = k` in the result file |

Infeasibility and the iteration limit are valid mathematical results, not errors. Upstream
exited 0 on *every* path -- including fatal errors -- so a crashed run was indistinguishable
from a solved one in any harness that checks exit codes.
