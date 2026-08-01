---
name: install-sdpa-omp
description: Build, verify and install the sdpa-dd-omp / sdpa-gmp-omp / sdpa-qd-omp SDP solvers from source on Linux or macOS. Use when the user asks to install, build, compile or set up this solver (or its sibling forks), or when a build of it has failed.
---

# Installing the sdpa-*-omp solvers

These are OpenMP-threaded forks of the SDPA multiprecision SDP solvers. Every instruction
here was verified by execution on Ubuntu x86-64 and macOS Apple Silicon; the traps listed
below are real failures that were hit, not hypotheticals.

## Standard procedure

From the repository root:

```bash
bash .claude/skills/install-sdpa-omp/scripts/install.sh            # auto-detects the solver
bash .claude/skills/install-sdpa-omp/scripts/install.sh --prefix ~/.local   # also installs
bash .claude/skills/install-sdpa-omp/scripts/install.sh --serial   # no-OpenMP build
```

(That path is where this skill is packaged inside each fork repository. If the skill was
loaded from somewhere else, run the `scripts/install.sh` that sits next to this SKILL.md.)

The script validates its arguments and required tools, auto-detects dd/gmp/qd from
`configure.ac` (and refuses a solver/repository mismatch), resolves the right compiler, applies
the SPOOLES rescue when needed, **verifies OpenMP linkage matches what was requested**, and
solves the bundled example, checking the objective (≈ −4.19e+01). Logs go to a temporary
directory (named in every FATAL message), the checkout stays clean, QD sources are
sha256-verified and cached under `~/.cache/sdpa-omp/`, and `DONE` is printed only after
every requested operation — including `--prefix` installation — has succeeded. Treat its
FATAL messages as authoritative — each names the trap it corresponds to.

To install a sibling solver, clone it first (they are separate repositories):
`git clone https://github.com/Canonical111/sdpa-{dd,gmp,qd}-omp.git`

## After installing

- Run with `OMP_NUM_THREADS=<physical cores>`. On Linux, pin: `taskset -c <cores>` plus
  `OMP_PROC_BIND=true OMP_PLACES=cores`. Hybrid CPUs: count P+E physical cores, not threads.
- Verify a solve prints identical `Iteration =` and `objValPrimal` at 1 thread and N threads
  — these forks are trajectory-stable by design; upstream is not.
- BENCHMARKS.md in the repo has measured expectations per machine class.

## Troubleshooting (symptom → cause → action)

| symptom | cause | action |
|---|---|---|
| macOS: configure succeeds but binary has no gomp | `/usr/bin/gcc` is Apple clang: **silent serial build** | `brew install gcc`; the script refuses this state — never bypass that check |
| `c99: illegal option -- f` inside SPOOLES | upstream's `.POSIX:` Make.inc defaults CC to `c99` | script's SPOOLES rescue; if manual, also `mkdir -p external/i/SPOOLES/lib` before the `cp` |
| `storage size of 'TZ' isn't known` | strict toolchain hides `struct timezone` even under gcc | `-D_GNU_SOURCE` in SPOOLES CFLAGS (script does this) |
| link: `Graph_free`/`IV_entries` undefined | incomplete `spooles.a` from an earlier failed pass | delete SPOOLES `.o`s and rebuild it (script's rescue does) |
| qd: `config.sub ... failed` on arm64 | 2009 `config.guess` predates Apple Silicon | copy `config.guess`/`config.sub` from automake's share dir |
| qd macOS link: `operator<<(ostream&, qd_real)` undefined | Homebrew qd bottle is clang/libc++; g++ uses libstdc++ | build QD from source with the same GCC (script does; `QD_DIR` overrides) |
| qd: fixes vanish on retry | failed make re-untars SPOOLES over your edits | fix only after first failure, then `cp spooles.a libspooles.a` to satisfy the target |
| gmp macOS: GMP test-suite `Bus error: 10` | bundled GMP 6.2.1 predates arm64 | `--with-system-gmp` against `brew install gmp` (script does) |
| threaded run not faster | binary is serial, or threads unpinned/oversubscribed | check gomp linkage (qd is static: use `nm \| grep GOMP`, not ldd); pin cores |
| qd timings look worse with more threads *in old builds* | upstream's clock reported process CPU time | this fork fixes it; distrust any historical sdpa-qd self-reported timing |

## Facts worth knowing while assisting

- dd/gmp need `autoreconf -fi` (upstream ships no `configure`); qd tracks its `configure`.
  autoreconf refreshes two *tracked* boilerplate files (`INSTALL`; on gmp also `aclocal.m4`),
  and on macOS the qd build replaces tracked `config.guess`/`config.sub`/`configure` with
  arm64-aware versions — the only working-tree changes a successful run leaves. Expected;
  `git checkout --` them if unwanted, but leave `aclocal.m4` and the qd config files alone
  before building again in the same tree.
- A serial build (`--serial`) is legitimate and CI-checked, not a failure mode.
- Expected example objective: dd/qd `-4.1899999999999999e+01`, gmp `-4.1900000000000000e+01`
  (gmp at its default 256-bit working precision).
- The three forks share design but not build systems; do not copy fixes across them blindly —
  SPOOLES lives at `external/spooles/work/internal` (dd, gmp) vs `spooles/build` (qd).
