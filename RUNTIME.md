# sdpa-dd-omp runtime reference

Every environment variable this fork adds, what it does, and whether it can change results.

Two rules hold for all of them, and are tested in CI rather than merely stated:

1. **A typo is refused, never treated as the default.** Setting `SDPA_SPCHOL_MODE=Auto` or
   `SDPA_DD_MIN_SPCHOL_WORK=-5` stops the run with a diagnostic. Silently falling back would let a
   mistyped variable select the very path the caller was trying to avoid — and `-5` in particular
   used to be *accepted*, wrapping to 18446744073709551611, i.e. "never parallel" while looking
   configured.
2. **Validation happens on every route.** A knob is parsed once per process at a point both the
   dense and sparse routes reach, so a dense-route problem refuses a malformed value it would
   never otherwise consult. This fork has been bitten five times by checks placed where only one
   route could reach them.

Anything marked **test hook** requires `-DSDPA_SPCHOL_TEST_HOOKS` at build time; a release build
refuses it rather than ignoring it, because a knob that silently does nothing misleads whoever set
it.

---

## Choosing how the Schur complement is stored and factored

| variable | values | changes results? |
|---|---|---|
| `SDPA_BMAT_MODE` | `auto` (default) · `fill` · `dense` · `sparse` | **yes** — see below |
| `SDPA_BMAT_MAX_GB` | positive number of GB | no — a safety cap; refuses to allocate beyond it |
| `SDPA_BMAT_LOG` | `1` on, `0` off | no — but see the advisory below |

`auto` is the released chooser, unchanged. `fill` uses the fill-derived rule, which routes some
problems to the sparse factorisation that `auto` sends dense.

**`SDPA_BMAT_MODE` is result-changing by design.** Dense and sparse are different factorisation
routes, so they follow different iterate trajectories: across a 48-cell harness 30 cells differ
between routes, while every cell that *converges* agrees on the objective. If you need
bit-reproducibility against an earlier run, keep the mode fixed.

`SDPA_BMAT_MAX_GB` fails **closed**: an empty value is an error rather than "no cap", because
`SDPA_BMAT_MAX_GB="$TYPO"` is the ordinary way an empty value reaches a job script and treating it
as no-cap would fail open on a safety limit.

**Available memory is reported, never acted on.** With `SDPA_BMAT_LOG=1`, and in the message a
`SDPA_BMAT_MAX_GB` refusal prints, the dense requirement is shown next to what the machine has
free, so "needs 21 GB" can be read against something. Nothing branches on it: a route chosen from
free memory would be a route that changes with the machine, which is the opposite of what this
fork guarantees. The figure is Linux-only (`MemAvailable` from `/proc/meminfo`; elsewhere the line
says it is not readable rather than guessing) and it reports the **host**, so inside a memory
cgroup it is an upper bound.

## Threading the sparse Cholesky

| variable | values | changes results? |
|---|---|---|
| `SDPA_SPCHOL_MODE` | `auto` (default) · `serial` · `force` · `legacy` | no |
| `SDPA_DD_MIN_SPCHOL_WORK` | integer, default **10000** | no |
| `SDPA_DD_MIN_SPCHOL_WIDTH` | integer, default **8** | no |
| `SDPA_DD_MIN_SPCHOL_TOTAL` | integer, default **0** (pass-through) | no |
| `SDPA_SPCHOL_LOG` | `1` on, `0` off | no |

The factorisation is **bit-identical at any thread count, by construction**: for a fixed pivot each
`k1` owns one destination row, every read comes from the pivot row, and a barrier per pivot
preserves pivot order. That is a property of the loop's dependence structure, not of the
arithmetic, so it does not depend on precision or on the gate values.

`legacy` runs an independent transcription of the pre-refactor expression. It exists as an
**oracle**: `serial` and the threaded path share their inner kernel, so comparing them cannot
detect a bug in that shared kernel, whereas `legacy` can.

`force` **forces or fails** — in a build without OpenMP it is refused rather than silently running
serial, so a benchmark cannot believe it measured the threaded path when it measured the other one.

The gates were calibrated on dd, on two architectures (i9-13900K/24 cores and Apple M1 Max/10
cores). They are runtime-overridable precisely so a machine that disagrees can be accommodated
without a rebuild.

## Threading the sparse bMat assembly

| variable | values | changes results? |
|---|---|---|
| `SDPA_BMAT_ASM_MODE` | `auto` (default) · `serial` · `parallel` | no |
| `SDPA_BMAT_ASM_MIN_PAIRS` | integer, default **8000** | no |
| `SDPA_BMAT_ASM_SCRATCH_MB` | integer MB, default **4096** | no — bounds memory, may reduce the team |
| `SDPA_BMAT_ASM_LOG` | `1` | no |
| `SDPA_BMAT_ASM_PROFILE` | `1` | no |
| `SDPA_BMAT_ASM_CENSUS` | `1` | no |

Parallel over i-groups, serial over blocks. Each worker gets private `work1`/`work2` matrices, so
concurrent groups cannot race through shared scratch; the accumulation into the Schur complement is
disjoint within a block. **Bit-identical**, and asserted by comparing the assembled matrix itself
across thread counts, not the printed solution.

`SDPA_BMAT_ASM_SCRATCH_MB` is a **production safety control**, not observability: each worker beyond
the first needs two dense matrices of the largest block's order. If the budget admits fewer workers
than requested, the team shrinks; it never silently exceeds the budget.

`SDPA_BMAT_ASM_MIN_PAIRS = 8000` was calibrated on two architectures. Below roughly 5,000 pairs
threading loses (the fork/join costs more than the work); above roughly 15,000 it wins by 8–10×.
The gate is a **proxy** — pair count does not fully determine the work, since a group's cost also
scales with its block dimension — so a problem near the boundary may benefit from overriding it.

## Reproducibility oracles

| variable | what it does |
|---|---|
| `SDPA_SPCHOL_DIGEST` | fingerprint of the finished factor: records, bytes, 64-bit FNV-1a |
| `SDPA_SPCHOL_DIGEST_DUMP=<file>` | the factor's canonical stream, appended — compare two runs with `cmp` |
| `SDPA_BMAT_ASM_DIGEST` | the same, for the assembled Schur complement |
| `SDPA_BMAT_ASM_DUMP=<file>` | the assembly's canonical stream |

The dumps are the proof-grade comparison; a fingerprint is strong evidence, not proof, and is
described that way wherever it is printed. Both streams carry a length frame before every
variable-length field, so two different structures cannot serialise to the same bytes, and they use
different tags so an assembly stream can never compare equal to a factor stream.

Use these rather than the printed solution when you need to know whether a change altered the
computation: a `dd_real` carries about 32 significant digits and the printed fields carry 17, so a
factor can change without any printed field moving.

## Test hooks

Require `-DSDPA_SPCHOL_TEST_HOOKS`; refused by a release build.

| variable | what it injects |
|---|---|
| `SDPA_SPCHOL_MUTATE=1` | one ulp into the finished factor — the digest's negative control |
| `SDPA_BMAT_ASM_MUTATE=1` | one ulp into the assembled matrix — the assembly oracle's negative control |
| `SDPA_BMAT_TEST_BREAK_INVARIANT=1` | a violation of the fill ≥ aggregate invariant |
| `SDPA_BMAT_TEST_F2_STALE=1` | the worst-case stale-`G` read a since-fixed defect could produce |
| `SDPA_SPCHOL_TEAM_OVERRIDE=n` | asks the runtime for `n` threads while leaving the gates' team alone. `1` is the only way to reach the "requested a team and received one thread" fallback |
| `SDPA_SPCHOL_FAIL_AT=i` | declares pivot `i` indefinite, so failure can be returned from a factorisation that has already run threaded pivots |

These exist because a comparison that cannot fail is worse than no comparison: it reports success.
Each has a CI assertion requiring it to actually change what it claims to change.

The last two are the other half of that idea: not negative controls for an oracle, but the only
way to execute two branches no input reaches on purpose. `TEAM_OVERRIDE=1` reaches the
one-member-team fallback — the dispatcher's own `team < 2` test exits earlier and by a different
path, so without the hook that branch is dead to every test — and CI additionally requires the
fallback's factor to be byte-identical to the threaded one. `FAIL_AT` makes the failure return
from a partly-threaded factorisation executable; CI requires threaded pivots to have run first,
or the test would prove nothing about that path. When it fires, the diagnostic says the failure
was **declared**, not measured.

## Profiling

`SDPA_BMAT_ASM_PROFILE=1` splits the `Make bMat` phase into zeroing, SDP assembly and LP assembly.
`SDPA_SOLVE_PROFILE=1` splits the triangular solve into its forward and backward passes. Both are
reported in the phase table and cost a few clock reads per iteration when off.

`SDPA_BMAT_ASM_CENSUS=1` reports, per SDP block: pair count, group count, largest group, the
formula mix, and how many pairs consume the F2/dense-`Aj` path. Useful for judging whether a
problem can balance across threads before trying.
