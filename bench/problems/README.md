# Benchmark problems that are not in SDPLIB

The tables in [`../../BENCHMARKS.md`](../../BENCHMARKS.md) mix classic SDPLIB problems
(`control1`, `theta1`, `gpp100`, `truss5`, `truss6`, `arch0` — get those from the
[SDPLIB collection](https://github.com/vsdp/SDPLIB)) with two large sparse problems from
conformal-bootstrap work that SDPLIB does not carry. **This directory ships those two**, so every
problem named in the tables is obtainable.

| problem | m | blocks | size (raw) | what it is used for |
|---|---:|---:|---:|---|
| `dE3` | 6067 | 17 | 7.5 MB | the route-promotion case: dense under upstream, sparse here since 2026-08-24 |
| `dE4` | 7401 | 17 | 7.5 MB | the headline: routes sparse under every policy, and upstream threads none of that path |

The sibling [`sdpa-gmp-omp`](https://github.com/Canonical111/sdpa-gmp-omp) publishes the same two
files plus a `min` series (`8_min`, `10_min`, `12_min`). **The `min` series is deliberately not
here**: no table in this repository measures it, and shipping inputs nothing cites would be
clutter rather than provenance. The two files here are byte-identical to that fork's copies —
the same `dE3.dat-s` and `dE4.dat-s` hashes appear in both repositories' checksum files.

## Use

```bash
xz -dk dE3.dat-s.xz dE4.dat-s.xz
shasum -a 256 -c SHA256SUMS          # sha256sum -c on Linux
```

Run exactly as benchmarked. **`-p` is mandatory** — without it the solver falls back to
compiled-in defaults and you are not running the published protocol:

```bash
OMP_PROC_BIND=true OMP_PLACES=cores OMP_MAX_ACTIVE_LEVELS=1 \
  taskset -c 0-23 OMP_NUM_THREADS=24 \
  ./sdpa_dd -ds dE4.dat-s -o out.result -p bench/problems/param_dd_4iter.sdpa
```

`param_dd_4iter.sdpa` is the repository's own `param.sdpa` with `maxIteration` set to **4**. Every
`dE3`/`dE4` figure published here is `main loop time` over that fixed four-iteration budget —
**per-iteration solver cost, not time to solution.** Neither problem converges at double-double
precision under either tolerance tested, so a four-iteration budget is what there is to compare.

## What is proven here, and what is not

Being exact, because a provenance directory that overstates itself is worse than none.

**The parameter file is provably the one the campaign used.** The metadata in
[`../dd-port3-2026-08-24/`](../dd-port3-2026-08-24/) records `param` as
`9197d394cb3907c83ac3c021af7eef9241952030d2647af2fb045eab8d308853`, and the file shipped here
hashes to exactly that. The documented recipe — take `param.sdpa`, set `maxIteration` to 4 —
reproduces it byte for byte, which is why the recipe can be trusted as written.

**The input files are *not* provably the ones the campaign measured.** No `.meta` in this archive
records an input hash — none of the fifteen — so nothing in the recorded provenance ties a
measured row to a specific `dE3.dat-s`. What can be said is what is true: these are the `dE3` and
`dE4` inputs from that campaign, byte-identical to the copies the sibling fork published and
measured, and their hashes are now recorded in `SHA256SUMS` so that every run from here forward
*can* be tied to them. An earlier version of `../README.md` claimed that every `.meta` carried the
input's sha256. It did not, and that claim has been removed rather than quietly softened.
