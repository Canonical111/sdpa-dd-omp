# Benchmark evidence

The **current headline figures** in [README.md](../README.md), [INSTALL.md](../INSTALL.md) and
[BENCHMARKS.md](../BENCHMARKS.md) are derived from the rows in
[dd-port3-2026-08-24/](dd-port3-2026-08-24/). Each `.tsv` has a `.tsv.meta` beside it describing
**the provenance that was captured for that experiment, and the known limits of it**. The two
paired campaigns — `dd_v3_scaling_1t_24t` and `dd_upstream_vs_v3` — carry the strongest record,
including exact binary, compiler, parameter, input, affinity, route and repeat information. Some
earlier diagnostic campaigns retain a narrower record; their metadata says what it has rather than
claiming a full tuple. Provenance is not manufactured after the fact.

**The dE4 headline is quoted as ~7.2×, not to three figures.** Nine runs of the same protocol on
the same host across three campaigns span 6.234–6.761 s — an **8.2%** spread — giving campaign
ratios of 7.13×, 7.24× and 7.25× and a pooled median of 7.21×. No single campaign's median deserves
three significant figures, and the number moved every time it was re-measured until this was said
out loud.

**[validate.py](validate.py) recomputes the listed current headline figures from these rows** and
fails if any one of them no longer follows. CI runs it on every push, so those tables and the
evidence cannot drift apart silently.

Read its total honestly: it is the number of **all** checks — structure, parsing, blank fields,
return codes, repeat coverage, determinism **and** numerical derivation — not a count of
independently recomputed figures. The derivation block covers the current headline, the scaling
table, every cell of the upstream table, all seven `fill` structures, and the memory claim. It does
**not** regenerate the older per-machine historical tables, which are backed by other public TSVs
outside this folder.

```bash
python3 bench/validate.py
```

## What is here

| file | backs |
|---|---|
| `dd_final_headline_postfix.tsv` | one of three campaigns behind the dE4 headline (its own median: 6.495 s, 7.13×) |
| `dd_v3_scaling_1t_24t.tsv` | 1→24 thread scaling, both ends from one `v7.1.3-omp.3` binary |
| `dd_current_scaling.tsv` | the same, re-measured on current `main` after the route promotion — INSTALL's table |
| `dd_upstream_vs_current.tsv` | the upstream comparison in README and BENCHMARKS — current `main` vs upstream, both on their default routes |
| `dd_upstream_vs_v3.tsv` | its predecessor, **superseded**: measured before the route promotion, so its dE3 fork row is a route the fork no longer defaults to |
| `dd_fill_seven_structures.tsv` | INSTALL's `fill` table — all seven route-switch structures |
| `dd_bitset_memory.tsv` | INSTALL's memory claim, before and after the two-bit overlap map |
| `dd_token_length.tsv` + `truncate_significant_digits.py` | the reader defect in BENCHMARKS, isolated at 45 significant digits |
| `dd_race_truss6.tsv`, `dd_pi_truss6.tsv` | the data race that `v7.1.3-omp.3` fixes, on aarch64 and x86-64 |
| `dd_de4_stream_census.tsv` | how often that race actually fired on dE4 — once in seventeen threaded runs |
| `dd_gate_separation_{m1,pi}_4iter.tsv` | why the assembly gate was **not** retuned |
| `dd_alloc_vs_barrier.tsv` | why the small-block loss is synchronisation and not allocation |
| `dd_historical_baseline.tsv` | the pre-threading baseline the dE4 headline is measured against — so both sides of that ratio come from published rows |

## Reproducing a row

The protocol is the same throughout: `main loop time` over a **fixed 4-iteration budget** — the
shipped `param.sdpa` with `maxIteration` set to 4 — pinned with `OMP_PROC_BIND=true
OMP_PLACES=cores` and `taskset`, medians of three repeats. Each `.meta` gives the exact build and
hashes for its own rows.

Generated inputs come from the scripts in [../tests/](../tests/):
`gen_crossblock_fixture.py`, `gen_spchol_fixture.py`, `gen_route_fixture.py`.

## What is not here, and why

**The `dE3` / `dE4` inputs and the `input/sdpaquestions/` problems are not published.** They are
research data belonging to a conformal-bootstrap campaign, not ours to redistribute. The
consequence is stated plainly rather than hidden: **the dE3/dE4 and user-problem rows cannot be
reproduced from this repository alone.** Their shapes are recorded — dE3 is m=6067 and dE4 m=7401,
both 17 blocks; the seven `fill` structures are m = 2439, 4489, 5278, 6067, 8359, 10614, 11227, all
17 blocks — and every `.meta` carries the input's sha256, so a holder of the same file can confirm
they are running what we ran.

Three different things get conflated by the word "reproducible", so to be exact about what this
archive gives you:

- **auditable** — you can recompute the published summaries from the published rows. True for
  everything here.
- **identifiable** — if you hold the same input, its sha256 in the metadata confirms it is the one
  we ran. True for the private-input campaigns.
- **reproducible** — a fresh public clone can rerun the experiment. True only for the public-input
  campaigns.

Everything measured on **public** inputs — SDPLIB (`truss6`, `arch0`, `gpp100`, `control1`) and the
generated fixtures — is fully reproducible here, and that includes the data race itself: `truss6`
at default settings is the clearest demonstration of the defect `v7.1.3-omp.3` fixes.
