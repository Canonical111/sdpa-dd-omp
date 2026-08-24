# Benchmark evidence

Every performance figure in [README.md](../README.md), [INSTALL.md](../INSTALL.md) and
[BENCHMARKS.md](../BENCHMARKS.md) is derived from the rows in
[dd-port3-2026-08-24/](dd-port3-2026-08-24/). Each `.tsv` has a `.tsv.meta` beside it recording
host, build commit and binary hash, compiler, parameter and input hashes, affinity, route and
overlap policy, iteration budget, and repeat count.

**[validate.py](validate.py) recomputes the published figures from these rows** and fails if any
one of them no longer follows. CI runs it on every push, so the tables and the evidence cannot
drift apart silently.

```bash
python3 bench/validate.py
```

## What is here

| file | backs |
|---|---|
| `dd_final_headline_postfix.tsv` | the README headline — dE4 **6.495 s**, **7.13×** over the historical baseline |
| `dd_v3_scaling_1t_24t.tsv` | INSTALL's 1→24 thread table, both ends from one `v7.1.3-omp.3` binary |
| `dd_upstream_vs_v3.tsv` | the upstream comparison in README and BENCHMARKS, 3 repeats, builds interleaved |
| `dd_fill_seven_structures.tsv` | INSTALL's `fill` table — all seven route-switch structures |
| `dd_bitset_memory.tsv` | INSTALL's memory claim, before and after the two-bit overlap map |
| `dd_token_length.tsv` + `truncate_significant_digits.py` | the reader defect in BENCHMARKS, isolated at 45 significant digits |
| `dd_race_truss6.tsv`, `dd_pi_truss6.tsv` | the data race that `v7.1.3-omp.3` fixes, on aarch64 and x86-64 |
| `dd_de4_stream_census.tsv` | how often that race actually fired on dE4 — once in seventeen threaded runs |
| `dd_gate_separation_{m1,pi}_4iter.tsv` | why the assembly gate was **not** retuned |
| `dd_alloc_vs_barrier.tsv` | why the small-block loss is synchronisation and not allocation |

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

Everything measured on **public** inputs — SDPLIB (`truss6`, `arch0`, `gpp100`, `control1`) and the
generated fixtures — is fully reproducible here, and that includes the data race itself: `truss6`
at default settings is the clearest demonstration of the defect `v7.1.3-omp.3` fixes.
