# Retired measurement probes

## `allocation-probe.patch`

**Question it answered:** the threaded sparse assembly loses badly on problems with many small
blocks (`truss6` runs its assembly ~19× slower threaded than serial). Two costs were confounded —
the per-block barrier, and `DenseMatrix::initialize` freeing and reallocating the per-worker
scratch whenever a block's area differs from the previous one. A gate recalibrated from a
confounded measurement would be worthless, so the two had to be separated before anything else.

**What the probe does.** Two independent compile flags:

- `SDPA_PROBE_NOREALLOC` — reserve per-worker capacity once at `max_dim` and reuse it, while
  keeping each block's **logical** dimensions and zeroing exactly the region the shipped path
  zeroes. It deliberately does **not** pad the computation to `max_dim`: that would change the
  flop count and memory traffic, and would measure the padding instead of the allocation.
- `SDPA_PROBE_ALLOC_COUNT` — count allocations that actually happen, on its own flag so the probe
  can be built with *and* without the reuse path. Without that separation, "removing reallocation
  changed nothing" is indistinguishable from "the probe removed no reallocations".

```bash
git apply bench/probes/allocation-probe.patch
./configure --enable-openmp=yes CXXFLAGS="-O2 -ffp-contract=off -DSDPA_PROBE_ALLOC_COUNT -DSDPA_PROBE_NOREALLOC"
```

**Answer — allocation is not a material cost.** Rows in
[../dd-port3-2026-08-24/dd_alloc_vs_barrier.tsv](../dd-port3-2026-08-24/dd_alloc_vs_barrier.tsv):

| | allocations removed | assembly time |
|---|---:|---:|
| `12_min` (14 distinct block sizes) | **−91%** | −0.3% |
| `10_min` (9) | −89% | +0.7% |
| `8_min` (6) | −86% | −3.1% |

Removing 86–91% of every `DenseMatrix` allocation in the process moves the threaded assembly by at
most 3.1%. And `truss6` — the case that prompted the question — turns out to be 150 blocks of order
3 plus one of order 1, so consecutive blocks almost always share an area and `initialize()` never
takes its reallocating path. Its loss is synchronisation, full stop.

**Precondition, checked before any timing was compared:** the probe build and the shipped build
produce byte-identical assembly streams on all eight test problems. A probe that changes the answer
measures nothing.

**Why it is a patch and not `#ifdef`s in the solver.** The question is answered. Leaving disabled
experimental branches and a global counter in `sdpa_newton.cpp` and `sdpa_struct.cpp` means
maintaining and reviewing them as solver code forever, for an experiment that concluded. The patch
keeps it reproducible without that cost.
