# sdpa-dd-omp

OpenMP-threaded, reproducible fork of sdpa-dd. ("-omp" and not "parallel": SDPARA is the
SDPA group's own MPI solver, and "mp" already means multiprecision in this family.)

Fork of [nakatamaho/sdpa-dd](https://github.com/nakatamaho/sdpa-dd) at commit
`6eaad8d9`, carrying four patches that fix performance and reproducibility
regressions. The issues were reported to the upstream maintainer and have not
been adopted there; this fork exists so the fixes are usable.

## What changed (4 commits on top of upstream)

| commit | change |
|---|---|
| 1 | restore netlib dgemm's zero-skip in `Rgemm_NN/NT_omp` — dropped in mplapack 2.0.1, up to 2× serial regression |
| 2 | gate OpenMP on work / width / nesting; run `Rdot` serially — up to 13× on small problems, removes a nondeterminism |
| 3 | thread the Schur-complement (`bMat`) constraint loop — the dominant gain |
| 4 | report the reduced per-formula timers as worker-seconds, not elapsed time |

Every modified file carries an in-file, dated change notice (GPLv2 §2a):
`sdpa_newton.cpp`, `sdpa_parts.cpp`, `mplapack/Rgemm_{NN,NT,TN,TT}_omp.cpp`,
`mplapack/{Rdot,Raxpy,Rcopy}_omp.cpp`; `mplapack/mplapack_omp_tuning.h` is new. `git log`
has the full rationale per change.

## Measured results

20 SDPLIB problems (m = 21…1106), external wall clock, median of 3 pinned repeats:

| | i9-13900K (24 cores) | M1 Max (8P+2E) |
|---|---|---|
| sdpa-dd 7.1.2 (2009) | 309.6 s | 360.0 s |
| upstream master, threaded | 136.7 s | 112.4 s |
| **this fork** | **58.8 s** | **86.0 s** |
| vs upstream (identical iteration counts, 20/20) | **2.33×** | **1.31×** |

Results are **thread-count independent**: a checksum over every scalar of the final
`xMat`/`yVec`/`zMat` is identical at 1/2/4/8/24 threads. Unpatched upstream is *not*
reproducible when threaded (iteration counts and even objectives vary run to run).

Full methodology, raw per-repeat data, and the investigation history live in the
companion repository (benchmarks, harness, and the patch generators).

## Building

```bash
autoreconf -fi                       # upstream ships configure.ac only, no configure
./configure --enable-openmp=yes
make -j$(nproc)
```

Verified on a fresh Ubuntu 24.04 clone (gcc 13) and by CI on every push. **If** `make` stops
inside SPOOLES with a `struct timezone` error (seen with some toolchains — upstream's
`.POSIX:` Make.inc selects `c99` as the compiler):

```bash
sed -i 's|^# CC = gcc|  CC = gcc|' external/spooles/work/internal/Make.inc
( cd external/spooles/work/internal && find . -name '*.o' -delete && rm -f spooles.a \
  && make global -f makefile && cp spooles.a ../../../i/SPOOLES/lib/libspooles.a )
make -j$(nproc)
```

Run with `OMP_NUM_THREADS=<physical cores>`, pinned (`taskset`/`OMP_PLACES=cores`).
Serial builds (`--enable-openmp=no`) are supported and CI-checked.

## License

GPL v2, unchanged from upstream (`COPYING`). Copyright remains with the original
SDPA authors; the patches in this fork are contributed under the same license.
