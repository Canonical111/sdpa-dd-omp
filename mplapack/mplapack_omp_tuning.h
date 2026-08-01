#ifndef MPLAPACK_OMP_TUNING_H
#define MPLAPACK_OMP_TUNING_H

/* Bump when macros are added or removed; the generator refuses a stale header. */
#define MPLAPACK_OMP_TUNING_VERSION 2

/* Tuning knobs for the mplapack OpenMP kernels. Override at compile time with -D. */

/* Minimum gemm work (m*n*k dd multiply-adds) before an OpenMP fork/join pays for itself.
   A dd multiply-add is ~25-50 ns; a fork/join across 24 threads costs ~5-20 us. Requiring
   ~1 ms of work keeps the overhead below a few percent. */
#ifndef MPLAPACK_OMP_MIN_GEMM_WORK
#define MPLAPACK_OMP_MIN_GEMM_WORK 20000.0
#endif

/* Minimum width of the parallelised gemm loop (over j, so n iterations). Work alone is not
   sufficient: a tall thin gemm with n=1 can clear MIN_GEMM_WORK yet offer one iteration to
   share out, paying team-creation cost for no parallelism. */
#ifndef MPLAPACK_OMP_MIN_GEMM_WIDTH
#define MPLAPACK_OMP_MIN_GEMM_WIDTH 2
#endif

/* Minimum vector length before parallelising axpy/copy. */
#ifndef MPLAPACK_OMP_MIN_VECTOR_N
#define MPLAPACK_OMP_MIN_VECTOR_N 20000
#endif

#endif /* MPLAPACK_OMP_TUNING_H */
