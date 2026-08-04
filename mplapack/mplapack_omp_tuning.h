#ifndef MPLAPACK_OMP_TUNING_H
#define MPLAPACK_OMP_TUNING_H

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: added the Left-side triangular-kernel
   work gates used by Rtrsm_omp/Rtrmm_omp. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: the two triangular work gates were
   calibrated on hardware and the thread-count scaling replaced; the constants that were here
   before had never been measured. See the block below and git log. */

/* Bump when macros are added or removed; the generator refuses a stale header.
   4: MIN_TRSM_WORK_8T/MIN_TRMM_WORK_8T renamed to MIN_TRSM_WORK/MIN_TRMM_WORK -- the
      measurement showed the break-even does not track the team size, so the "_8T" in the
      name (meaning "the value at 8 threads, to be scaled") no longer describes them. */
#define MPLAPACK_OMP_TUNING_VERSION 4

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

/* ---------------- Left-side triangular kernels: Rtrsm_omp / Rtrmm_omp ----------------

   These two are the Cholesky-inverse phase -- Rtrsm against an identity RHS inside
   Lal::getInvLowTriangularMatrix, and Rtrmm forming Z^-1 = L**T * L inside
   Jal::getInvCholAndInv. Both are parallelised over the n COLUMNS of B (see the header
   comment in either file for why it is columns and not rows).

   A work gate here is mandatory, not a nicety: below the crossover the fork/join is not
   repaid. The gate is on WORK, m*m*n (the dd multiply-add count to within a factor of two),
   rather than on one dimension, so a short wide B is judged by what it actually costs.

   CALIBRATION, 2026-08-04. Measured by calling these two kernels directly, with their gate
   removed, on the operands the solver actually hands them (Rtrsm: A = the lower Cholesky
   factor, B = identity, alpha = 1; Rtrmm: A = B = L^-1 lower triangular, alpha = 1). The
   serial baseline was taken in a separate process with OMP_NUM_THREADS=1: with a team alive
   its idle workers spin in the GOMP barrier and get co-scheduled onto the master's SMT
   sibling, which doubled the measured serial time at 8 threads and manufactured impossible
   12-14x speedups. Threads were bound one per physical core, and each size was warmed for
   0.4 s first -- thanos runs the schedutil governor with a 1500 MHz idle floor against a
   3100 MHz max, and without the warmup two identical runs disagreed by 30%. Statistic is
   the min over repetitions.

   Speedup (t_serial / t_threaded), thanos = AMD EPYC 7232P (8 physical cores),
   pi = i9-13900K (8 P-cores + 16 E-cores = 24 physical):

                       Rtrsm                                 Rtrmm
       n    thanos 2t/4t/8t   pi 2t/4t/8t/24t    thanos 2t/4t/8t   pi 2t/4t/8t/24t
       8     0.89 0.64 0.42   0.75 0.51 0.35 0.32  0.95 0.78 0.54   0.59 0.42 0.30 0.29
      10     1.03 0.84 0.52   0.98 0.76 0.49 0.42  1.12 0.99 0.71   0.80 0.71 0.50 0.44
      12     1.20 1.19 0.91   0.89 0.96 0.75 0.55  1.22 1.41 1.34   0.96 1.02 0.81 0.60
      16     1.68 1.78 1.67   1.14 1.37 1.21 0.83  1.73 2.12 2.30   1.10 1.51 1.52 1.05
      20     1.82 2.43 2.54   1.30 1.69 1.76 1.16  1.78 2.64 3.50   1.28 1.90 2.17 1.50
      25     1.90 2.84 3.37   1.40 2.04 2.48 1.66  1.80 2.98 4.70   1.37 2.22 3.14 2.21
      35     1.96 3.29 4.80   1.62 2.42 3.44 3.15  1.84 3.33 5.79   1.53 2.65 4.37 4.40
      50     1.98 3.64 6.33   1.95 3.30 4.59 6.25  1.86 3.52 6.58   1.83 3.20 5.06 7.74
     100     2.00 3.90 7.58   1.91 3.82 6.48 10.26 1.86 3.68 7.19   1.90 3.72 6.35 10.25
     161     2.00 3.97 7.83   2.00 3.99 6.95 13.54 1.86 3.69 7.32   1.94 3.85 6.68 13.50

   The break-even is at n = 10-13 on thanos and n = 12-18 on pi -- FAR below the n >= 64 /
   n >= 45 that stood here before, which had never been measured. Those old values left the
   entire Rtrsm of theta1 (whose one and only block is n=50, measured 6.3x at 8 threads)
   running serially.

   The thresholds below are therefore NOT the break-even; they are the break-even plus a
   wide margin, placed at n = 25. Two reasons that specific number:

     - margin: n=16 still loses at 24 threads (0.83x). n=20 is the smallest tabulated size
       that wins at EVERY (machine, thread-count) pair, but only by 1.16x in its worst cell
       -- too close to the crossover to sit a gate on. n=25 is the first size with real
       headroom: its worst cell over all fourteen (machine, thread-count, kernel)
       combinations above is 1.37x.
     - it costs nothing: the SDPLIB blocks the solver actually sees are 1, 5, 8, 10, 19,
       35, 50, 70, 100, 124, 150, 161, 294, 800, 801, 1600. Nothing lies between 19 and 35,
       so any threshold in [20, 35) admits and rejects exactly the same set of real blocks.
       n=25 sits in the middle of that dead zone.

   End-to-end on pi (idle, min of 5 runs, whole-solver wall time relative to a build whose
   gate always rejects), for this gate ("cand"), the old n>=64/n>=45 ("brief"), and no gate
   at all ("none"):

                        8 threads                 24 threads
                  cand   brief   none        cand   brief   none
       hinf10     0.99    1.00   0.55        1.00    1.00   0.52
       truss5     0.99    1.01   0.95        0.99    0.97   0.86
       truss8     0.99    1.00   1.02        1.00    1.02   0.99
       theta1     1.09    1.02   1.09        1.12    1.05   1.07
       theta3     1.04    1.03   1.05        1.02    1.03   1.04
       gpp124-1   1.17    1.17   1.17        1.27    1.31   1.23
       arch0      1.26    1.35   1.33        1.32    1.32   1.35

   That is the whole case in one table. Ungating costs 2x on hinf10 and up to 14% on
   truss5, so the gate is required. The old thresholds cost theta1 half its gain, because
   n >= 64 rejected its only block (n=50) for Rtrsm -- 1.02x where this gate gets 1.09x at
   8 threads and 1.05x vs 1.12x at 24. Everywhere both gates admit the same blocks
   (gpp124-1, arch0) they agree to within run-to-run noise, which on these problems is
   several percent -- do not read the arch0 8-thread column as a real difference.

   NOT scaled by team size within 2..24 threads. The previous linear scaling assumed the
   break-even work grows with nt; the table shows it does not usefully do so, and worse, no
   monotone function of nt fits -- the tightest case is 24 threads at small n but 2 threads
   at large n (pi at 2 threads only reaches 2.0x at n=161). A flat threshold validated at
   2, 4, 8 and 24 is the honest reading. Beyond 24 threads there is no measurement, so the
   threshold is grown linearly there as a guard rather than extrapolated silently.

   The two kernels get the same number. Rtrmm does cross slightly earlier than Rtrsm (Rtrsm
   carries a division by the diagonal and a longer dependent chain per column), but the gap
   between them is smaller than the thanos-to-pi spread, and no benchmark block lands
   between the two crossovers, so splitting them would be false precision. */

/* Rtrsm, Left/Lower/NoTranspose: minimum work to thread. n = 25 cubed. */
#ifndef MPLAPACK_OMP_MIN_TRSM_WORK
#define MPLAPACK_OMP_MIN_TRSM_WORK 15625.0
#endif

/* Rtrmm, Left/Lower/Transpose: minimum work to thread. n = 25 cubed. */
#ifndef MPLAPACK_OMP_MIN_TRMM_WORK
#define MPLAPACK_OMP_MIN_TRMM_WORK 15625.0
#endif

/* Flat across the measured range (2..24 threads); grown linearly past 24, where nothing
   has been measured, so an unvalidated team size errs towards the serial kernel. */
#ifndef MPLAPACK_OMP_TRI_WORK
#define MPLAPACK_OMP_TRI_WORK(w, nt) ((nt) <= 24 ? (double)(w) : (double)(w) * (double)(nt) / 24.0)
#endif

/* Minimum number of columns to share out. Work alone is not sufficient: a tall thin solve
   with n=1 can clear the work gate yet offer a single iteration to the team. */
#ifndef MPLAPACK_OMP_MIN_TRI_WIDTH
#define MPLAPACK_OMP_MIN_TRI_WIDTH 2
#endif

#endif /* MPLAPACK_OMP_TUNING_H */
