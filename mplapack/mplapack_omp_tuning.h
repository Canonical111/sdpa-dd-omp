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
      name (meaning "the value at 8 threads, to be scaled") no longer describes them.
   5: added MIN_POTRF_SYRK_WORK / MIN_POTRF_TRSM_WORK / POTRF_WORK / MIN_POTRF_WIDTH for the
      Rpotrf lower-path panel kernels (Rsyrk_omp, Rtrsm_omp's Right case). */
#define MPLAPACK_OMP_TUNING_VERSION 5

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

/* ------------- Rpotrf lower-path panel kernels: Rsyrk_omp / Rtrsm_omp (Right) -------------

   CALIBRATED, 2026-08-04/05, on both machines. Method as for the Left gates above: the two
   kernels were called directly with their gate removed, on the shapes Rpotrf actually
   produces (taken from a per-call census of 15 SDPLIB problems: 10263 Rpotrf calls, 8337
   panel iterations, 22888 logged lines), serial baseline in a separate OMP_NUM_THREADS=1
   process, threads one per physical core, each size warmed ~0.4 s first, statistic = min
   over repetitions. thanos = EPYC 7232P (8 physical), pi = i9-13900K (8 P + 16 E = 24).
   Full tables and the raw data: patches/b1_notes/b1_gate_calibration{,_pi}.md.

   STRUCTURE FIRST -- three facts that decide more than the constants do.

   * These kernels are reached only when the matrix order exceeds nb. iMlaenv returns nb = 64
     for POTRF (iMlaenv.cpp:222-225), and Rpotrf takes the unblocked Rpotrf2 recursion
     whenever nb >= n. So every Cholesky of order <= 64 bypasses both panel kernels entirely,
     whatever these constants say: 8728 of 10263 census calls (85%), including every Cholesky
     of hinf10, control1, truss5's blocks and theta1's n = 50 block. That fallback, not this
     gate, is what excludes the small per-SDP-block factorisations.

   * For Rsyrk the gate must be on WORK, not on a dimension: the panel Rsyrk is short and fat,
     n = jb <= 64 always while k grows to the whole factored prefix (up to 1536 in a 1600
     Cholesky), so a gate on n cannot tell the first panel from the last.

   * For the Right-side Rtrsm that argument is VACUOUS, and saying otherwise misled the first
     draft of this block. jb == 64 in every Right-side call without exception -- structurally,
     not as a workload accident: the call sits inside `if (j + jb <= n)` and jb < nb only on
     the last panel, so jb == nb == 64 whenever it runs (6807 of 6807 census calls). The work
     expression m*n*n is therefore identically 4096*m, and MIN_POTRF_TRSM_WORK IS a minimum
     row count however it is written. 131072 == m >= 32. Say "m >= 32" when explaining it.

   Work convention, matching the Left gate above: the expression compared against these
   constants is twice the true dd multiply-add count.
     Rsyrk_omp        k*n*n   vs true  k*n*(n+1)/2       (in Rpotrf terms: (j-1)*jb*jb)
     Rtrsm_omp Right  m*n*n   vs true  m*n*(n+1)/2       (in Rpotrf terms: (n-j-jb+1)*jb*jb)

   RSYRK -- the break-even is below every shape this solver produces.
   Speedup vs serial, thanos worst cell over four allocation offsets, k = 64 unless noted:

       jb   gate work      thanos 2t/4t/8t        pi (crossing 1000-2000)
        1          64      0.78 0.79 0.77   <- num_threads clamps to n: one thread, pure fork
        2         256      0.97 0.96 0.95
        3         576      1.44 1.24 1.24
        4        1024      1.70 1.70 1.61   0.96 worst at 2t, 1.32-1.42 at 24t
        6        2304      1.84 2.37 2.29   1.43 at 2t .. 1.93 at 4t, 1.68 at 24t
        8        4096      1.90 2.84 3.08
       18       20736      1.99 3.62 5.27
       32       65536      2.01 3.85 7.06
       40      102400      2.02 3.91 7.34
       46 (k=128) 270848   2.02 3.97 7.76

   thanos crosses between jb = 2 and 3 at k = 64 and below jb = 2 at k = 1024; pi crosses at
   gate work 1000-2000. 2000 is the conservative end of that: on pi every shape at or above
   2048 wins by >= 1.30x at every thread count, and the largest shape below it (1024) is
   0.96x worst-case. It admits 100% of the census (the smallest census shape is 2304,
   control7's jb = 6 k = 64, worth 1.43-2.41x). Its job is to guard shapes SDPLIB does not
   contain, not to reject anything SDPLIB has.

   RTRSM (Right) -- there IS a losing regime, and it is not fork/join.
   Speedup vs serial, worst cell; m is the row count, jb = 64 fixed:

       m      thanos 2t/4t/8t          pi
        6     0.95  0.94  0.94         1.07 at 4t, 0.88 at 24t
        8     1.10  0.87  0.84         1.55 at 4t, 0.93 at 8t
       16     1.44  0.98  1.37         first pi size winning everywhere, >= 1.20x
       18     1.53  0.97  1.44
       22     1.62  1.69  1.75
       26     1.72  1.36  1.98
       32+    1.83+ 1.90+ 1.51+        rises to 5.9-6.7x at m ~ 1000
   The m = 8 shape takes 167 us serially -- twenty fork/joins' worth -- and still loses at 8
   threads, and a 250 us serial gap before every call (cold team) does not move the crossing
   on pi. It is a memory/efficiency regime, not an overhead regime. Contributing causes,
   both located: the row split (t*m)/nblk is not cache-line aligned (4 dd_real rows per 64 B
   line, and the same shape reads 1.37x or 4.95x depending only on the address the allocator
   returned), and the n(n+1)/2 = 2080 scalar zero-tests/temp loads of the (k,j) loops are
   executed IN FULL by every thread, so that overhead grows with the team while the useful
   work per thread shrinks. Both are fixable; see patches/b1_notes/ for the follow-up.

   131072 (m >= 32) is the cross-machine conservative choice: thanos alone would take 131072,
   pi alone 65536 (m >= 16); the larger of the two is the only setting under which no admitted
   shape loses on either machine. It rejects census m in {6, 16, 18, 22, 26} -- 487 of 10482
   Right-side calls (4.6%), carrying 0.12% of the panel arithmetic -- giving up two mild wins
   (m = 22, 26) to remove three mild losses (m = 6, 8, 16/18 at 4t).

   THREAD-COUNT SCALING. Flat over 2..24 is kept, but it is NOT uniformly confirmed and the
   next person should know which half is which:
     - Rsyrk: flat is measured. Near the crossing it is worst at 2 threads on BOTH machines
       (more threads helps), so a flat gate is if anything conservative.
     - Rtrsm Right: thanos CONTRADICTS flat -- its break-even moves from m ~ 4 at 2 threads to
       m ~ 10-22 at 4 and 8, and m = 16/18 fall back to 0.97-0.98x at 4 threads after m = 10
       had already crossed. pi found flat adequate over 2..24. The flat form is made safe here
       by the CONSTANT, not by the form: m >= 32 sits above the worst break-even seen at any
       measured thread count on either machine. If either constant is ever lowered towards its
       break-even, this macro must be re-derived first.
   Beyond 24 threads nothing is measured on either machine, so the threshold is grown linearly
   there as a guard rather than extrapolated silently.

   DO NOT reintroduce a thread-count-scaled guard on the WIDTH (the "n < 4*nthreads" shape).
   With n = jb = 64 it rejects every panel at 24 threads, and at 24 threads those very shapes
   measure 9.6-12.5x on pi. (On the m axis of the Rtrsm gate, 4*nt at 8 threads happens to
   land on the same m >= 32 chosen here; that is a coincidence, and the constant above is
   flat, not scaled.)

   END TO END, whole-solver wall clock, ratio to a build whose panel gates always reject
   (= pre-B1). thanos 8 threads on cores 0-7, idle box, min of 3, interleaved, warm-up
   discarded, foreign-load-verified (2026-08-05); pi 8/24 threads, idle, min of 5:

       problem     thanos 8t     pi 8t    pi 24t
       theta2         2.07        1.90      2.20
       qap9           2.24        1.77      2.16
       theta3         2.41 [*]     --        --
       truss5         1.36        1.43      1.34
       control7       1.24        1.18      1.33
       theta1         1.12        1.16      1.18
       gpp124-1       1.10        1.04      1.13
       arch0          1.09        1.08      1.17
       hinf10         0.99-1.03   0.99      0.93   <- null control: order <= 21, never reaches
                                                      the blocked path; this is its noise floor
       [*] measured in the calibration leg (72.39 s -> 30.06 s), not in the verdict run.

   Two warnings about that table. (1) The gate constant is NOT what pays: builds with gates
   0 (always thread), 20000, and the values below agree to within 1-3% on every problem
   above, on both machines. Quote the B1 improvement -- the gate-off column against any of
   the others -- never a "gate improvement". The constants are justified by the KERNEL
   measurement. (2) Any such figure taken on a machine that is not exclusively yours is
   unproven: B1 adds ~2 OpenMP parallel regions per Cholesky panel (~1022 per theta2 solve),
   each ending in a team barrier, so one preempted worker stalls the whole team at a cost the
   serial baseline never pays. An earlier thanos run under a competing sweep read 0.68x on
   theta2 where the idle box reads 2.07x, with the baseline's own time inflated 2.2x.

   Finally, one number quoted above in MIN_GEMM_WORK's justification is wrong and should not
   be reused: a dd multiply-add is 2.6-4.1 ns on pi, not the 25-50 ns claimed there. The
   fork/join is ~2 us. Any "of order 1 ms of work" reasoning built on 25-50 ns is off by an
   order of magnitude. */

/* Rsyrk, Lower/NoTranspose: minimum work (k*jb*jb) to thread. Break-even is 256-1024 on
   thanos and 1000-2000 on pi; this is the conservative end of pi's crossing. Admits every
   panel shape SDPLIB produces. */
#ifndef MPLAPACK_OMP_MIN_POTRF_SYRK_WORK
#define MPLAPACK_OMP_MIN_POTRF_SYRK_WORK 2000.0
#endif

/* Rtrsm, Right/Lower/Transpose: minimum work (m*jb*jb) to thread. jb == 64 structurally, so
   this is exactly "m >= 32 rows to share out". thanos wants 131072, pi 65536; the larger is
   the only value under which no admitted shape loses on either machine. */
#ifndef MPLAPACK_OMP_MIN_POTRF_TRSM_WORK
#define MPLAPACK_OMP_MIN_POTRF_TRSM_WORK 131072.0
#endif

/* Flat across the measured range (2..24 threads); grown linearly past 24, where nothing has
   been measured on either machine, so an unvalidated team size errs towards the serial
   kernel. Flat is measured for Rsyrk and contradicted for Rtrsm Right on thanos between 2 and
   4 threads -- see the "THREAD-COUNT SCALING" paragraph above before changing either
   constant. */
#ifndef MPLAPACK_OMP_POTRF_WORK
#define MPLAPACK_OMP_POTRF_WORK(w, nt) ((nt) <= 24 ? (double)(w) : (double)(w) * (double)(nt) / 24.0)
#endif

/* Minimum number of iterations to share out, on whichever axis is being split: columns (n)
   for Rsyrk_omp, ROWS (m) for Rtrsm_omp's Right case. Work alone is not sufficient -- the
   last panel of a factorisation can have m as low as 1 and still carry real work. For Rsyrk
   this floor is the ONLY thing that can reject jb = 1 (measured 0.77-0.79x, flat in thread
   count, because num_threads is clamped to n and the "parallel" region runs on one thread):
   no work gate can, since k is unbounded. It is a floor on the parallel axis only, and it is
   NOT scaled by thread count, deliberately -- see the "n < 4*nthreads" note above. */
#ifndef MPLAPACK_OMP_MIN_POTRF_WIDTH
#define MPLAPACK_OMP_MIN_POTRF_WIDTH 2
#endif

#endif /* MPLAPACK_OMP_TUNING_H */
