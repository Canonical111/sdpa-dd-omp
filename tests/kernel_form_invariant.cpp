/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: new file, not part of the solver.
   Regression cell for the numerical invariants the Rpotrf lower-path threading rests on.
   See git log and the CI step "kernel accumulation-form invariant".

   Three things are pinned here. Each of them is silent when it breaks -- the solver still
   compiles, still converges, and quietly returns different bits or wrong answers -- which is
   why they get a test rather than a comment.

   (1) Rgemm_NT_omp's ACCUMULATION FORM.
       Rpotrf's lower path calls Rgemm("No transpose","Transpose",...) between the two panel
       kernels, so the factorisation's last bits depend on how that kernel accumulates. It
       presently uses
              C[i + j*ldc] += temp * A[i + l*lda];   temp = alpha * B[j + l*ldb];
       which is elementwise identical to Rsyrk's Lower/NoTranspose inner loop when B == A.
       Netlib's other gemm form -- accumulate TEMP from zero over l, then write
       C = alpha*TEMP + beta*C -- computes the same value in a different order, and in dd_real
       a different order is different bits. This fork has already edited that kernel once (A2,
       "restored netlib dgemm zero-skip"), so the risk is live. The cell asserts
       Rgemm("N","T") and Rsyrk("L","N") agree BIT for BIT on the lower triangle, and carries a
       positive control proving the comparison can actually see the TEMP-from-zero form.

   (2) Rsyrk_omp == Rsyrk, bit for bit, with a team up.
       Rsyrk_omp splits the output columns. If someone ever splits k -- the reduction axis --
       instead, the sum reorders and this fails.

   (3) Rtrsm_omp(Right/Lower/Transpose) == Rtrsm, bit for bit, with a team up.
       Rtrsm_omp splits Left over COLUMNS and Right over ROWS, because those are the independent
       axes of two different loop nests. Splitting Right over columns compiles and converges and
       is wrong; this is the check that catches it.

   Bit comparison is memcmp on the dd_real objects, not operator==, because operator== is a
   value comparison and the whole point is the representation. */

#include <mpblas_dd.h>

#include <cstdio>
#include <cstring>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

static int failures = 0;

/* Deterministic and self-contained: no rand(), whose sequence is not fixed across libcs. */
static double next_unit(unsigned long long &s) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((s >> 11) & 0x1FFFFFFFFFFFFFULL) / (double)0x20000000000000ULL * 2.0 - 1.0;
}

static bool same_bits(dd_real const &x, dd_real const &y) { return memcmp(&x, &y, sizeof(dd_real)) == 0; }

/* Returns the number of differing entries in the lower triangle of two n x n matrices. */
static long lower_diff(dd_real *p, dd_real *q, long n, long ld) {
    long bad = 0;
    for (long j = 0; j < n; j++)
        for (long i = j; i < n; i++)
            if (!same_bits(p[i + j * ld], q[i + j * ld]))
                bad++;
    return bad;
}

static long full_diff(dd_real *p, dd_real *q, long m, long n, long ld) {
    long bad = 0;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            if (!same_bits(p[i + j * ld], q[i + j * ld]))
                bad++;
    return bad;
}

static void check(const char *what, bool ok) {
    printf("%-58s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok)
        failures++;
}

/* The form we must NOT drift to: netlib's TEMP-from-zero accumulation. Used only as a positive
   control, to prove the bit comparison above is not vacuous. */
static void gemm_nt_tempform(long m, long n, long k, dd_real alpha, dd_real *A, long lda, dd_real *B, long ldb, dd_real beta, dd_real *C, long ldc) {
    for (long j = 0; j < n; j++) {
        for (long i = 0; i < m; i++) {
            dd_real temp = 0.0;
            for (long l = 0; l < k; l++)
                temp = temp + A[i + l * lda] * B[j + l * ldb];
            C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
        }
    }
}

/* Control for "never split the reduction axis k". A team given a slice of l each, accumulating
   into its own partial C and summing at the end, is the obvious way to parallelise Rsyrk over k
   and it is bit-wrong. Modelled deterministically here (two halves, then one add) so the
   control cannot itself be flaky. */
static void syrk_ln_ksplit(long n, long k, dd_real alpha, dd_real *A, long lda, dd_real beta, dd_real *C, long ldc) {
    long half = k / 2;
    dd_real *p0 = new dd_real[n * n];
    dd_real *p1 = new dd_real[n * n];
    for (long t = 0; t < n * n; t++) {
        p0[t] = 0.0;
        p1[t] = 0.0;
    }
    for (long pass = 0; pass < 2; pass++) {
        dd_real *P = pass ? p1 : p0;
        long l0 = pass ? half : 0;
        long l1 = pass ? k : half;
        for (long j = 0; j < n; j++)
            for (long l = l0; l < l1; l++)
                if (A[j + l * lda] != dd_real(0.0)) {
                    dd_real temp = alpha * A[j + l * lda];
                    for (long i = j; i < n; i++)
                        P[i + j * n] += temp * A[i + l * lda];
                }
    }
    for (long j = 0; j < n; j++)
        for (long i = j; i < n; i++)
            C[i + j * ldc] = beta * C[i + j * ldc] + (p0[i + j * n] + p1[i + j * n]);
    delete[] p0;
    delete[] p1;
}

/* Control for "Right side splits over ROWS, never over columns". Splitting the Right/Lower/
   Transpose branch over its outer loop -- which runs over the COLUMNS k -- is a data race, so
   it cannot be exhibited reproducibly by just adding a pragma. What is exhibited here instead
   is one legal execution of that racy loop: an omp parallel-for gives no ordering guarantee
   between threads, so the block of k assigned to the last thread may complete before the first
   thread starts. Running the blocks in reverse order is therefore a schedule the buggy code
   could actually produce -- and it produces a different answer, which is the point. (Running
   them in forward order would coincidentally reproduce the serial order and prove nothing;
   that near-miss is exactly why this control is written out rather than assumed.) */
static void trsm_right_colsplit_model(long nblk, long m, long n, dd_real alpha, dd_real *a, long lda, dd_real *b, long ldb) {
    for (long t = nblk - 1; t >= 0; t--) {
        long k0 = 1 + (t * n) / nblk;
        long k1 = ((t + 1) * n) / nblk;
        for (long k = k0; k <= k1; k++) {
            dd_real temp = dd_real(1.0) / a[(k - 1) + (k - 1) * lda];
            for (long i = 1; i <= m; i++)
                b[(i - 1) + (k - 1) * ldb] = temp * b[(i - 1) + (k - 1) * ldb];
            for (long j = k + 1; j <= n; j++)
                if (a[(j - 1) + (k - 1) * lda] != dd_real(0.0)) {
                    temp = a[(j - 1) + (k - 1) * lda];
                    for (long i = 1; i <= m; i++)
                        b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - temp * b[(i - 1) + (k - 1) * ldb];
                }
            if (alpha != dd_real(1.0))
                for (long i = 1; i <= m; i++)
                    b[(i - 1) + (k - 1) * ldb] = alpha * b[(i - 1) + (k - 1) * ldb];
        }
    }
}

int main() {
    const long n = 80;  /* order of the symmetric result   */
    const long k = 60;  /* reduction length                */
    const long mr = 200; /* rows of B for the Right-side trsm */

    unsigned long long seed = 0x9E3779B97F4A7C15ULL;

    dd_real *A = new dd_real[n * k];
    dd_real *C0 = new dd_real[n * n];
    dd_real *C1 = new dd_real[n * n];
    dd_real *C2 = new dd_real[n * n];
    dd_real *C3 = new dd_real[n * n];

    /* Every seventh entry of A is an exact zero, so the zero-skip branches of BOTH kernels are
       exercised. Rsyrk skips on a[j + l*lda] and Rgemm_NT skips on B[j + l*ldb] == a[j + l*lda];
       that they test the same element is itself part of the invariant. */
    for (long t = 0; t < n * k; t++)
        A[t] = (t % 7 == 3) ? dd_real(0.0) : dd_real(next_unit(seed));
    for (long t = 0; t < n * n; t++) {
        C0[t] = dd_real(next_unit(seed));
        C1[t] = C0[t];
        C2[t] = C0[t];
        C3[t] = C0[t];
    }

    const dd_real alpha = -1.0;
    const dd_real beta = 1.0;

#ifdef _OPENMP
    int nth = omp_get_max_threads();
    printf("openmp: yes, omp_get_max_threads() = %d\n", nth);
    if (nth < 2) {
        printf("FAIL: this cell must run with a team up (set OMP_NUM_THREADS >= 2), else the\n"
               "      threading invariants (2) and (3) are vacuous.\n");
        return 1;
    }
#else
    printf("openmp: no -- invariants (2) and (3) degenerate to serial-vs-serial and are skipped\n");
#endif

    /* ---- (1) Rsyrk("L","N") vs Rgemm("N","T") on the lower triangle ---- */
    Rsyrk("Lower", "No transpose", n, k, alpha, A, n, beta, C1, n);
    Rgemm("No transpose", "Transpose", n, n, k, alpha, A, n, A, n, beta, C2, n);
    check("(1) Rgemm(N,T) lower triangle == Rsyrk(L,N), bit for bit", lower_diff(C1, C2, n, n) == 0);

    /* positive control: the form we must not drift to has to be visibly different, or the
       comparison above proves nothing. */
    gemm_nt_tempform(n, n, k, alpha, A, n, A, n, beta, C3, n);
    check("(1c) control: TEMP-from-zero form IS detected as different", lower_diff(C1, C3, n, n) != 0);

    /* ---- (2) Rsyrk_omp vs Rsyrk ---- */
    for (long t = 0; t < n * n; t++)
        C2[t] = C0[t];
    Rsyrk_omp("Lower", "No transpose", n, k, alpha, A, n, beta, C2, n);
    check("(2) Rsyrk_omp == Rsyrk, bit for bit", lower_diff(C1, C2, n, n) == 0);

    for (long t = 0; t < n * n; t++)
        C3[t] = C0[t];
    syrk_ln_ksplit(n, k, alpha, A, n, beta, C3, n);
    check("(2c) control: splitting the k reduction IS detected", lower_diff(C1, C3, n, n) != 0);

    /* ---- (3) Rtrsm_omp Right/Lower/Transpose vs Rtrsm ---- */
    {
        const long nb = 64; /* columns of B, = Rpotrf's jb */
        dd_real *L = new dd_real[nb * nb];
        dd_real *rhs = new dd_real[mr * nb];
        dd_real *B1 = new dd_real[mr * nb];
        dd_real *B2 = new dd_real[mr * nb];
        dd_real *B3 = new dd_real[mr * nb];
        for (long j = 0; j < nb; j++)
            for (long i = 0; i < nb; i++)
                L[i + j * nb] = (i == j) ? dd_real(4.0 + next_unit(seed)) : (i > j ? dd_real(next_unit(seed)) : dd_real(0.0));
        for (long t = 0; t < mr * nb; t++) {
            rhs[t] = dd_real(next_unit(seed));
            B1[t] = rhs[t];
            B2[t] = rhs[t];
            B3[t] = rhs[t];
        }
        Rtrsm("Right", "Lower", "Transpose", "Non-unit", mr, nb, dd_real(1.0), L, nb, B1, mr);
        Rtrsm_omp("Right", "Lower", "Transpose", "Non-unit", mr, nb, dd_real(1.0), L, nb, B2, mr);
        check("(3) Rtrsm_omp(Right,Lower,T) == Rtrsm, bit for bit", full_diff(B1, B2, mr, nb, mr) == 0);

        trsm_right_colsplit_model(4, mr, nb, dd_real(1.0), L, nb, B3, mr);
        check("(3c) control: splitting Right over COLUMNS IS detected", full_diff(B1, B3, mr, nb, mr) != 0);

        delete[] L;
        delete[] rhs;
        delete[] B1;
        delete[] B2;
        delete[] B3;
    }

    delete[] A;
    delete[] C0;
    delete[] C1;
    delete[] C2;
    delete[] C3;

    if (failures) {
        printf("\n%d invariant(s) broken.\n", failures);
        return 1;
    }
    printf("\nall kernel form invariants hold.\n");
    return 0;
}
