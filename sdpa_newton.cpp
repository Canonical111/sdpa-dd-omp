/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

/* MODIFIED from upstream (GPLv2 2a notice), 2026-07-31: Schur-complement (bMat) construction threaded. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-05: the dense Schur complement bMat is built in its LOWER TRIANGLE ONLY; the strict upper half was accumulated every iteration and never read. See git log. */
#include <sdpa_newton.h>
#include <sdpa_parts.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

// review2 dimension edge 2: m and SDP_nBlock are each bounded by the reader,
// but their PRODUCT was formed in signed int at the allocation sites below.
// Bounding the factors by file size does not prove the product fits.
// Outside any _OPENMP guard: the call sites are unconditional, and a serial
// or flag-overridden build (CI's sanitizer/warnings jobs) needs this too.
static int checkedProductInt(int a, int b, const char *what) {
    const long long p = static_cast<long long>(a) * static_cast<long long>(b);
    if (a < 0 || b < 0 || p > INT_MAX) {
        std::cerr << "allocation size overflow: " << what << " = " << a << " * " << b << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return static_cast<int>(p);
}

// Strict parse of the census knob. Same contract as every other SDPA_* knob in this fork: an
// unrecognised value is REFUSED, never treated as the default. A typo that silently disabled a
// census would read as "the census found nothing", which is the worst possible failure mode for
// a measurement whose whole purpose is to establish reachability.
/* TEST-ONLY negative control for the hasF2Gcal fix, behind the spchol test gate.

   Reproduces the WORST-CASE indeterminate read the group-scope fix eliminated: forces the flag
   true on every pair that is not its group's first, which is exactly what the old per-pair
   declaration could deliver. calF2 then skips computing G = X*F and contracts Aj against
   whatever work2 still holds -- so if this changes the answer, the branch is consequential and
   a fixture that reaches it is non-vacuous.

   This hook exists because the OBVIOUS test does not work. `-ftrivial-auto-var-init=pattern`
   fills an unsigned char with 0xFE and an int with 0xFEFEFEFE, but a `bool` with ZERO -- a bool
   has no invalid representation to poison with. So zero-init and pattern-init both yield
   `false` here, and a differential build of the two cannot distinguish the readings at all. An
   earlier round of this port ran exactly that comparison across seven problems and reported "no
   dependence"; the comparison could not have failed. See git log. */
static bool bmat_test_f2_stale() {
    const char *e = getenv("SDPA_BMAT_TEST_F2_STALE");
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
#ifndef SDPA_SPCHOL_TEST_HOOKS
    rError("SDPA_BMAT_TEST_F2_STALE is a test hook and this binary was not built with"
           " -DSDPA_SPCHOL_TEST_HOOKS, so the hook does not exist here");
    return false;
#else
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_TEST_F2_STALE must be 0 or 1 (got \"" << e << "\")");
    return false;
#endif
}

/* TEST-ONLY negative control for the assembly oracle, behind the spchol test gate.

   Perturbs ONE assembled element by one ulp of its low limb. Without it, "every thread count
   produced the same stream" is equally consistent with a stream that cannot tell anything
   apart -- and this fork has now shipped two comparisons that could not fail. The low limb is
   chosen because it is invisible to every printed field and visible only to a comparison over
   the actual bits, which is the property being claimed. */
// unset/auto -> gated; serial -> never thread; parallel -> thread regardless of the gate.
enum BmatAsmMode { BMAT_ASM_AUTO, BMAT_ASM_SERIAL, BMAT_ASM_PARALLEL };
static BmatAsmMode bmat_asm_mode() {
    const char *e = getenv("SDPA_BMAT_ASM_MODE");
    if (e == NULL || e[0] == '\0' || strcmp(e, "auto") == 0) {
        return BMAT_ASM_AUTO;
    }
    if (strcmp(e, "serial") == 0) {
        return BMAT_ASM_SERIAL;
    }
    if (strcmp(e, "parallel") == 0) {
#ifndef _OPENMP
        // Force or fail, never silently downgrade -- the lesson from SDPA_SPCHOL_MODE=force.
        rError("SDPA_BMAT_ASM_MODE=parallel was requested but this binary was built without"
               " OpenMP, so there is no parallel path to take");
#endif
        return BMAT_ASM_PARALLEL;
    }
    rError("SDPA_BMAT_ASM_MODE must be auto, serial or parallel (got \"" << e << "\")");
    return BMAT_ASM_AUTO;
}

static bool bmat_asm_log() {
    const char *e = getenv("SDPA_BMAT_ASM_LOG");
    return e != NULL && e[0] != '\0' && strcmp(e, "0") != 0;
}

// Bounded unsigned knob, same strict contract as spchol_gate: empty refused, signs refused,
// overflow refused. A silently-wrapped negative here would read as an enormous budget.
static unsigned long long bmat_asm_u64(const char *name, unsigned long long dflt) {
    const char *e = getenv(name);
    if (e == NULL) {
        return dflt;
    }
    if (e[0] == '\0') {
        rError(name << " is set but empty; unset it to use the default");
    }
    const char *q = e;
    while (*q == ' ' || *q == '\t') {
        q++;
    }
    if (*q == '-' || *q == '+') {
        rError(name << " must be a non-negative integer without a sign (got \"" << e << "\")");
    }
    errno = 0;
    char *end = NULL;
    const unsigned long long v = strtoull(q, &end, 10);
    if (end == q || *end != '\0' || errno == ERANGE) {
        rError(name << " must be a non-negative integer (got \"" << e << "\")");
    }
    return v;
}

static bool bmat_asm_mutate() {
    const char *e = getenv("SDPA_BMAT_ASM_MUTATE");
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
#ifndef SDPA_SPCHOL_TEST_HOOKS
    rError("SDPA_BMAT_ASM_MUTATE is a test hook and this binary was not built with"
           " -DSDPA_SPCHOL_TEST_HOOKS, so the hook does not exist here");
    return false;
#else
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_ASM_MUTATE must be 0 or 1 (got \"" << e << "\")");
    return false;
#endif
}

static bool bmat_asm_digest_wanted() {
    const char *e = getenv("SDPA_BMAT_ASM_DIGEST");
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_ASM_DIGEST must be 0 or 1 (got \"" << e << "\")");
    return false;
}

// Exact-comparison sink, APPENDED to, so one file holds every assembly of a solve in order and
// two runs compare with cmp(1) -- byte identity of the whole stream, not equality of a hash.
// Openability is validated at configuration time (Make_bMat), not here, so that a dense-route
// problem cannot accept an unwritable path in silence.
static const char *bmat_asm_dump_path() {
    const char *e = getenv("SDPA_BMAT_ASM_DUMP");
    if (e == NULL || e[0] == '\0') {
        return NULL;
    }
    return e;
}

/* CALIBRATED on dd, 2026-08-23, on TWO architectures -- not inherited from gmp, whose 4000 was
   never measured here. Threading the assembly is forced either way and the whole problem set is
   swept; the crossover is what the gate encodes.

     problem     pairs        i9-13900K/24thr   M1 Max/10thr
     example1          6         0.003x            -
     control1        351         0.427x            -
     8_min           833         0.813x            -
     gpp100        5,151         0.133x          0.292x     <- largest problem that LOSES
     theta1        5,460         1.966x            -        <- marginal win
     arch0        15,225         9.650x          7.887x     <- smallest CLEAR win
     10_min       16,219         9.315x          9.698x
     truss5      124,368         7.208x            -
     12_min      397,430         6.440x            -
     dE3       8,113,226        10.398x         10.120x
     dE4      11,287,191         9.672x            -

   The two architectures agree on the boundary despite differing in core count (24 vs 10), core
   design (8P+16E heterogeneous vs homogeneous), ISA (x86-64 vs aarch64) and memory system --
   which is the point of requiring two.

   ADOPTED: 8000, between gpp100's 5,151 (loses on both) and arch0's 15,225 (wins ~8-10x on
   both). Deliberately NOT at the crossover: theta1 wins only 1.97x at 5,460, so the region
   just above gpp100 is where the gain is small and the risk of a different machine tipping
   negative is highest. 8000 keeps every clear win and refuses every measured loss.

   Note the gate is a proxy: pairs alone does not determine the work, since a group's cost also
   scales with its block dimension (8_min loses at 833 pairs while theta1 wins at 5,460). A
   work-based gate would be sharper. This one is what gmp uses, it is measured rather than
   inherited, and it is overridable at runtime -- which is how a machine that disagrees is
   accommodated without a rebuild. */
#ifndef SDPA_BMAT_ASM_MIN_PAIRS_DEFAULT
#define SDPA_BMAT_ASM_MIN_PAIRS_DEFAULT 8000ULL
#endif
#ifndef SDPA_BMAT_ASM_SCRATCH_MB_DEFAULT
#define SDPA_BMAT_ASM_SCRATCH_MB_DEFAULT 4096ULL
#endif

static bool bmat_asm_profile_wanted() {
    const char *e = getenv("SDPA_BMAT_ASM_PROFILE");
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_ASM_PROFILE must be 0 or 1 (got \"" << e << "\")");
    return false;
}

static bool bmat_asm_census_wanted() {
    const char *e = getenv("SDPA_BMAT_ASM_CENSUS");
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_ASM_CENSUS must be 0 or 1 (got \"" << e << "\")");
    return false;
}


// ---------------------------------------------------------------------------
// Thresholds for threading the Schur-complement (bMat) construction. Override with -D.
// The k1 x k2 loop costs roughly nConstraint^2 * blockDim; below this there is not
// enough work to amortise an OpenMP fork/join.
// ---------------------------------------------------------------------------
#ifndef SDPA_OMP_MIN_CONSTRAINTS
#define SDPA_OMP_MIN_CONSTRAINTS 8
#endif
#ifndef SDPA_OMP_MIN_BMAT_WORK
#define SDPA_OMP_MIN_BMAT_WORK 20000.0
#endif
// Hard ceiling on the extra memory used to privatise work1/work2 across threads.
// Cost is 2 * blockDim^2 * bytes-per-element per extra thread; on a large block that would
// otherwise grow without bound. If the full thread count would exceed this, the thread
// count for the block is reduced rather than the memory.
#ifndef SDPA_OMP_MAX_PRIV_MB
#define SDPA_OMP_MAX_PRIV_MB 256.0
#endif
// Bytes actually occupied by one scalar. For dd_real/qd_real the mantissa is stored inline,
// so sizeof() is exact. For mpf_class it is NOT: the object is a 24-byte descriptor whose
// limbs are allocated separately, so sizeof() undercounts by ~3x at 256-bit precision and
// the memory cap above would admit several times its nominal budget.
static inline double sdpa_omp_bytes_per_elem() {
    return (double)sizeof(dd_real);
}
// Choosing the parallel axis. For an F1/F2-dominated block the per-constraint setup is a
// blockDim^3 dense gemm, which Rgemm already threads well on its own; threading k1 instead
// makes each of those gemms serial and gains nothing. Below this gemm size Rgemm cannot
// parallelise effectively (blocks of 25-80 in the control* family) and threading k1 wins
// several-fold. Measured crossover on an i9-13900K lies between 80^3 and 100^3.
#ifndef SDPA_OMP_RGEMM_OWNS_BLOCK
#define SDPA_OMP_RGEMM_OWNS_BLOCK 700000.0
#endif

namespace sdpa {

Newton::Newton() {
    useFormula = NULL;

    bMat_type = DENSE;

    // Caution: if SDPA doesn't use sparse bMat,
    //          following variables are indefinite.
    this->SDP_nBlock = -1;
    SDP_blocks_disjoint = true;
    SDP_number = NULL;
    SDP_location_sparse_bMat = NULL;
    SDP_constraint1 = NULL;
    SDP_constraint2 = NULL;
    SDP_blockIndex1 = NULL;
    SDP_blockIndex2 = NULL;
    this->SOCP_nBlock = -1;
    SOCP_number = NULL;
    SOCP_location_sparse_bMat = NULL;
    SOCP_constraint1 = NULL;
    SOCP_constraint2 = NULL;
    SOCP_blockIndex1 = NULL;
    SOCP_blockIndex2 = NULL;
    this->LP_nBlock = -1;
    LP_number = NULL;
    LP_location_sparse_bMat = NULL;
    LP_constraint1 = NULL;
    LP_constraint2 = NULL;
    LP_blockIndex1 = NULL;
    LP_blockIndex2 = NULL;

    ordering = NULL;
    reverse_ordering = NULL;
    diagonalIndex = NULL;
}

Newton::Newton(int m, int SDP_nBlock, int *SDP_blockStruct, int SOCP_nBlock, int *SOCP_blockStruct, int LP_nBlock) { initialize(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock); }

Newton::~Newton() { terminate(); }

void Newton::initialize(int m, int SDP_nBlock, int *SDP_blockStruct, int SOCP_nBlock, int *SOCP_blockStruct, int LP_nBlock) {
    gVec.initialize(m);

    DxMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    DyVec.initialize(m);
    DzMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    r_zinvMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    x_rd_zinvMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

    rNewCheck();
    useFormula = new FormulaType[checkedProductInt(m, SDP_nBlock, "useFormula")];
    if (useFormula == NULL) {
        rError("Newton:: memory exhausted ");
    }

    bMat_type = DENSE;

    // Caution: if SDPA doesn't use sparse bMat,
    //          following variables are indefinite.
    this->SDP_nBlock = -1;
    SDP_blocks_disjoint = true;
    SDP_number = NULL;
    SDP_location_sparse_bMat = NULL;
    SDP_constraint1 = NULL;
    SDP_constraint2 = NULL;
    SDP_blockIndex1 = NULL;
    SDP_blockIndex2 = NULL;
    this->SOCP_nBlock = -1;
    SOCP_number = NULL;
    SOCP_location_sparse_bMat = NULL;
    SOCP_constraint1 = NULL;
    SOCP_constraint2 = NULL;
    SOCP_blockIndex1 = NULL;
    SOCP_blockIndex2 = NULL;
    this->LP_nBlock = -1;
    LP_number = NULL;
    LP_location_sparse_bMat = NULL;
    LP_constraint1 = NULL;
    LP_constraint2 = NULL;
    LP_blockIndex1 = NULL;
    LP_blockIndex2 = NULL;

    ordering = NULL;
    reverse_ordering = NULL;
    diagonalIndex = NULL;
}

void Newton::terminate() {

    if (bMat_type == SPARSE) {

        if (SDP_location_sparse_bMat && SDP_constraint1 && SDP_constraint2 && SDP_blockIndex1 && SDP_blockIndex2) {
            for (int k = 0; k < SDP_nBlock; ++k) {
                delete[] SDP_location_sparse_bMat[k];
                delete[] SDP_constraint1[k];
                delete[] SDP_constraint2[k];
                delete[] SDP_blockIndex1[k];
                delete[] SDP_blockIndex2[k];
                SDP_location_sparse_bMat[k] = NULL;
                SDP_constraint1[k] = NULL;
                SDP_constraint2[k] = NULL;
                SDP_blockIndex1[k] = NULL;
                SDP_blockIndex2[k] = NULL;
            }
            delete[] SDP_number;
            delete[] SDP_location_sparse_bMat;
            delete[] SDP_constraint1;
            delete[] SDP_constraint2;
            delete[] SDP_blockIndex1;
            delete[] SDP_blockIndex2;
            SDP_number = NULL;
            SDP_location_sparse_bMat = NULL;
            SDP_constraint1 = NULL;
            SDP_constraint2 = NULL;
            SDP_blockIndex1 = NULL;
            SDP_blockIndex2 = NULL;
        }
#if 0
    if (SOCP_location_sparse_bMat && SOCP_constraint1 && SOCP_constraint2
	&& SOCP_blockIndex1 && SOCP_blockIndex2) {
      for (int k=0; k<SOCP_nBlock; ++k) {
	delete[] SOCP_location_sparse_bMat[k];
	delete[] SOCP_constraint1[k];    delete[] SOCP_constraint2[k];
	delete[] SOCP_blockIndex1[k];    delete[] SOCP_blockIndex2[k];
	SOCP_location_sparse_bMat[k] = NULL;
	SOCP_constraint1[k] = NULL;   SOCP_constraint2[k] = NULL;
	SOCP_blockIndex1[k] = NULL;   SOCP_blockIndex2[k] = NULL;
      }
      delete[] SOCP_number;  delete[] SOCP_location_sparse_bMat;
      delete[] SOCP_constraint1;  delete[] SOCP_constraint2;
      delete[] SOCP_blockIndex1;  delete[] SOCP_blockIndex2;
      SOCP_number = NULL;  SOCP_location_sparse_bMat = NULL;
      SOCP_constraint1 = NULL;  SOCP_constraint2 = NULL;
      SOCP_blockIndex1 = NULL;  SOCP_blockIndex2 =NULL;
    }
#endif
        if (LP_location_sparse_bMat && LP_constraint1 && LP_constraint2 && LP_blockIndex1 && LP_blockIndex2) {
            for (int k = 0; k < LP_nBlock; ++k) {
                delete[] LP_location_sparse_bMat[k];
                delete[] LP_constraint1[k];
                delete[] LP_constraint2[k];
                delete[] LP_blockIndex1[k];
                delete[] LP_blockIndex2[k];
                LP_location_sparse_bMat[k] = NULL;
                LP_constraint1[k] = NULL;
                LP_constraint2[k] = NULL;
                LP_blockIndex1[k] = NULL;
                LP_blockIndex2[k] = NULL;
            }
            delete[] LP_number;
            delete[] LP_location_sparse_bMat;
            delete[] LP_constraint1;
            delete[] LP_constraint2;
            delete[] LP_blockIndex1;
            delete[] LP_blockIndex2;
            LP_number = NULL;
            LP_location_sparse_bMat = NULL;
            LP_constraint1 = NULL;
            LP_constraint2 = NULL;
            LP_blockIndex1 = NULL;
            LP_blockIndex2 = NULL;
        }

        if (ordering) {
            delete[] ordering;
            ordering = NULL;
        }
        if (reverse_ordering) {
            delete[] reverse_ordering;
            reverse_ordering = NULL;
        }
        if (diagonalIndex) {
            delete[] diagonalIndex;
            diagonalIndex = NULL;
        }
        sparse_bMat.terminate();

    } else { // bMat_type == DENSE
        bMat.terminate();
    }

    gVec.terminate();
    DxMat.terminate();
    DyVec.terminate();
    DzMat.terminate();
    r_zinvMat.terminate();
    x_rd_zinvMat.terminate();

    if (useFormula != NULL) {
        delete[] useFormula;
    }
    useFormula = NULL;
}

void Newton::initialize_dense_bMat(int m) {
    //  bMat_type = DENSE;
    //  printf("DENSE computations\n");
    bMat.initialize(m, m, DenseMatrix::DENSE);
}

// 2008/03/12 kazuhide nakata
void Newton::initialize_sparse_bMat(int m, IV *newToOldIV, IVL *symbfacIVL) {

    //  bMat_type = SPARSE;
    //  printf("SPARSE computation\n");

    int i, j, k;
    int *newToOld;

    newToOld = IV_entries(newToOldIV);

    rNewCheck();
    ordering = new int[m];
    if (ordering == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }
    for (i = 0; i < m; i++) {
        ordering[i] = newToOld[i];
    }

    rNewCheck();
    reverse_ordering = new int[m];
    if (reverse_ordering == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }
    for (i = 0; i < m; i++) {
        reverse_ordering[ordering[i]] = i;
    }

    // separate front or back node
    int *counter;
    int nClique = IVL_nlist(symbfacIVL);
    int psize;
    int *pivec;
    bool *bnode;
    int *nFront;

    rNewCheck();
    counter = new int[m];
    bnode = new bool[m];
    nFront = new int[nClique];

    if ((counter == NULL) || (bnode == NULL) || (nFront == NULL)) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }

    for (i = 0; i < m; i++) {
        bnode[i] = false;
        counter[i] = -1;
    }

    // search number of front
    for (int l = nClique - 1; l >= 0; l--) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (i = 0; i < psize; i++) {
            int ii = reverse_ordering[pivec[i]];
            if (bnode[ii] == false) {
                counter[ii] = psize - i;
                bnode[ii] = true;
            } else {
                nFront[l] = i;
                break;
            }
        }
        if (i == psize) {
            nFront[l] = psize;
        }
    }

    // error check
    for (i = 0; i < m; i++) {
        if (counter[i] == -1) {
            rError("Newton::initialize_sparse_bMat: program bug");
        }
    }

    // make index of diagonal
    rNewCheck();
    diagonalIndex = new int[m + 1];
    if (diagonalIndex == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }

    // Checked prefix sum. diagonalIndex[m] is the sparse bMat's stored-element count and is
    // accumulated in int; the per-row counters are each bounded, but their SUM is not. Formed
    // in 64-bit and refused before narrowing, naming the row where it overflowed.
    diagonalIndex[0] = 0;
    {
        long long running = 0;
        for (i = 1; i < m + 1; i++) {
            running += (long long)counter[i - 1];
            if (!sdpaFitsInt(running)) {
                rError("Newton::initialize_sparse_bMat: the sparse bMat needs " << running
                       << " stored elements by row " << i << ", past INT_MAX=" << INT_MAX
                       << ". The dense representation is also unavailable above m=46340, so"
                       << " this problem is too large for this build's index type.");
            }
            diagonalIndex[i] = (int)running;
        }
    }

    // initialize sparse_bMat
    sparse_bMat.initialize(m, m, SparseMatrix::SPARSE, diagonalIndex[m]);

    // initialize index of sparse_bmat
    int nonzeros = 0;
    for (int l = 0; l < nClique; l++) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (i = 0; i < nFront[l]; i++) {
            int ii = reverse_ordering[pivec[i]];
            for (j = i; j < psize; j++) {
                int jj = reverse_ordering[pivec[j]];
                int index = diagonalIndex[ii] + j - i;
                sparse_bMat.row_index[index] = ii;
                sparse_bMat.column_index[index] = jj;
                nonzeros++;
            }
        }
    }
    // error check
    if (nonzeros != sparse_bMat.NonZeroNumber) {
        rError("Newton::initialize_sparse_bMat  probram bug");
    }
    sparse_bMat.NonZeroCount = nonzeros;
    //  sparse_bMat.display();

    delete[] counter;
    delete[] bnode;
    delete[] nFront;
}

// 2008/03/12 kazuhide nakata
void Newton::initialize_bMat(int m, Chordal &chordal, InputData &inputData, FILE *fpOut) {
    /* Create clique tree */

    switch (chordal.best) {
    case -1: {
        bMat_type = DENSE;
        printf("DENSE computations\n");
        fprintf(fpOut, "DENSE computation\n");
        initialize_dense_bMat(m);
        break;
    }
    case 0: {
        rError("no support for METIS");
        break;
    }
    case 1: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_MMD, chordal.symbfacIVL_MMD);
        make_aggrigateIndex(inputData);
        break;
    }
    case 2: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_ND, chordal.symbfacIVL_ND);
        make_aggrigateIndex(inputData);
        break;
    }
    case 3: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_MS, chordal.symbfacIVL_MS);
        make_aggrigateIndex(inputData);
        break;
    }
    case 4: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_NDMS, chordal.symbfacIVL_NDMS);
        make_aggrigateIndex(inputData);
        break;
    }
    }
}

void Newton::make_aggrigateIndex_SDP(InputData &inputData) {
    int t, ii, jj;

    SDP_nBlock = inputData.SDP_nBlock;
    rNewCheck();
    SDP_number = new int[SDP_nBlock];
    if (SDP_number == NULL) {
        rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
    }

    // memory allocate for aggrigateIndex
    rNewCheck();
    SDP_constraint1 = new int *[SDP_nBlock];
    SDP_constraint2 = new int *[SDP_nBlock];
    SDP_blockIndex1 = new int *[SDP_nBlock];
    SDP_blockIndex2 = new int *[SDP_nBlock];
    SDP_location_sparse_bMat = new int *[SDP_nBlock];
    if ((SDP_constraint1 == NULL) || (SDP_constraint2 == NULL) || (SDP_blockIndex1 == NULL) || (SDP_blockIndex2 == NULL) || (SDP_location_sparse_bMat == NULL)) {
        rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
    }

    for (int l = 0; l < SDP_nBlock; l++) {
        // Checked: n(n+1)/2 overflows int at n >= 65536, and the reader's bounds do not
        // forbid that count for a block. Formed in 64-bit, refused before narrowing.
        const int n_sdp = inputData.SDP_nConstraint[l];
        const long long tmp64_sdp = sdpaTriangularCount(n_sdp);
        if (n_sdp < 0 || !sdpaFitsInt(tmp64_sdp)) {
            rError("Newton::make_aggrigateIndex_SDP: block " << l << " with "
                   << n_sdp << " constraints needs " << tmp64_sdp
                   << " pairs, which is not representable as an int");
        }
        int tmp = (int)tmp64_sdp;
        rNewCheck();
        SDP_number[l] = tmp;
        SDP_constraint1[l] = new int[tmp];
        SDP_constraint2[l] = new int[tmp];
        SDP_blockIndex1[l] = new int[tmp];
        SDP_blockIndex2[l] = new int[tmp];
        SDP_location_sparse_bMat[l] = new int[tmp];
        if ((SDP_constraint1[l] == NULL) || (SDP_constraint2[l] == NULL) || (SDP_blockIndex1[l] == NULL) || (SDP_blockIndex2[l] == NULL) || (SDP_location_sparse_bMat[l] == NULL)) {
            rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
        }
    }

    // ---- validate the reader-supplied source indices BEFORE they become array subscripts ----
    //
    // Everything this function stores is later used as an INDEX by the assembly:
    // inputData.A[constraint1[k]], A[i].SDP_sp_block[blockIndex1[k]] and
    // sparse_bMat.sp_ele[location[k]]. Since the sparse SDP assembly is threaded, an
    // out-of-range entry is an out-of-bounds WRITE from inside a parallel region -- which
    // surfaces as a corrupted factor a long way from its cause. The size guard above proves
    // the ALLOCATION is representable; it says nothing about the values, so this is the other
    // half of "the aggregate map is validated".
    //
    // The duplicate check is not decoration. The pair enumeration below assumes each
    // constraint appears at most once in a block's constraint list: if one repeated, the
    // i == j ties would be kept twice and NonZeroCount would no longer equal the triangular
    // count the array was sized to -- which is exactly the mismatch the post-check catches.
    {
        const int m_bmat = sparse_bMat.nRow;
        int *seen = new int[m_bmat > 0 ? m_bmat : 1];
        if (seen == NULL) {
            rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
        }
        for (int i = 0; i < m_bmat; i++) {
            seen[i] = -1;
        }
        for (int l = 0; l < SDP_nBlock; l++) {
            for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
                const int i = inputData.SDP_constraint[l][k1];
                if (i < 0 || i >= m_bmat) {
                    rError("Newton::make_aggrigateIndex_SDP: block " << l << " entry " << k1
                           << " names constraint " << i << ", outside [0," << m_bmat << ")");
                }
                const int ib = inputData.SDP_blockIndex[l][k1];
                if (ib < 0 || ib >= inputData.A[i].SDP_sp_nBlock) {
                    rError("Newton::make_aggrigateIndex_SDP: block " << l << " entry " << k1
                           << " names sub-block " << ib << " of constraint " << i
                           << ", which has " << inputData.A[i].SDP_sp_nBlock);
                }
                if (seen[i] == l) {
                    rError("Newton::make_aggrigateIndex_SDP: constraint " << i
                           << " appears more than once in block " << l
                           << "; the pair enumeration requires it to appear at most once");
                }
                seen[i] = l;
                const int ri = reverse_ordering[i];
                if (ri < 0 || ri >= m_bmat) {
                    rError("Newton::make_aggrigateIndex_SDP: reverse_ordering[" << i << "] = "
                           << ri << ", outside [0," << m_bmat << ")");
                }
            }
        }
        delete[] seen;
    }

    for (int l = 0; l < SDP_nBlock; l++) {
        int NonZeroCount = 0;

        for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
            int i = inputData.SDP_constraint[l][k1];
            int ib = inputData.SDP_blockIndex[l][k1];
            int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

            for (int k2 = 0; k2 < inputData.SDP_nConstraint[l]; k2++) {
                int j = inputData.SDP_constraint[l][k2];
                int jb = inputData.SDP_blockIndex[l][k2];
                int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

                if ((inz < jnz) || ((inz == jnz) && (i < j))) {
                    continue;
                }

                // Refuse to write past the triangular allocation rather than detect it after
                // the fact -- NonZeroCount is the subscript of the five stores below.
                if (NonZeroCount >= SDP_number[l]) {
                    rError("Newton::make_aggrigateIndex_SDP: block " << l
                           << " produced more than its " << SDP_number[l] << " allocated pairs");
                }
                // set index which A_i and A_j are not zero matrix
                SDP_constraint1[l][NonZeroCount] = i;
                SDP_constraint2[l][NonZeroCount] = j;
                SDP_blockIndex1[l][NonZeroCount] = ib;
                SDP_blockIndex2[l][NonZeroCount] = jb;
                if (reverse_ordering[i] < reverse_ordering[j]) {
                    ii = reverse_ordering[i];
                    jj = reverse_ordering[j];
                } else {
                    jj = reverse_ordering[i];
                    ii = reverse_ordering[j];
                }

                // binary search for index of sparse_bMat
                t = -1;
                int begin = diagonalIndex[ii];
                int end = diagonalIndex[ii + 1] - 1;
                int target = (begin + end) / 2;
                while (end - begin > 1) {
                    if (sparse_bMat.column_index[target] < jj) {
                        begin = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] > jj) {
                        end = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] == jj) {
                        t = target;
                        break;
                    }
                }
                if (t == -1) {
                    if (sparse_bMat.column_index[begin] == jj) {
                        t = begin;
                    } else if (sparse_bMat.column_index[end] == jj) {
                        t = end;
                    } else {
                        rError("Newton::make_aggrigateIndex_SDP  program bug");
                    }
                }

                if (t < 0 || t >= sparse_bMat.NonZeroCount) {
                    rError("Newton::make_aggrigateIndex_SDP: block " << l << " pair "
                           << NonZeroCount << " maps to sparse bMat element " << t
                           << ", outside [0," << sparse_bMat.NonZeroCount << ")");
                }
                SDP_location_sparse_bMat[l][NonZeroCount] = t;
                NonZeroCount++;
            }
        } // for k1

        // THE check, not a formality: both assembly routines loop to SDP_number[l], NOT to
        // NonZeroCount. A short fill leaves the tail of these five arrays as whatever new[]
        // returned, and the assembly then subscripts inputData.A[] with it.
        if (NonZeroCount != SDP_number[l]) {
            rError("Newton::make_aggrigateIndex_SDP: block " << l << " filled " << NonZeroCount
                   << " of its " << SDP_number[l] << " allocated pairs; the assembly reads all "
                   << SDP_number[l]);
        }
    }     // for k  kth block

    // ---- the property the threaded assembly's `nowait` depends on, MEASURED ----
    //
    // Two questions in one pass over the finished map, both about sparse_bMat slots:
    //
    //   within a block, are the slots pairwise distinct? The assembly writes `+=` into
    //     SDP_location_sparse_bMat[l][iter] from several threads at once, and it is only safe
    //     because distinct iters of one block hit distinct slots. That follows from the
    //     constraint indices being distinct, which the pre-pass checked -- but the derivation
    //     is one step long and the cost of confirming it is one comparison per pair, so it is
    //     confirmed rather than argued.
    //
    //   across blocks, is any slot written by more than one block? If not, two blocks may be in
    //     flight at once and the worksharing loop can carry `nowait`; if so, they may not.
    //     The released v7.1.3-omp.2 assumed the answer was always "no" and was wrong on any
    //     problem with a constraint pair nonzero in two blocks -- SDPLIB truss6 among them.
    //
    // TWO BITS PER STORED bMat ELEMENT, held only for this pass.
    //
    // The first version used one *int* per element, which cost dE4 about 24 MB -- a 6% rise in
    // peak RSS on a 382 MB run, visible in dd_final_headline_postfix.tsv and traced there to
    // exactly this array. Two bitsets are 16x smaller and answer the same two questions.
    //
    // The trick that removes the need for a per-block scratch list: `seen_blk` is cleared by
    // RE-WALKING the block's own locations, which is O(pairs) time and zero extra memory. That
    // second walk also sets `seen_any`, so the two passes cost one combined sweep, not two.
    SDP_blocks_disjoint = true;
    if (sparse_bMat.NonZeroCount > 0 && SDP_nBlock > 0) {
        const size_t nbytes = ((size_t)sparse_bMat.NonZeroCount + 7) / 8;
        std::vector<unsigned char> seen_any(nbytes, 0); // written by SOME earlier block
        std::vector<unsigned char> seen_blk(nbytes, 0); // written by THIS block; cleared per block
        // The whole map is walked even after sharing is found: stopping early would skip the
        // within-block duplicate check on every later block, and that check is the one that
        // guards a genuine memory error rather than a scheduling choice.
        for (int l = 0; l < SDP_nBlock; l++) {
            for (int iter = 0; iter < SDP_number[l]; ++iter) {
                const int t = SDP_location_sparse_bMat[l][iter];
                const size_t byte = (size_t)t >> 3;
                const unsigned char bit = (unsigned char)(1u << ((unsigned)t & 7u));
                if (seen_blk[byte] & bit) {
                    rError("Newton::make_aggrigateIndex_SDP: block " << l
                           << " writes sparse bMat element " << t
                           << " twice; the threaded assembly requires one writer per element"
                              " within a block");
                }
                seen_blk[byte] |= bit;
                if (seen_any[byte] & bit) {
                    // Shared with an earlier block. Nothing is wrong -- it is what an SDP with
                    // a constraint pair nonzero in two blocks looks like -- but those blocks
                    // may no longer overlap in time.
                    SDP_blocks_disjoint = false;
                }
            }
            // Publish this block into seen_any and clear seen_blk, in one walk over the same
            // locations. Clearing only what this block set is what keeps the cost O(pairs)
            // rather than O(NonZeroCount) per block.
            for (int iter = 0; iter < SDP_number[l]; ++iter) {
                const int t = SDP_location_sparse_bMat[l][iter];
                const size_t byte = (size_t)t >> 3;
                const unsigned char bit = (unsigned char)(1u << ((unsigned)t & 7u));
                seen_any[byte] |= bit;
                seen_blk[byte] = (unsigned char)(seen_blk[byte] & ~bit);
            }
        }
    }
}

// DISPOSITION, recorded rather than left ambiguous: this function is DEAD CODE. Its only call
// site is commented out in make_aggrigateIndex() below, and dd ships no SOCP reader, so nothing
// reaches it. It was size-hardened alongside its two live siblings (an unguarded n(n+1)/2 in
// dead code is still a landmine for whoever revives it), but its aggregate map is deliberately
// NOT validated: the checks in the SDP and LP versions are exercised on every sparse-route run
// and by the CI fixtures, whereas the same code here could never be executed and so could never
// be trusted. If SOCP is ever revived, port the two validation blocks with it -- that is a
// smaller and more honest debt than shipping an unexecutable validator.
void Newton::make_aggrigateIndex_SOCP(InputData &inputData) {
    int t, ii, jj;

    SOCP_nBlock = inputData.SOCP_nBlock;
    rNewCheck();
    SOCP_number = new int[SOCP_nBlock];
    if (SOCP_number == NULL) {
        rError("Newton::make_aggrigateIndex_SOCP memory exhausted ");
    }

    // memory allocate for aggrigateIndex
    rNewCheck();
    SOCP_constraint1 = new int *[SOCP_nBlock];
    SOCP_constraint2 = new int *[SOCP_nBlock];
    SOCP_blockIndex1 = new int *[SOCP_nBlock];
    SOCP_blockIndex2 = new int *[SOCP_nBlock];
    SOCP_location_sparse_bMat = new int *[SOCP_nBlock];
    if ((SOCP_constraint1 == NULL) || (SOCP_constraint2 == NULL) || (SOCP_blockIndex1 == NULL) || (SOCP_blockIndex2 == NULL) || (SOCP_location_sparse_bMat == NULL)) {
        rError("Newton::make_aggrigateIndex_SOCP memory exhausted ");
    }

    for (int l = 0; l < SOCP_nBlock; l++) {
        // Checked: n(n+1)/2 overflows int at n >= 65536, and the reader's bounds do not
        // forbid that count for a block. Formed in 64-bit, refused before narrowing.
        const int n_socp = inputData.SOCP_nConstraint[l];
        const long long tmp64_socp = sdpaTriangularCount(n_socp);
        if (n_socp < 0 || !sdpaFitsInt(tmp64_socp)) {
            rError("Newton::make_aggrigateIndex_SOCP: block " << l << " with "
                   << n_socp << " constraints needs " << tmp64_socp
                   << " pairs, which is not representable as an int");
        }
        int tmp = (int)tmp64_socp;
        rNewCheck();
        SOCP_number[l] = tmp;
        SOCP_constraint1[l] = new int[tmp];
        SOCP_constraint2[l] = new int[tmp];
        SOCP_blockIndex1[l] = new int[tmp];
        SOCP_blockIndex2[l] = new int[tmp];
        SOCP_location_sparse_bMat[l] = new int[tmp];
        if ((SOCP_constraint1[l] == NULL) || (SOCP_constraint2[l] == NULL) || (SOCP_blockIndex1[l] == NULL) || (SOCP_blockIndex2[l] == NULL) || (SOCP_location_sparse_bMat[l] == NULL)) {
            rError("Newton::make_aggrigateIndex_SOCP memory exhausted ");
        }
    }

    for (int l = 0; l < SOCP_nBlock; l++) {
        int NonZeroCount = 0;

        for (int k1 = 0; k1 < inputData.SOCP_nConstraint[l]; k1++) {
            int i = inputData.SOCP_constraint[l][k1];
            int ib = inputData.SOCP_blockIndex[l][k1];
            int inz = inputData.A[i].SOCP_sp_block[ib].NonZeroEffect;

            for (int k2 = 0; k2 < inputData.SOCP_nConstraint[l]; k2++) {
                int j = inputData.SOCP_constraint[l][k2];
                int jb = inputData.SOCP_blockIndex[l][k2];
                int jnz = inputData.A[j].SOCP_sp_block[jb].NonZeroEffect;

                if ((inz < jnz) || ((inz == jnz) && (i < j))) {
                    continue;
                }

                // set index which A_i and A_j are not zero matrix
                SOCP_constraint1[l][NonZeroCount] = i;
                SOCP_constraint2[l][NonZeroCount] = j;
                SOCP_blockIndex1[l][NonZeroCount] = ib;
                SOCP_blockIndex2[l][NonZeroCount] = jb;
                if (reverse_ordering[i] < reverse_ordering[j]) {
                    ii = reverse_ordering[i];
                    jj = reverse_ordering[j];
                } else {
                    jj = reverse_ordering[i];
                    ii = reverse_ordering[j];
                }

                // binary search for index of sparse_bMat
                t = -1;
                int begin = diagonalIndex[ii];
                int end = diagonalIndex[ii + 1] - 1;
                int target = (begin + end) / 2;
                while (end - begin > 1) {
                    if (sparse_bMat.column_index[target] < jj) {
                        begin = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] > jj) {
                        end = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] == jj) {
                        t = target;
                        break;
                    }
                }
                if (t == -1) {
                    if (sparse_bMat.column_index[begin] == jj) {
                        t = begin;
                    } else if (sparse_bMat.column_index[end] == jj) {
                        t = end;
                    } else {
                        rError("Newton::make_aggrigateIndex_SDP  program bug");
                    }
                }

                SOCP_location_sparse_bMat[l][NonZeroCount] = t;
                NonZeroCount++;
            }
        } // for k1
    }     // for k  kth block
}

void Newton::make_aggrigateIndex_LP(InputData &inputData) {
    int t, ii, jj;

    LP_nBlock = inputData.LP_nBlock;
    rNewCheck();
    LP_number = new int[LP_nBlock];
    if (LP_number == NULL) {
        rError("Newton::make_aggrigateIndex_LP memory exhausted ");
    }

    // memory allocate for aggrigateIndex
    rNewCheck();
    LP_constraint1 = new int *[LP_nBlock];
    LP_constraint2 = new int *[LP_nBlock];
    LP_blockIndex1 = new int *[LP_nBlock];
    LP_blockIndex2 = new int *[LP_nBlock];
    LP_location_sparse_bMat = new int *[LP_nBlock];
    if ((LP_constraint1 == NULL) || (LP_constraint2 == NULL) || (LP_blockIndex1 == NULL) || (LP_blockIndex2 == NULL) || (LP_location_sparse_bMat == NULL)) {
        rError("Newton::make_aggrigateIndex_LP memory exhausted ");
    }

    for (int l = 0; l < LP_nBlock; l++) {
        // Checked: n(n+1)/2 overflows int at n >= 65536, and the reader's bounds do not
        // forbid that count for a block. Formed in 64-bit, refused before narrowing.
        const int n_lp = inputData.LP_nConstraint[l];
        const long long tmp64_lp = sdpaTriangularCount(n_lp);
        if (n_lp < 0 || !sdpaFitsInt(tmp64_lp)) {
            rError("Newton::make_aggrigateIndex_LP: block " << l << " with "
                   << n_lp << " constraints needs " << tmp64_lp
                   << " pairs, which is not representable as an int");
        }
        int tmp = (int)tmp64_lp;
        rNewCheck();
        LP_number[l] = tmp;
        LP_constraint1[l] = new int[tmp];
        LP_constraint2[l] = new int[tmp];
        LP_blockIndex1[l] = new int[tmp];
        LP_blockIndex2[l] = new int[tmp];
        LP_location_sparse_bMat[l] = new int[tmp];
        if ((LP_constraint1[l] == NULL) || (LP_constraint2[l] == NULL) || (LP_blockIndex1[l] == NULL) || (LP_blockIndex2[l] == NULL) || (LP_location_sparse_bMat[l] == NULL)) {
            rError("Newton::make_aggrigateIndex_LP memory exhausted ");
        }
    }

    // Same validation as the SDP map above, and for the same reason: compute_bMat_sparse_LP
    // uses these entries as subscripts into inputData.A[] and sparse_bMat.sp_ele[], and loops
    // to LP_number[l] rather than to NonZeroCount. See the comment there.
    {
        const int m_bmat = sparse_bMat.nRow;
        int *seen = new int[m_bmat > 0 ? m_bmat : 1];
        if (seen == NULL) {
            rError("Newton::make_aggrigateIndex_LP memory exhausted ");
        }
        for (int i = 0; i < m_bmat; i++) {
            seen[i] = -1;
        }
        for (int l = 0; l < LP_nBlock; l++) {
            for (int k1 = 0; k1 < inputData.LP_nConstraint[l]; k1++) {
                const int i = inputData.LP_constraint[l][k1];
                if (i < 0 || i >= m_bmat) {
                    rError("Newton::make_aggrigateIndex_LP: block " << l << " entry " << k1
                           << " names constraint " << i << ", outside [0," << m_bmat << ")");
                }
                const int ib = inputData.LP_blockIndex[l][k1];
                if (ib < 0 || ib >= inputData.A[i].LP_sp_nBlock) {
                    rError("Newton::make_aggrigateIndex_LP: block " << l << " entry " << k1
                           << " names sub-block " << ib << " of constraint " << i
                           << ", which has " << inputData.A[i].LP_sp_nBlock);
                }
                if (seen[i] == l) {
                    rError("Newton::make_aggrigateIndex_LP: constraint " << i
                           << " appears more than once in block " << l
                           << "; the pair enumeration requires it to appear at most once");
                }
                seen[i] = l;
                const int ri = reverse_ordering[i];
                if (ri < 0 || ri >= m_bmat) {
                    rError("Newton::make_aggrigateIndex_LP: reverse_ordering[" << i << "] = "
                           << ri << ", outside [0," << m_bmat << ")");
                }
            }
        }
        delete[] seen;
    }

    for (int l = 0; l < LP_nBlock; l++) {
        int NonZeroCount = 0;

        for (int k1 = 0; k1 < inputData.LP_nConstraint[l]; k1++) {
            int i = inputData.LP_constraint[l][k1];
            int ib = inputData.LP_blockIndex[l][k1];

            for (int k2 = 0; k2 < inputData.LP_nConstraint[l]; k2++) {
                int j = inputData.LP_constraint[l][k2];
                int jb = inputData.LP_blockIndex[l][k2];

                if (i < j) {
                    continue;
                }

                if (NonZeroCount >= LP_number[l]) {
                    rError("Newton::make_aggrigateIndex_LP: block " << l
                           << " produced more than its " << LP_number[l] << " allocated pairs");
                }
                // set index which A_i and A_j are not zero matrix
                LP_constraint1[l][NonZeroCount] = i;
                LP_constraint2[l][NonZeroCount] = j;
                LP_blockIndex1[l][NonZeroCount] = ib;
                LP_blockIndex2[l][NonZeroCount] = jb;
                if (reverse_ordering[i] < reverse_ordering[j]) {
                    ii = reverse_ordering[i];
                    jj = reverse_ordering[j];
                } else {
                    jj = reverse_ordering[i];
                    ii = reverse_ordering[j];
                }

                // binary search for index of sparse_bMat
                t = -1;
                int begin = diagonalIndex[ii];
                int end = diagonalIndex[ii + 1] - 1;
                int target = (begin + end) / 2;
                while (end - begin > 1) {
                    if (sparse_bMat.column_index[target] < jj) {
                        begin = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] > jj) {
                        end = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] == jj) {
                        t = target;
                        break;
                    }
                }
                if (t == -1) {
                    if (sparse_bMat.column_index[begin] == jj) {
                        t = begin;
                    } else if (sparse_bMat.column_index[end] == jj) {
                        t = end;
                    } else {
                        rError("Newton::make_aggrigateIndex_SDP  program bug");
                    }
                }

                if (t < 0 || t >= sparse_bMat.NonZeroCount) {
                    rError("Newton::make_aggrigateIndex_LP: block " << l << " pair "
                           << NonZeroCount << " maps to sparse bMat element " << t
                           << ", outside [0," << sparse_bMat.NonZeroCount << ")");
                }
                LP_location_sparse_bMat[l][NonZeroCount] = t;
                NonZeroCount++;
            }
        } // for k1

        if (NonZeroCount != LP_number[l]) {
            rError("Newton::make_aggrigateIndex_LP: block " << l << " filled " << NonZeroCount
                   << " of its " << LP_number[l] << " allocated pairs; the assembly reads all "
                   << LP_number[l]);
        }
    }     // for k  kth block
}

void Newton::make_aggrigateIndex(InputData &inputData) {
    make_aggrigateIndex_SDP(inputData);
    //  make_aggrigateIndex_SOCP(inputData);
    make_aggrigateIndex_LP(inputData);
}

void Newton::computeFormula_SDP(InputData &inputData, dd_real DenseRatio, dd_real Kappa) {
    int m = inputData.b.nDim;
    int SDP_nBlock = inputData.SDP_nBlock;

    int *upNonZeroCount;
    rNewCheck();
    upNonZeroCount = new int[checkedProductInt(m, SDP_nBlock, "upNonZeroCount")];
    if (upNonZeroCount == NULL) {
        rError("Newton:: memory exhausted ");
    }

    // We have no chance to use DenseRatio
    if (upNonZeroCount == NULL || useFormula == NULL) {
        rError("Newton:: failed initialization");
    }

    SparseLinearSpace *A = inputData.A;

#if 0
  for (int k=0; k<m; ++k) {
    for (int l=0; l<inputData.A[0].nBlock; ++l) {
      rMessage("A[" << k << "].ele[" << l << "] ="
	       << inputData.A[k].ele[l].NonZeroEffect);
    }
  }
#endif

    // Count sum of number of elements
    // that each number of elements are less than own.

    for (int iter = 0; iter < m * SDP_nBlock; iter++) {
        upNonZeroCount[iter] = 0;
    }

    for (int l = 0; l < SDP_nBlock; ++l) {
        for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
            int i = inputData.SDP_constraint[l][k1];
            int ib = inputData.SDP_blockIndex[l][k1];
            int inz = A[i].SDP_sp_block[ib].NonZeroEffect;
            int up = inz;
            // rMessage("up = " << up);

            for (int k2 = 0; k2 < inputData.SDP_nConstraint[l]; k2++) {
                int j = inputData.SDP_constraint[l][k2];
                int jb = inputData.SDP_blockIndex[l][k2];
                int jnz = A[j].SDP_sp_block[jb].NonZeroEffect;
                //	printf("%d %d %d %d %d %d\n",i,ib,inz, j, jb,jnz);
                if (jnz < inz) {
                    up += jnz;
                }
#if 1
                else if ((jnz == inz) && (j < i)) {
                    up += jnz;
                }
#endif
            }
            upNonZeroCount[i * SDP_nBlock + l] = up;
            // rMessage("up = " << up);
        }
    }

    // Determine which formula
    for (int l = 0; l < SDP_nBlock; ++l) {
        int countf1, countf2, countf3;
        countf1 = countf2 = countf3 = 0;
        for (int k = 0; k < inputData.SDP_nConstraint[l]; k++) {
            int i = inputData.SDP_constraint[l][k];
            int ib = inputData.SDP_blockIndex[l][k];
            dd_real inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

            dd_real f1, f2, f3;
            dd_real n = inputData.A[i].SDP_sp_block[ib].nRow;
            dd_real up = upNonZeroCount[i * SDP_nBlock + l];

            f1 = Kappa * n * inz + n * n * n + Kappa * up;
            f2 = Kappa * n * inz + Kappa * (n + 1) * up;
#if 1
            f3 = Kappa * (2 * Kappa * inz + 1) * up / Kappa;
#else
            f3 = Kappa * (2 * Kappa * inz + 1) * up;
#endif
            // rMessage("up = " << up << " nonzero = " << nonzero);
            // rMessage("f1=" << f1 << " f2=" << f2 << " f3=" << f3);
            // printf("%d %d %lf %lf %lf %lf\n",k,l,nonzero,f1,f2,f3);
            if (inputData.A[i].SDP_sp_block[ib].type == SparseMatrix::DENSE) {
                // if DENSE, we use only F1 or F2,
                // that is we don't use F3
                if (f1 < f2) {
                    useFormula[i * SDP_nBlock + l] = F1;
                    countf1++;
                } else {
                    useFormula[i * SDP_nBlock + l] = F2;
                    countf2++;
                }
            } else {
                // this case is SPARSE
                if (f1 < f2 && f1 < f3) {
                    //	   rMessage("line " << k << " is F1");
                    useFormula[i * SDP_nBlock + l] = F1;
                    countf1++;
                } else if (f2 < f3) {
                    //	   rMessage("line " << k << " is F2");
                    useFormula[i * SDP_nBlock + l] = F2;
                    countf2++;
                } else {
                    //	   rMessage("line " << k << " is F3");
                    useFormula[i * SDP_nBlock + l] = F3;
                    countf3++;
                }
            }
        }
// rMessage("Kappa = " << Kappa);
#if 0
    rMessage("count f1 = " << countf1
	     << ":: count f2 = " << countf2
	     << ":: count f3 = " << countf3);
#endif
    } // end of 'for (int l)'

    if (upNonZeroCount != NULL) {
        delete[] upNonZeroCount;
    }
    upNonZeroCount = NULL;

    return;
}

void Newton::compute_rMat(Newton::WHICH_DIRECTION direction, AverageComplementarity &mu, DirectionParameter &beta, Solutions &currentPt, WorkVariables &work) {

    //     CORRECTOR ::  r_zinv = (-XZ -dXdZ + mu I)Z^{-1}
    // not CORRECTOR ::  r_zinv = (-XZ + mu I)Z^{-1}
    dd_real target = beta.value * mu.current;
    Lal::let(r_zinvMat, '=', currentPt.invzMat, '*', &target);
    Lal::let(r_zinvMat, '=', r_zinvMat, '+', currentPt.xMat, &MMONE);

    if (direction == CORRECTOR) {
        // work.DLS1 = Dx Dz Z^{-1}
        Jal::ns_jordan_triple_product(work.DLS1, DxMat, DzMat, currentPt.invzMat, work.DLS2);
        Lal::let(r_zinvMat, '=', r_zinvMat, '+', work.DLS1, &MMONE);
    }

    //  rMessage("r_zinvMat = ");
    //  r_zinvMat.display();
}

void Newton::Make_gVec(Newton::WHICH_DIRECTION direction, InputData &inputData, Solutions &currentPt, Residuals &currentRes, AverageComplementarity &mu, DirectionParameter &beta, Phase &phase, WorkVariables &work, ComputeTime &com) {
    TimeStart(START1);
    // rMessage("mu = " << mu.current);
    // rMessage("beta = " << beta.value);
    compute_rMat(direction, mu, beta, currentPt, work);

    TimeEnd(END1);

    com.makerMat += TimeCal(START1, END1);

    TimeStart(START2);
    TimeStart(START_GVEC_MUL);

    // work.DLS1 = R Z^{-1} - X D Z^{-1} = r_zinv - X D Z^{-1}
    if (phase.value == SolveInfo::pFEAS || phase.value == SolveInfo::noINFO) {

        if (direction == CORRECTOR) {
            // x_rd_zinvMat is computed in PREDICTOR step
            Lal::let(work.DLS1, '=', r_zinvMat, '+', x_rd_zinvMat, &MMONE);
        } else {
            // currentPt is infeasilbe, that is the residual
            // dualMat is not 0.
            //      x_rd_zinvMat = X D Z^{-1}
            Jal::ns_jordan_triple_product(x_rd_zinvMat, currentPt.xMat, currentRes.dualMat, currentPt.invzMat, work.DLS2);
            Lal::let(work.DLS1, '=', r_zinvMat, '+', x_rd_zinvMat, &MMONE);
        } // if (direction == CORRECTOR)

    } else {
        // dualMat == 0
        work.DLS1.copyFrom(r_zinvMat);
    }

    //  rMessage("work.DLS1");
    //  work.DLS1.display();

    TimeEnd(END_GVEC_MUL);
    com.makegVecMul += TimeCal(START_GVEC_MUL, END_GVEC_MUL);

    inputData.multi_InnerProductToA(work.DLS1, gVec);
    Lal::let(gVec, '=', gVec, '*', &MMONE);
    // rMessage("gVec =  ");
    // gVec.display();

#if 0
  if (phase.value == SolveInfo:: dFEAS
      || phase.value == SolveInfo::noINFO) {
#endif
    Lal::let(gVec, '=', gVec, '+', currentRes.primalVec);
#if 0
  }
#endif

    TimeEnd(END2);
    com.makegVec += TimeCal(START2, END2);
}

void Newton::calF1(dd_real &ret, DenseMatrix &G, SparseMatrix &Aj) { Lal::let(ret, '=', Aj, '.', G); }

void Newton::calF2(dd_real &ret, DenseMatrix &F, DenseMatrix &G, DenseMatrix &X, SparseMatrix &Aj, bool &hasF2Gcal) {
    int alpha, beta;
    dd_real value1, value2;

    int n = Aj.nRow;
    // rMessage(" using F2 ");
    switch (Aj.type) {
    case SparseMatrix::SPARSE:
        // rMessage("F2::SPARSE  " << Aj.NonZeroCount);
        ret = 0.0;
        for (int index = 0; index < Aj.NonZeroCount; ++index) {
            alpha = Aj.row_index[index];
            beta = Aj.column_index[index];
            value1 = Aj.sp_ele[index];

            // value2 = F77_FUNC (ddot, DDOT)(&n, &X.de_ele[alpha+n*0], &n,
            //	     &F.de_ele[0+n*beta], &IONE);
            value2 = Rdot(n, X.de_ele + alpha, n, F.de_ele + (n * beta), 1);
            ret += value1 * value2;
            if (alpha != beta) {
                // value2 = F77_FUNC (ddot, DDOT)(&n, &X.de_ele[beta+n*0], &n,
                //        &F.de_ele[0+n*alpha], &IONE);
                value2 = Rdot(n, X.de_ele + beta, n, F.de_ele + (n * alpha), 1);
                ret += value1 * value2;
            }
        }
        break;
    case SparseMatrix::DENSE:
        // G is temporary matrix
        // rMessage("F2::DENSE");
        if (hasF2Gcal == false) {
            // rMessage(" using F2 changing to F1");
            Lal::let(G, '=', X, '*', F);
            hasF2Gcal = true;
        }
        Lal::let(ret, '=', Aj, '.', G);
        break;
    } // end of switch
}

void Newton::calF3(dd_real &ret, DenseMatrix &F, DenseMatrix &G, DenseMatrix &X, DenseMatrix &invZ, SparseMatrix &Ai, SparseMatrix &Aj) {
    // Ai and Aj are SPARSE
    ret = 0.0;
    dd_real sum;
    // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
    for (int index1 = 0; index1 < Aj.NonZeroCount; ++index1) {
        int alpha = Aj.row_index[index1];
        int beta = Aj.column_index[index1];
        dd_real value1 = Aj.sp_ele[index1];
        sum = 0.0;
        for (int index2 = 0; index2 < Ai.NonZeroCount; ++index2) {
            int gamma = Ai.row_index[index2];
            int delta = Ai.column_index[index2];
            dd_real value2 = Ai.sp_ele[index2];
            dd_real plu = value2 * invZ.de_ele[delta + invZ.nCol * beta] * X.de_ele[alpha + X.nCol * gamma];
            sum += plu;
            if (gamma != delta) {
                dd_real plu2 = value2 * invZ.de_ele[gamma + invZ.nCol * beta] * X.de_ele[alpha + X.nCol * delta];
                sum += plu2;
            }
        }
        ret += value1 * sum;
        if (alpha == beta) {
            continue;
        }
        sum = 0.0;
        for (int index2 = 0; index2 < Ai.NonZeroCount; ++index2) {
            int gamma = Ai.row_index[index2];
            int delta = Ai.column_index[index2];
            dd_real value2 = Ai.sp_ele[index2];
            dd_real plu = value2 * invZ.de_ele[delta + invZ.nCol * alpha] * X.de_ele[beta + X.nCol * gamma];
            sum += plu;
            if (gamma != delta) {
                dd_real plu2 = value2 * invZ.de_ele[gamma + invZ.nCol * alpha] * X.de_ele[beta + X.nCol * delta];
                sum += plu2;
            }
        }
        ret += value1 * sum;
    } // end of 'for (index1)'
    return;
}

void Newton::compute_bMat_dense_SDP(InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    int m = currentPt.mDim;
    int SDP_nBlock = inputData.SDP_nBlock;

    for (int l = 0; l < SDP_nBlock; ++l) {
        DenseMatrix &xMat = currentPt.xMat.SDP_block[l];
        DenseMatrix &invzMat = currentPt.invzMat.SDP_block[l];
        DenseMatrix &work1_master = work.DLS1.SDP_block[l];
        DenseMatrix &work2_master = work.DLS2.SDP_block[l];
        const int nConstraint = inputData.SDP_nConstraint[l];

        // ------------------------------------------------------------------
        // Decide whether to run the k1 loop in parallel.
        //
        // Threading is only safe if k1 -> i is injective within this block, because the
        // proof that bMat writes are disjoint relies on each constraint index having a
        // single owner. SDP_constraint[l] is built by appending i once per sub-block of
        // A_i that lands in block l, so a repeat is possible in principle. Check it
        // rather than assume it, and fall back to serial if it does not hold.
        //
        // Also require enough work to be worth a fork/join, in the spirit of SDPB's
        // minimal_split_factor.
        // ------------------------------------------------------------------
        bool injective = true;
#ifdef _OPENMP
        {
            std::vector<char> seen(m, 0);
            for (int k = 0; k < nConstraint; ++k) {
                const int ii = inputData.SDP_constraint[l][k];
                if (ii < 0 || ii >= m || seen[ii]) {
                    injective = false;
                    break;
                }
                seen[ii] = 1;
            }
        }
        bool anyF12 = false;
        for (int k = 0; k < nConstraint; ++k) {
            const FormulaType f = useFormula[inputData.SDP_constraint[l][k] * SDP_nBlock + l];
            if (f == F1 || f == F2) {
                anyF12 = true;
                break;
            }
        }
        // If the block has any F1/F2 constraint AND its setup gemm is big enough for
        // Rgemm to thread well, leave k1 serial and let Rgemm own the parallelism.
        // Threading k1 there makes each blockDim^3 gemm serial inside a thread and gains
        // nothing -- measured 1.35x SLOWER than upstream on gpp124-1 (blockDim 124).
        // The test is on PRESENCE, not count: a single F1 constraint can dominate the
        // block's cost, so a count-based majority test misses it (gpp124-1 has few F1
        // constraints but they carry ~97% of the bMat time). Blocks with no F1/F2 at all
        // (arch0, truss5) have no setup gemm and keep their large k1-threading win.
        const double setup_gemm =
            (double)xMat.nRow * (double)xMat.nRow * (double)xMat.nRow;
        // Let Rgemm own the block, rather than threading k1, under two conditions:
        //   (a) a single setup gemm is big enough for Rgemm to thread well at all, and
        //   (b) there is NOT abundant k1 work relative to the block dimension.
        //
        // (b) is an EMPIRICALLY CALIBRATED PROXY, not a cost model. It compares constraint
        // count against block dimension (nConstraint >= 2*blockDim) and nothing else -- it
        // does not total the setup gemms and weigh them against the k1 x k2 pair work,
        // which is what a real cost comparison would do. It is kept because it is the only
        // rule of four tried that got both ends right, not because it is derived:
        //   theta3    1106 constraints over a 150 block (7.4x) -- wants k1 threading
        //   gpp124-1   125 constraints over a 124 block (1.0x) -- wants Rgemm
        // Measured over 4 runs each, both effects consistent and outside run-to-run noise:
        //   gpp124-1  upstream ~0.542  k1-threaded ~0.706  Rgemm-owned ~0.448
        //   theta3    upstream ~12.6   k1-threaded ~12.2   Rgemm-owned ~13.9
        // Without (b) theta3 would be handed to Rgemm and lose the k1 parallelism it has in
        // abundance (1.18x slower than upstream); without (a) the control* family
        // (blockDim 25-80) would be handed blocks Rgemm cannot thread, forfeiting ~3x.
        // Do not re-derive this from a constraint-type count without new benchmarks.
        const bool enough_k1_work =
            (double)nConstraint >= 2.0 * (double)xMat.nRow;
        const bool rgemm_owns_block = anyF12 && (setup_gemm >= SDPA_OMP_RGEMM_OWNS_BLOCK) &&
                                      !enough_k1_work;
        const bool par = injective && !rgemm_owns_block && nConstraint >= SDPA_OMP_MIN_CONSTRAINTS &&
                         (double)nConstraint * (double)nConstraint * (double)xMat.nRow >= SDPA_OMP_MIN_BMAT_WORK;

        // Cap the team by the work actually available. There are exactly nConstraint k1
        // tasks, so a larger team only pays fork/join cost and privatises scratch for
        // workers that will never be handed a task. This is reachable in practice: the
        // work threshold admits blocks with as few as SDPA_OMP_MIN_CONSTRAINTS (8)
        // constraints, which on a 24-core machine would otherwise start 24 threads.
        int max_threads = omp_get_max_threads();
        if (max_threads > nConstraint)
            max_threads = nConstraint < 1 ? 1 : nConstraint;

        // Then cap so privatising work1/work2 cannot exceed the memory budget.
        // Only relevant when the block actually has F1/F2 constraints; F3 needs no scratch.
        if (anyF12 && max_threads > 1) {
            const double per_thread_mb =
                2.0 * (double)work1_master.nRow * (double)work1_master.nCol *
                sdpa_omp_bytes_per_elem() / 1048576.0;
            if (per_thread_mb > 0.0) {
                const int allowed = 1 + (int)(SDPA_OMP_MAX_PRIV_MB / per_thread_mb);
                if (allowed < max_threads)
                    max_threads = allowed < 1 ? 1 : allowed;
            }
        }
#else
        const bool par = false;
        const bool anyF12 = true;
        (void)injective;
#endif

        double acc_pre = 0.0, acc_f1 = 0.0, acc_f2 = 0.0, acc_f3 = 0.0;

        // The body is a lambda so that the SERIAL path can run without entering any
        // OpenMP construct at all. This matters: "#pragma omp parallel if(false)" still
        // creates a parallel region (a team of one), which makes every inner Rgemm call
        // *nested* -- and nested parallelism is off by default, so Rgemm's own threading
        // would be silently disabled. On gpp124-1 that cost 7.7x in bMat (0.035s -> 0.269s)
        // versus upstream, because the k1 loop did not engage while Rgemm's threading was
        // lost anyway.
        auto run_k1 = [&](int k1, DenseMatrix *w1, DenseMatrix *w2,
                          double &a_pre, double &a_f1, double &a_f2, double &a_f3,
                          bool may_need_priv, bool &owns_priv,
                          DenseMatrix &priv1, DenseMatrix &priv2) {
            // Per-thread scratch. Thread 0 (and the serial case) reuses the existing
            // per-block work matrices, so only the extra threads allocate.
                int i = inputData.SDP_constraint[l][k1];
                int ib = inputData.SDP_blockIndex[l][k1];
                int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
                SparseMatrix &Ai = inputData.A[i].SDP_sp_block[ib];

                FormulaType formula = useFormula[i * SDP_nBlock + l];

                // Plain locals, not TimeStart/TimeEnd: those macros declare `static
                // double`, which would be shared across threads.
                const double t_start1 = Time::rGetUseTime();
                const double t_start2 = t_start1;

                if (may_need_priv && !owns_priv && (formula == F1 || formula == F2)) {
                    priv1.initialize(work1_master.nRow, work1_master.nCol, work1_master.type);
                    priv2.initialize(work2_master.nRow, work2_master.nCol, work2_master.type);
                    owns_priv = true;
                }
                DenseMatrix &work1 = owns_priv ? priv1 : *w1;
                DenseMatrix &work2 = owns_priv ? priv2 : *w2;

                bool hasF2Gcal = false;
                if (formula == F1) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    Lal::let(work2, '=', xMat, '*', work1);
                } else if (formula == F2) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    hasF2Gcal = false;
                }
                a_pre += Time::rGetUseTime() - t_start2;

                for (int k2 = 0; k2 < nConstraint; k2++) {
                    int j = inputData.SDP_constraint[l][k2];
                    int jb = inputData.SDP_blockIndex[l][k2];
                    int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
                    SparseMatrix &Aj = inputData.A[j].SDP_sp_block[jb];

                    // Select the formula A[i] or the formula A[j].
                    // Use formula that has more NonZeroEffects than others.
                    // We must calculate i==j.
                    // This test is also what makes the bMat writes below disjoint across
                    // k1: it gives each unordered pair {i,j} exactly one owner.
                    if ((inz < jnz) || ((inz == jnz) && (i < j))) {
                        continue;
                    }

                    dd_real value;
                    switch (formula) {
                    case F1:
                        calF1(value, work2, Aj);
                        break;
                    case F2:
                        calF2(value, work1, work2, xMat, Aj, hasF2Gcal);
                        break;
                    case F3:
                        calF3(value, work1, work2, xMat, invzMat, Ai, Aj);
                        break;
                    } // end of switch
                    // Write the LOWER triangle only (row >= col).
                    //
                    // Every consumer of the dense bMat is Lower-only: Rpotrf("Lower") in
                    // Lal::choleskyFactorWithAdjust, and the two Rtrsv("Lower") in
                    // Lal::solveSystems that the '/' operator dispatches to. The strict
                    // upper half was therefore accumulated on every iteration and never
                    // read. The one routine that would read it, Newton::permuteMat, has
                    // no call sites -- see the note on its definition below.
                    //
                    // Disjointness across k1 is unchanged: the (inz, i) vs (jnz, j) test
                    // above gives each unordered pair {i, j} exactly one owner, and this
                    // writes a strict subset of what that owner wrote before.
                    const int brow = (i > j) ? i : j;
                    const int bcol = (i > j) ? j : i;
                    bMat.de_ele[brow + m * bcol] += value;
                } // end of 'for (int j)'

                const double t = Time::rGetUseTime() - t_start1;
                switch (formula) {
                case F1:
                    a_f1 += t;
                    break;
                case F2:
                    a_f2 += t;
                    break;
                case F3:
                    a_f3 += t;
                    break;
                }
        }; // end of run_k1 lambda

        // Decide AFTER every cap, not before. `par` is computed from the work thresholds,
        // but max_threads is then reduced by the constraint count and the scratch-memory
        // budget, and either can bring it to 1. Entering `omp parallel num_threads(1)`
        // creates a team of one, which is exactly the case the serial path below exists to
        // avoid: it makes any inner Rgemm call nested, and nested parallelism is off by
        // default, so Rgemm's own threading is silently lost. That is most likely to bite
        // large GMP blocks, where the memory cap does reduce the team.
        // !omp_in_parallel() additionally keeps this correct if the routine is ever reached
        // from an enclosing parallel region.
        // The WHOLE decision is inside the guard: max_threads exists only when _OPENMP is
        // defined, so referencing it outside fails to compile in a serial build.
#ifdef _OPENMP
        const bool use_parallel = par && max_threads > 1 && !omp_in_parallel();
#else
        const bool use_parallel = false;
#endif

        if (use_parallel) {
#ifdef _OPENMP
#pragma omp parallel num_threads(max_threads) reduction(+ : acc_pre, acc_f1, acc_f2, acc_f3)
#endif
            {
                // Scratch is only needed by threads other than 0, and only for F1/F2, so
                // allocate LAZILY on the first F1/F2 constraint a thread actually reaches.
                // A block may hold a handful of F1 constraints whose cost rounds to zero;
                // allocating eagerly for every thread then wastes 2*blockDim^2*32 bytes
                // each. On theta3 (blockDim 150, 24 threads) that was +15 MB for nothing.
                DenseMatrix *w1 = &work1_master;
                DenseMatrix *w2 = &work2_master;
                DenseMatrix priv1, priv2;
                bool owns_priv = false;
                bool may_need_priv = false;
#ifdef _OPENMP
                may_need_priv = (omp_get_num_threads() > 1 && omp_get_thread_num() > 0);
#endif
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (int k1 = 0; k1 < nConstraint; k1++)
                    run_k1(k1, w1, w2, acc_pre, acc_f1, acc_f2, acc_f3,
                           may_need_priv, owns_priv, priv1, priv2);
                if (owns_priv) {
                    priv1.terminate();
                    priv2.terminate();
                }
            }
        } else {
            // No OpenMP construct at all here, so inner Rgemm/Rdot keep their own threading.
            DenseMatrix priv1, priv2;
            bool owns_priv = false;
            for (int k1 = 0; k1 < nConstraint; k1++)
                run_k1(k1, &work1_master, &work2_master, acc_pre, acc_f1, acc_f2, acc_f3,
                       false, owns_priv, priv1, priv2);
        }

        com.B_PRE += acc_pre;
        com.B_F1 += acc_f1;
        com.B_F2 += acc_f2;
        com.B_F3 += acc_f3;
    }     // end of 'for (int l)'
}

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23: added. Ported from gmp.

   Counts what the sparse assembly's pair loop actually contains, per SDP block. The column
   that motivates it is `f2dense_after_first`: the number of pairs that CONSUME hasF2Gcal
   without being the first such pair in their group -- i.e. exactly the pairs whose behaviour
   was undefined before the group-scope fix above. Reasoning about reachability from the source
   got its answer wrong in both forks, so this prints what the data says instead.

   Also the eligibility census the assembly-threading work needs: group count, largest group,
   and the formula mix decide whether parallelising over i-groups can balance at all. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23: the ASSEMBLY's own canonical stream.

   The assembled sparse bMat, serialised with the same framing the factor stream uses (shared
   writer in sdpa_linear.cpp) but its OWN TAG, so an assembly stream can never compare equal to
   a factor stream by accident. This is the oracle the Phase-4 threading will be judged by, and
   it exists BEFORE that threading deliberately: the previous port shipped its parallelism
   first and its factor oracle months later, which is how "byte-identical" came to mean one
   printed line for a whole release.

   Emitting the assembly's own stream rather than inferring identity from the factor's is not
   fussiness: the factor is computed FROM the assembled matrix, so a Cholesky that happened to
   be insensitive to a small assembly error would mask it. See git log. */
void Newton::emit_bMat_stream(FILE *dump, bool want_fingerprint) {
    CanonicalStream d;
    canonicalInit(d, dump);
    const char *tag = "DDBMATASMv1"; // never equal to the factor stream's DDSPCHOLv1
    canonicalBytes(d, tag, strlen(tag));
    canonicalI64(d, (long long)bMat_type);
    if (bMat_type == SPARSE) {
        canonicalI64(d, sparse_bMat.nRow);
        canonicalI64(d, sparse_bMat.nCol);
        canonicalI64(d, sparse_bMat.NonZeroNumber);
        canonicalI64(d, sparse_bMat.NonZeroCount);
        canonicalU64(d, (uint64_t)sparse_bMat.NonZeroCount);
        for (int k = 0; k < sparse_bMat.NonZeroCount; ++k) {
            canonicalByte(d, 'E'); // record tag
            canonicalI64(d, k);
            canonicalI64(d, sparse_bMat.row_index[k]);
            canonicalI64(d, sparse_bMat.column_index[k]);
            // Both IEEE limbs: a dd_real IS its 128 bits, so no formatting step is involved.
            canonicalDouble(d, sparse_bMat.sp_ele[k].x[0]);
            canonicalDouble(d, sparse_bMat.sp_ele[k].x[1]);
            d.records++;
        }
    } else {
        canonicalI64(d, bMat.nRow);
        canonicalI64(d, bMat.nCol);
        const long long n = sdpaProduct(bMat.nRow, bMat.nCol);
        canonicalU64(d, (uint64_t)n);
        for (long long k = 0; k < n; ++k) {
            canonicalByte(d, 'D');
            canonicalI64(d, k);
            canonicalDouble(d, bMat.de_ele[k].x[0]);
            canonicalDouble(d, bMat.de_ele[k].x[1]);
            d.records++;
        }
    }
    canonicalByte(d, '.'); // terminator, so a truncated stream cannot look complete
    if (d.io_error) {
        rError("SDPA_BMAT_ASM_DUMP: writing failed after " << (unsigned long long)d.bytes
               << " bytes. The dump is truncated, so any comparison against it would be"
               << " meaningless; failing rather than leaving a file that looks complete.");
    }
    if (want_fingerprint) {
        // Counts printed alongside the fingerprint deliberately: equality of a 64-bit hash is
        // strong evidence, not proof, and two streams of different length or record count are
        // not the same matrix whatever their hashes do.
        printf("bmat assembly    : %llu records, %llu stream bytes, fingerprint %016llx\n",
               (unsigned long long)d.records, (unsigned long long)d.bytes,
               (unsigned long long)d.fnv);
        fflush(stdout);
    }
}

void Newton::census_bMat_sparse_SDP(InputData &inputData, FILE *fp) {
    const int SDP_nBlock_local = inputData.SDP_nBlock;
    long long tot_pairs = 0, tot_groups = 0, tot_f2dense_after = 0, tot_f2dense_total = 0;
    fprintf(fp, "bmat census: block  pairs groups maxgrp   F1    F2    F3  denseAj  f2dense_after_first  f2dense_total\n");
    for (int l = 0; l < SDP_nBlock_local; ++l) {
        long long pairs = 0, groups = 0, maxgrp = 0, nf[3] = {0, 0, 0};
        long long denseAj = 0, f2dense_after = 0, f2dense_total = 0;
        int previous_i = -1;
        long long cur = 0;
        long long seen_in_group_f2dense = 0;
        for (int iter = 0; iter < SDP_number[l]; ++iter) {
            const int i = SDP_constraint1[l][iter];
            const int j = SDP_constraint2[l][iter];
            const int jb = SDP_blockIndex2[l][iter];
            const FormulaType formula = useFormula[i * SDP_nBlock + l];
            if (i != previous_i) {
                if (cur > maxgrp) {
                    maxgrp = cur;
                }
                groups++;
                cur = 0;
                seen_in_group_f2dense = 0;
            }
            cur++;
            pairs++;
            nf[(int)formula]++;
            const bool aj_dense = (inputData.A[j].SDP_sp_block[jb].type == SparseMatrix::DENSE);
            if (aj_dense) {
                denseAj++;
            }
            if (formula == F2 && aj_dense) {
                f2dense_total++;
                // The first such pair in a group SETS hasF2Gcal; any after it READ it. Before
                // the group-scope fix, those reads were of indeterminate storage.
                if (seen_in_group_f2dense > 0) {
                    f2dense_after++;
                }
                seen_in_group_f2dense++;
            }
            previous_i = i;
        }
        if (cur > maxgrp) {
            maxgrp = cur;
        }
        fprintf(fp, "bmat census: %5d %6lld %6lld %6lld %4lld %5lld %5lld %8lld %10lld %13lld\n",
                l, pairs, groups, maxgrp, nf[0], nf[1], nf[2], denseAj, f2dense_after, f2dense_total);
        tot_pairs += pairs;
        tot_groups += groups;
        tot_f2dense_after += f2dense_after;
        tot_f2dense_total += f2dense_total;
    }
    fprintf(fp, "bmat census: TOTAL pairs=%lld groups=%lld f2dense_total=%lld f2dense_after_first=%lld\n",
            tot_pairs, tot_groups, tot_f2dense_total, tot_f2dense_after);
    fflush(fp);
}

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23: threaded sparse bMat assembly.

   PARALLEL OVER i-GROUPS. WHETHER BLOCKS MAY OVERLAP IN TIME IS DECIDED PER PROBLEM, from the
   finished map -- see SDP_blocks_disjoint. THE INVARIANT, in the order it must be read:

     1. WITHIN a block the destination slots are INJECTIVE: distinct iters of one block write
        distinct sparse_bMat entries. This is CHECKED, not argued --  make_aggrigateIndex_SDP
        walks the finished map and refuses a block that writes one entry twice.
     2. Two different blocks may overlap in time ONLY when the finished map proves their
        destination sets are DISJOINT. That is what SDP_blocks_disjoint records; it is computed
        in the same pass as (1).
     3. Otherwise each block's worksharing loop supplies the BARRIER, so all threads are in the
        same block at once. dE3 and dE4 are examples of THIS arm, not the fast one: both report
        blocks_disjoint=0. Take no example from the CI fixture here -- it is block-diagonal by
        construction and takes the fast arm, which is exactly why it could not catch the race.
     4. work1/work2 are per-GROUP scratch: written once when a new i begins, read by every pair
        in that group. Concurrent groups would overwrite each other's, so each worker gets its
        OWN pair of matrices. This is the one piece of gmp's scratch that dd genuinely needs --
        what dd does NOT need is gmp's per-flop scalar temporary, because a dd_real allocates
        nothing.
     5. hasF2Gcal is group state, and lives in assemble_group, so it is per worker and per group
        without anyone having to remember to make it so.

   An earlier version wrote `nowait` unconditionally, on the strength of (1) alone. (1) is true
   and does not imply (2). See the note at the loop.

   Bit-identity therefore holds by construction, as for the Cholesky: every element receives
   exactly the same contributions computed by the same expressions, and no sum is reordered.
   Asserted rather than assumed -- the assembly stream (emit_bMat_stream) compares the whole
   matrix across thread counts, and its mutation control proves that comparison can fail.

   Blocks are processed collectively: every thread walks the same block sequence so the
   worksharing constructs are encountered by all of them in the same order. See git log. */
#ifdef SDPA_PROBE_NOREALLOC
/* PROBE ONLY. Resize a private scratch matrix within capacity already held, without freeing or
   reallocating, zeroing exactly the logical region the shipped path zeroes.

   Exists to answer one question Phase 4 of the plan asks: on problems with many differently
   shaped blocks, how much of the threaded assembly's loss is the BARRIER and how much is
   DenseMatrix::initialize reallocating once per block per worker per iteration? Those two costs
   are currently confounded, and a gate recalibrated from a confounded measurement would be
   worthless.

   Falls back to initialize() if the capacity is somehow insufficient, so the probe can never
   silently read past the buffer. */
static void probe_resize_reusing(DenseMatrix &m, int nRow, int nCol, DenseMatrix::Type type,
                                 int cap_dim) {
    const long long want = (long long)nRow * (long long)nCol;
    const long long have = (long long)cap_dim * (long long)cap_dim;
    if (m.de_ele == NULL || type != DenseMatrix::DENSE || want > have) {
        m.initialize(nRow, nCol, type);
        return;
    }
    m.nRow = nRow;
    m.nCol = nCol;
    m.type = type;
    for (long long k = 0; k < want; ++k) {
        m.de_ele[k] = 0.0;
    }
}
#endif

/* One i-group's contribution to the Schur complement, in one place.

   Extracted when the assembly grew two worksharing forms -- with and without `nowait` -- so
   that the choice of form cannot change the arithmetic. Two copies of this body would be two
   things to keep in step, and the stream oracle compares the assembled matrix against the
   SERIAL routine, which is a third copy again. See git log. */
void Newton::assemble_group(InputData &inputData, int l, int lo, int hi, DenseMatrix &xMat,
                            DenseMatrix &invzMat, DenseMatrix &work1, DenseMatrix &work2) {
    // Group state: one worker owns this group entirely, so the flag is local to the call. Its
    // being function-local rather than hoisted is finding A's fix, kept structural.
    bool hasF2Gcal = false;
    const int i = SDP_constraint1[l][lo];
    const int ib = SDP_blockIndex1[l][lo];
    SparseMatrix &Ai = inputData.A[i].SDP_sp_block[ib];
    const FormulaType formula = useFormula[i * SDP_nBlock + l];

    if (formula == F1) {
        Lal::let(work1, '=', Ai, '*', invzMat);
        Lal::let(work2, '=', xMat, '*', work1);
    } else if (formula == F2) {
        Lal::let(work1, '=', Ai, '*', invzMat);
    }

    for (int iter = lo; iter < hi; ++iter) {
        const int j = SDP_constraint2[l][iter];
        const int jb = SDP_blockIndex2[l][iter];
        SparseMatrix &Aj = inputData.A[j].SDP_sp_block[jb];
        dd_real value;
        switch (formula) {
        case F1:
            calF1(value, work2, Aj);
            break;
        case F2:
            calF2(value, work1, work2, xMat, Aj, hasF2Gcal);
            break;
        case F3:
            calF3(value, work1, work2, xMat, invzMat, Ai, Aj);
            break;
        }
        // Disjoint within the block, and make_aggrigateIndex_SDP proves it rather than
        // asserting it: distinct iter -> distinct location.
        sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += value;
    }
}

bool Newton::compute_bMat_sparse_SDP_parallel(InputData &inputData, Solutions &currentPt,
                                              WorkVariables &work, ComputeTime &com,
                                              int team, const char *why_serial_out) {
    (void)why_serial_out;
#ifndef _OPENMP
    (void)inputData; (void)currentPt; (void)work; (void)com; (void)team;
    return false;
#else
    // Group boundaries per block, computed once. A group is a maximal run of pairs sharing i.
    // Building the index first keeps the parallel loop a flat, statically-sized iteration space,
    // which is what lets `omp for` distribute it without any thread scanning for boundaries.
    std::vector<std::vector<int> > gstart(SDP_nBlock);
    long long total_groups = 0;
    int max_dim = 0;
    for (int l = 0; l < SDP_nBlock; ++l) {
        int previous_i = -1;
        for (int iter = 0; iter < SDP_number[l]; ++iter) {
            const int i = SDP_constraint1[l][iter];
            if (i != previous_i) {
                gstart[l].push_back(iter);
                previous_i = i;
            }
        }
        gstart[l].push_back(SDP_number[l]); // sentinel: group g spans [gstart[g], gstart[g+1])
        total_groups += (long long)gstart[l].size() - 1;
        if (currentPt.xMat.SDP_block[l].nRow > max_dim) {
            max_dim = currentPt.xMat.SDP_block[l].nRow;
        }
    }
    if (total_groups <= 0) {
        return false;
    }

    bool ok = true;
    int actual_team = 1;
    long long groups_done = 0;
    // Read once, outside the region, so every thread branches on the same value.
    const bool blocks_disjoint = SDP_blocks_disjoint;
#pragma omp parallel num_threads(team) reduction(+ : groups_done) firstprivate(blocks_disjoint)
    {
        // THE team, read inside the region that owns it: num_threads is a REQUEST, and acting
        // on the requested size when the runtime gave fewer is the error that reopened gate 6
        // in the Cholesky review.
        const int here = omp_get_num_threads();
#pragma omp master
        actual_team = here;

        // Private scratch, one pair per worker, first sized to the largest block so the peak
        // allocation happens here rather than mid-loop. Thread 0 keeps using the shared work
        // matrices: with a team of one this is then exactly the serial path.
        const int me = omp_get_thread_num();
        DenseMatrix priv1, priv2;
        bool owns_priv = false;
        if (me != 0 && max_dim > 0) {
            priv1.initialize(max_dim, max_dim, DenseMatrix::DENSE);
            priv2.initialize(max_dim, max_dim, DenseMatrix::DENSE);
            owns_priv = true;
        }

        for (int l = 0; l < SDP_nBlock; ++l) {
            DenseMatrix &xMat = currentPt.xMat.SDP_block[l];
            DenseMatrix &invzMat = currentPt.invzMat.SDP_block[l];
            DenseMatrix &shared1 = work.DLS1.SDP_block[l];
            DenseMatrix &shared2 = work.DLS2.SDP_block[l];
            const int ngroups = (int)gstart[l].size() - 1;
            // Resized to THIS block. Note what initialize() actually does: it frees and
            // reallocates whenever nRow*nCol differs from the previous dimensions -- it does
            // NOT keep the largest buffer and use a leading submatrix, which is what an earlier
            // version of this comment claimed.
            //
            // MEASURED, so this is no longer a worry to carry: a probe build that reuses the
            // capacity instead (see probe_resize_reusing) removes 50-91% of ALL DenseMatrix
            // allocations on problems with heterogeneous blocks -- 12_min 8665 -> 748, 10_min
            // 5166 -> 587, 8_min 4079 -> 553 -- and moves the threaded assembly time by at most
            // 3.1%, and by 0.3% on the largest reduction. truss6, the loss case that prompted
            // the question, is 150 blocks of order 3 and one of order 1, so it barely
            // reallocates at all (3129 -> 3124). Allocation is therefore NOT a material part of
            // the small-block loss; the barrier and thread overhead are essentially all of it.
            if (owns_priv) {
#ifdef SDPA_PROBE_NOREALLOC
                /* MEASUREMENT PROBE, never in a release build (-DSDPA_PROBE_NOREALLOC).
                   Reuses the capacity already reserved at max_dim instead of reallocating, while
                   keeping this block's LOGICAL dimensions and zeroing exactly the same region the
                   shipped path zeroes.

                   It must not become a padding probe. The matrices are column-major with leading
                   dimension nRow -- every kernel indexes de_ele[row + nRow*col] -- so the live
                   region is exactly [0, nRow*nCol) and a larger buffer beyond it is never read
                   or written. Same arithmetic, same memory traffic, same results; the only thing
                   removed is the allocator call. Padding the COMPUTATION to max_dim would change
                   the flop count and would measure the padding, not the allocation. */
                probe_resize_reusing(priv1, shared1.nRow, shared1.nCol, shared1.type, max_dim);
                probe_resize_reusing(priv2, shared2.nRow, shared2.nCol, shared2.type, max_dim);
#else
                priv1.initialize(shared1.nRow, shared1.nCol, shared1.type);
                priv2.initialize(shared2.nRow, shared2.nCol, shared2.type);
#endif
            }
            DenseMatrix &work1 = owns_priv ? priv1 : shared1;
            DenseMatrix &work2 = owns_priv ? priv2 : shared2;

            // ONE loop body, TWO worksharing forms. The branch is on a firstprivate constant,
            // so it is collective: every thread takes the same arm and the constructs are
            // encountered by all of them in the same order.
            if (blocks_disjoint) {
                // No sparse_bMat entry is written by more than one block, so letting a worker
                // roll into block l+1 while another finishes block l cannot make two threads
                // meet on one element. This is the fast path and it is taken because the map
                // was CHECKED, not because the shape was assumed.
#pragma omp for schedule(dynamic) nowait
                for (int g = 0; g < ngroups; ++g) {
                    const int lo = gstart[l][g], hi = gstart[l][g + 1];
                    if (lo < hi) {
                        assemble_group(inputData, l, lo, hi, xMat, invzMat, work1, work2);
                        groups_done++;
                    }
                }
            } else {
#pragma omp for schedule(dynamic)
                for (int g = 0; g < ngroups; ++g) {
                    const int lo = gstart[l][g], hi = gstart[l][g + 1];
                    if (lo < hi) {
                        assemble_group(inputData, l, lo, hi, xMat, invzMat, work1, work2);
                        groups_done++;
                    }
                }
            }
            // WHY THE SECOND ARM ABOVE HAS NO `nowait`.
            //
            // The first version of this function had one, justified by "the bMat writes are +=
            // into distinct slots within a block". That is true and irrelevant: `nowait` lets
            // one worker start block l+1's groups while another is still in block l, and ACROSS
            // blocks the slots are NOT distinct -- a constraint pair nonzero in both blocks
            // accumulates into the same sparse_bMat entry from both. Two threads then did
            // read-modify-write on the same dd_real.
            //
            // It was invisible for two reasons, and the second is the uncomfortable one.
            //
            // The CI fixture (tests/gen_spchol_fixture.py) is block-diagonal by construction --
            // each constraint touches exactly one block -- so no slot is ever shared and the
            // race has nothing to land on. Every stream-identity assertion ran on it. A fixture
            // that cannot expose a defect is not evidence that the defect is absent.
            //
            // dE3 and dE4 are NOT immune: both report blocks_disjoint=0. They escaped notice
            // because at that size the race fires RARELY -- few blocks, long groups, a narrow
            // window -- and one dE4 run in SEVENTEEN on pi produced a differing assembled matrix.
            // A rare race compared once is a race not tested for.
            //
            // It is NOT confined to forced modes. SDPLIB truss6.dat-s routes sparse under the
            // default `auto` policy, and the released v7.1.3-omp.2 binary returned five
            // different objectives in five runs with no environment variables set at all:
            // 2.16e+06, 9.72e+04, -9.0093934e+02, 1.66e+04, 9.82e+05, against the correct
            // -9.0100140477088098e+02. 10_min.dat-s under SDPA_BMAT_MODE=sparse shows the
            // mechanism most plainly: 16,219 pairs are written into at most m(m+1)/2 = 2,775
            // slots, so cross-block sharing is forced by pigeonhole.
            //
            // The barrier restores the property the region is supposed to have: all threads are
            // in the same block at the same time. It is not free -- on an M1 Max it costs 1.8x
            // to 3.5x of the main loop on 10_min and 12_min, and 8.8x of the assembly kernel on
            // a 24-block fixture, because it is paid once per SDP block per iteration. So it is
            // paid only where it is needed, and where that is gets COMPUTED: see the pass at the
            // end of make_aggrigateIndex_SDP, which walks the finished map and asks whether any
            // sparse_bMat entry is written by two different blocks.
            //
            // Measured: arch0 (2 blocks) and the block-diagonal CI fixture take `nowait`;
            // truss6, 10_min, 12_min, the cross-block fixture, AND dE3 and dE4 take the barrier.
            // The two headline problems being on the barrier arm is the point worth remembering
            // -- it costs them +0.8% and +0.1%, and an earlier version of this comment asserted
            // the opposite from their block structure without measuring it. The distinction is
            // the one the released code got wrong by assuming rather than checking; do not
            // reintroduce it by assuming which arm a problem takes.
        }
        if (owns_priv) {
            priv1.terminate();
            priv2.terminate();
        }
    }
    if (bmat_asm_log()) {
#ifdef SDPA_PROBE_ALLOC_COUNT
        {
            extern unsigned long long g_probe_dense_allocs;
            rMessage("bmat asm probe: dense_allocs_so_far=" << g_probe_dense_allocs);
        }
#endif
        rMessage("bmat asm: team=" << actual_team << " groups=" << total_groups
                                   << " groups_executed=" << groups_done
                                   << " blocks_disjoint=" << (SDP_blocks_disjoint ? 1 : 0)
                                   << " overlap=" << (SDP_blocks_disjoint ? "nowait" : "barrier"));
    }
    // A team of one means the region did the serial thing; report it so the caller's log does
    // not claim concurrency that never happened.
    ok = ok && (actual_team >= 1);
    return ok;
#endif
}

void Newton::compute_bMat_sparse_SDP(InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    TimeStart(B_NDIAG_START1);
    TimeStart(B_NDIAG_START2);

    // Once per solve, and only when asked. Parsed strictly like every other knob: an
    // unrecognised value is refused rather than silently treated as "off", because a typo that
    // silently disables a census reads as "the census found nothing".
    {
        static bool census_done = false;
        // Already validated in Make_bMat; this only reads the answer.
        if (!census_done && bmat_asm_census_wanted()) {
            census_done = true;
            census_bMat_sparse_SDP(inputData, stdout);
        }
    }

#ifdef _OPENMP
    // ---- admission for the threaded assembly ------------------------------------------
    // Gates mirror the Cholesky's: a mode, a work floor, a nesting guard, and -- unique to the
    // assembly -- a MEMORY budget, because each extra worker needs two private dense matrices.
    {
        const BmatAsmMode amode = bmat_asm_mode();
        long long pairs = 0;
        int max_dim = 0;
        for (int l = 0; l < SDP_nBlock; ++l) {
            pairs += (long long)SDP_number[l];
            if (currentPt.xMat.SDP_block[l].nRow > max_dim) {
                max_dim = currentPt.xMat.SDP_block[l].nRow;
            }
        }
        int team = (amode == BMAT_ASM_SERIAL) ? 1 : omp_get_max_threads();
        const char *why = NULL;
        if (amode == BMAT_ASM_SERIAL) {
            why = "SDPA_BMAT_ASM_MODE=serial";
        } else if (omp_get_level() != 0) {
            why = "nested inside another OpenMP region";
            team = 1;
        } else {
            const int tl = omp_get_thread_limit();
            if (tl > 0 && tl < team) {
                team = tl;
            }
            // MEMORY. Every worker beyond thread 0 needs work1 and work2 privately:
            // 2 * max_dim^2 * sizeof(dd_real) bytes each. Computed in 64-bit and compared
            // before any narrowing, so a large but valid budget cannot wrap.
            const unsigned long long budget_mb =
                bmat_asm_u64("SDPA_BMAT_ASM_SCRATCH_MB", SDPA_BMAT_ASM_SCRATCH_MB_DEFAULT);
            const long long elems = sdpaProduct(max_dim, max_dim);
            if (elems > 0 && budget_mb > 0) {
                const unsigned long long per_worker =
                    2ULL * (unsigned long long)elems * (unsigned long long)sizeof(dd_real);
                const unsigned long long budget = budget_mb * 1048576ULL;
                unsigned long long affordable = 1ULL + budget / (per_worker ? per_worker : 1ULL);
                if (affordable > (unsigned long long)INT_MAX) {
                    affordable = (unsigned long long)INT_MAX;
                }
                if ((int)affordable < team) {
                    team = (int)affordable;
                    if (team < 2) {
                        why = "SDPA_BMAT_ASM_SCRATCH_MB budget admits only one worker";
                    }
                }
            }
            const unsigned long long min_pairs =
                bmat_asm_u64("SDPA_BMAT_ASM_MIN_PAIRS", SDPA_BMAT_ASM_MIN_PAIRS_DEFAULT);
            if (amode != BMAT_ASM_PARALLEL && (unsigned long long)pairs < min_pairs) {
                why = "below SDPA_BMAT_ASM_MIN_PAIRS";
                team = 1;
            }
        }
        if (team >= 2) {
            if (compute_bMat_sparse_SDP_parallel(inputData, currentPt, work, com, team, why)) {
                TimeEnd(B_NDIAG_END1);
                return;
            }
        } else if (bmat_asm_log()) {
            rMessage("bmat asm: serial (" << (why ? why : "team would be 1") << ")");
        }
    }
#endif

    for (int l = 0; l < SDP_nBlock; ++l) {
        DenseMatrix &xMat = currentPt.xMat.SDP_block[l];
        DenseMatrix &invzMat = currentPt.invzMat.SDP_block[l];
        DenseMatrix &work1 = work.DLS1.SDP_block[l];
        DenseMatrix &work2 = work.DLS2.SDP_block[l];
        int previous_i = -1;
        /* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23: hasF2Gcal is GROUP state, not
           per-pair state, and declaring it inside the pair loop made every pair after a group's
           first one read indeterminate storage.

           It is passed to calF2 BY REFERENCE, and calF2 both reads and writes it:

               if (hasF2Gcal == false) { Lal::let(G, '=', X, '*', F); hasF2Gcal = true; }
               Lal::let(ret, '=', Aj, '.', G);

           so an indeterminate `true` does not merely invoke undefined behaviour -- it SKIPS
           computing G = X*F and contracts Aj against a stale G, producing a silently wrong
           Schur-complement entry. The value is consumed only when the formula is F2 and Aj is
           dense, so reaching it needs a group with at least two such pairs; SDPA_BMAT_ASM_CENSUS
           counts exactly those pairs (f2dense_after_first) so reachability is measurable rather
           than argued. Ported from gmp, which fixed the same defect after a review supplied a
           -k counterexample to an earlier claim of unreachability. See git log. */
        bool hasF2Gcal = false;
        const bool f2_stale = bmat_test_f2_stale();

        for (int iter = 0; iter < SDP_number[l]; iter++) {
            //      TimeStart(B_NDIAG_START1);
            int i = SDP_constraint1[l][iter];
            int ib = SDP_blockIndex1[l][iter];
            SparseMatrix &Ai = inputData.A[i].SDP_sp_block[ib];
            FormulaType formula = useFormula[i * SDP_nBlock + l];

            if (i != previous_i) {
                // ---------------------------------------------------
                // formula = F3; // this is force change
                // ---------------------------------------------------
                TimeStart(B_NDIAG_START2);

                hasF2Gcal = false;
                if (formula == F1) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    Lal::let(work2, '=', xMat, '*', work1);
                } else if (formula == F2) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    hasF2Gcal = false;
                    // Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
                }
                TimeEnd(B_NDIAG_END2);
                com.B_PRE += TimeCal(B_NDIAG_START2, B_NDIAG_END2);
            } else if (f2_stale) {
                // NEGATIVE CONTROL: the worst-case value the old per-pair declaration could
                // hand to calF2 on a pair after its group's first.
                hasF2Gcal = true;
            }

            int j = SDP_constraint2[l][iter];
            int jb = SDP_blockIndex2[l][iter];
            SparseMatrix &Aj = inputData.A[j].SDP_sp_block[jb];

            dd_real value;
            switch (formula) {
            case F1:
                // rMessage("calF1");
                calF1(value, work2, Aj);
                break;
            case F2:
                // rMessage("calF2 ");
                calF2(value, work1, work2, xMat, Aj, hasF2Gcal);
                // calF1(value2,gMat.ele[l],A[j].ele[l]);
                // rMessage("calF2:  " << (value-value2));
                break;
            case F3:
                // rMessage("calF3");
                calF3(value, work1, work2, xMat, invzMat, Ai, Aj);
                break;
            } // end of switch
            sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += value;
            previous_i = i;
        } // end of 'for (int index)'
#if 0
    TimeEnd(B_NDIAG_END1);
    dd_real t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    switch (formula) {
    case F1: com.B_F1 += t; break;
    case F2: com.B_F2 += t; break;
    case F3: com.B_F3 += t; break;
    }
#endif
    } // end of 'for (int l)'
}

#if 0
void Newton::compute_bMat_dense_SCOP(InputData& inputData,
				     Solutions& currentPt,
				     WorkVariables& work,
				     ComputeTime& com)
{
    rError("current version does not support SOCP");
}

void Newton::compute_bMat_sparse_SOCP(InputData& inputData,
				      Solutions& currentPt,
				      WorkVariables& work,
				      ComputeTime& com)
{
    rError("current version does not support SOCP");
}
#endif

void Newton::compute_bMat_dense_LP(InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    int m = currentPt.mDim;
    int LP_nBlock = inputData.LP_nBlock;

    TimeEnd(B_DIAG_START1);
    for (int l = 0; l < LP_nBlock; ++l) {
        dd_real xMat = currentPt.xMat.LP_block[l];
        dd_real invzMat = currentPt.invzMat.LP_block[l];

        for (int k1 = 0; k1 < inputData.LP_nConstraint[l]; k1++) {
            int i = inputData.LP_constraint[l][k1];
            int ib = inputData.LP_blockIndex[l][k1];
            //	int inz = inputData.A[i].LP_sp_block[ib].NonZeroEffect;
            dd_real Ai = inputData.A[i].LP_sp_block[ib];

            for (int k2 = k1; k2 < inputData.LP_nConstraint[l]; k2++) {
                int j = inputData.LP_constraint[l][k2];
                int jb = inputData.LP_blockIndex[l][k2];
                //	  int jnz = inputData.A[j].LP_sp_block[jb].NonZeroEffect;
                dd_real Aj = inputData.A[j].LP_sp_block[jb];

                dd_real value;
                value = xMat * invzMat * Ai * Aj;

                // Lower triangle only -- see compute_bMat_dense_SDP above.
                const int brow = (i > j) ? i : j;
                const int bcol = (i > j) ? j : i;
                bMat.de_ele[brow + m * bcol] += value;
            } // end of 'for (int j)'
        }     // end of 'for (int i)'
    }         // end of 'for (int l)'
    TimeEnd(B_DIAG_END1);
    com.B_DIAG += TimeCal(B_DIAG_START1, B_DIAG_END1);
}

void Newton::compute_bMat_sparse_LP(InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    TimeEnd(B_DIAG_START1);
    for (int l = 0; l < LP_nBlock; ++l) {
        dd_real xMat = currentPt.xMat.LP_block[l];
        dd_real invzMat = currentPt.invzMat.LP_block[l];

        for (int iter = 0; iter < LP_number[l]; iter++) {
            int i = LP_constraint1[l][iter];
            int ib = LP_blockIndex1[l][iter];
            dd_real Ai = inputData.A[i].LP_sp_block[ib];

            int j = LP_constraint2[l][iter];
            int jb = LP_blockIndex2[l][iter];
            dd_real Aj = inputData.A[j].LP_sp_block[jb];

            dd_real value;
            value = xMat * invzMat * Ai * Aj;
            sparse_bMat.sp_ele[LP_location_sparse_bMat[l][iter]] += value;
        } // end of 'for (int iter)
    }     // end of 'for (int l)'
    TimeEnd(B_DIAG_END1);
    com.B_DIAG += TimeCal(B_DIAG_START1, B_DIAG_END1);
}

void Newton::Make_bMat(InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    TimeStart(START3);
    // Parse-and-validate the assembly knobs HERE, not inside the sparse routine: Make_bMat is
    // reached on BOTH routes, the sparse assembly is not. A validator only one route can reach
    // is not a validator -- this fork has now been bitten by that exact shape three times
    // (SDPA_SPCHOL_MODE, SDPA_SPCHOL_DIGEST_DUMP, and this).
    {
        static bool asm_cfg_done = false;
        if (!asm_cfg_done) {
            asm_cfg_done = true;
            (void)bmat_asm_census_wanted();
            (void)bmat_asm_profile_wanted();
            (void)bmat_asm_digest_wanted();
            (void)bmat_asm_mutate();
            (void)bmat_asm_mode();
            (void)bmat_asm_u64("SDPA_BMAT_ASM_MIN_PAIRS", SDPA_BMAT_ASM_MIN_PAIRS_DEFAULT);
            (void)bmat_asm_u64("SDPA_BMAT_ASM_SCRATCH_MB", SDPA_BMAT_ASM_SCRATCH_MB_DEFAULT);
            {
                const char *dp = bmat_asm_dump_path();
                if (dp != NULL) {
                    FILE *probe = fopen(dp, "ab");
                    if (probe == NULL) {
                        rError("SDPA_BMAT_ASM_DUMP: cannot open \"" << dp << "\" for append");
                    }
                    fclose(probe);
                }
            }
            (void)bmat_test_f2_stale();
        }
    }
    if (bMat_type == SPARSE) {
        /* SUBPHASE PROFILE (SDPA_BMAT_ASM_PROFILE=1). `Make bMat` is 36.8% of dE4's main loop
           and serial, but that figure covers the WHOLE phase: the sparse zeroing, the SDP
           assembly and the LP assembly. Attributing all of it to compute_bMat_sparse_SDP was
           an inference from the source, and this fork has now been wrong twice reasoning that
           way. So measure the three parts separately BEFORE anyone parallelises one of them.
           Costs three clock reads per iteration when off. */
        const bool prof = bmat_asm_profile_wanted();
        double t_zero = 0.0, t_sdp = 0.0, t_lp = 0.0;
        double t0 = prof ? Time::rGetUseTime() : 0.0;
        // set sparse_bMat zero
        for (int iter = 0; iter < sparse_bMat.NonZeroCount; ++iter) {
            sparse_bMat.sp_ele[iter] = 0.0;
        }
        if (prof) {
            const double t1 = Time::rGetUseTime();
            t_zero = t1 - t0;
            t0 = t1;
        }
        compute_bMat_sparse_SDP(inputData, currentPt, work, com);
        if (prof) {
            const double t1 = Time::rGetUseTime();
            t_sdp = t1 - t0;
            t0 = t1;
        }
        //   compute_bMat_sparse_SOCP(inputData,currentPt,work,com);
        compute_bMat_sparse_LP(inputData, currentPt, work, com);
        if (prof) {
            t_lp = Time::rGetUseTime() - t0;
            com.bmat_asm_zero += t_zero;
            com.bmat_asm_sdp += t_sdp;
            com.bmat_asm_lp += t_lp;
        }
    } else {
        // Keep this a FULL-matrix zero. Only the lower triangle is written below,
        // but leaving the strict upper half uninitialised would put indeterminate
        // values in a live allocation for no measurable gain.
        bMat.setZero();
        compute_bMat_dense_SDP(inputData, currentPt, work, com);
        //    compute_bMat_dense_SOCP(inputData,currentPt,work,com);
        compute_bMat_dense_LP(inputData, currentPt, work, com);
    }
    // rMessage("bMat =  ");
    // bMat.display();
    // sparse_bMat.display();
    TimeEnd(END3);
    com.makebMat += TimeCal(START3, END3);

    // The assembly oracle. Off unless asked for: O(nnz) work is fine for a fixture and not
    // something to pay for in a solve.
    {
        const char *dumpf = bmat_asm_dump_path();
        const bool want_fp = bmat_asm_digest_wanted();
        if (bmat_asm_mutate()) {
            // One ulp on the low limb of the middle stored element.
            if (bMat_type == SPARSE && sparse_bMat.NonZeroCount > 0) {
                double &lo = sparse_bMat.sp_ele[sparse_bMat.NonZeroCount / 2].x[1];
                lo = nextafter(lo, HUGE_VAL);
            } else if (bMat_type != SPARSE && bMat.nRow > 0) {
                double &lo = bMat.de_ele[(sdpaProduct(bMat.nRow, bMat.nCol)) / 2].x[1];
                lo = nextafter(lo, HUGE_VAL);
            }
        }
        if (dumpf != NULL || want_fp) {
            FILE *dump = NULL;
            if (dumpf != NULL) {
                dump = fopen(dumpf, "ab");
                if (dump == NULL) {
                    rError("SDPA_BMAT_ASM_DUMP: cannot open \"" << dumpf << "\" for append");
                }
            }
            emit_bMat_stream(dump, want_fp);
            if (dump != NULL) {
                if (fflush(dump) != 0 || ferror(dump) != 0 || fclose(dump) != 0) {
                    rError("SDPA_BMAT_ASM_DUMP: writing \"" << dumpf << "\" failed");
                }
            }
        }
    }
}

// nakata 2004/12/01
// WARNING: permuteMat is NOT CALLED anywhere in this tree. A whole-tree grep finds
// only this definition and the declaration in sdpa_newton.h, and the method is not
// virtual, so it cannot be reached indirectly either. It is also the only routine that
// would read the dense bMat's strict UPPER triangle: it copies arbitrary (i, j) chosen
// by ordering[]. Since 2026-08-05 the dense bMat is accumulated in its lower triangle
// only, so the strict upper half holds whatever bMat.setZero() left there, i.e. zero.
// A future caller must either mirror the lower half up first, or index with row >= col.
void Newton::permuteMat(DenseMatrix &bMat, SparseMatrix &sparse_bMat) {
    int i, j, k;
    int mDIM = bMat.nRow;

    for (k = 0; k < sparse_bMat.NonZeroCount; k++) {
        i = ordering[sparse_bMat.row_index[k]];
        j = ordering[sparse_bMat.column_index[k]];
        sparse_bMat.sp_ele[k] = bMat.de_ele[i + j * mDIM];
    }
}

// nakata 2004/12/01
void Newton::permuteVec(Vector &gVec, Vector &gVec2) {
    int i, k;
    int mDIM = gVec2.nDim;

    for (k = 0; k < mDIM; k++) {
        i = ordering[k];
        gVec2.ele[k] = gVec.ele[i];
    }
}

// nakata 2004/12/01
void Newton::reverse_permuteVec(Vector &DyVec2, Vector &DyVec) {
    int i, k;
    int mDIM = DyVec.nDim;

    for (k = 0; k < mDIM; k++) {
        i = ordering[k];
        DyVec.ele[i] = DyVec2.ele[k];
    }
}

bool Newton::compute_DyVec(Newton::WHICH_DIRECTION direction, InputData &inputData, Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    if (direction == PREDICTOR) {
        TimeStart(START3_2);

        if (bMat_type == SPARSE) {
            bool ret = Lal::getCholesky(sparse_bMat, diagonalIndex);
            if (ret == FAILURE) {
                return FAILURE;
            }
        } else {
            bool ret = Lal::choleskyFactorWithAdjust(bMat);
            if (ret == FAILURE) {
                return FAILURE;
            }
        }
        // rMessage("Cholesky of bMat =  ");
        // bMat.display();
        // sparse_bMat.display();
        TimeEnd(END3_2);
        com.choleskybMat += TimeCal(START3_2, END3_2);
    }
    // bMat is already cholesky factorized.

    TimeStart(START4);
    if (bMat_type == SPARSE) {
        permuteVec(gVec, work.DV1);
        Lal::let(work.DV2, '=', sparse_bMat, '/', work.DV1);
        reverse_permuteVec(work.DV2, DyVec);
    } else {
        Lal::let(DyVec, '=', bMat, '/', gVec);
    }
    TimeEnd(END4);
    com.solve += TimeCal(START4, END4);
    // rMessage("DyVec =  ");
    // DyVec.display();
    return _SUCCESS;
}

void Newton::compute_DzMat(InputData &inputData, Residuals &currentRes, Phase &phase, ComputeTime &com) {
    TimeStart(START_SUMDZ);
    inputData.multi_plusToA(DyVec, DzMat);
    Lal::let(DzMat, '=', DzMat, '*', &MMONE);
    if (phase.value == SolveInfo::pFEAS || phase.value == SolveInfo::noINFO) {
        Lal::let(DzMat, '=', DzMat, '+', currentRes.dualMat);
    }
    TimeEnd(END_SUMDZ);
    com.sumDz += TimeCal(START_SUMDZ, END_SUMDZ);
}

void Newton::compute_DxMat(Solutions &currentPt, WorkVariables &work, ComputeTime &com) {
    TimeStart(START_DX);
    // work.DLS1 = dX dZ Z^{-1}
    Jal::ns_jordan_triple_product(work.DLS1, currentPt.xMat, DzMat, currentPt.invzMat, work.DLS2);
    // dX = R Z^{-1} - dX dZ Z^{-1}
    Lal::let(DxMat, '=', r_zinvMat, '+', work.DLS1, &MMONE);
    TimeEnd(END_DX);
    TimeStart(START_SYMM);
    Lal::getSymmetrize(DxMat);
    TimeEnd(END_SYMM);
    // rMessage("DxMat =  ");
    // DxMat.display();
    com.makedX += TimeCal(START_DX, END_DX);
    com.symmetriseDx += TimeCal(START_SYMM, END_SYMM);
}

bool Newton::Mehrotra(Newton::WHICH_DIRECTION direction, InputData &inputData, Solutions &currentPt, Residuals &currentRes, AverageComplementarity &mu, DirectionParameter &beta, Switch &reduction, Phase &phase, WorkVariables &work, ComputeTime &com) {
    //   rMessage("xMat, yVec, zMat =  ");
    //   currentPt.xMat.display();
    //   currentPt.yVec.display();
    //   currentPt.zMat.display();

    Make_gVec(direction, inputData, currentPt, currentRes, mu, beta, phase, work, com);

    if (direction == PREDICTOR) {
        Make_bMat(inputData, currentPt, work, com);
    }

    // rMessage("gVec, bMat =  ");
    //   gVec.display();
    //   bMat.display();
    //   sparse_bMat.display();  //
    //   display_sparse_bMat();  // with reverse ordering

    bool ret = compute_DyVec(direction, inputData, currentPt, work, com);
    if (ret == FAILURE) {
        return FAILURE;
    }
    //  rMessage("cholesky factorization =  ");
    //  sparse_bMat.display();

    TimeStart(START5);

    compute_DzMat(inputData, currentRes, phase, com);
    compute_DxMat(currentPt, work, com);

    TimeEnd(END5);
    com.makedXdZ += TimeCal(START5, END5);

    // rMessage("DxMat, DyVec, DzMat =  ");
    //   DxMat.display();
    //   DyVec.display();
    //   DzMat.display();

    return true;
}

void Newton::display(FILE *fpout) {
    if (fpout == NULL) {
        return;
    }

    fprintf(fpout, "rNewton.DxMat = \n");
    DxMat.display(fpout);
    fprintf(fpout, "rNewton.DyVec = \n");
    DyVec.display(fpout);
    fprintf(fpout, "rNewton.DzMat = \n");
    DzMat.display(fpout);
}

void Newton::display_index(FILE *fpout) {
    if (fpout == NULL) {
        return;
    }
    printf("display_index: %d %d %d\n", SDP_nBlock, SOCP_nBlock, LP_nBlock);

    for (int b = 0; b < SDP_nBlock; b++) {
        printf("SDP:%dth block\n", b);
        for (int i = 0; i < SDP_number[b]; i++) {
            printf("cons1:%d const2:%d block1:%d block2:%d sp_bMat:%d \n", SDP_constraint1[b][i], SDP_constraint2[b][i], SDP_blockIndex1[b][i], SDP_blockIndex2[b][i], SDP_location_sparse_bMat[b][i]);
        }
    }

    for (int b = 0; b < SOCP_nBlock; b++) {
        printf("SOCP:%dth block\n", b);
        for (int i = 0; i < SOCP_number[b]; i++) {
            printf("cons1:%d const2:%d block1:%d block2:%d sp_bMat:%d \n", SOCP_constraint1[b][i], SOCP_constraint2[b][i], SOCP_blockIndex1[b][i], SOCP_blockIndex2[b][i], SOCP_location_sparse_bMat[b][i]);
        }
    }

    for (int b = 0; b < LP_nBlock; b++) {
        printf("LP:%dth block\n", b);
        for (int i = 0; i < LP_number[b]; i++) {
            printf("cons1:%d const2:%d block1:%d block2:%d sp_bMat:%d \n", LP_constraint1[b][i], LP_constraint2[b][i], LP_blockIndex1[b][i], LP_blockIndex2[b][i], LP_location_sparse_bMat[b][i]);
        }
    }
}

void Newton::display_sparse_bMat(FILE *fpout) {
    if (fpout == NULL) {
        return;
    }
    fprintf(fpout, "{");
    for (int index = 0; index < sparse_bMat.NonZeroCount; ++index) {
        int i = sparse_bMat.row_index[index];
        int j = sparse_bMat.column_index[index];
        dd_real value = sparse_bMat.sp_ele[index];
        int ii = ordering[i];
        int jj = ordering[j];
        fprintf(fpout, "val[%d,%d] = %e\n", ii, jj, value.x[0]);
    }
    fprintf(fpout, "}\n");
}

} // namespace sdpa
