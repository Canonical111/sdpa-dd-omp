/* -------------------------------------------------------------

This file is a component of SDPA-C
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

#define UseMETIS 0
#define PrintSparsity 0
#define OrderOnlyByMDO 1

#include <sdpa_chordal.h>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>

namespace sdpa {

Chordal::Chordal() { initialize(); }

Chordal::~Chordal() {
    //
}

void Chordal::initialize() {
    // condition of sparse computation
    // m_threshold < mDim,
    // b_threshold < nBlock,
    // aggregate_threshold >= aggrigated sparsity ratio
    // extend_threshold    >= extended sparsity ratio
    m_threshold = 100;
    b_threshold = 5;
    aggregate_threshold = 0.25;
    extend_threshold = 0.4;

#if 0 // DENSE computation for debugging
  m_threshold = 10000000;
  b_threshold = 1000000; 
  aggregate_threshold = 0.0; 
  extend_threshold = 0.0;
#endif
#if 0 // SPARSE computation for debugging
  m_threshold = 0;
  b_threshold = 0; 
  aggregate_threshold = 2.0; 
  extend_threshold = 2.0;
#endif

    /* indicates the used ordering method */
    /* 0: METIS 4.0.1 - nested dissection <--- not support */
    /* 1: Spooles 2.2 - mininum degree */
    /* 2: Spooles 2.2 - generalized nested dissection */
    /* 3: Spooles 2.2 - multisection */
    /* 4: Spooles 2.2 - better of 2 and 3 */
#if OrderOnlyByMDO
    Method[0] = 0;
    Method[1] = 1;
    Method[2] = 0;
    Method[3] = 0;
    Method[4] = 0;
#else
    Method[0] = 0;
    Method[1] = 1;
    Method[2] = 1;
    Method[3] = 1;
    Method[4] = 1;
#endif
    best = -1;
}

void Chordal::terminate() {
    if (Method[0]) {
        rError("no support for METIS");
    }
    if (Method[1] > 1) {
        IV_free(newToOldIV_MMD);
        IVL_free(symbfacIVL_MMD);
    }
    if (Method[2] > 1) {
        IV_free(newToOldIV_ND);
        IVL_free(symbfacIVL_ND);
    }
    if (Method[3] > 1) {
        IV_free(newToOldIV_MS);
        IVL_free(symbfacIVL_MS);
    }
    if (Method[4] > 1) {
        IV_free(newToOldIV_NDMS);
        IVL_free(symbfacIVL_NDMS);
    }
}

// marge array1 to array2
void Chordal::margeArray(int na1, int *array1, int na2, int *array2) {

    int ptr = na1 + na2 - 1;
    int ptr1 = na1 - 1;
    int ptr2 = na2 - 1;
    int idx1, idx2;

    while ((ptr1 >= 0) || (ptr2 >= 0)) {

        if (ptr1 >= 0) {
            idx1 = array1[ptr1];
        } else {
            idx1 = -1;
        }
        if (ptr2 >= 0) {
            idx2 = array2[ptr2];
        } else {
            idx2 = -1;
        }
        if (idx1 > idx2) {
            array2[ptr] = idx1;
            ptr1--;
        } else {
            array2[ptr] = idx2;
            ptr2--;
        }
        ptr--;
    }

    // error check
    if (ptr != -1) {
        rMessage("Chordal::margeArray:: program bug");
    }
}

// make aggrigate sparsity pattern
void Chordal::makeGraph(InputData &inputData, int m) {

    int i, j, k, l;
    int SDP_nBlock = inputData.SDP_nBlock;
    int SOCP_nBlock = inputData.SOCP_nBlock;
    int LP_nBlock = inputData.LP_nBlock;

    int *counter;
    counter = new int[m];
    for (int i = 0; i < m; i++) {
        counter[i] = 0;
    }

    // count maximum mumber of index
    for (l = 0; l < SDP_nBlock; l++) {
        int SDP_nConstraint = inputData.SDP_nConstraint[l];
        for (k = 0; k < SDP_nConstraint; k++) {
            i = inputData.SDP_constraint[l][k];
            counter[i] += SDP_nConstraint;
        }
    }
    for (l = 0; l < SOCP_nBlock; l++) {
        int SOCP_nConstraint = inputData.SOCP_nConstraint[l];
        for (k = 0; k < SOCP_nConstraint; k++) {
            i = inputData.SOCP_constraint[l][k];
            counter[i] += SOCP_nConstraint;
        }
    }
    for (l = 0; l < LP_nBlock; l++) {
        int LP_nConstraint = inputData.LP_nConstraint[l];
        for (k = 0; k < LP_nConstraint; k++) {
            i = inputData.LP_constraint[l][k];
            counter[i] += LP_nConstraint;
        }
    }

    // allocate temporaly workspace
    int **tmp;
    tmp = new int *[m];
    for (i = 0; i < m; i++) {
        tmp[i] = new int[counter[i]];
    }

    // merge index
    for (int i = 0; i < m; i++) {
        counter[i] = 0;
    }
    // marge index of for SDP
    for (l = 0; l < SDP_nBlock; l++) {
        for (k = 0; k < inputData.SDP_nConstraint[l]; k++) {
            i = inputData.SDP_constraint[l][k];
            margeArray(inputData.SDP_nConstraint[l], inputData.SDP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.SDP_nConstraint[l];
        }
    }
    // marge index of for SOCP
    for (l = 0; l < SOCP_nBlock; l++) {
        for (k = 0; k < inputData.SOCP_nConstraint[l]; k++) {
            i = inputData.SOCP_constraint[l][k];
            margeArray(inputData.SOCP_nConstraint[l], inputData.SOCP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.SOCP_nConstraint[l];
        }
    }
    // marge index of for LP
    for (l = 0; l < LP_nBlock; l++) {
        for (k = 0; k < inputData.LP_nConstraint[l]; k++) {
            i = inputData.LP_constraint[l][k];
            margeArray(inputData.LP_nConstraint[l], inputData.LP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.LP_nConstraint[l];
        }
    }

    // construct adjacency list of SPOOLES
    IVL_init1(adjIVL, IVL_CHUNKED, m);
    int isize, previous;
    int *ivec;
    ivec = new int[m];
    for (i = 0; i < m; i++) {
        isize = 0;
        previous = -1;
        for (j = 0; j < counter[i]; j++) {
            if (tmp[i][j] != previous) {
                ivec[isize] = tmp[i][j];
                previous = ivec[isize];
                isize++;
            }
        }
        IVL_setList(adjIVL, i, isize, ivec);
    }

    // constract graph of SPOOLES
    Graph_init2(graph, 0, m, 0, IVL_tsize(adjIVL), m, IVL_tsize(adjIVL), adjIVL, NULL, NULL);

    delete[] counter;
    for (int i = 0; i < m; i++) {
        delete[] tmp[i];
    }
    delete[] tmp;
    delete[] ivec;
}

int Chordal::countNonZero(int m, IVL *symbfacIVL) {
    int nonzeros = 0;
    bool *bnode;

    // count non-zero element
    rNewCheck();
    bnode = new bool[m];
    if (bnode == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }
    for (int i = 0; i < m; i++) {
        bnode[i] = false;
    }

    int nClique = IVL_nlist(symbfacIVL);
    int psize;
    int *pivec;
    for (int l = nClique - 1; l >= 0; l--) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (int i = 0; i < psize; i++) {
            int ii = pivec[i];
            if (bnode[ii] == false) {
                nonzeros += psize - i;
                bnode[ii] = true;
            }
        }
    }

    delete[] bnode;
    // Representability guard. Every caller returns `nonzeros * 2 - m`, so the DOUBLED value is
    // what has to fit in the int they return, not the count itself. Validate in 64 bits and
    // then narrow: the multiplication would otherwise overflow silently, and a negative or
    // wrapped fill count would be compared against the gate-4 cutoff as if it were small,
    // routing a problem too large for either representation straight into the sparse path.
    if ((long long)nonzeros * 2LL - (long long)m > (long long)INT_MAX ||
        (long long)nonzeros > (long long)INT_MAX) {
        rError("Chordal::countNonZero: the symbolic factor has " << nonzeros
               << " stored entries, so 2*nnz-m exceeds INT_MAX=" << INT_MAX
               << ". This problem is too large for either bMat representation in this build.");
    }
    return nonzeros;
}

int Chordal::Spooles_MMD(int m) {
    int seed = 0, msglvl = 0;
    FILE *fp = NULL;

    //  rMessage("orderViaMMD:start");
    etree = orderViaMMD(graph, seed, msglvl, fp);
    //  rMessage("orderViaMMD:end");
    newToOldIV_MMD = ETree_newToOldVtxPerm(etree);
    symbfacIVL_MMD = SymbFac_initFromGraph(etree, graph);
    //  IVL_writeForHumanEye(symbfacIVL_MMD,stdout);

    int nonzeros = countNonZero(m, symbfacIVL_MMD);

    return nonzeros * 2 - m;
}

int Chordal::Spooles_NDMS(int m) {
    int seed = 0, msglvl = 0;
    FILE *fp = NULL;

    int maxdomainsize = m / 16 + 1;
    int maxzeros = m / 10 + 1;
    int maxsize = 64;

    //  rMessage("orderViaMMD:start");
    etree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros, maxsize, seed, msglvl, fp);
    //  rMessage("orderViaMMD:end");
    newToOldIV_NDMS = ETree_newToOldVtxPerm(etree);
    symbfacIVL_NDMS = SymbFac_initFromGraph(etree, graph);
    //  IVL_writeForHumanEye(symbfacIVL_NDMS,stdout);

    int nonzeros = countNonZero(m, symbfacIVL_NDMS);

    return nonzeros * 2 - m;
}

int Chordal::Spooles_ND(int m) {
    int seed = 0, msglvl = 0;
    FILE *fp = NULL;
    bool *bnode;

    int maxdomainsize = m / 16 + 1;

    //  rMessage("orderViaMMD:start");
    etree = orderViaND(graph, maxdomainsize, seed, msglvl, fp);
    //  rMessage("orderViaMMD:end");
    newToOldIV_ND = ETree_newToOldVtxPerm(etree);
    symbfacIVL_ND = SymbFac_initFromGraph(etree, graph);
    //  IVL_writeForHumanEye(symbfacIVL_ND,stdout);

    int nonzeros = countNonZero(m, symbfacIVL_ND);

    return nonzeros * 2 - m;
}

int Chordal::Spooles_MS(int m) {
    int seed = 0, msglvl = 0;
    FILE *fp = NULL;
    bool *bnode;

    int maxdomainsize = m / 16 + 1;

    //  rMessage("orderViaMMD:start");
    etree = orderViaMS(graph, maxdomainsize, seed, msglvl, fp);
    //  rMessage("orderViaMMD:end");
    newToOldIV_MS = ETree_newToOldVtxPerm(etree);
    symbfacIVL_MS = SymbFac_initFromGraph(etree, graph);
    //  IVL_writeForHumanEye(symbfacIVL_MS,stdout);

    int nonzeros = countNonZero(m, symbfacIVL_MS);

    return nonzeros * 2 - m;
}

/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-23, revised 2026-08-24: adds SDPA_BMAT_MODE,
   an explicit and strictly-parsed route selector, and decouples gate 2 from gate 3's constant.
   The fill-derived policy was opt-in when this was first written; it is now the DEFAULT, with
   SDPA_BMAT_MODE=legacy restoring the previous chooser -- see the block below. Ported from
   sdpa-gmp-omp; see review/GATE3-DECISION-RULE-FINAL-REPORT.md and
   review/DD-PORT-PLAN-2026-08-23.md in the recipe repo.

   THE DEFAULT IS THE FILL-DERIVED POLICY (2026-08-24). Gate 2 keeps a frozen cutoff of 0.5*m
   through a constant of its own; gate 3 is demoted from an independent 0.25 policy to the early
   implication of gate 4's own threshold; gate 4 decides the rest. SDPA_BMAT_MODE=legacy restores
   the released chooser exactly.

   WHY GATE 3's 0.25 IS WRONG, AND WHY DEMOTING IT IS SOUND. Under the legacy rule gate 3 rejects
   on aggregate density > 0.25 while gate 4 rejects on ordered fill > 0.40, so gate 3 is the
   STRICTER test and a problem with aggregate in (0.25, 0.40] never reaches gate 4 at all -- the
   cheap pre-screen overrules the real one. Demoting gate 3 to gate 4's own constant is sound
   because symbolic factorisation only ADDS entries: aggregate > F proves fill > F, so the early
   exit can only skip work, never change gate 4's verdict. F = 0.40 is therefore the LARGEST
   sound cutoff for gate 3, and 0.25 was conservatism with nothing behind it.

   ROUTE DATA TRANSFERS ACROSS PRECISIONS because these gates read only structure -- the same
   census governs dd and gmp:
     SDPLIB     92 instances:   0 change route (84 decided at gate 1, 6 at gate 2, and the one
                                that reaches gate 3, truss5, has aggregate = fill = 1.0 and goes
                                dense under BOTH rules)
     bootstrap 221 instances: 167 change DENSE -> SPARSE, across 7 structures. Every one of the
                                167 has aggregate in 0.257-0.362 and ordered fill in 0.259-0.367:
                                above the legacy 0.25, and NOT ONE above 0.40. The real test says
                                sparse on all of them and the pre-screen was overriding it.

   MEASURED ON dd, all seven switch structures, 24 threads, medians of 3
   (bench/dd-port3-2026-08-24/dd_fill_seven_structures.tsv):
     m= 2439  1.92 s/116 MB -> 0.61 s/ 50 MB   3.17x, 2.30x less memory
     m= 4489  9.22 s/392 MB -> 3.59 s/213 MB   2.57x, 1.84x
     m= 5278 14.50 s/519 MB -> 3.88 s/234 MB   3.74x, 2.22x
     m= 6067 22.97 s/671 MB -> 4.50 s/277 MB   5.11x, 2.42x
     m= 8359 66.20 s/1.3 GB -> 15.13 s/587 MB  4.38x, 2.20x
     m=10614  150.5 s/2.0 GB -> 22.36 s/831 MB 6.73x, 2.45x
     m=11227  175.8 s/2.4 GB -> 47.94 s/1.2 GB 3.67x, 1.99x
   Seven of seven, no reversal.

   AND THE SIBLING FORK ALREADY DID THIS, on 2026-08-18 after five review rounds, with evidence dd
   cannot generate for itself: sdpa-gmp-omp promoted the same policy on the SAME seven structures,
   measured AT FULL CONVERGENCE (m=6067: dense 2708 s against sparse 310 s), on two architectures,
   at 256 and 512 bits, at nine thread points, with no reversal -- and at T=1, where sparse still
   won 2.60-3.06x, proving the 0.25 was already mis-routing in the SERIAL regime it was calibrated
   for. Where the dense arm converged, phases and iteration counts were EQUAL per pair and the
   solution vectors agreed to 6.3e-38. dd cannot reproduce that half: these instances do not
   converge at double-double precision under any tolerance tested, so dd's own timings are
   per-iteration only. See review/GATE3-DECISION-RULE-FINAL-REPORT.md in the recipe repo.

   NOT "SPARSE ALWAYS WINS", and the counterexample is in our own corpus: SDPLIB truss5 has an
   ordered fill of 1.0 -- a completely dense factor -- and dense IS faster there. Both rules route
   it dense. That is why the rule tests fill rather than assuming an answer.

   WHAT THE CHANGE COSTS. Selecting a different factorisation is RESULT-CHANGING in the last
   digits: the two routes follow different iterate trajectories. SDPA_BMAT_MODE=legacy restores
   the pre-promotion chooser exactly -- the released gate expressions, including gate 2 derived
   from m*sqrt(aggregate_threshold) and gate 3 at aggregate_threshold. Anyone reproducing a
   pre-promotion dd result should set it. `fill` is retained as an explicit synonym for the new
   default so scripts written during the opt-in phase keep working. */

namespace {

/* NOTE THE MISSING BMAT_AUTO. `auto` is still the default mode NAME, but its policy is now
   BMAT_FILL, so an enumerator by that name would be a trap: `mode == BMAT_AUTO` reads as "is this
   the default?" and would be false for every default run. The name still exists and is still the
   default; it simply is not its own policy any more. Same choice, and same reasoning, as the gmp
   fork. */
enum BMatMode { BMAT_DENSE, BMAT_SPARSE, BMAT_FILL, BMAT_LEGACY };

// Gate 2's own constant. Under `legacy` the cutoff stays coupled to aggregate_threshold via
// sqrt(), exactly as upstream wrote it; under the default it is this number instead, which has
// the SAME value (0.5) but is no longer moved by retuning gate 3. That coupling is why gate 3
// could not be retuned without silently changing gate 2.
const double BMAT_FILL_BLOCK_FRACTION = 0.5;

// Cap on the DENSE bMat allocation, in GB. Absent means no cap requested; present-but-empty is
// an error rather than "unset", because SDPA_BMAT_MAX_GB="$TYPO" is the case worth catching.
// dd-SPECIFIC BYTE ACCOUNTING: gmp estimates from limb counts at a runtime precision. A dd_real
// is a fixed pair of doubles, so the dense bMat is exactly m*m*sizeof(dd_real) -- no precision
// term, and nothing to look up at runtime.
double bmat_max_gb() {
    const char *e = getenv("SDPA_BMAT_MAX_GB");
    if (e == NULL) {
        return -1.0;
    }
    if (e[0] == '\0') {
        rError("SDPA_BMAT_MAX_GB is set but empty; unset it to request no cap");
    }
    errno = 0;
    char *end = NULL;
    const double v = strtod(e, &end);
    if (end == e || *end != '\0' || errno == ERANGE || !(v > 0.0) || !std::isfinite(v)) {
        rError("SDPA_BMAT_MAX_GB must be a positive finite number of GB (got \"" << e << "\")");
    }
    return v;
}

/* ADVISORY available memory, in GB, or -1 when it cannot be known. Finding F.

   This is a DIAGNOSTIC and never an admission decision, which is the whole design. A route
   chosen from what the machine happens to have free is a route that changes with the machine,
   and this fork's central claim is that the answer does not. So nothing branches on this
   number: it is printed alongside the requirement so that "needs 21 GB" can be read next to
   "about 6 GB available" instead of leaving the user to find that out from the OOM killer.

   Linux only, and deliberately only MemAvailable -- the kernel's own estimate of what a new
   allocation can get, which is the question being asked. MemFree is the wrong field (it
   excludes reclaimable page cache and reads far too low). On every other platform, and on a
   Linux without /proc mounted, it returns -1 and the caller says nothing rather than guessing.

   CAVEAT, printed with the number rather than buried here: /proc/meminfo reports the HOST. In a
   container with a memory cgroup limit below the host's free memory, this over-reports -- the
   cgroup files that would answer correctly differ between cgroup v1 and v2 and between
   runtimes, and an advisory line does not justify that. Treat it as an upper bound. */
// "0" means OFF. Until now this knob was tested for PRESENCE only, so the documented way to
// turn the route log off turned it on -- the identical defect that was found and fixed in
// SDPA_SPCHOL_LOG, still present in its sibling because the fix was applied to one call site
// rather than to the pattern. The runtime reference already documented `1` here; the code now
// agrees. (That reference was RUNTIME.md until 2026-08-25; it is now doc/technical.tex, sec. 9.)
bool bmat_log_wanted() {
    const char *e = getenv("SDPA_BMAT_LOG");
    return (e != NULL && e[0] != '\0' && strcmp(e, "0") != 0);
}

double bmat_available_gb() {
#ifdef __linux__
    FILE *fp = fopen("/proc/meminfo", "r");
    if (fp == NULL) {
        return -1.0;
    }
    char line[256];
    double gb = -1.0;
    while (fgets(line, (int)sizeof(line), fp) != NULL) {
        unsigned long long kb = 0ULL;
        // The field is documented in kB regardless of the kernel's page size.
        if (sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
            gb = (double)kb / (1024.0 * 1024.0);
            break;
        }
    }
    fclose(fp);
    return gb;
#else
    return -1.0;
#endif
}

// Enforced on EVERY route to a dense bMat -- gate 1, gate 2, gate 3, gate 4 and forced dense --
// rather than at one of them, because a cap only some paths consult is not a cap.
void bmat_dense_cap_check(int m, const char *why) {
    const double cap = bmat_max_gb();
    if (cap < 0.0) {
        // No cap requested. The advisory still runs when the caller asked to see the route
        // decision, because "which route, and could this machine have held the other one" is
        // one question, and answering only half of it is what sends people to the OOM killer.
        if (bmat_log_wanted()) {
            const double need = (double)m * (double)m * (double)sizeof(dd_real)
                                / (1024.0 * 1024.0 * 1024.0);
            const double avail = bmat_available_gb();
            if (avail >= 0.0) {
                rMessage("bmat: dense bMat for m=" << m << " needs " << need
                         << " GB; about " << avail
                         << " GB available (MemAvailable, HOST-wide -- an upper bound inside a"
                            " memory cgroup). Advisory only; no route depends on it");
            } else {
                rMessage("bmat: dense bMat for m=" << m << " needs " << need
                         << " GB; available memory not readable on this platform");
            }
        }
        return;
    }
    // m*m in 64 bits: the product overflows int for m > 46341, and this check exists precisely
    // for the large problems where that matters.
    const double bytes = (double)m * (double)m * (double)sizeof(dd_real);
    const double gb = bytes / (1024.0 * 1024.0 * 1024.0);
    if (gb > cap) {
        // The advisory belongs HERE above all: the message a user actually reads is the one
        // that stopped their job. Reporting the requirement without the machine's side of it
        // leaves "is 21 GB a lot?" unanswered at the only moment it matters.
        const double avail = bmat_available_gb();
        std::ostringstream note;
        if (avail >= 0.0) {
            note << " For reference this host reports about " << avail
                 << " GB available (MemAvailable; HOST-wide, so an upper bound inside a memory"
                    " cgroup). It is advisory: no route was chosen from it.";
        } else {
            note << " Available memory is not readable on this platform, so this message cannot"
                    " say whether the machine could have held it.";
        }
        rError("dense bMat for m=" << m << " needs " << gb << " GB (" << sizeof(dd_real)
               << " B/element), over the SDPA_BMAT_MAX_GB cap of " << cap
               << " GB; route chosen by " << why
               << ". Raise the cap, or use SDPA_BMAT_MODE=fill/sparse if the problem admits it."
               << note.str());
    }
}

/* TEST-ONLY, behind the same compile gate as the spchol hooks: perturb the value COMPARED
   against the aggregate count so the fill-policy invariant's error branch is executable in CI
   rather than merely source-audited. Without it, "the invariant never fired" is equally
   consistent with an invariant that CANNOT fire -- the same vacuity trap SDPA_SPCHOL_MUTATE
   closes for the factor oracle.

   Absent unless -DSDPA_SPCHOL_TEST_HOOKS, and a build without the gate REFUSES the variable
   rather than ignoring it: a knob that silently does nothing misleads whoever set it. Ported
   from gmp. See git log. */
static bool bmat_test_break_invariant() {
    const char *e = getenv("SDPA_BMAT_TEST_BREAK_INVARIANT");
    // "0" means "I am not asking for the hook", which every build can satisfy, so it is
    // accepted before the compile gate. This follows SDPA_SPCHOL_MUTATE's convention in this
    // fork rather than gmp's, which refuses even "0" in a release build -- being consistent
    // with dd's own neighbouring knob matters more here than matching gmp.
    if (e == NULL || e[0] == '\0' || strcmp(e, "0") == 0) {
        return false;
    }
#ifndef SDPA_SPCHOL_TEST_HOOKS
    rError("SDPA_BMAT_TEST_BREAK_INVARIANT is a test hook and this binary was not built with"
           " -DSDPA_SPCHOL_TEST_HOOKS, so the hook does not exist here");
    return false;
#else
    if (strcmp(e, "1") == 0) {
        return true;
    }
    rError("SDPA_BMAT_TEST_BREAK_INVARIANT must be 0 or 1 (got \"" << e << "\")");
    return false;
#endif
}

BMatMode bmat_mode() {
    const char *e = getenv("SDPA_BMAT_MODE");
    // Unset, "auto" and "fill" all mean the current default policy, which is the fill-derived
    // one. "fill" is kept as an explicit synonym so opt-in-era scripts keep working.
    if (e == NULL || e[0] == '\0' || strcmp(e, "auto") == 0 || strcmp(e, "fill") == 0) {
        return BMAT_FILL;
    }
    if (strcmp(e, "dense") == 0) {
        return BMAT_DENSE;
    }
    if (strcmp(e, "sparse") == 0) {
        return BMAT_SPARSE;
    }
    // The pre-promotion chooser, exactly.
    if (strcmp(e, "legacy") == 0) {
        return BMAT_LEGACY;
    }
    // Strict: a typo must not silently select the default, which is the very route a caller
    // setting this variable is usually trying to select against.
    rError("SDPA_BMAT_MODE must be auto, fill, dense, sparse or legacy (got \"" << e
           << "\"). auto is the current policy; legacy is the pre-2026-08-24 chooser;"
           << " fill is an explicit synonym for auto.");
    return BMAT_FILL;
}

} // namespace

void Chordal::ordering_bMat(int m, int nBlock, InputData &inputData, FILE *fpOut) {
    const BMatMode mode = bmat_mode();
    // Parsed HERE, at the top, not at the invariant it feeds: gates 1-3 return early, so a
    // dense-route problem would otherwise never validate the hook variable and a build without
    // the test gate would silently accept it. Same lesson as bmat_mode() above.
    const bool break_inv = bmat_test_break_invariant();
    const bool want_log = bmat_log_wanted();
    // Gate 2's cutoff: auto's expression verbatim, or fill's own decoupled constant.
    const double g2 = (mode == BMAT_LEGACY) ? (m * sqrt(aggregate_threshold))
                                            : (BMAT_FILL_BLOCK_FRACTION * m);
    // Gate 3's cutoff: the legacy 0.25, or the default's demotion to gate 4's own F.
    const double g3 = (mode == BMAT_LEGACY) ? (aggregate_threshold * m * (double)m)
                                            : (extend_threshold * m * (double)m);

    if (mode == BMAT_DENSE) {
        if (want_log)
            rMessage("bmat: mode=dense (forced) -> DENSE");
        bmat_dense_cap_check(m, "forced dense");
        best = -1;
        return;
    }
    if (mode != BMAT_SPARSE) {
        if ((m <= m_threshold) || (nBlock <= b_threshold)) {
            if (want_log)
                rMessage("bmat: gate1 m=" << m << " nBlock=" << nBlock << " -> DENSE");
            bmat_dense_cap_check(m, "gate 1");
            best = -1;
            return;
        }
        for (int b = 0; b < inputData.SDP_nBlock; b++) {
            if (inputData.SDP_nConstraint[b] > g2) {
                if (want_log)
                    rMessage("bmat: gate2 SDP block " << b << " -> DENSE");
                bmat_dense_cap_check(m, "gate 2");
                best = -1;
                return;
            }
        }
        for (int b = 0; b < inputData.SOCP_nBlock; b++) {
            if (inputData.SOCP_nConstraint[b] > g2) {
                if (want_log)
                    rMessage("bmat: gate2 SOCP block " << b << " -> DENSE");
                bmat_dense_cap_check(m, "gate 2");
                best = -1;
                return;
            }
        }
        for (int b = 0; b < inputData.LP_nBlock; b++) {
            if (inputData.LP_nConstraint[b] > g2) {
                if (want_log)
                    rMessage("bmat: gate2 LP block " << b << " -> DENSE");
                bmat_dense_cap_check(m, "gate 2");
                best = -1;
                return;
            }
        }
    }

    adjIVL = IVL_new();
    graph = Graph_new();

    makeGraph(inputData, m);

    // Capture the aggregate count NOW, while adjIVL is still live. Graph_init2 hands adjIVL to
    // the graph, so the Graph_free below frees it -- and the fill-vs-aggregate invariant is
    // checked AFTER that free. Reading it there was a heap-use-after-free, which the sanitizer
    // job caught on the first push. A local costs nothing and cannot dangle.
    const int aggregate_nnz = IVL_tsize(adjIVL);

    if (mode != BMAT_SPARSE && aggregate_nnz > g3) {
        if (want_log)
            rMessage("bmat: gate3 aggregate=" << aggregate_nnz << " cutoff=" << g3
                                              << (mode == BMAT_LEGACY ? " [legacy policy]" : "")
                                              << " -> DENSE");
        bmat_dense_cap_check(m, "gate 3");
        best = -1;
        Graph_free(graph);
        return;
    }
#if PrintSparsity
    /* print sparsity information */
    printf("dense matrix               :\t\t\t%14d elements\n", m * m);
    fprintf(fpOut, "dense matrix               :\t\t\t%14d elements\n", m * m);
    printf("aggregate sparsity pattern :\t\t\t%14d elements\n", IVL_tsize(adjIVL));
    fprintf(fpOut, "aggregate sparsity pattern :\t\t\t%14d elements\n", IVL_tsize(adjIVL));
#endif

    /* Uses METIS */
    if (Method[0]) {
        rError("no support for METIS");
    }

    /* Uses Spooles */
    if (Method[1]) { /* Spooles 2.2 - minimum degree */
        Method[1] = Spooles_MMD(m);
        ETree_free(etree);
#if PrintSparsity
        printf("\tSpooles2.2 (minimum degree)\t\t%14d elements\n", Method[1]);
        fprintf(fpOut, "\tSpooles2.2 (minimum degree)\t\t%14d elements\n", Method[1]);
#endif
    }
    if (Method[2]) { /* Spooles 2.2 - generalized nested dissection */
        Method[2] = Spooles_ND(m);
        ETree_free(etree);
#if PrintSparsity
        printf("\tSpooles2.2 (generalized nested dissection)%12d elements\n", Method[2]);
        fprintf(fpOut, "\tSpooles2.2 (generalized nested dissection)%12d elements\n", Method[2]);
#endif
    }
    if (Method[3]) { /* Spooles 2.2 - multisection */
        Method[3] = Spooles_MS(m);
        ETree_free(etree);
#if PrintSparsity
        printf("\tSpooles2.2 (multisection)\t\t%14d elements\n", Method[3]);
        fprintf(fpOut, "\tSpooles2.2 (multisection)\t\t%14d elements\n", Method[3]);
#endif
    }
    if (Method[4]) { /* Spooles 2.2 - best between nested
                        dissection and multisection */
        Method[4] = Spooles_NDMS(m);
        ETree_free(etree);
#if PrintSparsity
        printf("\tSpooles2.2 (best of ND and MS)\t\t%14d elements\n", Method[4]);
        fprintf(fpOut, "\tSpooles2.2 (best of ND and MS)\t\t%14d elements\n", Method[4]);
#endif
    }
    /* Select the best ordering */

    Graph_free(graph);

    best = Best_Ordering(Method);

    // The invariant the fill policy's early exit rests on: symbolic factorisation only ADDS
    // entries, so the ordered fill can never be below the aggregate pattern in the same counting
    // convention. A violation means the counts have left those units -- a build defect, not a
    // data condition -- and it must not be papered over, because the gate-3 skip would then be
    // unsound.
    const int fill_for_invariant = break_inv ? (aggregate_nnz - 1) : Method[best];
    if (fill_for_invariant < aggregate_nnz) {
        // The diagnostic prints the value that was COMPARED, and says so when that value was
        // injected -- otherwise it would report a real fill count as "below" an aggregate it
        // actually exceeds, which is a lying diagnostic (caught in gmp's review round).
        rError("Chordal::ordering_bMat: ordered fill " << fill_for_invariant
               << (break_inv ? " (TEST INJECTION ACTIVE: real count was " : " (real count ")
               << Method[best] << ")"
               << " is below the aggregate pattern " << aggregate_nnz
               << "; the fill-policy early exit's invariant is broken. This is a solver build"
               << " defect, not a data condition: run with SDPA_BMAT_MODE=dense, which does not"
               << " consult this count, and report it.");
    }
    if (mode == BMAT_SPARSE) {
        if (want_log)
            rMessage("bmat: mode=sparse (forced) fill=" << Method[best] << " -> SPARSE");
        return;
    }
    if (Method[best] > extend_threshold * m * (double)m) {
        if (want_log)
            rMessage("bmat: gate4 fill=" << Method[best] << " cutoff="
                                         << extend_threshold * m * (double)m << " -> DENSE");
        bmat_dense_cap_check(m, "gate 4");
        best = -1;
    } else if (want_log) {
        rMessage("bmat: gate4 fill=" << Method[best] << " -> SPARSE");
    }
}

int Chordal::Best_Ordering(int *Method)
/************************************************************************
        Determine the best ordering.
************************************************************************/
{
    int i, best;

    for (i = 0; Method[i] == 0; i++)
        ;
    best = i++;
    while (i < 5) {
        for (; i < 5; i++) {
            if (Method[i] != 0)
                break;
        }
        if (i < 5) {
            if (Method[i] < Method[best])
                best = i;
            i++;
        }
    }
    return best;
}

} // namespace sdpa
