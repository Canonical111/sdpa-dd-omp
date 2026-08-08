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
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-03: fatal parameter-file error exits non-zero; solver status propagated to exit code. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-05: the residual is built once into currentRes; initRes keeps only the scalars the solve actually reads, so the second full copy is never allocated. See git log. */
/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-06: numerical solve failure is tracked, suppresses the solution section, is diagnosed on stderr and exits non-zero; exit-status policy documented in this file; fpu_fix_end() moved before the return. See git log. */

#ifndef _MAIN_
#define _MAIN_
#endif

#define LengthOfBuffer 1024
static double KAPPA = 1.2;

#include <sdpa_io.h>

#if UseMETIS
#ifdef __cplusplus
extern "C" {
#endif
#include "metis.h"
#ifdef __cplusplus
}
#endif
#endif

namespace sdpa {

#define SDPA_PROGRAM_NAME "sdpa_dd"

// ---------------------------------------------------------------------------
// EXIT-STATUS POLICY  (2026-08-06, review1 blocker 3)
//
// The exit status answers exactly one question:
//
//     "Did this run compute an answer at all?"
//
// It deliberately does NOT answer "is the answer optimal?", and -- this is the
// part a first attempt at this fix got wrong -- it does NOT answer "did every
// factorisation succeed?" either.
//
//   0  SDPA_EXIT_OK
//        The solver completed at least one iteration and printed the iterate it
//        ended on.  This covers every valid mathematical outcome, including the
//        many that are not pdOPT:
//          - primal and/or dual infeasibility detected   (pdINF, pINF_dFEAS, ...)
//          - unboundedness detected                      (pUNBD, dUNBD)
//          - the iteration limit was reached             (noINFO, pdFEAS, ...)
//          - the "objValPrimal < objValDual" stopping criterion fired
//          - the step length collapsed while the iterate was still a valid
//            interior point ("Step length is too small." / "cannot move")
//        Judge the quality of such an answer from phase.value, the two
//        objectives and the gap.  Not from the exit status.
//
//        MEASURED, and this is why the policy is shaped this way.  Two tempting
//        rules were tried against the ten-problem regression set at DD precision
//        and the measurements killed both:
//
//        (a) "nonzero unless phase == pdOPT" -- ZERO of the ten problems reach
//            pdOPT (noINFO x4, pFEAS x4, pdINF x2), yet control1 agrees with its
//            dual objective to 17 digits and hinf1 reaches its published optimum
//            to seven.  This rule fails all ten.
//
//        (b) "nonzero whenever a Cholesky factorisation failed" -- the literal
//            reading of the review -- fails FOUR of the ten:
//              control1  iteration 71, mu = 3.9e-28, objP = objD = 1.78e+01;
//                        X/Z lost positive definiteness -- because it converged
//              theta1    iteration 50, mu = 1.4e-30, objP = objD = 23.0, which
//                        is the published SDPLIB optimum; same mechanism
//              arch0     iteration 72, Schur complement singular at pivot 90
//              qap5      iteration 21, Schur complement singular at pivot 86
//            Driving an interior-point method to the limit of its arithmetic is
//            precisely what makes X, Z and the Schur complement numerically
//            singular.  A failed factorisation at the END of a solve is how this
//            solver normally STOPS; it is not a sign that anything went wrong.
//
//        Since review2 finding 1 such an end-of-solve failure is no longer
//        folded into exit 0: it exits 3, with the last VALID iterate printed.
//        The measurements above are why it is exit 3 and not exit 2.
//
//   3  SDPA_EXIT_PARTIAL
//        A factorisation failed after k >= 1 completed iterations and a valid
//        iterate survives: either the Schur complement failed before X/Z were
//        touched, or the failed X/Z update was rolled back and refactored
//        (sdpa_dataset.cpp).  The output carries "solveStatus = PARTIAL",
//        "failureReason" and "failureIteration = k", followed by the last
//        valid iterate.  At this fork's default parameters the four cases
//        from measurement (b) -- control1, theta1, arch0, qap5 -- end here.
//
//   2  SDPA_EXIT_NUMERICAL_FAILURE
//        Nothing valid survives to print: a factorisation failed before a
//        single iteration completed (pIteration == 0), the supplied initial
//        point was not positive definite, or a failed X/Z update could not be
//        rolled back (the restored point does not refactor).  In the
//        zero-iteration case what upstream printed was the untouched starting
//        point -- xVec all zeros, xMat = zMat = lambdaStar*I, objValPrimal =
//        -0.0 -- dressed in the
//        upstream printed here was the untouched starting point -- xVec all
//        zeros, xMat = zMat = lambdaStar*I, objValPrimal = -0.0 -- dressed in the
//        normal solution format, with a phase, and exit 0.  That is the
//        wrong-answer bug this policy exists to kill.  Here the solution section
//        is SUPPRESSED, the output file carries "solveStatus = FAILURE", and the
//        exit status is nonzero.
//
//   1  SDPA_EXIT_INPUT
//        Unusable input or usage: unreadable file, malformed data, bad options.
//        The pre-existing rError() / exit(EXIT_FAILURE) paths, unchanged here;
//        named so that the whole policy is readable in one place.
//
//        2 and 1 are both nonzero, so a harness testing "rc != 0" is unaffected
//        by the split; the split only lets a harness that cares tell "your file
//        is broken" from "the arithmetic broke".
//
// INDEPENDENTLY of the exit status: whenever a factorisation fails, at ANY
// iteration, a line naming the stage and the iteration is written to STDERR.
// Before this change stderr was empty in every such run -- rMessage and rError
// write to cout -- so a harness watching stderr could not see a numerical
// breakdown at all.  stderr is now the channel for "something went numerically
// wrong"; the exit status says whether, and how much of, an answer exists.
// Runs that end on a mathematical stopping condition write nothing extra to
// the output file or to stdout, so their output stays byte-for-byte identical;
// PARTIAL runs additionally carry the solveStatus block described above.
// ---------------------------------------------------------------------------
enum SolveResult { SDPA_EXIT_OK = 0, SDPA_EXIT_INPUT = 1, SDPA_EXIT_NUMERICAL_FAILURE = 2, SDPA_EXIT_PARTIAL = 3 };

// Report a solve failure on stderr and, per severity, in the output file.
// One function so the wording a harness greps for cannot drift between sites.
/* Three failure severities, decided by what state survives (review2 finding 1):
     FATAL      nothing was computed (iteration 0)          -> FAILURE, no solution, exit 2
     CORRUPTED  a failed X/Z update could not be rolled     -> FAILURE, no solution, exit 2
                back: the restored point does not refactor,
                so no valid in-memory iterate survives
     PARTIAL    the Schur factorisation failed AFTER k good -> PARTIAL, last valid
                iterations (currentPt untouched), or the       iterate printed, exit 3
                failed X/Z update was rolled back and
                refactored (sdpa_dataset.cpp)
   Before this, every non-zero-iteration failure printed the ordinary final
   section and exited 0 -- an updated-X/Z breakdown was indistinguishable from
   success, and bench_v2.sh recorded it as an ok row. */
enum FailureKind { FAIL_FATAL, FAIL_CORRUPTED, FAIL_PARTIAL };
static void reportSolveFailure(FILE *fpOut, const char *outFile, const char *reason, int iteration, enum FailureKind kind) {
    fprintf(stderr, "%s: %s :: %s", SDPA_PROGRAM_NAME, kind == FAIL_PARTIAL ? "PARTIAL RESULT" : "SOLVE FAILED", reason);
    if (iteration >= 0) {
        fprintf(stderr, " at iteration %d", iteration);
    }
    fprintf(stderr, "\n");
    if (kind == FAIL_PARTIAL) {
        fprintf(stderr, "%s: %d iteration(s) completed before the failure; the LAST VALID iterate is reported below and labelled solveStatus = PARTIAL. Judge it from phase.value and the two objectives; exit status is 3, not 0.\n", SDPA_PROGRAM_NAME, iteration);
        if (fpOut) {
            fprintf(fpOut, "\nsolveStatus = PARTIAL\n");
            fprintf(fpOut, "failureReason = %s\n", reason);
            fprintf(fpOut, "failureIteration = %d\n", iteration);
            fprintf(fpOut, "The section below is the last iterate BEFORE the failure; it is a valid interior point but not a certified solution.\n");
        }
        return;
    }
    if (kind == FAIL_CORRUPTED) {
        fprintf(stderr, "%s: the update had already replaced X and Z when the factorisation failed, so the in-memory iterate is not positive definite; %s contains no solution section.\n", SDPA_PROGRAM_NAME, outFile ? outFile : "(no output file)");
    } else {
        fprintf(stderr, "%s: no iteration completed, so no solution was computed; %s contains no solution section.\n", SDPA_PROGRAM_NAME, outFile ? outFile : "(no output file)");
    }
    if (fpOut) {
        fprintf(fpOut, "\nsolveStatus = FAILURE\n");
        fprintf(fpOut, "failureReason = %s\n", reason);
        if (iteration >= 0) {
            fprintf(fpOut, "failureIteration = %d\n", iteration);
        }
        if (kind == FAIL_CORRUPTED) {
            fprintf(fpOut, "No solution is presented: the failed update left X or Z outside the positive definite cone, so the in-memory iterate is invalid (see sdpa_dataset.h).\n");
        } else {
            fprintf(fpOut, "No solution is presented: no iteration completed, so the last iterate is the untouched initial point.\n");
        }
    }
}

int pinpal(char *dataFile, char *initFile, char *outFile, char *paraFile, bool isInitFile, bool isInitSparse, bool isDataSparse, bool isParameter, Parameter::parameterType parameterType, FILE *Display) {

    TimeStart(TOTAL_TIME_START1);
    TimeStart(FILE_READ_START1);
    ComputeTime com;

    FILE *fpData = NULL;
    FILE *fpOut = NULL;

    if ((fpOut = fopen(outFile, "w")) == NULL) {
        rError("Cannot open out file " << outFile);
    }
    Parameter param;
    param.setDefaultParameter(parameterType);
    if (isParameter) {
        FILE *fpParameter = NULL;
        if ((fpParameter = fopen(paraFile, "r")) == NULL) {
            fprintf(Display, "Cannot open parameter file %s \n", paraFile);
            exit(EXIT_FAILURE);
        } else {
            param.readFile(fpParameter);
            fclose(fpParameter);
        }
    }
    // param.display(Display);

    if ((fpData = fopen(dataFile, "r")) == NULL) {
        rError("Cannot open data file " << dataFile);
    }
    char titleAndComment[LengthOfBuffer];
    int m;
    time_t ltime;
    time(&ltime);
    fprintf(fpOut, "SDPA-DD start at %s", ctime(&ltime));
    IO::read(fpData, fpOut, m, titleAndComment);
    fprintf(fpOut, "data      is %s\n", dataFile);
    if (paraFile) {
        fprintf(fpOut, "parameter is %s\n", paraFile);
    }
    if (initFile) {
        fprintf(fpOut, "initial   is %s\n", initFile);
    }
    fprintf(fpOut, "out       is %s\n", outFile);

#if 1 // 2007/11/28 nakata    for multi LP block

    int SDP_nBlock, SOCP_nBlock, LP_nBlock, nBlock;
    IO::read(fpData, nBlock);

    int *blockStruct = NULL;
    int *blockType = NULL;
    int *blockNumber = NULL;
    int *SDP_blockStruct = NULL;
    int *SOCP_blockStruct = NULL;
    blockStruct = new int[nBlock];
    if (blockStruct == NULL) {
        rError("Memory exhausted about blockStruct");
    }
    blockType = new int[nBlock];
    if (blockType == NULL) {
        rError("Memory exhausted about blockType");
    }
    blockNumber = new int[nBlock];
    if (blockNumber == NULL) {
        rError("Memory exhausted about blockNumber");
    }
    IO::read(fpData, nBlock, blockStruct);

    SDP_nBlock = 0;
    SOCP_nBlock = 0;
    LP_nBlock = 0;
    for (int i = 0; i < nBlock; i++) {
        if (blockStruct[i] >= 2) {
            blockType[i] = 1;
            blockNumber[i] = SDP_nBlock;
            SDP_nBlock++;
        } else if (blockStruct[i] < 0) {
            blockType[i] = 3;
            blockStruct[i] = -blockStruct[i];
            blockNumber[i] = LP_nBlock;
            LP_nBlock += blockStruct[i];
        } else if (blockStruct[i] == 1) {
            blockType[i] = 3;
            blockNumber[i] = LP_nBlock;
            LP_nBlock += blockStruct[i];
        } else {
            rError("block struct");
        }
    }

    SDP_blockStruct = new int[SDP_nBlock];
    if (SDP_blockStruct == NULL) {
        rError("Memory exhausted about SDP blockStruct");
    }
    SOCP_blockStruct = new int[SOCP_nBlock];
    if (SOCP_blockStruct == NULL) {
        rError("Memory exhausted about SOCP blockStruct");
    }

    SDP_nBlock = 0;
    for (int i = 0; i < nBlock; i++) {
        if (blockType[i] == 1) {
            SDP_blockStruct[SDP_nBlock] = blockStruct[i];
            SDP_nBlock++;
        }
    }

    InputData inputData;
    //  rMessage("read input data: start");
    IO::read(fpData, m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock, nBlock, blockStruct, blockType, blockNumber, inputData, isDataSparse);
    //  rMessage("read input data: end");
    inputData.initialize_index(SDP_nBlock, SOCP_nBlock, LP_nBlock, com);
#else

    int SDP_nBlock, SOCP_nBlock, LP_nBlock;
    IO::read(fpData, SDP_nBlock, SOCP_nBlock, LP_nBlock);
    int *SDP_blockStruct = NULL;
    int *SOCP_blockStruct = NULL;
    SDP_blockStruct = new int[SDP_nBlock];
    if (SDP_blockStruct == NULL) {
        rError("Memory exhausted about SDP blockStruct");
    }
    SOCP_blockStruct = new int[SOCP_nBlock];
    if (SOCP_blockStruct == NULL) {
        rError("Memory exhausted about SOCP blockStruct");
    }
    IO::read(fpData, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

    for (int i = 0; i < SDP_nBlock - 1; i++) {
        if (SDP_blockStruct[i] < 0) {
            rError("LP block must be in last block");
        }
    }
    // muriyari nyuuryoku saseru
    if (SDP_blockStruct[SDP_nBlock - 1] < 0) {
        LP_nBlock = -SDP_blockStruct[SDP_nBlock - 1];
        SDP_nBlock--;
    }

    InputData inputData;
    IO::read(fpData, m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock, inputData, isDataSparse);
    inputData.initialize_index(SDP_nBlock, SOCP_nBlock, LP_nBlock, com);

#endif

    fclose(fpData);

#if 0
  inputData.display();
#endif

#if 1
    TimeStart(FILE_CHANGE_START1);
    // if possible , change C and A to Dense
    inputData.C.changeToDense();
    for (int k = 0; k < m; ++k) {
        inputData.A[k].changeToDense();
    }
    TimeEnd(FILE_CHANGE_END1);
    com.FileChange += TimeCal(FILE_CHANGE_START1, FILE_CHANGE_END1);
#endif

    // rMessage("C = ");
    // inputData.C.display(Display);
    // for (int k=0; k<m; ++k) {
    //   rMessage("A["<<k<<"] = ");
    //   inputData.A[k].display(Display);
    //   }

    // the end of initialization of C and A

    Newton newton(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    int nBlock2 = SDP_nBlock + SOCP_nBlock + LP_nBlock;
    // 2008/03/12 kazuhide nakata
    Chordal chordal;
    // rMessage("ordering bMat: start");
    chordal.ordering_bMat(m, nBlock2, inputData, fpOut);
    // rMessage("ordering bMat: end");
    newton.initialize_bMat(m, chordal, inputData, fpOut);
    chordal.terminate();

    //  rMessage("newton.computeFormula_SDP: start");
    newton.computeFormula_SDP(inputData, 0.0, KAPPA);
    //  rMessage("newton.computeFormula_SDP: end");

    // set initial solutions.
    Solutions currentPt;
    WorkVariables work;
    DenseLinearSpace initPt_xMat;
    DenseLinearSpace initPt_zMat;

    currentPt.initialize(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock, param.lambdaStar, com);
    work.initialize(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

    if (isInitFile) {
        FILE *fpInit = NULL;
        if ((fpInit = fopen(initFile, "r")) == NULL) {
            rError("Cannot open init file " << initFile);
        }
        IO::read(fpInit, currentPt.xMat, currentPt.yVec, currentPt.zMat, nBlock, blockStruct, blockType, blockNumber, isInitSparse);
        fclose(fpInit);
        // This return value used to be discarded.  computeInverse() Cholesky-
        // factorises the supplied X and Z; if either is not positive definite
        // there is no invCholeskyX / invCholeskyZ, and every later step runs on
        // whatever was left in those arrays.  Stop instead of solving on garbage.
        // (No delete[] here for the same reason rError() has none: the process is
        // about to exit, and the OS reclaims the arrays.)
        if (currentPt.computeInverse(work, com) == false) {
            reportSolveFailure(fpOut, outFile, "the initial point is not positive definite (Cholesky factorisation of X or Z failed)", -1, FAIL_FATAL);
            fclose(fpOut);
            return SDPA_EXIT_NUMERICAL_FAILURE;
        }

        initPt_xMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

        initPt_zMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
        initPt_xMat.copyFrom(currentPt.xMat);
        initPt_zMat.copyFrom(currentPt.zMat);
    }
    //  rMessage("initial xMat = "); initPt_xMat.display(Display);
    //  rMessage("initial yVec = "); currentPt.yVec.display(Display);
    //  rMessage("initial zMat = "); initPt_zMat.display(Display);
    //  rMessage("current pt = "); currentPt.display(Display);

    TimeEnd(FILE_READ_END1);
    com.FileRead += TimeCal(FILE_READ_START1, FILE_READ_END1);
    // -------------------------------------------------------------
    // the end of file read
    // -------------------------------------------------------------

    // 2026-08-05 ("C2"): never allocate the second copy of the residual at all.
    //
    // Upstream built the residual into `initRes`, copied the whole thing into
    // `currentRes`, and then held BOTH for the entire solve. The arrays are
    // primalVec (m elements) + dualMat (sum_l n_l^2 elements), i.e. exactly
    // (m + sum_l n_l^2) * sizeof(dd_real) per copy -- 748 MB per copy on maxG60, measured.
    //
    // initRes is NOT dead after the copy: it is still read by the RatioInitResCurrentRes
    // ctor, by the Phase ctor, and by theta.update_exact() inside the main loop. But every
    // one of those reads touches ONLY the scalars normPrimalVec / normDualMat, so initRes
    // never needs its arrays at all. Build the residual straight into currentRes and keep
    // an initRes that carries the scalars and nothing else: Residuals' default ctor leaves
    // primalVec.ele and every dualMat block pointer NULL, so no array is ever allocated for
    // it and ~Residuals()'s terminate() is a no-op on it.
    //
    // Bit-identical: currentRes now runs the same initialize()+compute() that initRes used
    // to run on the same inputData/currentPt, instead of receiving an element-wise copy of
    // that same result, and the three scalars are carried across verbatim.
    //
    // MEASURED, and this is why the change is shaped this way rather than as an
    // `initRes.terminate()` after the copy: peak RSS is reached AT the copy and is never
    // approached again (mcp250-1 15.06 -> 17.09 MB across these two lines, then 17.09 MB
    // flat for the whole main loop, +0.016 MB by the end; maxG60 8980 -> 10476 MB, same
    // shape). Freeing after the copy therefore cannot move the high-water mark -- only not
    // allocating the copy can.
    Residuals currentRes(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock, inputData, currentPt);
    Residuals initRes;
    initRes.normPrimalVec = currentRes.normPrimalVec;
    initRes.normDualMat = currentRes.normDualMat;
    initRes.centerNorm = currentRes.centerNorm;
    // rMessage("initial currentRes = ");
    // currentRes.display(Display);

    StepLength alpha;
    DirectionParameter beta(param.betaStar);
    Switch reduction(Switch::ON);
    AverageComplementarity mu(param.lambdaStar);

    // rMessage("init mu"); mu.display();

    if (isInitFile) {
        mu.initialize(currentPt);
    }

    RatioInitResCurrentRes theta(param, initRes);
    SolveInfo solveInfo(inputData, currentPt, mu.initial, param.omegaStar);
    Phase phase(initRes, solveInfo, param, currentPt.nDim);

    int pIteration = 0;

    // Solve-failure state.  See the EXIT-STATUS POLICY at the top of this file.
    // A NULL reason means "no failure": every break out of the main loop that
    // leaves this NULL is a legitimate stopping condition and keeps exit 0.
    const char *failureReason = NULL;
    int failureIteration = -1;
    bool iterateCorrupted = false;

    IO::printHeader(fpOut, Display);
    // -----------------------------------------------------
    // Here is MAINLOOP
    // -----------------------------------------------------

    TimeStart(MAIN_LOOP_START1);

    // explicit maxIteration
    // param.maxIteration = 2;
    while (phase.updateCheck(currentRes, solveInfo, param) && pIteration < param.maxIteration) {
        // rMessage(" turn hajimari " << pIteration );
        // Mehrotra's Predictor
        TimeStart(MEHROTRA_PREDICTOR_START1);
        // set variable of Mehrotra
        reduction.MehrotraPredictor(phase);
        beta.MehrotraPredictor(phase, reduction, param);
        // rMessage("reduction = "); reduction.display();
        // rMessage("phase = "); phase.display();
        // rMessage("beta.predictor.value = " << beta.value);
        // rMessage(" mu = " << mu.current);
        // rMessage("currentPt = "); currentPt.display();
        bool isSuccessCholesky;
        isSuccessCholesky = newton.Mehrotra(Newton::PREDICTOR, inputData, currentPt, currentRes, mu, beta, reduction, phase, work, com);
        if (isSuccessCholesky == false) {
            // Before this change this break was the ONLY reaction to a failed
            // factorisation: the loop stopped, the last iterate was printed as a
            // solution, and the process exited 0.
            failureReason = "Cholesky factorisation of the Schur complement matrix failed in the Mehrotra predictor";
            failureIteration = pIteration;
            break;
        }
        // rMessage("newton predictor = "); newton.display();

        TimeEnd(MEHROTRA_PREDICTOR_END1);
        com.Predictor += TimeCal(MEHROTRA_PREDICTOR_START1, MEHROTRA_PREDICTOR_END1);
        TimeStart(STEP_PRE_START1);
        alpha.MehrotraPredictor(inputData, currentPt, phase, newton, work, com);
        // rMessage("alpha predictor = "); alpha.display();

        TimeStart(STEP_PRE_END1);
        com.StepPredictor += TimeCal(STEP_PRE_START1, STEP_PRE_END1);

        // rMessage("alphaStar = " << param.alphaStar);
        // Mehrotra's Corrector
        // rMessage(" Corrector ");

        TimeStart(CORRECTOR_START1);
        beta.MehrotraCorrector(phase, alpha, currentPt, newton, mu, param);
        // rMessage("beta corrector = " << beta.value);

#if 1 // 2007/08/29 kazuhide nakata
      // add stopping criteria: objValPrimal < ObjValDual
      //	if ((pIteration > 10) &&
        if ((phase.value == SolveInfo::pdFEAS) && ((beta.value > 5) || (solveInfo.objValPrimal < solveInfo.objValDual))) {
            break;
        }
#endif
        // The corrector's return value used to be discarded here.  MEASURED, and
        // this qualifies the review's claim that "Mehrotra's return value is not
        // checked at all": Newton::compute_DyVec() factorises only for the
        // PREDICTOR (sdpa_newton.cpp), so the CORRECTOR call cannot return false
        // today and this check is provably behaviour-neutral.  That is precisely
        // why it is both safe and necessary to add -- nothing but this line stops
        // a future factorisation on the corrector path from being silently
        // ignored the way the predictor's was.
        if (newton.Mehrotra(Newton::CORRECTOR, inputData, currentPt, currentRes, mu, beta, reduction, phase, work, com) == false) {
            failureReason = "Cholesky factorisation of the Schur complement matrix failed in the Mehrotra corrector";
            failureIteration = pIteration;
            break;
        }

        // rMessage("currentPt = "); currentPt.display();
        // rMessage("newton corrector = "); newton.display();

        TimeEnd(CORRECTOR_END1);
        com.Corrector += TimeCal(CORRECTOR_START1, CORRECTOR_END1);
        TimeStart(CORRECTOR_STEP_START1);
        alpha.MehrotraCorrector(inputData, currentPt, phase, reduction, newton, mu, theta, work, param, com);
        // rMessage("alpha corrector = "); alpha.display();
        TimeEnd(CORRECTOR_STEP_END1);
        com.StepCorrector += TimeCal(CORRECTOR_STEP_START1, CORRECTOR_STEP_END1);
        // the end of Corrector

        IO::printOneIteration(pIteration, mu, theta, solveInfo, alpha, beta, fpOut, Display);
        if (currentPt.update(alpha, newton, work, com) == false) {
            // update() fails for two different reasons and upstream treated both
            // as "cannot move".  Solutions::notPositiveDefinite now separates
            // them (see sdpa_dataset.h):
            //   - step length collapsed: the iterate is still a valid interior
            //     point, so this stays a legitimate stop with exit status 0 and
            //     the same "cannot move" message as before.
            //   - Cholesky of the new X or Z failed: the iterate has left the
            //     positive definite cone and is not a solution.
            if (currentPt.restoredToLastIterate) {
                // The update left the cone but was rolled back: currentPt again holds
                // the last valid iterate (to one rounding per entry). Report PARTIAL
                // and print that point, exit 3 -- not the old behaviour of printing a
                // corrupted point with exit 0, and not the needless one of printing
                // nothing: at this fork's default parameters real SDPLIB problems end
                // exactly here.
                failureReason = "the updated X or Z left the positive definite cone; the step was rolled back and the last valid iterate is reported";
                failureIteration = pIteration;
            } else if (currentPt.notPositiveDefinite) {
                failureReason = "Cholesky factorisation failed on the updated X or Z and the rollback could not be refactored: the iterate is corrupted";
                failureIteration = pIteration;
                iterateCorrupted = true;
            } else {
                // if step length is too short,
                // we finish algorithm
                rMessage("cannot move");
                //   memo by kazuhide nakata
                //   StepLength::MehrotraCorrector
                //   thetaMax*mu.initial -> thetamax*thetaMax*mu.initial
            }
            break;
        }

        // rMessage("currentPt = "); currentPt.display();
        // rMessage("updated");
        theta.update(reduction, alpha);
        mu.update(currentPt);
        currentRes.update(m, inputData, currentPt, com);
        theta.update_exact(initRes, currentRes);
        if (isInitFile) {
            solveInfo.update(inputData, initPt_xMat, initPt_zMat, currentPt, currentRes, mu, theta, param);
        } else {
            solveInfo.update(param.lambdaStar, inputData, currentPt, currentRes, mu, theta, param);
        }
        // 2007/09/18 kazuhide nakata
        // print information of ObjVal, residual, gap, complementarity
        //	solveInfo.check(inputData, currentPt, currentRes, mu, theta, param);
        pIteration++;
    } // end of MAIN_LOOP

    TimeEnd(MAIN_LOOP_END1);

    com.MainLoop = TimeCal(MAIN_LOOP_START1, MAIN_LOOP_END1);
    currentRes.compute(m, inputData, currentPt);
    TimeEnd(TOTAL_TIME_END1);

    com.TotalTime = TimeCal(TOTAL_TIME_START1, TOTAL_TIME_END1);
#if REVERSE_PRIMAL_DUAL
    phase.reverse();
#endif
    // failureIteration == 0 means NO iteration completed, so currentPt still holds
    // the untouched initial point (xVec zeros, xMat = zMat = lambdaStar*I) and
    // printLastInfo() would dress it up as a solution.  That, and only that, is
    // fatal.  A factorisation that failed after real iterations is how this solver
    // normally terminates -- see the EXIT-STATUS POLICY at the top of this file --
    // so its iterate is printed exactly as before and only stderr is told.
    const bool fatalFailure = (failureReason != NULL && failureIteration == 0);
    const bool suppressSolution = fatalFailure || iterateCorrupted;
    if (failureReason != NULL) {
        reportSolveFailure(fpOut, outFile, failureReason, failureIteration,
                           fatalFailure ? FAIL_FATAL : (iterateCorrupted ? FAIL_CORRUPTED : FAIL_PARTIAL));
    }
    if (suppressSolution) {
        // Deliberately no printLastInfo(): either nothing was computed, or the
        // update already destroyed the last valid iterate (sdpa_dataset.h says
        // that point "is not a solution and must not be printed as one").
    } else {
#if 1
        IO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta, currentRes, phase, currentPt, com.TotalTime, nBlock, blockStruct, blockType, blockNumber, inputData, work, com, param, fpOut, Display);
#else
        IO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta, currentRes, phase, currentPt, com.TotalTime, inputData, work, com, param, fpOut, Display);
#endif
    }
    // com.display(fpOut);

    if (SDP_blockStruct) {
        delete[] SDP_blockStruct;
        SDP_blockStruct = NULL;
    }
    if (SOCP_blockStruct) {
        delete[] SOCP_blockStruct;
        SOCP_blockStruct = NULL;
    }
    if (blockStruct) {
        delete[] blockStruct;
        blockStruct = NULL;
    }
    if (blockType) {
        delete[] blockType;
        blockType = NULL;
    }
    if (blockNumber) {
        delete[] blockNumber;
        blockNumber = NULL;
    }

    fprintf(Display, "  main loop time = %.6f\n", com.MainLoop);
    fprintf(fpOut, "    main loop time = %.6f\n", com.MainLoop);
    fprintf(Display, "      total time = %.6f\n", com.TotalTime);
    fprintf(fpOut, "        total time = %.6f\n", com.TotalTime);
#if 0
  fprintf(Display,   "file  check time = %.6f\n",com.FileCheck);
  fprintf(fpOut,   "  file  check time = %.6f\n",com.FileCheck);
  fprintf(Display,   "file change time = %.6f\n",com.FileChange);
  fprintf(fpOut,   "  file change time = %.6f\n",com.FileChange);
#endif
    fprintf(Display, "file   read time = %.6f\n", com.FileRead);
    fprintf(fpOut, "  file   read time = %.6f\n", com.FileRead);
    fclose(fpOut);

#if 0
  rMessage("memory release");
  currentRes.terminate();
  initRes.terminate();
  currentPt.terminate();
  initPt_xMat.terminate();
  initPt_zMat.terminate();
  newton.terminate();
  work.terminate();
  inputData.terminate();
  com.~ComputeTime();
  param.~Parameter();
  alpha.~StepLength();
  beta.~DirectionParameter();
  reduction.~Switch();
  mu.~AverageComplementarity();
  theta.~RatioInitResCurrentRes();
  solveInfo.~SolveInfo();
  phase.~Phase();
#endif

    if (fatalFailure || iterateCorrupted) {
        return SDPA_EXIT_NUMERICAL_FAILURE;
    }
    return failureReason != NULL ? SDPA_EXIT_PARTIAL : SDPA_EXIT_OK;
}

static void message(char *argv0) {
    cout << endl;
    cout << "*** Please assign data file and output file.***" << endl;
    cout << endl;
    cout << "---- option type 1 ------------" << endl;
    cout << argv0
         << " DataFile OutputFile [InitialPtFile]"
            " [-pt parameters]"
         << endl;
    cout << "parameters = 0 default, 1 fast (unstable),"
            " 2 slow (stable)"
         << endl;
    cout << "example1-1: " << argv0 << " example1.dat example1.result" << endl;
    cout << "example1-2: " << argv0 << " example1.dat-s example1.result" << endl;
    cout << "example1-3: " << argv0 << " example1.dat example1.result example1.ini" << endl;
    cout << "example1-4: " << argv0 << " example1.dat example1.result -pt 2" << endl;

    cout << endl;
    cout << "---- option type 2 ------------" << endl;
    cout << argv0 << " [option filename]+ " << endl;
    cout << "  -dd : data dense :: -ds : data sparse     " << endl;
    cout << "  -id : init dense :: -is : init sparse     " << endl;
    cout << "  -o  : output     :: -p  : parameter       " << endl;
    cout << "  -pt : parameters , 0 default, 1 fast (unstable)" << endl;
    cout << "                     2 slow (stable)         " << endl;
    // cout << "  -k  : Kappa(RealValue)" << endl;
    cout << "example2-1: " << argv0 << " -o example1.result -dd example1.dat" << endl;
    cout << "example2-2: " << argv0 << " -ds example1.dat-s -o example2.result "
         << "-p param.sdpa" << endl;
    cout << "example2-3: " << argv0 << " -ds example1.dat-s -o example3.result "
         << "-pt 2" << endl;
    exit(1);
}

} // namespace sdpa

using namespace sdpa;

int main(int argc, char **argv) {
    FILE *Display = stdout;
    setbuf(Display, NULL);

    time_t ltime;
    time(&ltime);

    unsigned int oldcw;
    fpu_fix_start(&oldcw);

    cout << "SDPA-DD start at    " << ctime(&ltime);
    // << "... (built at "<< __DATE__ << " " <<__TIME__ ")" << endl;
    // cout << "let me see your ..." << endl;

    bool isInitFile = false;
    bool isInitSparse = false;
    bool isOutFile = false;
    bool isDataSparse = false;
    bool isParameter = false;

    char *dataFile = NULL;
    char *initFile = NULL;
    char *outFile = NULL;
    char *paraFile = NULL;

    Parameter::parameterType parameterType = Parameter::PARAMETER_DEFAULT;

    if (argc == 1) {
        message(argv[0]);
    }
    if (argv[1][0] == '-') {
        // rsdpa argument

        for (int index = 0; index < argc; ++index) {
            char *target = argv[index];
            if (strcmp(target, "-dd") == 0 && index + 1 < argc) {
                dataFile = argv[index + 1];
                index++;
                continue;
            }
            if (strcmp(target, "-ds") == 0 && index + 1 < argc) {
                dataFile = argv[index + 1];
                index++;
                isDataSparse = true;
                continue;
            }
            if (strcmp(target, "-id") == 0 && index + 1 < argc) {
                initFile = argv[index + 1];
                index++;
                isInitFile = true;
                continue;
            }
            if (strcmp(target, "-is") == 0 && index + 1 < argc) {
                initFile = argv[index + 1];
                index++;
                isInitFile = true;
                isInitSparse = true;
                continue;
            }
            if (strcmp(target, "-o") == 0 && index + 1 < argc) {
                outFile = argv[index + 1];
                index++;
                isOutFile = true;
                continue;
            }
            if (strcmp(target, "-p") == 0 && index + 1 < argc) {
                paraFile = argv[index + 1];
                index++;
                isParameter = true;
                continue;
            }
            if (strcmp(target, "-k") == 0 && index + 1 < argc) {
                KAPPA = atof(argv[index + 1]);
                rMessage("Kappa = " << KAPPA);
                index++;
                continue;
            }
            if (strcmp(target, "-pt") == 0 && index + 1 < argc) {
                int tmp = atoi(argv[index + 1]);
                switch (tmp) {
                case 0:
                    parameterType = Parameter::PARAMETER_DEFAULT;
                    break;
                case 1:
                    parameterType = Parameter::PARAMETER_UNSTABLE_BUT_FAST;
                    break;
                case 2:
                    parameterType = Parameter::PARAMETER_STABLE_BUT_SLOW;
                    break;
                default:
                    parameterType = Parameter::PARAMETER_DEFAULT;
                }
                index++;
                paraFile = NULL;
                isParameter = false;
                continue;
            }
        }
    } else { // SDPA argument
        dataFile = argv[1];
        int len = strlen(dataFile);
        if (dataFile[len - 1] == 's' && dataFile[len - 2] == '-') {
            isDataSparse = true;
        }

        outFile = argv[2];

        paraFile = (char *)"./param.sdpa";
        isParameter = true;

        for (int index = 3; index < argc; ++index) {
            if (strcmp(argv[index], "-pt") == 0 && index + 1 < argc) {
                int tmp = atoi(argv[index + 1]);
                switch (tmp) {
                case 0:
                    parameterType = Parameter::PARAMETER_DEFAULT;
                    break;
                case 1:
                    parameterType = Parameter::PARAMETER_UNSTABLE_BUT_FAST;
                    break;
                case 2:
                    parameterType = Parameter::PARAMETER_STABLE_BUT_SLOW;
                    break;
                default:
                    parameterType = Parameter::PARAMETER_DEFAULT;
                }
                index++;
                paraFile = NULL;
                isParameter = false;
            } // end of "-pt"
            else {
                initFile = argv[index];
                isInitFile = true;
                int len = strlen(initFile);
                if (initFile[len - 1] == 's' && initFile[len - 2] == '-') {
                    isInitSparse = true;
                }
            }
        } // end of 'for'
    }

    if (dataFile == NULL || outFile == NULL) {
        message(argv[0]);
    }

    cout << "data      is " << dataFile;
    if (isDataSparse) {
        cout << " : sparse" << endl;
    } else {
        cout << " : dense" << endl;
    }
    if (paraFile) {
        cout << "parameter is " << paraFile << endl;
    }
    if (outFile) {
        cout << "out       is " << outFile << endl;
    }
    if (initFile) {
        cout << "initial   is " << initFile;
    }
    if (isInitFile) {
        if (isInitSparse) {
            cout << " : sparse" << endl;
        } else {
            cout << " : dense" << endl;
        }
    } else {
        cout << endl;
    }
    if (paraFile == NULL) {
        if (parameterType == Parameter::PARAMETER_DEFAULT) {
            cout << "set       is DEFAULT" << endl;
        } else if (parameterType == Parameter::PARAMETER_UNSTABLE_BUT_FAST) {
            cout << "set       is UNSTABLE_BUT_FAST" << endl;
        } else if (parameterType == Parameter::PARAMETER_STABLE_BUT_SLOW) {
            cout << "set       is STABLE_BUT_SLOW" << endl;
        }
    }
    const int status = pinpal(dataFile, initFile, outFile, paraFile, isInitFile, isInitSparse, isDataSparse, isParameter, parameterType, Display);
    // fpu_fix_end() used to sit AFTER the return and was therefore dead code: the
    // extended-precision FPU mode that fpu_fix_start() installed was never undone.
    // It has to run before the return, not after it.
    fpu_fix_end(&oldcw);
    return status;
}
