// Driver for the real solver, mirrors beta_AD_correctu/RunADRE.m (and, for
// Milestone 1, RunADRE_no_LS.m). Kept flat and thin on purpose -- the actual
// physics lives in the module files below as they get ported, matching the
// MATLAB folder layout (see CMakeLists.txt).
//
// Roadmap only below -- nothing implemented yet. Order matches the staged
// plan: Milestone 1 (flow, frozen geometry) first, validated against the
// MATLAB _no_LS reference trace, before Milestone 2 (+transport) and
// Milestone 3 (+level-set/moving geometry). See
// beta_AD_correctu/PARALLELIZATION_PLAN.md for the full rationale.
#include <petsc.h>
#include <mpi.h>
#include <cstdio>
#include <string>
#include <filesystem>
#include <iostream>
#include "ControlVar.h"
#include "VariableNonDim.h"
#include "Utilities/Coordinates.h"
#include "Utilities/DiffFlux.h"
#include "Utilities/Timer.h"
#include "Models/flowSStransportUS.h"
#include "IBM/LSIBMcoeffs.h"

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // argv[1] is the project directory (e.g. .../Example_Project), which
    // must contain config.json. Passed in by run.sbatch as $PROJECT_DIR.
    if (argc < 2) {
        if (rank == 0) fprintf(stderr, "usage: %s <project_directory>\n", argv[0]);
        PetscFinalize();
        return 1;
    }
    std::string projectDir = argv[1];
    std::string configPath = projectDir + "/config.json";

    // 1. Set up control variables
    //    - mirrors setUpControlVar(VARIABLES, DOMAIN)
    //    - tolerances, iteration caps, output settings, model selection
    //    - moved ahead of step 2 (vs. RunADRE.m's call order) so
    //      output_folder exists before we need it to save coordinates.json
    ControlVar controlVar(configPath);
    std::string outputDir = projectDir + "/" + controlVar.output_folder;
    std::filesystem::create_directories(outputDir);

    // 2. Load variables + coordinates
    //    - mirrors setUpVariablesNonDim.m's 7 return values: LSCase is a
    //      field already; the rest are VariableNonDim methods, called in
    //      dependency order (each depends on Domain, which getDomain()
    //      produces first).
    //    - getDomain() is real now (Utilities/Coordinates.cpp); StateVar/LS
    //      still need LSInitialize()/LSnormals() -- see VariableNonDim.cpp
    //      for exactly what's real vs. TODO in each method.
    VariableNonDim varSetup(configPath);
    LSCase lsCase = varSetup.getLSCase();
    // Note: named `domain`, not `DOMAIN` -- `DOMAIN` is a macro defined by
    // <cmath> (the old SVID matherr() constants: DOMAIN/SING/OVERFLOW/
    // UNDERFLOW/TLOSS/PLOSS), pulled in transitively via <petsc.h>. Using
    // it as an identifier silently preprocesses into a numeric literal.
    Domain domain = varSetup.getDomain();
    Variables VARIABLES = varSetup.getVariables(domain);
    IBM ibm = varSetup.getIBM(domain);
    LS ls = varSetup.getLS(domain, lsCase);                 // masks (psiU/psiV/psi)
    StateVar stateVar = varSetup.getStateVar(domain, ls);   // needs the masks
    BC bc = varSetup.getBC(domain, stateVar, ls);
    VARIABLES.defineReactivity(domain, "A_0.2.json");
    VARIABLES.defineLSvariables(domain);
    


    // Grid is identical on every rank until domain decomposition exists --
    // save once, from rank 0, so xu/yu/xv/yv/xp/yp can be reloaded for
    // plotting later (psi and phi live on the xp/yp nodes, same as p).
    if (rank == 0) {
        saveCoordinatesJson(domain, "coordinates.json");
    }

    if (rank == 0) {
        printf("LS_IBM solver -- no modules ported yet\n");
        printf("  lsCase.caseId=%d xc=%g yc=%g diamcyl=%g\n",
               lsCase.caseId, lsCase.xc, lsCase.yc, lsCase.diamcyl);
        printf("  domain: lx=%g ly=%g imax=%d jmax=%d \n",
               domain.lx, domain.ly, domain.imax, domain.jmax);
        printf("coordinates.json saved to %s", outputDir.c_str());
        printf("  variables: Re=%g Pe=%g D=%g dt=%g Np=%d\n",
               VARIABLES.Re, VARIABLES.Pe, VARIABLES.D, VARIABLES.dt, VARIABLES.Np);
        printf("  bc: BC_e_u=%d BC_w_u=%d P0_e=%g\n", bc.BC_e_u, bc.BC_w_u, bc.P0_e);
        printf("  ibm: q=%g alpha=%g beta=%g\n", ibm.q, ibm.alpha, ibm.beta);
        printf("  stateVar: U.size()=%zu V.size()=%zu P.size()=%zu\n",
               stateVar.U.data().size(), stateVar.V.data().size(), stateVar.P.data().size());
        printf("  ls: caseId=%d psi.size()=%zu\n", ls.caseId, ls.psi.data().size());
        // DIAGNOSTIC (temporary): raw P-grid psi at the two cells MATLAB
        // reports LS.psi(195,62)=1.5784e-04 for -- compare directly to see
        // whether the level-set field itself already disagrees here, vs.
        // agreeing on psi but disagreeing in the U-grid averaging/flagging
        // built on top of it.
        printf("  ls.psi(195,62)=%.10e  ls.psi(196,62)=%.10e\n", ls.psi(195, 62), ls.psi(196, 62));

        // reaction matrix A (Np x Np, row-major flat) from defineReactivity()
        printf("  VARIABLES.A (%dx%d):\n", VARIABLES.Np, VARIABLES.Np);
        for (int i = 0; i < VARIABLES.Np; ++i) {
            printf("   ");
            for (int j = 0; j < VARIABLES.Np; ++j)
                printf(" %12.6g", VARIABLES.A[i * VARIABLES.Np + j]);
            printf("\n");
        }
    }

    // 3. Compute diffusive flux setup
    //    - mirrors getDiffFlux(VARIABLES, DOMAIN, BC) -- solution-independent
    //      (mesh + Re/D only), so computed once here, not per iteration.
    Flux flux = computeDiffFlux(domain, VARIABLES, bc);
    if (rank == 0) {
        printf("  flux: Diffu_U.De.nx()=%d Diffu_U.De.ny()=%d Diffu_V.De.nx()=%d Diffu_V.De.ny()=%d\n",
               flux.Diffu_U.De.nx(), flux.Diffu_U.De.ny(), flux.Diffu_V.De.nx(), flux.Diffu_V.De.ny());
    }

    // 4. Run model -- inlined directly here rather than a separate
    //    Models/modelSimulation.cpp file. Mirrors
    //    Models/modelSimulation_no_LS.m, model==4 branch only (models
    //    1-3 exist in MATLAB but aren't used in practice, so not
    //    translated). LS.psi is frozen (geometry doesn't move yet),
    //    matching RunADRE_no_LS.m / Milestone 1's scope.
    //
    //    Only the driving structure (prints, CurrentStateVar save) is
    //    real for now -- the actual per-iteration solve
    //    (LSVelocityExtrapolation, LSIBMcoeffs, flowSStransportUS) isn't
    //    ported yet, so StateVar/LS never actually change across
    //    iterations until those land.

    // Only flow_steady=true && transport_steady=false ("model 4" in
    // MATLAB's numbering) is translated -- bail out on every rank
    // together, not just rank 0, so a differently-configured run doesn't
    // leave other ranks waiting on one that already exited.
    if (controlVar.flow_steady && !controlVar.transport_steady) {
        if (rank == 0) {
        printf("=========================================================================== \n");
        printf("       Model 4 : Simulation starts with flow steady state and transport and LS being unsteady:  (LS frozen -- no_LS variant)\n");
        printf("=========================================================================== \n");
        }
    } else {
        if (rank == 0) {
            fprintf(stderr,
                    "only Model 4 flow_steady=true && transport_steady=false is implemented "
                    "(got flow_steady=%d, transport_steady=%d)\n",
                    controlVar.flow_steady, controlVar.transport_steady);
        }
        PetscFinalize();
        return 1;
    }

    

    // DEBUG (temporary): dump the initial state (before any time-stepping,
    // in particular the initial psi) in the same format saveCurrentStateJson()
    // uses inside the loop below, so it can be compared directly against
    // MATLAB's own pre-modelSimulation_no_LS() snapshot -- see
    // debug_compare/compare_initial_state.m.
    if (rank == 0) {
        std::string initPath = "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/cpp_initial_state.json";
        varSetup.saveCurrentStateJson(stateVar, ls, controlVar.time, initPath);
        printf("  initial state saved to %s\n", initPath.c_str());
    }

    double endTime = controlVar.time + controlVar.noLStime * VARIABLES.dt;
    int iTime = controlVar.iStart;

    while (controlVar.time < endTime) {
        controlVar.time += VARIABLES.dt;
        iTime += 1;

        if (rank == 0) {
            printf("=========================================================================== \n");
            printf("                     updating geometry at time %d = %8.6f / %8.6f             \n",
                   iTime,controlVar.time, endTime);
            printf("=========================================================================== \n");
        }

        // ---- LEVEL SET EQUATION (FROZEN) ----
        // TODO: LSVelocityExtrapolation(variables, stateVar, domain, ls)
        // -- not ported yet; populates LS.u/LS.v/LS.q_out/LS.beta_out
        // (the latter two are what the P/scalar grid would need for a
        // real per-species Robin BC -- LSPointIdent() is still a
        // placeholder either way, see IBM/LSPointIdent.h).
        //
        // One call, both grids at once, matching LSIBMcoeffs.m -- index
        // 0=U, 1=V. No P slot -- see IBM/LSIBMcoeffs.h for why P's own
        // IBM_coeffP is never computed at all.
        std::vector<IBMCoeff> ibmCoeff = LSIBMcoeffs(ibm, domain, ls, stateVar.phi);
        const IBMCoeff &ibmCoeffU = ibmCoeff[0];
        const IBMCoeff &ibmCoeffV = ibmCoeff[1];
        if (rank == 0) {
            printf("  ibmCoeff numg: U=%d V=%d\n", ibmCoeffU.numg, ibmCoeffV.numg);
        }
        // DIAGNOSTIC (temporary): does any U ghost cell's own bilinear
        // stencil corner (I1..I4) land on a non-fluid cell (flag!=0 --
        // itself a ghost or solid cell, not a directly-known value)? Ported
        // faithfully from LSmirPointsBQ.m's own bilinear corner search,
        // which -- in both languages -- never checks the corner's flag at
        // all, so this can happen; checking whether it actually does here,
        // and where, for this geometry.
        if (rank == 0) {
            int nBad = 0;
            for (int k = 0; k < ibmCoeffU.numg; ++k) {
                int corners[4][2] = {{ibmCoeffU.I1[k], ibmCoeffU.J1[k]},
                                      {ibmCoeffU.I2[k], ibmCoeffU.J2[k]},
                                      {ibmCoeffU.I3[k], ibmCoeffU.J3[k]},
                                      {ibmCoeffU.I4[k], ibmCoeffU.J4[k]}};
                for (auto &c : corners) {
                    double f = ibmCoeffU.flag(c[0], c[1]);
                    if (f != 0.0) {
                        printf("  U ghost (%d,%d) stencil corner (%d,%d) has flag=%.0f (non-fluid)\n",
                               ibmCoeffU.I_g[k], ibmCoeffU.J_g[k], c[0], c[1], f);
                        ++nBad;
                    }
                }
            }
            printf("  U ghost-stencil non-fluid-corner count: %d\n", nBad);

            // DIAGNOSTIC (temporary): full stencil dump for specific U
            // ghost cells requested by (I_g,J_g) -- ghost cell location,
            // its 4 bilinear stencil corners, their lambda_g_k weights, and
            // A1_g, so this can be compared directly against MATLAB's own
            // IBM_coeffU fields for the same cells.
            int wanted[][2] = {{194, 61}, {195, 61}};
            bool found[2] = {false, false};
            for (int k = 0; k < ibmCoeffU.numg; ++k) {
                for (int w = 0; w < 2; ++w) {
                    if (ibmCoeffU.I_g[k] == wanted[w][0] && ibmCoeffU.J_g[k] == wanted[w][1]) {
                        found[w] = true;
                        printf("  U ghost k=%d (I_g=%d,J_g=%d): I1..I4=(%d,%d) (%d,%d) (%d,%d) (%d,%d) "
                               "lambda_g_1..4=%.6e %.6e %.6e %.6e A1_g=%.6e\n",
                               k, ibmCoeffU.I_g[k], ibmCoeffU.J_g[k], ibmCoeffU.I1[k], ibmCoeffU.J1[k],
                               ibmCoeffU.I2[k], ibmCoeffU.J2[k], ibmCoeffU.I3[k], ibmCoeffU.J3[k],
                               ibmCoeffU.I4[k], ibmCoeffU.J4[k], ibmCoeffU.lambda_g_1[k],
                               ibmCoeffU.lambda_g_2[k], ibmCoeffU.lambda_g_3[k], ibmCoeffU.lambda_g_4[k],
                               ibmCoeffU.A1_g[k]);
                    }
                }
            }
            // If not found as a ghost cell, report what flag it actually
            // has (0=fluid, 1=ghost, 2=solid), plus its 4 orthogonal
            // neighbors' flags, to see what our classification really did
            // there vs. what was expected.
            for (int w = 0; w < 2; ++w) {
                if (found[w]) continue;
                int i = wanted[w][0], j = wanted[w][1];
                printf("  U (%d,%d) is NOT a ghost cell -- flag=%.0f, neighbor flags "
                       "E=%.0f W=%.0f N=%.0f S=%.0f, ls.psi(%d,%d)=%.10e\n",
                       i, j, ibmCoeffU.flag(i, j), ibmCoeffU.flag(i + 1, j), ibmCoeffU.flag(i - 1, j),
                       ibmCoeffU.flag(i, j + 1), ibmCoeffU.flag(i, j - 1), i, j, ls.psi(i, j));
            }
        }

        // ---- model == 4: SS flow, unsteady transport, frozen LS ----
        // mirrors modelSimulation_no_LS.m's model==4 branch: split dt into
        // nLSupdate sub-steps, flag steady flow + unsteady transport, run
        // flowSStransportUS (which modifies stateVar in place), then
        // restore dt.
        controlVar.PISO = 0;
        controlVar.flow_steady = true;
        controlVar.transport_steady = false;
        controlVar.noTime = VARIABLES.nLSupdate;

        VARIABLES.dt = VARIABLES.dt / VARIABLES.nLSupdate;
        flowSStransportUS(stateVar, controlVar, domain, VARIABLES, ibm, ibmCoeffU, ibmCoeffV, bc,
                           ls, flux);
        VARIABLES.dt = VARIABLES.dt * VARIABLES.nLSupdate;

        // ----  save current state ----
        if (rank == 0 && iTime % controlVar.savedat == 0) {
            std::string flow_filename = "dataRDE" + std::to_string(iTime) + "dt.json";
            std::string output_filepath = outputDir + "/" + flow_filename;
            varSetup.saveCurrentStateJson(stateVar, ls, controlVar.time, output_filepath);
            printf("output saved to %s \n", output_filepath.c_str());
        }

        // TODO: NaN/blowup check (`if sum(sum(isnan(StateVar.U))) ...
        // break`) -- meaningless until the flow solve actually populates
        // StateVar.U.
    }

    // Mirrors RunADRE.m's TIMING SUMMARY block + timer_results.mat save --
    // see Utilities/Timer.h. Rank 0 only, same as every other summary/save.
    if (rank == 0) {
        Timer::printSummary();
        Timer::saveJson(outputDir + "/timer_results.json");
    }

    PetscFinalize();
    return 0;
}
