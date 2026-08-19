#include "flowSStransportUS.h"
#include "../SolveUVP/SolveUVP.h"
#include "../SolveTransport/SolveTransportADRE.h"
#include <cmath>
#include <cstdio>
#include <mpi.h>

namespace {

// MATLAB break condition: sum(sum(isnan(StateVar.U))) || max(...phi...) > 1e4
bool blewUp(const StateVar &stateVar) {
    for (double u : stateVar.U.data())
        if (std::isnan(u)) return true;
    for (const Field2D &ph : stateVar.phi)
        for (double v : ph.data())
            if (v > 1e4) return true;
    return false;
}

}  // namespace

void flowSStransportUS(StateVar &stateVar, ControlVar &controlVar,
                       const Domain &domain, const Variables &variables, const IBM &ibm,
                       const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV,
                       std::vector<IBMCoeff> &ibmCoeffPhi, const BC &bc,
                       const LS &ls, const Flux &diffu_flux, const VariableNonDim &varSetup,
                       const std::string &outputAdreDir) {
    // ---- SOLVE FOR FLOW (steady state), once ----
    solveUVP(controlVar, domain, variables, stateVar, ibm, ibmCoeffU, ibmCoeffV, bc, diffu_flux);
    printf("%s\n", controlVar.messageFlow.c_str());

    int rank = 0;
    if (!outputAdreDir.empty()) MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ---- SCALAR TRANSPORT: march noTime unsteady sub-steps ----
    // (MATLAB advances a local `time` here for logging; the outer model
    // loop owns controlVar.time, so this is display-only.)
    double time = controlVar.time;
    for (int iTime = 1; iTime <= controlVar.noTime; ++iTime) {
        time += variables.dt;
        printf("\n ~~~~~~~~~~~~~~~~~~~~ time = %8.6f ~~~~~~~~~~~~~~~~~~~~~ \n", time);

        solveTransportADRE(controlVar, domain, variables, stateVar, ibm, ibmCoeffPhi, bc, ls, diffu_flux, iTime,
                           rank == 0 ? outputAdreDir : "");

        // DEBUG (temporary): per-substep dump, see flowSStransportUS.h's
        // comment -- distinct filename from main.cpp's own once-per-outer-
        // iTime dataRDE<iTime>dt.json so the two don't collide.
        if (!outputAdreDir.empty() && rank == 0) {
            std::string path = outputAdreDir + "/dataRDE" + std::to_string(controlVar.iTime) + "_" +
                                std::to_string(iTime) + "dt.json";
            varSetup.saveCurrentStateJson(stateVar, ls, time, path);
            printf("  [ADRE debug] saved %s\n", path.c_str());
        }

        if (blewUp(stateVar)) break;
    }
}
