#include "flowSStransportUS.h"
#include "../SolveUVP/SolveUVP.h"
#include "../SolveTransport/SolveTransportADRE.h"
#include <cmath>
#include <cstdio>

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
                       const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV, const BC &bc,
                       const LS &ls, const Flux &flux) {
    // ---- SOLVE FOR FLOW (steady state), once ----
    solveUVP(controlVar, domain, variables, stateVar, ibm, ibmCoeffU, ibmCoeffV, bc, flux);
    printf("%s\n", controlVar.messageFlow.c_str());

    // ---- SCALAR TRANSPORT: march noTime unsteady sub-steps ----
    // (MATLAB advances a local `time` here for logging; the outer model
    // loop owns controlVar.time, so this is display-only.)
    double time = controlVar.time;
    for (int iTime = 1; iTime <= controlVar.noTime; ++iTime) {
        time += variables.dt;
        printf("\n ~~~~~~~~~~~~~~~~~~~~ time = %8.6f ~~~~~~~~~~~~~~~~~~~~~ \n", time);

        solveTransportADRE(controlVar, domain, variables, stateVar, ibm, bc, ls);

        if (blewUp(stateVar)) break;
    }
}
