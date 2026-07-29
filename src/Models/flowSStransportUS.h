#pragma once
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "../Utilities/DiffFlux.h"
#include "../IBM/IBMCoeff.h"

// Mirrors Models/flowSStransportUS.m (model 4): solve the flow to steady
// state once, then march the scalar transport equation for
// controlVar.noTime unsteady sub-steps, breaking early on blow-up.
// Modifies stateVar in place.
//
// SolveTransportADRE is NOT ported yet -- it's a no-op stub inside the
// .cpp. solveUVP() is real for the momentum solve (see SolveUVP.h/.cpp);
// its own pressure-correction step is still a placeholder.
//
// ibmCoeffU/V are LSIBMcoeffs()'s indices 0/1, threaded straight through
// to solveUVP() -> coeffU()/coeffV(). There's no ibmCoeffP -- see
// IBM/LSIBMcoeffs.h for why P's own IBM_coeffP is never computed at all.
void flowSStransportUS(StateVar &stateVar, ControlVar &controlVar,
                       const Domain &domain, const Variables &variables, const IBM &ibm,
                       const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV, const BC &bc,
                       const LS &ls, const Flux &flux);
