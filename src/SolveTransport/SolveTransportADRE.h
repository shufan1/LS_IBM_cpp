#pragma once
#include "../VariableNonDim.h"
#include "../ControlVar.h"

// Mirrors SolveTransport/SolveTransportADRE.m -- one unsteady
// advection-diffusion-reaction step for the scalar field phi. NOT ported
// yet: no-op stub (phi unchanged).
//
// When ported this also needs Flux (getDiffFlux) and IBM_coeffP
// (LSIBMcoeffs) -- add them to the signature then.
void solveTransportADRE(ControlVar &controlVar, const Domain &domain,
                        const Variables &variables, StateVar &stateVar,
                        const IBM &ibm, const BC &bc, const LS &ls);
