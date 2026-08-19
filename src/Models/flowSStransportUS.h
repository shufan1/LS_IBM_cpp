#pragma once
#include <vector>
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "../Utilities/DiffuFlux.h"
#include "../IBM/IBMCoeff.h"

// Mirrors Models/flowSStransportUS.m (model 4): solve the flow to steady
// state once, then march the scalar transport equation for
// controlVar.noTime unsteady sub-steps, breaking early on blow-up.
// Modifies stateVar in place.
//
// ibmCoeffU/V are LSIBMcoeffs()'s indices 0/1, threaded straight through
// to solveUVP() -> coeffU()/coeffV(). ibmCoeffPhi (one per species) is
// threaded straight through to solveTransportADRE().
//
// varSetup/outputAdreDir: DEBUG ONLY, for validating solveTransportADRE()
// call-by-call against MATLAB -- if outputAdreDir is non-empty, dumps
// stateVar/ls to <outputAdreDir>/dataRDE<controlVar.iTime>_<sub>dt.json
// after every one of the noTime sub-step calls below (sub = 1..noTime),
// not just the single post-loop save main.cpp's own outer iTime loop
// does. No MATLAB reference exists at this granularity yet -- pass ""
// to skip (varSetup is unused in that case).
void flowSStransportUS(StateVar &stateVar, ControlVar &controlVar,
                       const Domain &domain, const Variables &variables, const IBM &ibm,
                       const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV,
                       std::vector<IBMCoeff> &ibmCoeffPhi, const BC &bc,
                       const LS &ls, const Flux &flux, const VariableNonDim &varSetup,
                       const std::string &outputAdreDir);
