#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../Utilities/DiffuFlux.h"
#include "../Utilities/ConvFlux.h"
#include "../IBM/IBMCoeff.h"

// QUICK deferred-correction terms (S_e/S_w/S_n/S_s, from stateVar.phi[i_s]
// directly, using domain.g1c_*/g2c_* QUICK-interpolation weights from
// Coordinates.cpp) + S0 (old-timestep) + S_ur (under-relaxation, needs
// ap_p from coeffPhiADRE) + the ghost-cell forcing A1_g (this
// iteration's reactive boundary flux, from updateA1gPhi).
//
// rebuilding every QUICK iteration: A1_g changes with q_G every pass.

void rhsPhiADRE(const StateVar &stateVar, const BC &bc, const Flux &flux,
                 const ConvFluxCoeffs &convFluxPhi, const Domain &domain, const IBMCoeff &ibmCoeffPhi,
                 const IBM &ibm, const Variables &variables, const LS &ls, const Field2D &ap_p,
                 const std::vector<double> &A1_g, int speciesIndex, Vec RHS);
