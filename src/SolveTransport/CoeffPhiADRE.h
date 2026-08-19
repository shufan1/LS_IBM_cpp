#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../Utilities/DiffuFlux.h"
#include "../Utilities/ConvFlux.h"
#include "../IBM/IBMCoeff.h"

// Result of building one species' scalar-transport coefficient matrix.
// Mirrors COEFFPHIADRE.m's outputs -- only ap_p is kept (consumed by
// RHSPHIADRE's S_ur term); aw_p/ae_p/as_p/an_p are intermediate-only
// there (folded into CM, never read back out afterward).
struct TransportCoeffs {
    Field2D ap;
};

// Mirrors COEFFPHIADRE.m end-to-end: interior/boundary ae/aw/an/as, ap
// assembly, immersed-boundary solid/ghost diagonal pinning, AND the
// ghost rows' -lambda_g_k mirror-point coupling (MATLAB's A_g_sparse) --
// all built once and inserted into the caller-allocated/zeroed CM (same
// discipline as coeffU()/coeffV(), see COEFFU.h).
//
// Called once per species before the QUICK loop starts, not once per
// QUICK iteration -- diffusive/convective flux (and, under the linear
// reaction-rate model, lambda_g_k) are frozen for the whole transport
// solve (see SolveTransportADRE.cpp's comment). For the nonlinear-
// reaction-BC case, UpdateGhostReactionPhi.h's updateCoeffGhost() is the
// one that re-pushes a freshly-rederived lambda_g_k into CM's ghost rows
// every QUICK iteration instead -- this function itself is never rerun
// mid-solve.
TransportCoeffs coeffPhiADRE(const Domain &domain, const Flux &flux, const ConvFluxCoeffs &convFluxPhi,
                              const IBMCoeff &ibmCoeffPhi, const IBM &ibm, const Variables &variables,
                              Mat CM);
