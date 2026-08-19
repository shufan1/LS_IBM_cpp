#pragma once
#include <vector>
#include "../VariableNonDim.h"
#include "../IBM/IBMCoeff.h"

// Mirrors update_A1g.m's A1_g-only formula, under the linear-reaction
// model (beta_G fixed, so lambda_g_k never depends on q/phi -- see
// LSIBMcoeffsPhi()'s own comment). For the nonlinear case, where beta
// itself depends on phi, see UpdateGhostReactionPhi.h instead.

// Per-iteration (per species): reaction-rate-dependent ghost forcing
// A1_g only, from this iteration's q_G for this species (one value per
// ghost cell, parallel to ibmCoeffPhi.I_g/J_g). Reads ibmCoeffPhi.betaG/
// r_g/Delta (all cached once by LSmirPointsBQ() -- see IBMCoeff.h)
// rather than re-deriving them, since this can't redo that call's own
// mirror-point search/near-domain-edge safeguard itself.
std::vector<double> updateA1gPhi(const IBM &ibm, const IBMCoeff &ibmCoeffPhi,
                                  const std::vector<double> &q_G_species);
