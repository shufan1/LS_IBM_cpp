#pragma once
#include "../VariableNonDim.h"

// compute interfacial velocity from flux computed from interface concentration
// updates
//   ls.u, ls.v -> solveHJEq(), which advects psi
// MATLAB's LS.q_out/LS.beta_out are not produced, overwrite in calculateqG()
//
// ONE DELIBERATE NUMERICAL DIVERGENCE from the reference. The recession
// speed has to undo defineReactivity()'s `A = inv(Pd)*A` for the species
// that drives the interface -- the scaled A is right for the transport
// BC, but a recession speed needs the true, diffusivity-independent molar
// flux. MATLAB writes that undo as a bare `/1.261` literal at line 139,
// which is 1/Pd(2,2) = 1/0.793 rounded to four figures. Here it is
// `* variables.Pd[1]`, which cancels exactly and stays correct if Pd
// changes -- but 1/1.261 = 0.7930214 vs 0.793 is a 2.7e-5 relative
// difference in the interface speed, so this is NOT bit-comparable
// against MATLAB. To restore bit-exactness during validation, divide by
// a literal 1.261 instead.
void LSVelocityExtrapolation(LS &ls, const StateVar &stateVar, const Domain &domain,
                              const Variables &variables);
