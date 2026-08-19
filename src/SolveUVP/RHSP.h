#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "../IBM/IBMCoeff.h"

// Fills RHS_P2, the pressure-correction system's right-hand side: the
// net mass imbalance (divergence) of the current momentum-only velocity
// field (U_star, V_star) over each interior P cell -- west/south inflow
// minus east/north outflow. Zero means that cell already satisfies
// continuity; nonzero is what the pressure correction has to fix.
//

// Before computing the residual, RHSP.m forces U_star/V_star to a fixed
// placeholder value wherever flag_u/flag_v==1 (ghost cells --
// LSPointIdent.m's flag==1, not ==2/solid; see RHSP.cpp) -- ghost values
// are extrapolated to satisfy the immersed-boundary condition, not
// derived from a divergence-respecting formula, so using them as-is
// here would inject a non-physical contribution into a neighboring real
// fluid cell's continuity check. Applied inline when reading U_star/
// V_star below (not by mutating them -- they're const references here,
// same as MATLAB's own pass-by-value scoping keeps the substitution
// local to this function). Takes ibmCoeffU/ibmCoeffV -- the same fresh,
// per-iteration ghost/solid classification coeffU()/coeffV() already
// consume (IBMCoeff::flag) -- rather than the IBM struct's own
// flag_u/flag_v, which are one-time-zeroed in VariableNonDim.cpp and
// never updated (an earlier version of this function read those before
// LSPointIdent()/LSIBMcoeffs() were ported; it was left pointing at the
// permanently-zero placeholder after they were, silently turning this
// masking step into a no-op).
//
// RHS_P2 is caller-allocated (same size/lifetime discipline as
// coeffU()/coeffV()/coeffP()'s Mat/Vec -- see COEFFU.h) and filled in
// place via VecSetValue.
void rhsP(const Field2D &U_star, const Field2D &V_star, const Domain &domain,
          const ControlVar &controlVar, const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV,
          Vec RHS_P2);
