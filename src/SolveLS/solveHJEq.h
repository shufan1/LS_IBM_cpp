#pragma once
#include "../VariableNonDim.h"
#include "../Utilities/Field2D.h"

// Mirrors SolveLS/solveHJEq.m. Advances the level set one geometry step by
// solving the Hamilton-Jacobi advection equation
//
//     dpsi/dt + u*dpsi/dx + v*dpsi/dy = 0
//
// with the interface velocity from computeLSVelocityExtrapolation().
// Returns the advected psi; does not reinitialize it (that is
// LSreinitialization()'s job, called immediately after).
//
// Time step: MATLAB sets dtLS = VARIABLES.dt, the BIG dt. At the point
// LSeqSolve is called from modelSimulation.m:48, VARIABLES.dt has not yet
// been divided by nLSupdate -- that happens later at line 97 and is undone
// at line 102 -- so the level set advances by the full geometry-update
// interval in one shot. (The nLSupdate-scaled alternative is present but
// commented out at solveHJEq.m:3.)
//
// ---------------------------------------------------------------------
// DELIBERATE BUG-FOR-BUG PORT -- do not "fix" this without asking.
//
// MATLAB's RK3 third stage evaluates the derivative at the WRONG stage
// value. solveHJEq.m:43 reads
//
//     [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation,h,psi_o);
//
// which is byte-identical to the call at line 37 -- it recomputes a
// derivative already sitting in `psinp`. SSP-RK3's third stage needs
// L(psi_n12), the value formed on the line immediately above it:
//
//     psi^(2) = 3/4 psi^n + 1/4 psi^(1) + 1/4 dt L(psi^(1))   <- psi_n12
//     psi^n+1 = 1/3 psi^n + 2/3 psi^(2) + 2/3 dt L(psi^(2))   <- needs psi_n12
//
// The final combination `(psi_n + 2*psi_n32)/3` is right; only the
// argument is wrong. Effect: RK3 silently drops below third order and
// loses the SSP (non-oscillatory) guarantee that was the reason to use
// it. And it is live -- setUpVariablesNonDim.m:428 sets
// TimeSchemeLS = "RK3".
//
// That it is a typo rather than a variant is settled by the sibling:
// LSreinitialization.m:50 does the same third stage correctly, passing
// psi_m12. Same folder, same author, same pattern.
//
// Reproduced here anyway, on purpose. The C++ has to match MATLAB
// bit-for-bit first, so that a genuine porting error cannot hide behind
// a deliberate divergence. FIX ONLY AFTER validation passes -- pass the
// stage-2 value to the third derivative evaluation, and drop the now-
// redundant duplicate call. Fix the MATLAB at the same time, or the two
// stop being comparable.
// ---------------------------------------------------------------------
//
// Time scheme: variables.TimeSchemeLS, one of "RK1"/"RK2"/"RK3". Note the
// RK stages re-evaluate only the spatial derivatives -- u and v are read
// once and held frozen across all stages, exactly as MATLAB does. So a
// higher-order scheme buys advection accuracy, not coupling accuracy: the
// concentration->geometry coupling stays first-order regardless, because
// the velocity itself is lagged one big step.
//
// Spatial discretization: WENO5 with upwinding selected by the sign of
// u/v, via wenoDerivative(..., LSEquation::LevelSetEqn) in
// WENODerivative.h. That whole layer is implemented for this branch.
//
// psi_prev is the level set as it stood at the START of this geometry step,
// BEFORE any advection. It is used only to define the narrow tube
// |psi_prev| < LSgamma that the derivative stencils are restricted to, so
// that the tube does not move underneath the RK stages. Pass the same
// psi_prev to LSreinitialization().
Field2D solveHJEq(const Field2D &psi, const Field2D &u, const Field2D &v, const Field2D &psi_prev,
                  const Domain &domain, const Variables &variables);
