#pragma once
#include "../VariableNonDim.h"
#include "../Utilities/Field2D.h"

// Mirrors SolveLS/LSreinitialization.m. Restores psi to a true signed
// distance function after advection has distorted it, by iterating the
// reinitialization equation to steady state:
//
//     dpsi/dtau = sign(psi) * (1 - |grad psi|)
//
// which drives |grad psi| -> 1 while holding the zero contour (the actual
// interface) in place. Returns the corrected psi.
//
// This matters more than it looks. |grad psi| = 1 is what makes abs(psi) a
// real distance rather than an arbitrary scalar, and two things downstream
// depend on that directly: r_g = abs(psi) in LSmirPointsBQ (the ghost
// cell's distance to the surface, which sets the B/E coefficients), and
// the conditioning of the gradient that computeLSNormals() normalizes.
//
// It is also the single most expensive step in the level-set update -- in
// the MATLAB profile, 15.25 s over 51 calls, roughly 23x the advection it
// is correcting and about 80% of the total LS cost. If any part of this
// folder is worth parallelizing, it is this one.
//
// Iteration count is variables.n_iter_ReLS (4 by default), pseudo-time
// step variables.dtau, time scheme variables.TimeSchemeRLS ("RK3").
//
// Unlike solveHJEq(), this function's RK3 is CORRECT: line 50 evaluates
// the third stage at psi_m12, the stage-2 value, as SSP-RK3 requires.
// solveHJEq.m:43 passes psi_n1 there instead -- see the bug-for-bug note
// in solveHJEq.h. Port this one straight; there is nothing to reproduce.
//
// Two details easy to miss when porting:
//   - The sign function is smoothed, not exact:
//         SignPsi = psi_n / sqrt(psi_n^2 + epsSign^2),  epsSign = min(dxp)
//     which keeps the interface from chattering across zero.
//   - Every update is gated by `fl = double(G ~= 0)`, so nodes where the
//     Godunov flux is exactly zero (outside the tube, where all the
//     one-sided derivatives were left at zero) are frozen rather than
//     pushed by the (G-1) = -1 that would otherwise apply there.
//   - After the iteration loop, psi is clamped to +/-h (lines 63-64,
//     h = LSgamma). Values outside the tube are pinned to the band edge,
//     not left at whatever the advection produced.
//
// Spatial discretization differs from solveHJEq() in one important way.
// There is no velocity to upwind against, so the shared derivative layer
// (equation == "ReinitializationEqn") returns BOTH one-sided derivatives
// in each direction -- psi_xn/psi_xp/psi_yn/psi_yp -- and
// LSreinitilizationCoeff.m combines them with a Godunov flux
//
//     G = sqrt(max(aP^2, bM^2) + max(cP^2, dM^2))    where psi >= 0
//     G = sqrt(max(aM^2, bP^2) + max(cM^2, dP^2))    where psi <  0
//
// picking the entropy-correct branch pointwise by the sign of psi. That
// sign selection replaces the velocity upwinding.
//
// Note MATLAB's signature takes the whole LS struct and pulls u/v out of
// it, but the reinitialization branch never reads them -- confirmed by
// inspection of LSdirDerivates.m. They are dropped here.
//
// psi_prev is the pre-advection level set, same value passed to solveHJEq(),
// and again used only to fix the narrow tube the stencils operate in.
Field2D LSreinitialization(const Field2D &psi, const Field2D &psi_prev, const Domain &domain,
                           const Variables &variables);
