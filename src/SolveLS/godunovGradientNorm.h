#pragma once
#include "../Utilities/Field2D.h"
#include "WENODerivative.h"

// Mirrors SolveLS/LSreinitilizationCoeff.m.
//
// Renamed on the way across. The MATLAB name is a misnomer twice over: it
// computes no coefficient of anything -- it computes ||grad psi||, the
// Godunov numerical Hamiltonian -- and it misspells "reinitialization"
// (missing the second 'a'), unlike LSreinitialization.m next to it.
//
// in : grad  -- the four one-sided derivatives at the current iterate,
//               as returned by
//               wenoDerivative(..., LSEquation::ReinitializationEqn)
//      psi_n -- the level set FROZEN at entry to reinitialization, NOT
//               the live iterate. Only its sign is used, to pick which
//               way the characteristics run; freezing it stops nodes near
//               the zero contour flipping branch mid-sweep.
// out: G = ||grad psi||, same shape as psi_n. Exactly 0 outside the tube
//      (all four derivatives zero) and at local minima of psi -- which is
//      what the caller's `fl = (G != 0)` gate keys on.
Field2D godunovGradientNorm(const LSGradient &grad, const Field2D &psi_n);
