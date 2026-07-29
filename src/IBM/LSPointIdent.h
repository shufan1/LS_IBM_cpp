#pragma once
#include <vector>
#include "IBMCoeff.h"
#include "../VariableNonDim.h"  // Domain, LS

// Classifies every node of one grid (U, V, or P/scalar -- selected by
// UVP) into fluid/ghost/solid against the level-set field ls.psi, then
// (once implemented) builds each ghost cell's mirror-point
// interpolation stencil via LSmirPointsBQ/LSmirPointsBQnew, enforcing
// the Robin boundary condition -alpha*dphi/dn - beta*phi = q at the
// true interface.
//
// UVP selects which grid: -1 = U, 0 = V, 1 = P/scalar -- mirrors
// LSPointIdent.m's own UVP convention exactly. LSPointIdentnew.m's
// separate P/scalar MATLAB file collapses into this same UVP==1 branch
// here (see LSIBMcoeffs.h for why).
//
// Returns one grid's worth of coefficients as a single IBMCoeff.
// LSPointIdentnew.m's real per-species output (one stencil-weight set
// per chemical species, since each has its own Robin-BC strength
// sourced from ls.q_out/ls.beta_out) isn't represented here -- that
// distinction is deferred until this is actually implemented; for now
// LSIBMcoeffs() calls this once per grid (see LSIBMcoeffs.h), not once
// per species.
//
// alpha/beta/q/BQ: that grid's own Robin-BC coefficients and
// bilinear(0)/biquadratic(1) stencil choice -- LSIBMcoeffs() picks
// these per grid before calling here.
// treshold: LSPointIdent.m's own numerical tolerance parameter
// (IBM.treshold).
// phi: species concentration fields (StateVar.phi) -- LSPointIdentnew.m
// needs these (their count and values) for the P/scalar grid's
// per-species Robin BC; unused by the current placeholder (see
// LSPointIdent.cpp).
//
// Classification (fluid/ghost/solid) and the ghost-cell coordinate
// extraction are real; the ghost-cell mirror-point stencil itself
// (LSmirPointsBQ(), called at the end) is still a placeholder -- see
// LSmirPointsBQ.h/.cpp. LSmirPointsBQnew.m's separate P/scalar variant
// isn't ported either.
IBMCoeff LSPointIdent(const Domain &domain, double alpha, double beta, double q, int BQ,
                       const LS &ls, int UVP, const std::vector<Field2D> &phi, double treshold);
