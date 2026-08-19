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
// bilinear(0)/biquadratic(1) stencil choice -- callers pick these per
// grid (LSIBMcoeffsUV()) or per species (LSIBMcoeffsPhi(), which is
// where the per-species beta actually comes from -- this function
// itself has no notion of species count, since each call is already
// scoped to one grid/species via its own alpha/beta) before calling
// here.
// q is used to build rhs, A1_g, not really used and comptued for coeffPhi 
// because q dependes on concentration, this will actually be computed in
// quick iteration and update A1_g

// treshold: LSPointIdent.m's own numerical tolerance parameter
// (IBM.treshold).
//

//
// Classification (fluid/ghost/solid), ghost-cell coordinate extraction,
// and the mirror-point stencil itself (LSmirPointsBQ(), called at the
// end) are all real. LSmirPointsBQnew.m's separate P/scalar variant
// isn't ported as its own function -- its one real behavioral
// difference (the near-domain-edge beta=0 demotion safeguard) isn't
// replicated here yet either; everything else about it is either
// identical math or dead/unreachable code in one file or the other.
//
// computeA1g forwards straight to LSmirPointsBQ() -- see its own
// comment for when to pass false.
IBMCoeff LSPointIdent(const Domain &domain, double alpha, double beta, double q, int BQ,
                       const LS &ls, int UVP, double treshold, bool computeA1g = true);
