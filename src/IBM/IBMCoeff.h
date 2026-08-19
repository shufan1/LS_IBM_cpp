#pragma once
#include <string>
#include <vector>
#include "../Utilities/Field2D.h"

// One grid's worth of immersed-boundary coefficients -- what
// LSPointIdent()/LSIBMcoeffs.m compute for a single grid (U, V, or
// P/scalar). Mirrors IBM_coeffU/IBM_coeffV/IBM_coeffP's fields in
// LSIBMcoeffs.m, with two exceptions: the raw `landa` matrix (numg x
// landanu, LSmirPointsBQ.m's own internal combined form) isn't kept --
// only its already-split lambda_g_1..6 columns are, since nothing
// downstream reads the combined matrix (confirmed in COEFFU.m: it
// destructures IBM_coeffU.landa_g_1_u.. directly, never
// IBM_coeffU.landa_u); and LSPointIdentnew.m's per-species beta_G/q_G
// sampling isn't kept either, since it's only an intermediate used to
// build landa_g_1..6/A1_g for each species, not read again afterward.
struct IBMCoeff {
    // Every node's classification, matching LSPointIdent.m's own
    // convention: 0=fluid, 1=ghost, 2=solid. Shape matches whichever
    // grid (U/V/P) this instance belongs to. This is what
    // coeffU()/coeffV()/rhsP() actually branch on to decide whether a
    // cell needs any special treatment at all.

    Field2D flag;

    // ---- Solid cells (flag==2) -- cells deep enough inside the grain ----
    // record the row and column idx of solid cells
    // coeffU()/coeffV() would loop these and zero the off-diagonal coefficients,
    // pin ap=1, and set the source to u_inside_psi/phi_inside_psi -- a plain Dirichlet pin,
    std::vector<int> I_solid, J_solid;

    // ---- Ghost cells (flag==1) -- the fluid cells right next to the
    // boundary, whose normal 5-point stencil would reach across it.
    // ecord the row and column idx of  ghost cells
    std::vector<int> I_g, J_g;

    // --- Mirror cells
    // row and index of the mirror point of the ghost point. where the 4/6 neighbors based off
    std::vector<int> I_m, J_m;

    // The 4 grid neighbors of the *mirror point* (the point projected
    // 2*Delta outward along the normal, across the boundary) used to
    // build the bilinear interpolation stencil, when BQ==0.
    std::vector<int> I1, J1, I2, J2, I3, J3, I4, J4;

    // Two extra neighbors used only when BQ==1 (biquadratic, 6-point
    // stencil instead of bilinear 4-point).
    std::vector<int> I5, J5, I6, J6;

    // The interpolation/extrapolation weight for each of I1..I4's
    // neighbor -- the coefficient that neighbor's field value gets
    // multiplied by when reconstructing the ghost cell's value. These
    // are exactly what get spliced into the coefficient matrix as
    // A_g_sparse: the ghost row's equation becomes (roughly)
    // phi_ghost = lambda_g_1*phi(I1,J1) + ... + lambda_g_4*phi(I4,J4) + A1_g,
    // instead of the usual 5-point stencil.
    std::vector<double> lambda_g_1, lambda_g_2, lambda_g_3, lambda_g_4;

    // The matching weights for I5/J5, I6/J6 -- only populated when BQ==1.
    std::vector<double> lambda_g_5, lambda_g_6;

    // The ghost cell's RHS contribution (the q-derived forcing term
    // baked into the extrapolation formula), parallel to I_g/J_g.
    std::vector<double> A1_g;

    // The Robin-BC beta actually used for this ghost cell's B/E/A1_g,
    // parallel to I_g/J_g -- normally just the caller's own beta
    // (unchanged), but forced to 0 for ghost cells whose mirror point
    // landed too close to the domain's outer edge (LSmirPointsBQnew.m's
    // near-edge safeguard, phi/scalar grid only -- see LSmirPointsBQ.h).
    // Written once here, alongside lambda_g_k/A1_g; computeA1gPhi()
    // reads it back per QUICK iteration instead of ibm.beta_phi[i_s]
    // directly, since it can't redo the mirror-point search itself to
    // re-derive this.
    std::vector<double> betaG;

    // Each ghost cell's own distance to the true interface (abs(psi) at
    // its own (I_g,J_g)), parallel to I_g/J_g -- static per-ghost-cell
    // geometry, feeds A1_g's formula alongside betaG. Cached here for
    // the same reason as betaG: computeA1gPhi() can't redo any of
    // LSmirPointsBQ()'s own per-ghost-cell work.
    std::vector<double> r_g;

    // Delta = sqrt(2)*dx for this grid, the same single constant for
    // every ghost cell in this IBMCoeff -- cached here (rather than
    // recomputed from Domain) purely so computeA1gPhi() doesn't need a
    // Domain reference just for this one scalar.
    double Delta = 0.0;

    // == I_g.size(), kept as its own field only because LSPointIdent.m
    // returns it as an explicit output too, rather than making every
    // caller call .size().
    int numg = 0;

    // A human-readable status string (min/max of the ghost-cell-
    // extrapolated concentration, phi_IB) -- only ever populated for
    // the P/scalar grid call (LSIBMcoeffs.m disp()'s it right after
    // computing it). Empty for U/V.
    std::string message;
};
