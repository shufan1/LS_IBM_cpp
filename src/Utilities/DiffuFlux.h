#pragma once
#include <vector>
#include "../VariableNonDim.h"

// Mirrors Utilities/DiffFlux.m + getDiffFlux.m: precomputed diffusive
// conductance coefficients (face area / (diffusivity * distance)) for
// each of the three equations (U-momentum, V-momentum, scalar
// transport). These depend only on Domain (mesh) and a diffusivity
// constant (Re for momentum, D=1/Pe for the scalar) -- not on the
// current solution -- so they're computed once before the time loop,
// not recomputed per iteration (see main.cpp step 3).
//
// One struct is reused for all three equations rather than three
// separate types: De/Dw/Dn/Ds are always populated; D1_a/D2_a are only
// meaningful for the scalar-transport variant (computeDiffFluxPhi) and
// stay empty otherwise.
struct DiffFluxCoeffs {
    // East/west/north/south diffusive coefficients. Shape matches the
    // owning field: (imax, jmax+1) for U, (imax+1, jmax) for V,
    // (imax+1, jmax+1) for the scalar. Zero outside the interior region
    // where the formula applies (see computeDiffFluxU's comment for
    // exactly which rows/columns and why).
    Field2D De, Dw, Dn, Ds;

    // Scalar-transport only: a one-sided 3-point (non-uniform-grid)
    // correction for the diffusive flux at the west/inlet boundary --
    // needed for 2nd-order accuracy since that cell sits directly
    // against a prescribed Dirichlet face. D1_a augments the coefficient
    // linking to the 2nd interior node; D2_a folds into the boundary
    // source term (see COEFFPHIADRE.m/RHSPHIADRE.m for exactly how).
    // Length jmax+1 (a west-boundary strip, not a full 2D field).
    //
    // NOT included here: MATLAB's DiffFlux.m also computes an
    // `A0_p_p = dV/dt` field in this same branch, but it's dead code for
    // the ADRE pipeline this project follows -- SolveTransportADRE.m
    // reads the real (per-iteration, transport_steady-gated) transient
    // term from a separate Flux.ConvF_phi.A0_p_p instead (computed in
    // ConvFlux.m, not here). Deliberately omitted.
    std::vector<double> D1_a, D2_a;
};


struct Flux {
    DiffFluxCoeffs Diffu_U, Diffu_V, Diffu_Phi;
};

// U-momentum diffusive coefficients (DiffFlux.m's flag==1 branch).
// For interior i in [1, imax-2], j in [1, jmax-1] (0-indexed):
//   De(i,j) = dyv[j-1] / (Re*dxu[i]),   Dw(i,j) = dyv[j-1] / (Re*dxu[i-1])
//   Dn(i,j) = dxv[i]   / (Re*dyu[j]),   Ds(i,j) = dxv[i]   / (Re*dyu[j-1])
// Zero everywhere else:
//   i=0 (west/inlet row): prescribed Dirichlet value (BC.U_a), not solved
//     via this stencil at all.
//   i=imax-1 (east/outlet row): filled by extrapolation from the last
//     interior row (+ mass-flux rescaling), not simultaneously coupled
//     into the diffusion stencil.
//   j=0 / j=jmax (south/north ghost columns): U's y-axis borrows P's
//     padded axis, so these are ghost points filled by an antisymmetric
//     wall-mirror (U_ghost = 2*U_c - U_interior), not real interior
//     locations needing a diffusion coefficient.
// If bc.BC_e_p==1, Dw/Dn/Ds additionally get real values at the east row
// (i=imax-1) -- De's analogous correction is commented out in the
// original MATLAB, so De stays zero there even in that case.
DiffFluxCoeffs computeDiffFluxU(const Domain &domain, const Variables &variables, const BC &bc);

// V-momentum (DiffFlux.m's flag==-1 branch). Same structure as
// computeDiffFluxU with x/y roles swapped -- De/Dw/Dn/Ds sized
// (imax+1, jmax) instead of (imax, jmax+1), and the boundary-exclusion
// reasoning mirrors accordingly (V's north/south walls are literal
// domain-boundary Dirichlet rows needing no interior stencil, while its
// east/west columns are the ones borrowing P's padded axis).
DiffFluxCoeffs computeDiffFluxV(const Domain &domain, const Variables &variables, const BC &bc);

// Scalar transport (DiffFlux.m's flag==0 branch). Like computeDiffFluxU
// but using D=1/Pe instead of Re, sized (imax+1, jmax+1), plus the real
// D1_a/D2_a west-boundary correction (A0_p_p deliberately NOT computed
// -- see DiffFluxCoeffs's comment). One shared diffusivity for every
// species (variables.D), matching IBM.alpha_phi/q_phi's own
// shared-across-species convention -- unlike beta_phi, this genuinely
// isn't per-species here.
DiffFluxCoeffs computeDiffFluxPhi(const Domain &domain, const Variables &variables, const BC &bc);

// Mirrors getDiffFlux.m: computes all three at once.
Flux computeDiffFlux(const Domain &domain, const Variables &variables, const BC &bc);
