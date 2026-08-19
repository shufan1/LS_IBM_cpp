#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../Utilities/DiffFlux.h"
#include "../Utilities/ConvFlux.h"
#include "../IBM/IBMCoeff.h"

// Result of building the U-momentum (or, via COEFFV.h's reuse of this
// same type, V-momentum) coefficient matrix. Field names are generic
// (not "_u"/"_v"-suffixed) since the same struct holds either equation's
// values.
//
// Two of MATLAB's COEFFU.m outputs are deliberately NOT here: `M` and
// `A_g_sparse` are computed and returned there, but nothing anywhere in
// the codebase ever reads them back out afterward (confirmed by tracing
// every consumer) -- dead outputs, not ported.
struct MomentumCoeffs {
    // Per-direction coefficients for every cell (west/east/south/north/
    // diagonal), before assembly into the sparse matrix. Needed both to
    // build CM below and, separately, by the PISO second-corrector's own
    // right-hand side (RHSP_PISO.m reads these directly, not through CM).
    Field2D aw, ae, as, an, ap;

    // Right-hand-side source term (S0 + S_Pres + S_ur + S_u in COEFFU.m:
    // old-timestep contribution, pressure gradient, under-relaxation, and
    // boundary-value contributions). Kept as a full-sized Field2D rather
    // than MATLAB's reshaped-and-truncated vector -- only the interior
    // cells (matching aw/ae/as/an/ap's own range) are ever populated.
    Field2D S;

    // Velocity-correction sensitivity coefficient for the SIMPLEC
    // corrector: dyv / (ap+aw+ae+as+an) -- used by newUVP().
    Field2D d_SIMPLC; //SIMPLEC A_ij/(aij-Sum a_nb), AI_ij = dyv for U, AI_ij = dxu for V,

    // Velocity-correction sensitivity coefficient for PISO's second
    // corrector: dyv / ap alone (uses only the diagonal term, not the
    // full coefficient sum) -- used by the PISO second-corrector's RHS,
    // not by newUVP()'s standard correction.
    Field2D d_piso; // PISO  A_ij/(aij), AI_ij = dyv for U, AI_ij = dxu for V, 
};

// Builds the U-momentum coefficients and RHS source term: combines the
// diffusive coefficients (computeDiffFluxU) with the convective fluxes
// (computeConvFluxU) for every interior U cell, applies the boundary
// corrections, and forms ap/d_SIMPLC/d_piso plus the source term.
//
// Only disc_scheme==2 (central difference) is implemented, matching the
// only branch that's live in COEFFU.m for BC_e_p!=1 -- disc_scheme==2's
// own BC_e_p==1 extension is commented-out dead code in the MATLAB
// source itself, so it was never ported here either. disc_scheme==1 or
// 3 (upwind/QUICK) throws std::runtime_error; nothing in this project
// currently sets disc_scheme to anything but 2 (see config.json's
// disc_scheme_vel).
//
// The immersed-boundary treatment (COEFFU.m's "Immersed Boundary
// Treating" section) is implemented for both cell types, BC_e_p!=1
// only (BC_e_p==1 uses a differently-sized linear system this port's
// own k(i,j) bijection doesn't support, matching every other BC_e_p==1
// omission in this project): the I_solid loop zeroes ae/aw/an/as, pins
// ap=1/d_SIMPLC=0/d_piso=0, and sets S to ibm.u_inside_psi at every
// fully-interior solid cell; the I_g (ghost-cell) loop does the same
// zeroing/pinning plus sets S to ibmCoeffU.A1_g, and couples each ghost
// row to its (up to 6) mirror-point stencil neighbors via extra
// -lambda_g_k off-band CM entries (MATLAB's A_g_sparse), added directly
// once k(i,j) exists, since a stencil neighbor can be 2 grid cells away
// -- outside what ae/aw/an/as (the standard 4 orthogonal neighbors)
// can represent.
//
// CM/RHS are caller-allocated (one row per interior U cell, bijective
// index k(i,j)=(i-1)+(j-1)*(imax-2); PETSC_COMM_SELF since there's no
// domain decomposition yet) and filled here in place (MatZeroEntries +
// MatSetValue, not MatCreateSeqAIJ) -- this runs every SIMPLE iteration
// with genuinely different coefficient values each time (Fe/Fw/Fn/Fs
// depend on the current velocity field), but the sparsity pattern never
// changes, so solveUVP() allocates CM/RHS once before its SIMPLE loop
// and destroys them once after, instead of paying malloc/free every
// iteration for the same-shaped memory.
MomentumCoeffs coeffU(const StateVar &stateVar, const Field2D &U_star_old,
                       const Domain &domain, const DiffFluxCoeffs &diffFluxU,
                       const ConvFluxCoeffs &convFluxU, const IBM &ibm,
                       const IBMCoeff &ibmCoeffU, const Variables &variables, const BC &bc,
                       int disc_scheme, Mat CM, Vec RHS);
