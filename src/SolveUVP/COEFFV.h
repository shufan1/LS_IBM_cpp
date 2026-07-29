#pragma once
#include "COEFFU.h"  // reuses MomentumCoeffs -- same shape, V's own values
#include "../IBM/IBMCoeff.h"

// Builds the V-momentum coefficients and RHS source term -- same overall
// structure as coeffU() (see COEFFU.h): combines computeDiffFluxV with
// computeConvFluxV for every interior V cell, applies boundary
// corrections, forms ap/d_SIMPLC/d_piso and the source term, then assembles
// CM/RHS.
//
// NOT a simple x/y-relabeling of coeffU(): COEFFV.m's own formulas
// differ from COEFFU.m's in several specific ways (ported faithfully,
// not "symmetrized"):
//   - aw/ae use CoNSu (not CoEWv) -- moving east/west from a V-node
//     lands between U-grid columns, so the interpolation weight has to
//     come from U's own spacing, not V's. Symmetric to how coeffU's
//     an/as use CoEWv (moving north/south from a U-node lands between
//     V-grid rows).
//   - unlike coeffU's aw (which uses CoEWu[i-1], one less than its own
//     row), coeffV's aw uses CoNSu[i] -- its own row, not shifted.
//   - west/east boundary corrections carry a factor of 2 (matching
//     coeffU's south/north); south/north don't (matching coeffU's
//     west/east) -- because for V, south/north are its own real
//     boundary rows, while west/east border the ghost-padded U-axis
//     (the axis needing the factor-2 mirror correction), the opposite
//     of U's situation.
//
// Only disc_scheme==2 (central difference) is implemented; disc_scheme
// ==1 or 3 throws std::runtime_error -- same rationale as coeffU() (the
// upwind branch is dead code in practice, since disc_scheme_vel is
// hardcoded to 2 in setUpControlVar.m and never overridden).
//
// The immersed-boundary treatment (COEFFV.m's "Immersed Boundary
// Treating" section) is implemented for both cell types, BC_e_p!=1 only
// -- same rationale as coeffU() (see COEFFU.h): the I_solid loop zeroes
// ae/aw/an/as, pins ap=1/d_SIMPLC=0/d_piso=0, and sets S to
// ibm.u_inside_psi (COEFFV.m reuses this same constant, not a separate
// V-specific one) at every fully-interior solid cell; the I_g
// (ghost-cell) loop does the same zeroing/pinning plus sets S to
// ibmCoeffV.A1_g, and couples each ghost row to its (up to 6)
// mirror-point stencil neighbors via extra -lambda_g_k off-band CM
// entries (MATLAB's A_g_sparse), added directly once k(i,j) exists.
//
// CM/RHS are caller-allocated and filled here in place (MatZeroEntries +
// MatSetValue) -- see COEFFU.h's comment on coeffU() for why.
MomentumCoeffs coeffV(const StateVar &stateVar, const Field2D &V_star_old,
                       const Domain &domain, const DiffFluxCoeffs &diffFluxV,
                       const ConvFluxCoeffs &convFluxV, const IBM &ibm,
                       const IBMCoeff &ibmCoeffV, const Variables &variables, const BC &bc,
                       int disc_scheme, Mat CM, Vec RHS);
