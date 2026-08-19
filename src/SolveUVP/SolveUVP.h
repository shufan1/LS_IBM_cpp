#pragma once
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "../Utilities/DiffFlux.h"
#include "COEFFU.h"  // MomentumCoeffs

// Runs the SIMPLE (and, if ControlVar.PISO==1, PISO) pressure-velocity
// coupling loop for one timestep: each pass computes the convective
// fluxes, builds and solves the U and V momentum systems (coeffU/coeffV,
// KSPSolve), applies boundary conditions (formUV), builds and solves the
// pressure-correction system, corrects U/V/P (newUVP), and repeats until
// the residual drops below tolerance or the iteration cap (100 passes)
// is hit. PISO mode always runs exactly two pressure-correction passes
// per call, then stops, instead of iterating to convergence.
//
// The momentum and pressure-correction systems are both assembled and
// solved every iteration now (coeffU/coeffV/coeffP, rhsP/rhsPPiso,
// KSPSolve), and convergenceResiduals() drives the loop's own resi>tol
// exit test every pass. Still missing: the immersed-boundary ghost-cell
// treatment inside rhsP/rhsPPiso, and any at all inside coeffP (COEFFP.m
// has no "Immersed Boundary Treating" section to port) -- coeffU()/
// coeffV()'s own solid- and ghost-cell treatment is real, see COEFFU.h/
// COEFFV.h -- see SolveUVP.cpp for exactly which pieces of each loop
// pass are real.
//
// `flux` holds the precomputed (solution-independent) diffusive
// coefficients from computeDiffFlux(), called once before the time loop.
// `ibmCoeffU`/`ibmCoeffV` are U's/V's own IBMCoeff (LSIBMcoeffs()'s
// indices 0/1), threaded straight through to coeffU()/coeffV(). There's
// no ibmCoeffP parameter -- P's own IBM_coeffP is never computed by
// LSIBMcoeffs() at all (see IBM/LSIBMcoeffs.h for why), and nothing in
// this project's pressure path (coeffP(), rhsP()) reads it.
void solveUVP(ControlVar &controlVar, const Domain &domain,
              const Variables &variables, StateVar &stateVar, const IBM &ibm,
              const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV, const BC &bc,
              const Flux &flux);

// Applies boundary conditions to U_star and V_star, in place, for both
// variables in one call:
//   U: west/inlet is set to the prescribed inflow value; east/outlet is
//      either a zero-gradient copy of the last interior column (rescaled
//      to match the inlet's total mass flow once iter>1) or left
//      untouched if BC_e_p==1; north/south walls are either an
//      antisymmetric mirror (enforcing a prescribed wall velocity) or a
//      zero-gradient copy, depending on BC_n_u/BC_s_u; two corner points
//      are set last, overriding whatever the wall formulas wrote there.
//   V: north/south walls are either a prescribed value or a zero-gradient
//      copy, depending on BC_n_v/BC_s_v; west/east edges are always an
//      antisymmetric mirror / zero-gradient copy, with no BC-type check.
//
// This function does not compute U_star/V_star's interior values itself
// -- callers (solveUVP()) fill those in first, from the momentum KSP
// solve. Missing: a multiplicative mask on U's west boundary value that
// should come from the immersed-boundary setup (not implemented).
void formUV(Field2D &U_star, Field2D &V_star, const BC &bc, int iter);

// Corrects U_star/V_star toward a divergence-free velocity field using
// the pressure-correction field PCOR's gradient (scaled by the
// velocity-correction sensitivity coefficients d_u/d_v), then applies
// U/V/P's boundary conditions against the just-corrected interior
// (NEWUVP.m:63-153 for U, 155-195 for V, 202-213 for P -- structurally
// similar to formUV()'s own U/V boundary treatment and a few genuinely
// different ranges are called out inline, not symmetrized), and updates
// the pressure field by adding the under-relaxed correction alpha_p*PCOR
// plus its own boundary. Overwrites U_star/V_star/P_star in place -- the
// caller is responsible for syncing the result into stateVar.U/V
// (P_star is typically stateVar.P itself already, passed in directly,
// so no separate sync is needed for it). `bc`/`iter` are only needed for
// this boundary block (`iter` gates the outlet mass-rescale exactly like
// formUV()'s own `iter>1` check).
//
// d_u/d_v come from coeffU()/coeffV()'s own d_SIMPLC field. Missing: the
// IBM mask on U's west boundary value, and the outlet-row correction
// that applies when BC_e_p==1 (currently that row is only ever touched
// by formUV()'s extrapolation, never by this function) -- neither is
// implemented (IBM isn't ported; BC_e_p==1 is dead for the active
// config).
void newUVP(Field2D &U_star, Field2D &V_star, Field2D &P_star,
            const Field2D &PCOR, const Field2D &d_u, const Field2D &d_v,
            double alpha_p, const BC &bc, int iter);

// PISO's second-corrector version of newUVP(): applies the exact same
// correction newUVP() does (PCOR gradient for U/V, U/V/P boundary
// conditions, alpha_p*PCOR for P -- implemented by calling newUVP()
// directly), then layers on the extra neighbor-coupling term NEWUVP.m's
// own `if PISO==1` branch adds on top:
//   U_star(i,j) -= (sysU.aw(i,j)*deltaU(i-1,j) + sysU.ae(i,j)*deltaU(i+1,j) +
//                   sysU.as(i,j)*deltaU(i,j-1) + sysU.an(i,j)*deltaU(i,j+1)) / sysU.ap(i,j);
// (deltaV/V mirrors this with sysV). Divides by the raw diagonal ap, not
// d_SIMPLC/d_piso -- this is the same neighbor-coupling quantity
// rhsPPiso() uses for the pressure RHS, but without the dyv/dxu area
// scaling that's needed there and not here (see RHSP_PISO.h).
//
// deltaU/deltaV are the velocity change from this iteration's first
// correction alone (this iteration's stateVar.U/V minus their value
// right before that correction) -- the same deltaU/deltaV passed to
// rhsPPiso(). `bc`/`iter` are threaded straight through to the internal
// newUVP() call for its boundary block.
void newUVPPiso(Field2D &U_star, Field2D &V_star, Field2D &P_star, const Field2D &PCOR,
                const MomentumCoeffs &sysU, const MomentumCoeffs &sysV,
                const Field2D &deltaU, const Field2D &deltaV, double alpha_p,
                const BC &bc, int iter);

// Pins one interior grid point of the pressure-correction field PCOR to
// zero, and applies zero-gradient boundary conditions on all four domain
// edges (unconditionally -- no BC-type check on any edge). The pin
// matters because the pressure-correction equation has Neumann
// conditions on every edge, which leaves the system singular -- solvable
// only up to an arbitrary additive constant -- unless one point's value
// is fixed. The pinned point is always the interior column closest to
// the outlet, at a row chosen by ControlVar.imposePresBC ("top"/
// "bottom"/"middle").
//
// This function does not compute PCOR's interior values itself -- the
// caller (solveUVP(), via its solvePressureSystem() helper) fills those
// in first, from the pressure KSP solve.
void formPCor(Field2D &PCOR, const Domain &domain, const ControlVar &controlVar);

// The reference pressure-correction cell (0-indexed): i0 is always the
// interior column closest to the outlet (imax-1); j0 depends on
// imposePresBC ("top" = closest to the north wall, "bottom" = closest to
// the south wall, "middle" = vertically centered). Shared by formPCor()
// (pins PCOR there directly, after solving the full system) and
// coeffP() (excludes this cell from the linear system entirely instead)
// -- both enforce the same physical reference point via different
// mechanisms, so they must agree on which cell that is.
struct PressurePin {
    int i0, j0;
};
PressurePin pinnedPressureCell(int imax, int jmax, const std::string &imposePresBC);
