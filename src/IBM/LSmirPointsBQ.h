#pragma once
#include <vector>
#include "IBMCoeff.h"
#include "../Utilities/Field2D.h"

// Computes each ghost cell's mirror-point interpolation stencil weights
// (lambda_g_1..6), RHS contribution (A1_g), and stencil neighbor indices
// (I1..I6/J1..J6, I_e/J_e), enforcing the Robin boundary condition
// -alpha*dphi/dn - beta*phi = q at the true interface.

// x/y: this grid's own 1D axis arrays. 
// alpha/beta/q: this grid's
// Robin-BC coefficients. BQ: bilinear(0)/biquadratic(1) stencil choice.
// dx: this grid's own minimum spacing used to build Delta=sqrt(2)*dx, distance to get mirror point. 
// X_g/Y_g: each ghost cell's physical coordinates (from LSPointIdent()).
// I_g/J_g: ghost cell row and colunmn idx (from LSPointIdent()). 
// psi/nx/ny: this grid's own level-set field and surface normal (from LSPointIdent()).
//
// Mutates ibm_coeff in place: update ibm_coeff.lambda_g_1, lambda_g_2, lambda_g_3,...
//                                    ibm_coeff.I1,J1, I2,J2, ...
//                                    ibm_coeff.A1_g (unless computeA1g is false)
//
// computeA1g: set false when the caller knows this call's A1_g will never
// be read -- e.g. LSIBMcoeffsPhi(), called once per level-set update to
// get the (fixed, under the linear-reaction model) geometry/lambda_g_k
// only; the real, current-concentration-dependent A1_g is instead
// rebuilt every QUICK iteration by computeA1gPhi(), which this call's q
// argument has no bearing on either way. `q` itself is otherwise unused
// when this is false -- callers don't need a real q_G value at all (and,
// for the ghost-cell-classification call in particular, can't produce
// one anyway: computing a real q_G needs I_g/J_g, which don't exist
// until this same call has run). Defaults to true so U/V's own call
// (LSIBMcoeffsUV(), whose A1_g IS the real, final value COEFFU.cpp/
// COEFFV.cpp read directly -- no per-iteration recompute exists for
// them) needs no signature change.
//
// applyEdgeSafeguard/x_0/lx: LSmirPointsBQnew.m's own near-domain-edge
// check (phi/scalar grid only -- U/V's LSmirPointsBQ.m has no
// equivalent) -- forces this ghost cell's *effective* beta to 0 when its
// mirror point lands within 5 cells of the y-domain edges or within x_0
// of the x-domain edges, before B/E/A1_g are computed, so the safeguard
// actually reaches lambda_g_k (not just A1_g). Written to
// ibm_coeff.betaG (see its own comment) so a later, lighter-weight A1_g-
// only recompute (computeA1gPhi(), which can't redo this ghost cell's
// mirror-point search itself) can reuse the same decision. Defaults to
// disabled (x_0/lx unused) so U/V's call needs no signature change.
void LSmirPointsBQ(const std::vector<double> &x, const std::vector<double> &y, double alpha,
                    double beta, double q, const std::vector<double> &X_g,
                    const std::vector<double> &Y_g, int BQ, double dx, const Field2D &psi,
                    const Field2D &nx, const Field2D &ny, IBMCoeff &ibm_coeff,
                    bool computeA1g = true, bool applyEdgeSafeguard = false, double x_0 = 0.0,
                    double lx = 0.0);
