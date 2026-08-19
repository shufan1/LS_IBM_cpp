#pragma once
#include <vector>
#include "IBMCoeff.h"
#include "../VariableNonDim.h"  // IBM, Domain, LS

// Computes U's and V's IBM ghost/solid-cell coefficients in one call --
// returns a 2-element std::vector<IBMCoeff> (index 0=U, 1=V) instead of
// two separate named structs/output params.
//
// Internally calls LSPointIdent() twice, once per grid, with that
// grid's own Robin-BC coefficients: U and V both use ibm.alpha/beta/q
// (with ibm.BQu/BQv respectively). No phi parameter -- confirmed dead
// for UVP=-1/0 by reading LSmirPointsBQ.m directly: its two
// "check if phi_IB" sections are empty comment stubs, and the one place
// phi-derived data would matter (the phiIB summary message) is gated by
// `if UVP==1`, which this function is never called with (that's
// LSmirPointsBQnew.m's job, reached only through LSIBMcoeffsPhi() below).
//
// Renamed from LSIBMcoeffs() -- see LSIBMcoeffsPhi() below for the
// P/scalar-grid counterpart, now that the transport solver needs it.
std::vector<IBMCoeff> LSIBMcoeffsUV(const IBM &ibm, const Domain &domain, const LS &ls);

// Mirrors LSIBMcoeffs.m's IBM_coeffP output (built via LSPointIdentnew.m
// in MATLAB, collapsed into LSPointIdent.cpp's own UVP==1 branch here --
// see LSPointIdent.h), but NOT a literal mirror of ITS OWN inputs:
// LSIBMcoeffs.m seeds this with a dead single-scalar beta_phi and relies
// on SolveTransportADRE.m's update_A1g to overwrite every species'
// result with the real beta_G before first use -- this instead reads
// ibm.beta_phi[i_s] (the real per-species beta_G, populated by main.cpp
// from -diag(variables.A) once A is loaded) directly, skipping that
// redundant round-trip (see LSIBMcoeffs.cpp's comment and
// IBM::beta_phi's comment in VariableNonDim.h).
//
// One IBMCoeff per species (size Np): I_g/J_g/flag geometry is
// identical across species (purely psi-based), but each has its own
// beta_G, so lambda_g_k differs per species even though the stencil
// corners don't. Classification, ghost-coordinate extraction, and the
// mirror-point stencil (LSPointIdent()/LSmirPointsBQ()) are all real.
// Calls LSPointIdent() with computeA1g=false -- this function is meant
// to be called once per level-set update (not once per QUICK iteration),
// so whatever A1_g it would produce is always stale by the time
// anything reads it; the real, current-concentration-dependent A1_g
// gets rebuilt every QUICK iteration instead (SolveTransportADRE.cpp's
// computeA1gPhi() calls, using a fresh q_G each time). No phi parameter
// needed for the species count (see LSPointIdent.h's own comment for
// why phi isn't threaded through at all right now).
//
// Takes Np directly rather than the whole Variables struct -- that's
// the only field besides ibm/domain/ls this function actually needs
// (ibm.beta_phi already carries what A would have provided).
std::vector<IBMCoeff> LSIBMcoeffsPhi(const IBM &ibm, const Domain &domain, const LS &ls, int Np);
