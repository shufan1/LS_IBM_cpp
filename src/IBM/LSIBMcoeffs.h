#pragma once
#include <vector>
#include "IBMCoeff.h"
#include "../VariableNonDim.h"  // IBM, Domain, LS

// Computes U's and V's IBM ghost/solid-cell coefficients in one call --
// returns a 2-element std::vector<IBMCoeff> (index 0=U, 1=V) instead of
// two separate named structs/output params. P's own (IBM_coeffP in
// LSIBMcoeffs.m) isn't computed here at all: COEFFP.m never reads it
// (confirmed -- no "Immersed Boundary Treating" section there), and
// RHSP.m's own IBM awareness only ever reads flag_u/flag_v (velocity's
// own classification, not a separate flag_p), so nothing in this port
// has a use for it.
//
// Internally calls LSPointIdent() twice, once per grid, with that
// grid's own Robin-BC coefficients: U and V both use ibm.alpha/beta/q
// (with ibm.BQu/BQv respectively).
std::vector<IBMCoeff> LSIBMcoeffs(const IBM &ibm, const Domain &domain, const LS &ls,
                                   const std::vector<Field2D> &phi);
