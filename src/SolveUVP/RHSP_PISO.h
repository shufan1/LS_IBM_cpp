#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "COEFFU.h"  // MomentumCoeffs

// Fills RHS_P2 for PISO's second pressure-correction pass: not the raw
// mass imbalance (rhsP(), RHSP.h) -- that's already ~zero after the
// first corrector -- but the divergence of the "neighbor coupling" term
// the first corrector's diagonal-only velocity correction ignored.
//
// For each momentum cell, this term is:
//   corrU(i,j) = -sysU.d_piso(i,j) * ( sysU.aw(i,j)*deltaU(i-1,j) +
//                sysU.ae(i,j)*deltaU(i+1,j) + sysU.as(i,j)*deltaU(i,j-1) +
//                sysU.an(i,j)*deltaU(i,j+1) )
//  The pressure cell's RHS is then this term's own divergence, same face convention as rhsP():
//   S0(i,j) = corrU(i-1,j) - corrU(i,j) + corrV(i,j-1) - corrV(i,j)

// corrU/corrV are zero outside sysU.d_piso/sysV.d_piso's own valid range
// (matching coeffU()/coeffV()'s interior), which naturally reproduces
// RHSP_PISO.m's own boundary handling (its Su1/Su2/Sv1/Sv2 insertion
// ranges miss the same cells, defaulting to the same zero) without
// needing to special-case the P-domain edges here.
//
//
// Only RHS_P2 (reduced, same pinned-cell exclusion as coeffP()'s CM2/
// rhsP()'s RHS_P2) is computed -- RHSP_PISO.m's RHSP (full) is dead for
// the same reason RHSP.m's own RHSP/coeffP()'s CM are.
//
// RHS_P2 is caller-allocated and filled in place via VecSetValue -- same
// discipline as rhsP()/coeffP() (see COEFFU.h).
void rhsPPiso(const MomentumCoeffs &sysU, const MomentumCoeffs &sysV,
              const Field2D &deltaU, const Field2D &deltaV, const Domain &domain,
              const ControlVar &controlVar, Vec RHS_P2);
