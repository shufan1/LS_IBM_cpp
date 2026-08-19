#pragma once
#include "../VariableNonDim.h"

// Mirrors SolveLS/LSeqSolve.m (single-mineral `else` branch, lines 71-90).
// Called once per geometry update from main.cpp's time loop -- BEFORE the
// flow/transport solve, so the StateVar.phi it reads is the field left at
// the end of the PREVIOUS big step. The concentration -> geometry coupling
// is first-order explicit operator splitting.
//
// Updates `ls` in place and returns nothing, matching MATLAB's
// `[LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN)`. Everything the step
// produces lands on the LS struct:
//
//   ls.psi            advected, then reinitialized to a signed distance
//   ls.nx, ls.ny      recomputed on the new psi
//   ls.u, ls.v        the interface velocity that did the advecting, kept
//                     for diagnostics and the interface CFL check
//
// MATLAB's LSeqSolve also leaves LS.q_out/LS.beta_out behind for the
// scalar IBM. Those are not ported -- SolveTransportADRE recomputes both
// itself and overwrites the coefficients they seed before they are used.
// See the note in struct LS (VariableNonDim.h).
//
// Four stages, strictly sequential -- each one's output is the next
// one's input:
//
//   1. computeLSVelocityExtrapolation  chemistry -> u/v on the band
//   2. solveHJEq                       advect psi with u/v
//   3. LSreinitialization              restore |grad psi| = 1
//   4. computeLSNormals                recompute nx/ny on the new psi
//
// MATLAB times these four separately (Timer.ibm_velocity_extrapolation /
// ls_hj_solve / ls_reinitialization / ls_normals); keep that split here so
// the C++ profile stays comparable against the MATLAB baseline of
// 2.23 / 0.65 / 15.25 / 0.10 seconds over 51 calls. Stage 3 is ~80% of it.
//
// NOT YET WIRED INTO CMakeLists.txt -- stages 1-3 are declarations only,
// so adding this to the ls_ibm target would fail to link.
void LSeqSolve(LS &ls, const StateVar &stateVar, const Domain &domain,
               const Variables &variables);
