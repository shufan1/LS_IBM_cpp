#pragma once
#include <vector>
#include "../VariableNonDim.h"

// TODO: mirrors calculate_qG.m -- the reactive boundary flux at each
// ghost cell, for every species at once. Genuinely cross-species: for
// each ghost cell, interpolates every species' phi out to two
// extrapolation points along the level-set normal, then applies
// variables.A/inv_A (Np x Np reaction-rate matrix) to get each
// species' q_G -- NOT splittable into a per-species call, since A's
// off-diagonal terms mix species. Not ported yet.
//
// Called once per QUICK iteration (SolveTransportADRE.cpp), reading
// stateVar.phi as it stood at the end of the previous iteration (or the
// initial phi, on the first pass) -- see SolveTransportADRE.cpp's
// comment for why this placement is equivalent to MATLAB's own
// before-the-loop-plus-end-of-iteration calls.
//
// Returns one q_G value per (ghost cell, species): outer index =
// species, inner index parallels ibmCoeffPhi[0].I_g/J_g (ghost-cell
// geometry is identical across species, so any one species' I_g/J_g
// works as the shared index).
std::vector<std::vector<double>> calculateQG(const Variables &variables,
                                              const std::vector<Field2D> &phi, const Domain &domain,
                                              const LS &ls, const std::vector<int> &I_g,
                                              const std::vector<int> &J_g);
