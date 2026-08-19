#pragma once
#include <string>
#include <vector>
#include "../VariableNonDim.h"
#include "../ControlVar.h"
#include "../Utilities/DiffuFlux.h"
#include "../IBM/IBMCoeff.h"

// Mirrors SolveTransport/SolveTransportADRE.m -- one unsteady
// advection-diffusion-reaction step for the scalar field phi (all
// Np species, coupled only through the reactive ghost-cell boundary
// flux -- see CalculateQG.h).
//
// ibmCoeffPhi: one IBMCoeff per species (size variables.Np), mirroring
// LSIBMcoeffs()'s ibmCoeffU/V pattern -- I_g/J_g/flag geometry is
// identical across species (purely psi-based), but lambda_g_k differs
// per species since each has its own Robin-BC beta.
//
// subIdx/outputAdreDir: DEBUG ONLY, for validating the ghost-cell
// reactive boundary machinery (calculateQG()/updateA1gPhi()/lambda_g_k)
// directly against MATLAB, rather than only the final interior phi --
// see compare_output_adre.m. If outputAdreDir is non-empty, dumps
// per-ghost-cell q_G/A1_g/lambda_g_1..6 (keyed by I_g/J_g) from the
// QUICK loop's first pass to
// <outputAdreDir>/dataRDE<controlVar.iTime>_<subIdx>_ghost.json. Pass ""
// to skip (subIdx is unused in that case).
void solveTransportADRE(ControlVar &controlVar, const Domain &domain,
                        const Variables &variables, StateVar &stateVar,
                        const IBM &ibm, std::vector<IBMCoeff> &ibmCoeffPhi,
                        const BC &bc, const LS &ls, const Flux &flux,
                        int subIdx = 0, const std::string &outputAdreDir = "");
