#pragma once
#include "../VariableNonDim.h"
#include "../Utilities/Field2D.h"

// Mirrors SolveLS/LSInitialize.m: builds the initial level-set field psi
// on the scalar (xp,yp) grid -- signed distance to the immersed
// boundary, negative inside the solid, positive in the fluid (matches
// the psi>=0 fluid-region convention used elsewhere, e.g. masking U/V/phi).
//
// Only LSCase.caseId==1 (grain/circle) is implemented, matching the rest
// of this project's "only geometry=='grain' is real" scope. Throws
// std::runtime_error for any other case.
Field2D computeLSInitialize(const Domain &domain, const LSCase &lsCase);
