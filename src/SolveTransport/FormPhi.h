#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"

// TODO: unpacks one species' solved KSP vector into a full
// (imax+1, jmax+1) field + boundary conditions -- west/east/north/south
// per bc.BC_w_phi etc, then overwrites solid cells (flag==2 only, NOT
// ghost -- different masking convention than SolveUVP.m's fU/fV) to
// phiInside. Not ported yet.
Field2D formPhi(const Domain &domain, Vec phiVec, const BC &bc, const Field2D &phiA,
                const Field2D &flag, double phiInside);
