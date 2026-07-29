#pragma once
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../ControlVar.h"

// Fills CM2, the pressure-correction coefficient matrix, built from the
// momentum solves' sensitivity coefficients (d_u from coeffU(), d_v from
// coeffV()). The right-hand side is computed separately, by rhsP(

// d_u = A_iJ / a_iJ, a_iJ coeff in u-momentum
// d_v = A_Ij / a_Ij, a_Ij coeff in v-momentum

//  
// CM2 is caller-allocated, sized (imax-1)*(jmax-1)-1 (excludes the
// reference cell -- pinnedPressureCell(), SolveUVP.h -- entirely, rather
// than solving the full singular system and fixing it up afterward), and
// filled here in place (MatZeroEntries + MatSetValue) -- same reasoning
// as coeffU()/coeffV() (COEFFU.h): the coefficients (via d_u/d_v, which
// trace back through coeffU()/coeffV()'s own convection-dependent ap)
// change every SIMPLE iteration, but CM2's sparsity pattern doesn't, so
// solveUVP() should allocate it once and reuse it, not recreate it every
// iteration.
void coeffP(const Field2D &d_u, const Field2D &d_v, const Domain &domain,
            const ControlVar &controlVar, Mat CM2);
