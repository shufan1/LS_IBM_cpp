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
//                                    ibm_coeff.A1_g

void LSmirPointsBQ(const std::vector<double> &x, const std::vector<double> &y, double alpha,
                    double beta, double q, const std::vector<double> &X_g,
                    const std::vector<double> &Y_g, int BQ, double dx, const Field2D &psi,
                    const Field2D &nx, const Field2D &ny, IBMCoeff &ibm_coeff);
