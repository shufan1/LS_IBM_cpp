#pragma once
#include <vector>
#include "Field2D.h"

// C++ analog of MATLAB's interp2 (linear), for the specific way it's used
// in setUpVariablesNonDim.m: interpolate a field V sampled on a
// tensor-product grid onto an arbitrary query point.
//
// V is a Field2D where V(i,j) is the sample at physical location
// (x[i], y[j]); x (length V.nx()) and y (length V.ny()) are the 1D axis
// vectors, each strictly ascending. Returns the bilinearly-interpolated
// value at (xq, yq).
//
// Note the argument order is physical/intuitive (x first, then y) --
// unlike MATLAB interp2(X,Y,V,Xq,Yq), whose first coordinate argument is
// the *column* axis. setUpVariablesNonDim.m's psiU/psiV calls feed the
// yp-grid as interp2's X and the xp-grid as its Y to compensate; here you
// just pass (xp, yp, psi, xu[i], yu[j]) directly.
//
// Out-of-range queries are clamped to the nearest edge cell (constant
// extrapolation), NOT NaN like MATLAB's default -- for the grain/uniform
// case every psiU/psiV query lands inside the psi grid, so this only
// guards against floating-point overshoot at the boundary.
double bilinearInterp(const std::vector<double> &x, const std::vector<double> &y,
                      const Field2D &V, double xq, double yq);
