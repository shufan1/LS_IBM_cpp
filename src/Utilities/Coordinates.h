#pragma once
#include <string>
#include "../VariableNonDim.h"

// Mirrors beta_AD_correctu/Utilities/Coordinates.m: builds the staggered
// MAC grid point coordinates (xu/yu/xv/yv/xp/yp) and cell spacings
// (dxu/dyu/dxv/dyv/dxp/dyp), for either a uniform or a graded/stretched
// mesh (refined near the immersed boundary), matching Grid.uniform.
//
// NOT ported here (only needed once the modules that consume them are
// ported): the QUICK advection-interpolation coefficients (g1/g2, 48
// arrays for u/v/transport x {e,w,n,s} x {positive,negative flow}), the
// central-difference interpolation weights (CoEWu/CoNSu/CoEWv/CoNSv/
// CoEWp/CoNSp), and the finite-volume cell volumes (dV_u/dV_v/dV_p).
// Those are consumed only by SolveUVP/COEFFU/COEFFV/COEFFP and
// SolveTransportADRE, none of which are ported yet -- add them there
// when that happens, following this same file's translation pattern.
// The full-resolution meshgrids (Xu/Yu/Xv/Yv/Xp/Yp) aren't stored either
// -- they're just the tensor product of the 1D axis vectors below, so
// any plotting tool can reconstruct them from xu/yu etc. directly.
Domain computeCoordinates(double lx, double ly, double diamcyl, const Grid &grid);

// Dumps domain's coordinate axis vectors (xu/yu/xv/yv/xp/yp, plus
// imax/jmax/lx/ly) to a JSON file for later plotting. psi and phi live on
// the same nodes as p (xp/yp), so no separate coordinate arrays are
// needed for them -- just read xp/yp when plotting those fields.
// Call once, from rank 0 only (the grid is identical on every rank until
// domain decomposition exists). Throws std::runtime_error if the file
// can't be written.
void saveCoordinatesJson(const Domain &domain, const std::string &path);
