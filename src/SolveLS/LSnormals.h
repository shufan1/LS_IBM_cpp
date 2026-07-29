#pragma once
#include "../VariableNonDim.h"
#include "../Utilities/Field2D.h"

// Result of computeLSNormals(): the (nx,ny) unit outward normal to the
// immersed boundary at every node of whatever grid `psi` lives on, same
// shape as `psi`. Zero away from the interface (only computed within a
// narrow band, see computeLSNormals()'s own comment).
struct LSNormals {
    Field2D nx, ny;
};

// Mirrors SolveLS/LSnormals.m. Computes a central-difference gradient of
// `psi` at every strictly-interior node (row/col 0 and the last row/col
// are skipped -- no one-sided stencil at the array's own edge), but
// only where |psi(i,j)| < epsi = 20*dx (a narrow band around the
// interface -- points far from it don't need a normal, since they'll
// never be a ghost cell). The gradient is then normalized to unit
// length; nodes outside the narrow band, or where the raw gradient
// magnitude is numerically zero, are left at Field2D's zero default.
//
// dx/dy (the narrow-band width and derivative step) always come from
// domain.dxp/dyp -- LSnormals.m uses DOMAIN.dxp/dyp regardless of which
// grid `psi` itself is defined on (LSPointIdent.m calls this with `psi`
// already averaged onto the U/V grid, but always passes the original,
// un-staggered DOMAIN through for this).
LSNormals computeLSNormals(const Field2D &psi, const Domain &domain);
