#include "LSeqSolve.h"
#include "LSVelocityExtrapolation.h"
#include "solveHJEq.h"
#include "LSreinitialization.h"
#include "LSnormals.h"
#include "../Utilities/Interp.h"
#include "../Utilities/Timer.h"
#include <petsc.h>
#include <algorithm>
#include <cmath>
#include <utility>

// Category names match MATLAB's Timer fields in LSeqSolve.m exactly
// (ibm_velocity_extrapolation / ls_hj_solve / ls_reinitialization /
// ls_normals) so the two profiles can be diffed directly. MATLAB
// baseline over 51 calls: 2.23 / 0.65 / 15.25 / 0.10 seconds.

void LSeqSolve(LS &ls, const StateVar &stateVar, const Domain &domain, const Variables &variables ){
    // psi as it stands BEFORE this step's advection. Stages 2 and 3 both
    // take it to pin down the narrow tube their WENO stencils operate in,
    // so the tube can't drift underneath the RK stages.
    const Field2D psi_prev = ls.psi;

    // 1.compute interface velocity u_gamma from reactive flux
    {
        ScopedTimer t("ibm_velocity_extrapolation");
        LSVelocityExtrapolation(ls, stateVar, domain, variables);
    }

    // Interface CFL check. MATLAB has this (LSVelocityExtrapolation.m:152)
    // but it is dead there -- it reads uI/vI, which are allocated as zeros
    // and never written, so it always reports 0. Restored here against the
    // real ls.u/ls.v, because it is the only guard against the interface
    // moving more than a cell per level-set step, which is what breaks the
    // Hamilton-Jacobi scheme.
    const double dxMin = *std::min_element(domain.dxp.begin(), domain.dxp.end());
    const double dyMin = *std::min_element(domain.dyp.begin(), domain.dyp.end());
    double maxU = 0.0, maxV = 0.0;
    int nBand = 0;
    for (size_t n = 0; n < ls.u.data().size(); ++n) {
        const double a = std::abs(ls.u.data()[n]), b = std::abs(ls.v.data()[n]);
        maxU = std::max(maxU, a);
        maxV = std::max(maxV, b);
        if (a > 0.0 || b > 0.0) ++nBand;
    }
    const double travelU = maxU * variables.dt / dxMin;
    const double travelV = maxV * variables.dt / dyMin;
    const double travel = std::max(travelU, travelV);

    if (variables.verbose) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "  [LS] band=%d  max|u_gamma|=%.4e  travel=%.4f cells\n", nBand,
                    std::max(maxU, maxV), travel);
    }
    // The CFL warning is NOT gated -- a run that silently violates it
    // produces garbage geometry, so it has to survive verbose=false.
    if (travel > 1.0) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "  [LS] WARNING: interface moves %.2f cells this step (>1) -- the HJ "
                    "scheme's tube and stencils assume <1; reduce dt or raise nLSupdate\n",
                    travel);
    }

    // 2.advect levelset function psi by interface velocity u_gamma
    {
        ScopedTimer t("ls_hj_solve");
        ls.psi = solveHJEq(ls.psi, ls.u, ls.v, psi_prev, domain, variables);
    }

    // 3.reinitialize LS to a signed distance function
    {
        ScopedTimer t("ls_reinitialization");
        ls.psi = LSreinitialization(ls.psi, psi_prev, domain, variables);
    }

    // 4.compute new normal vector on the p grid
    // Returns rather than writing into ls, unlike the stages above:
    // computeLSNormals() is shared with LSPointIdent.cpp, which calls it on
    // a psi averaged onto the U/V grid and has no LS to write back to.
    {
        ScopedTimer t("ls_normals");
        LSNormals normals = computeLSNormals(ls.psi, domain);
        ls.nx = std::move(normals.nx);
        ls.ny = std::move(normals.ny);
    }

    // How far the interface actually moved, and how much solid is left.
    // psi<0 cell count is the discrete grain area in cells -- the whole
    // point of the simulation, and the cheapest possible check that
    // dissolution is proceeding at a sane rate.
    int nSolid = 0;
    double maxShift = 0.0;
    for (size_t n = 0; n < ls.psi.data().size(); ++n) {
        if (ls.psi.data()[n] < 0.0) ++nSolid;
        maxShift = std::max(maxShift, std::abs(ls.psi.data()[n] - psi_prev.data()[n]));
    }
    if (variables.verbose) {
        PetscPrintf(PETSC_COMM_WORLD, "  [LS] solid cells=%d  max|dpsi|=%.4e (%.4f cells)\n",
                    nSolid, maxShift, maxShift / dxMin);
    }

    // 5.refresh the U/V-grid fluid masks on the moved geometry
    // Same definition getLS() uses -- psi bilinearly interpolated onto
    // (xu,yu) and (xv,yv). NOT the 2-point average LSPointIdent computes
    // internally for the IBM; those stay two distinct quantities, as in
    // MATLAB (setUpVariablesNonDim.m:367 interp2 vs LSPointIdent.m:78).
    ls.psiU = Field2D(domain.imax, domain.jmax + 1);
    for (int i = 0; i < domain.imax; ++i)
        for (int j = 0; j < domain.jmax + 1; ++j)
            ls.psiU(i, j) = bilinearInterp(domain.xp, domain.yp, ls.psi, domain.xu[i], domain.yu[j]);

    ls.psiV = Field2D(domain.imax + 1, domain.jmax);
    for (int i = 0; i < domain.imax + 1; ++i)
        for (int j = 0; j < domain.jmax; ++j)
            ls.psiV(i, j) = bilinearInterp(domain.xp, domain.yp, ls.psi, domain.xv[i], domain.yv[j]);
}
