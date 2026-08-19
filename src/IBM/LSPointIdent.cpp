#include "LSPointIdent.h"
#include "LSmirPointsBQ.h"
#include "../SolveLS/LSnormals.h"
#include <algorithm>

namespace {
// psi averaged onto the U grid: LSPointIdent.m's own local average
// (PSI(1:end-1,:)+PSI(2:end,:))/2 across the x-direction -- NOT the
// same computation as ls.psiU (which comes from a full bilinear
// interp2 in getLS()/setUpVariablesNonDim.m). MATLAB keeps these as two
// independently-computed approximations of "psi on the U grid" for two
// different purposes; ported faithfully as two separate quantities
// here too, not merged.
Field2D psiOnU(const Field2D &psi, int imax, int jmax) {
    Field2D result(imax, jmax + 1);
    for (int i = 0; i < imax; ++i)
        for (int j = 0; j <= jmax; ++j) result(i, j) = (psi(i, j) + psi(i + 1, j)) / 2.0;
    return result;
}

// Same idea for V, averaging across the y-direction instead.
Field2D psiOnV(const Field2D &psi, int imax, int jmax) {
    Field2D result(imax + 1, jmax);
    for (int i = 0; i <= imax; ++i)
        for (int j = 0; j < jmax; ++j) result(i, j) = (psi(i, j) + psi(i, j + 1)) / 2.0;
    return result;
}
}  // namespace

IBMCoeff LSPointIdent(const Domain &domain, double alpha, double beta, double q, int BQ,
                       const LS &ls, int UVP, double /*treshold*/, bool computeA1g) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // ---- Put psi onto this grid's own node locations, and get its
    // surface normal there (LSPointIdent.m:34-194, single-mineral
    // `else` branch only -- no LS.LS1 bi-mineral handling). ----
    Field2D psi;
    Field2D nx, ny;
    std::vector<double> x, y;  // this grid's own 1D axis arrays -- LSPointIdent.m's `x`/`y`
    double dx = 0.0;           // this grid's own min spacing -- LSPointIdent.m's `dx`
    if (UVP == -1) {  // u-velocity
        x = domain.xu;
        y = domain.yu;
        // put psi on p-grid to u grid:
        // averaging the two psi values at the P/scalar nodes immediately to its west and eas
        psi = psiOnU(ls.psi, imax, jmax);
        LSNormals normals = computeLSNormals(psi, domain);
        nx = normals.nx;
        ny = normals.ny;
        dx = *std::min_element(domain.dxu.begin(), domain.dxu.end());
    } else if (UVP == 0) {  // v-velocity
        x = domain.xv;
        y = domain.yv;
        psi = psiOnV(ls.psi, imax, jmax);
        LSNormals normals = computeLSNormals(psi, domain);
        nx = normals.nx;
        ny = normals.ny;
        dx = *std::min_element(domain.dxv.begin(), domain.dxv.end());
    } else {  // scalar / pressure grid
        x = domain.xp;
        y = domain.yp;
        //P/scalar dont need to project to a different grid, reuses LS.nx/LS.ny as already computed
        psi = ls.psi;
        nx = ls.nx;
        ny = ls.ny;
        dx = *std::min_element(domain.dxp.begin(), domain.dxp.end());
    }

    // ---- Classification on cell center/point (LSPointIdent.m:196-205): 
    // 1) start all-fluid, 0 = fluid
    // 2) mark solid wherever psi<0. 2= solid (the cell center is on the solid side from the boundary)
    // 3) mark ghost point, 1 = ghost
    //      fine the solid cell near solid-fluid interface. if all 4 neighbors are solid, the sum of flag = 8
    //      first fine in solid cell, which has fluid neighbors
    IBMCoeff ibm_coeff;
    // mark 0 or 2
    ibm_coeff.flag = Field2D(psi.nx(), psi.ny());
    for (int i = 0; i < psi.nx(); ++i)
        for (int j = 0; j < psi.ny(); ++j) ibm_coeff.flag(i, j) = (psi(i, j) < 0.0) ? 2.0 : 0.0;

    // find ghost cell by checking neighbors
    const int Nx = ibm_coeff.flag.nx();
    const int Ny = ibm_coeff.flag.ny();
    Field2D neighbor_flag_sum(Nx, Ny);
    for (int i = 1; i <= Nx - 2; ++i)
        for (int j = 1; j <= Ny - 2; ++j)
            neighbor_flag_sum(i, j) = ibm_coeff.flag(i + 1, j) + ibm_coeff.flag(i - 1, j) +
                                       ibm_coeff.flag(i, j + 1) + ibm_coeff.flag(i, j - 1);
    // update solid cell but has >= 1 fluid neighbor to ghost cell/point
    for (int i = 1; i <= Nx - 2; ++i)
        for (int j = 1; j <= Ny - 2; ++j)
            if (ibm_coeff.flag(i, j) == 2.0 && neighbor_flag_sum(i, j) != 8.0) ibm_coeff.flag(i, j) -= 1.0;

    // ---- Ghost-cell coordinates (LSPointIdent.m:240-244): collect
    // every (i,j) with flag==1, 
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            if (ibm_coeff.flag(i, j) == 1.0) {
                ibm_coeff.I_g.push_back(i);
                ibm_coeff.J_g.push_back(j);
            }
    // then look up each one's physical (x,y) position 
    // used only to feed LSmirPointsBQ ( ).        
    std::vector<double> X_g(ibm_coeff.I_g.size()), Y_g(ibm_coeff.J_g.size());
    for (size_t k = 0; k < ibm_coeff.I_g.size(); ++k) {
        X_g[k] = x[ibm_coeff.I_g[k]];
        Y_g[k] = y[ibm_coeff.J_g[k]];
    }

    // ---- Solid-cell indices (LSPointIdent.m:274-277: `[J_solid,I_solid]
    // = find(flag' == 2);` -- MATLAB transposes flag and swaps the two
    // outputs purely to control scan order; functionally identical to
    // this plain double loop, since nothing downstream cares what order
    // these come in). flag==2 here means "fully-interior solid" (all 4
    // neighbors also solid), post-demotion. ----
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            if (ibm_coeff.flag(i, j) == 2.0) {
                ibm_coeff.I_solid.push_back(i);
                ibm_coeff.J_solid.push_back(j);
            }

    // ---- Ghost-cell mirror-point stencil. Fills ibm_coeff.I_m/
    // J_m, I1..I6/J1..J6, lambda_g_1..6, A1_g in place, reading back the
    // I_g/J_g just collected above. The near-domain-edge safeguard
    // (LSmirPointsBQnew.m's own, see LSmirPointsBQ.h) only ever applies
    // to the P/scalar grid (UVP==1) in MATLAB -- x_0 mirrors
    // LSPointIdent.m's own `else: x_0=dx` branch (LS.case==3||4's
    // DOMAIN.x_0 alternative isn't ported, matching every other
    // case-3/4 omission in this project). ----
    bool applyEdgeSafeguard = (UVP == 1);
    LSmirPointsBQ(x, y, alpha, beta, q, X_g, Y_g, BQ, dx, psi, nx, ny, ibm_coeff, computeA1g,
                  applyEdgeSafeguard, /*x_0=*/dx, /*lx=*/domain.lx);

    return ibm_coeff;
}
