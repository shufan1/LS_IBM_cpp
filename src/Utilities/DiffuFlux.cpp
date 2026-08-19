#include "DiffuFlux.h"
#include "Timer.h"
#include <cstdio>

DiffFluxCoeffs computeDiffFluxU(const Domain &domain, const Variables &variables, const BC &bc) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const double Re = variables.Re;

    DiffFluxCoeffs coeffs;
    coeffs.De = Field2D(imax, jmax + 1);
    coeffs.Dw = Field2D(imax, jmax + 1);
    coeffs.Dn = Field2D(imax, jmax + 1);
    coeffs.Ds = Field2D(imax, jmax + 1);

    // Interior only -- i in [1, imax-2], j in [1, jmax-1] (0-indexed).
    // Everything outside this range stays at Field2D's zero default; see
    // DiffuFlux.h's comment on computeDiffFluxU for why each excluded
    // row/column doesn't need a diffusion coefficient here.
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            coeffs.De(i, j) = domain.dyv[j - 1] / (Re * domain.dxu[i]);
            coeffs.Dw(i, j) = domain.dyv[j - 1] / (Re * domain.dxu[i - 1]);
            coeffs.Dn(i, j) = domain.dxv[i] / (Re * domain.dyu[j]);
            coeffs.Ds(i, j) = domain.dxv[i] / (Re * domain.dyu[j - 1]);
        }
    }

    // MATLAB: `if BC_e_p == 1 ... end` -- Dw/Dn/Ds get real values at the
    // east row (i=imax-1) too, in that case. De's analogous line is
    // commented out in the original source, so De stays zero there
    // regardless. Dead for the active config (BC_e_p=3), kept faithful
    // to the reference for when it isn't.
    if (bc.BC_e_p == 1) {
        const int i = imax - 1;
        for (int j = 1; j <= jmax - 1; ++j) {
            coeffs.Dw(i, j) = domain.dyv[j - 1] / (Re * domain.dxu[imax - 2]);
            coeffs.Dn(i, j) = domain.dxv[imax - 1] / (Re * domain.dyu[j]);
            coeffs.Ds(i, j) = domain.dxv[imax - 1] / (Re * domain.dyu[j - 1]);
        }
    }

    // DEBUG (temporary): De/Dw/Dn/Ds at the specific U cells being traced
    // through the coeffU pipeline (see COEFFU.cpp/ConvFlux.cpp for the
    // matching prints) -- computeDiffFluxU() only runs once (solution-
    // independent), so this fires exactly once per program run, not once
    // per SIMPLE iteration.
    {
        // (i,j) here are C++ 0-based -- MATLAB-native equivalent is (i+1,j+1).
        int wanted[][2] = {{214, 64}, {186, 137}};
        for (auto &w : wanted) {
            int i = w[0], j = w[1];
            printf("  [DiffFluxU] (%d,%d) [0-based]: De=%.8e Dw=%.8e Dn=%.8e Ds=%.8e\n", i, j,
                   coeffs.De(i, j), coeffs.Dw(i, j), coeffs.Dn(i, j), coeffs.Ds(i, j));
        }
    }

    return coeffs;
}

DiffFluxCoeffs computeDiffFluxV(const Domain &domain, const Variables &variables, const BC & /*bc*/) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const double Re = variables.Re;

    DiffFluxCoeffs coeffs;
    coeffs.De = Field2D(imax + 1, jmax);
    coeffs.Dw = Field2D(imax + 1, jmax);
    coeffs.Dn = Field2D(imax + 1, jmax);
    coeffs.Ds = Field2D(imax + 1, jmax);

    // Interior only -- i in [1, imax-1], j in [1, jmax-2] (0-indexed) --
    // the x/y-mirror of computeDiffFluxU's range, matching V's own
    // mirrored (imax+1, jmax) shape and BCs (V's north/south walls are
    // literal domain-boundary Dirichlet rows needing no stencil here;
    // its east/west columns borrow P's padded axis, same role U's
    // north/south ghost columns played). No BC_e_p==1 special case --
    // DiffFlux.m's flag==-1 branch doesn't have one.
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            coeffs.De(i, j) = domain.dyu[j] / (Re * domain.dxv[i]);
            coeffs.Dw(i, j) = domain.dyu[j] / (Re * domain.dxv[i - 1]);
            coeffs.Dn(i, j) = domain.dxu[i - 1] / (Re * domain.dyv[j]);
            coeffs.Ds(i, j) = domain.dxu[i - 1] / (Re * domain.dyv[j - 1]);
        }
    }

    return coeffs;
}

DiffFluxCoeffs computeDiffFluxPhi(const Domain &domain, const Variables &variables, const BC & /*bc*/) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const double Gamma = variables.D;  // shared by every species D/Pe

    DiffFluxCoeffs coeffs;
    coeffs.De = Field2D(imax + 1, jmax + 1);
    coeffs.Dw = Field2D(imax + 1, jmax + 1);
    coeffs.Dn = Field2D(imax + 1, jmax + 1);
    coeffs.Ds = Field2D(imax + 1, jmax + 1);

    // Interior only -- i in [1, imax-1], j in [1, jmax-1] (0-indexed):
    // De/Dw/Dn/Ds all occupy the identical raw block here (unlike the
    // momentum branches, which split De/Dw from Dn/Ds by one row/column
    // -- DiffFlux.m's flag==0 branch pads all four into the same
    // (2:imax,2:jmax) 1-based slot). Everything outside this range stays
    // at Field2D's zero default.
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            coeffs.De(i, j) = Gamma * domain.dyv[j - 1] / domain.dxp[i];
            coeffs.Dw(i, j) = Gamma * domain.dyv[j - 1] / domain.dxp[i - 1];
            coeffs.Dn(i, j) = Gamma * domain.dxu[i - 1] / domain.dyp[j];
            coeffs.Ds(i, j) = Gamma * domain.dxu[i - 1] / domain.dyp[j - 1];
        }
    }

    // West/inlet boundary correction (DiffFlux.m's D1_a/D2_a): a one-
    // sided 3-point formula using only the fixed x-positions xu[0]/
    // xp[0]/xp[1] (MATLAB's xu(1)/xp(2)/xp(3)) -- no j-dependence at all
    // (D1/D2 are scalars), and MATLAB's own dxv(1,1:jmax-1) term is
    // uniform across j too (dxv only varies along i in this codebase's
    // 1D storage), so both D1_a/D2_a end up as a single constant value
    // repeated across the interior j range, zero at j=0 and j=jmax.
    const double D1 = Gamma * (domain.xu[0] - domain.xp[1]) /
                       ((domain.xp[2] - domain.xp[1]) * (domain.xp[2] - domain.xu[0]));
    const double D2 = Gamma * (domain.xp[2] - 2.0 * domain.xu[0] + domain.xp[1]) /
                       ((domain.xp[1] - domain.xu[0]) * (domain.xp[2] - domain.xu[0]));
    coeffs.D1_a.assign(jmax + 1, 0.0);
    coeffs.D2_a.assign(jmax + 1, 0.0);
    for (int j = 1; j <= jmax - 1; ++j) {
        coeffs.D1_a[j] = D1 * domain.dxv[0];
        coeffs.D2_a[j] = D2 * domain.dxv[0];
    }

    return coeffs;
}

Flux computeDiffFlux(const Domain &domain, const Variables &variables, const BC &bc) {
    Flux flux;
    {
        ScopedTimer t("diff_flux_u");
        flux.Diffu_U = computeDiffFluxU(domain, variables, bc);
    }
    {
        ScopedTimer t("diff_flux_v");
        flux.Diffu_V = computeDiffFluxV(domain, variables, bc);
    }
    {
        ScopedTimer t("diff_flux_p");
        flux.Diffu_Phi = computeDiffFluxPhi(domain, variables, bc);
    }
    return flux;
}
