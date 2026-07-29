#include "DiffFlux.h"
#include "Timer.h"

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
    // DiffFlux.h's comment on computeDiffFluxU for why each excluded
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

DiffFluxCoeffs computeDiffFluxPhi(const Domain & /*domain*/, const Variables & /*variables*/, const BC & /*bc*/) {
    // TODO: mirrors DiffFlux.m's flag==0 branch (De/Dw/Dn/Ds using D=1/Pe,
    // plus the real D1_a/D2_a west-boundary correction) -- not ported yet.
    return DiffFluxCoeffs{};
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
    // Diffu_P: computeDiffFluxPhi() isn't implemented yet (returns an
    // all-empty DiffFluxCoeffs{}), so no timer here yet either -- add one
    // when that lands, using category "diff_flux_p" to match MATLAB's
    // getDiffFlux.m instrumentation.
    flux.Diffu_P = computeDiffFluxPhi(domain, variables, bc);
    return flux;
}
