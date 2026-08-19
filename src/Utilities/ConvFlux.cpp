#include "ConvFlux.h"
#include <cstdio>

namespace {
double upwind(double f) { return f > 0.0 ? 1.0 : 0.0; }
}  // namespace

ConvFluxCoeffs computeConvFluxU(const StateVar &stateVar, const ControlVar &controlVar,
                                 const Domain &domain, const Variables &variables, const BC &bc) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const Field2D &U = stateVar.U;
    const Field2D &V = stateVar.V;

    ConvFluxCoeffs cf;
    cf.Fe = Field2D(imax, jmax + 1);
    cf.Fw = Field2D(imax, jmax + 1);
    cf.Fn = Field2D(imax, jmax + 1);
    cf.Fs = Field2D(imax, jmax + 1);
    cf.alphae = Field2D(imax, jmax + 1);
    cf.alphaw = Field2D(imax, jmax + 1);
    cf.alphan = Field2D(imax, jmax + 1);
    cf.alphas = Field2D(imax, jmax + 1);
    cf.dF = Field2D(imax, jmax + 1);
    cf.A0_p = Field2D(imax, jmax + 1);

    // Interior: i in [1, imax-2], j in [1, jmax-1] -- same range as
    // computeDiffFluxU's interior.
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            cf.Fe(i, j) = (domain.CoEWu[i] * U(i + 1, j) + (1.0 - domain.CoEWu[i]) * U(i, j)) * domain.dyv[j - 1];
            cf.Fw(i, j) = (domain.CoEWu[i - 1] * U(i, j) + (1.0 - domain.CoEWu[i - 1]) * U(i - 1, j)) * domain.dyv[j - 1];
            cf.Fn(i, j) = (domain.CoNSu[i] * V(i + 1, j) + (1.0 - domain.CoNSu[i]) * V(i, j)) * domain.dxv[i];
            cf.Fs(i, j) = (domain.CoNSu[i] * V(i + 1, j - 1) + (1.0 - domain.CoNSu[i]) * V(i, j - 1)) * domain.dxv[i];
        }
    }

    // Outlet row (i=imax-1), only when the pressure BC there is
    // Dirichlet (BC_e_p==1), meaning U(imax-1,:) is a genuinely solved
    // unknown. Fw/Fn/Fs extend the same interior formula to this row;
    // Fe has no further east neighbor to interpolate against, so it
    // falls back to U's own value at the face.
    if (bc.BC_e_p == 1) {
        const int i = imax - 1;
        for (int j = 1; j <= jmax - 1; ++j) {
            cf.Fe(i, j) = U(i, j) * domain.dyv[j - 1];
            cf.Fw(i, j) = (domain.CoEWu[i - 1] * U(i, j) + (1.0 - domain.CoEWu[i - 1]) * U(i - 1, j)) * domain.dyv[j - 1];
            cf.Fn(i, j) = (domain.CoNSu[i] * V(i + 1, j) + (1.0 - domain.CoNSu[i]) * V(i, j)) * domain.dxv[i];
            cf.Fs(i, j) = (domain.CoNSu[i] * V(i + 1, j - 1) + (1.0 - domain.CoNSu[i]) * V(i, j - 1)) * domain.dxv[i];
        }
    }

    // Upwind switches and net flux imbalance over the whole array --
    // wherever Fe/Fw/Fn/Fs were never set, they read as 0, giving
    // alpha=0 there too, so this doesn't need to be restricted to the
    // interior.
    for (int i = 0; i < imax; ++i) {
        for (int j = 0; j <= jmax; ++j) {
            cf.alphae(i, j) = upwind(cf.Fe(i, j));
            cf.alphaw(i, j) = upwind(cf.Fw(i, j));
            cf.alphan(i, j) = upwind(cf.Fn(i, j));
            cf.alphas(i, j) = upwind(cf.Fs(i, j));
            cf.dF(i, j) = cf.Fe(i, j) - cf.Fw(i, j) + cf.Fn(i, j) - cf.Fs(i, j);
            cf.A0_p(i, j) = controlVar.flow_steady ? 0.0 : domain.dV_u(i, j) / variables.dt;
        }
    }

    // DEBUG (temporary): Fe/Fw/Fn/Fs/dF/A0_p at the specific U cells being
    // traced through the coeffU pipeline (see DiffuFlux.cpp/COEFFU.cpp for
    // the matching prints) -- this runs every SIMPLE iteration, since
    // Fe/Fw/Fn/Fs are relinearized around the current velocity field.
    // {
    //     // (i,j) here are C++ 0-based -- MATLAB-native equivalent is (i+1,j+1).
    //     int wanted[][2] = {{214, 64}, {186, 137}};
    //     for (auto &w : wanted) {
    //         int i = w[0], j = w[1];
    //         printf("  [ConvFluxU] (%d,%d) [0-based]: Fe=%.8e Fw=%.8e Fn=%.8e Fs=%.8e dF=%.8e A0_p=%.8e\n",
    //                i, j, cf.Fe(i, j), cf.Fw(i, j), cf.Fn(i, j), cf.Fs(i, j), cf.dF(i, j), cf.A0_p(i, j));
    //     }
    // }

    return cf;
}

ConvFluxCoeffs computeConvFluxV(const StateVar &stateVar, const ControlVar &controlVar,
                                 const Domain &domain, const Variables &variables, const BC & /*bc*/) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const Field2D &U = stateVar.U;
    const Field2D &V = stateVar.V;

    ConvFluxCoeffs cf;
    cf.Fe = Field2D(imax + 1, jmax);
    cf.Fw = Field2D(imax + 1, jmax);
    cf.Fn = Field2D(imax + 1, jmax);
    cf.Fs = Field2D(imax + 1, jmax);
    cf.alphae = Field2D(imax + 1, jmax);
    cf.alphaw = Field2D(imax + 1, jmax);
    cf.alphan = Field2D(imax + 1, jmax);
    cf.alphas = Field2D(imax + 1, jmax);
    cf.dF = Field2D(imax + 1, jmax);
    cf.A0_p = Field2D(imax + 1, jmax);

    // Interior: i in [1, imax-1], j in [1, jmax-2] -- same range as
    // computeDiffFluxV's interior. No outlet-row extension exists for V
    // (ConvFlux.m's V-branch has no boundary-override block), and the
    // upwind switches are only ever populated here too (left at 0
    // elsewhere), so both are computed together in this one loop.
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            cf.Fe(i, j) = (domain.CoEWv[j] * U(i, j + 1) + (1.0 - domain.CoEWv[j]) * U(i, j)) * domain.dyu[j];
            cf.Fw(i, j) = (domain.CoEWv[j] * U(i - 1, j + 1) + (1.0 - domain.CoEWv[j]) * U(i - 1, j)) * domain.dyu[j];
            cf.Fn(i, j) = (domain.CoNSv[j] * V(i, j + 1) + (1.0 - domain.CoNSv[j]) * V(i, j)) * domain.dxu[i - 1];
            cf.Fs(i, j) = (domain.CoNSv[j - 1] * V(i, j) + (1.0 - domain.CoNSv[j - 1]) * V(i, j - 1)) * domain.dxu[i - 1];

            cf.alphae(i, j) = upwind(cf.Fe(i, j));
            cf.alphaw(i, j) = upwind(cf.Fw(i, j));
            cf.alphan(i, j) = upwind(cf.Fn(i, j));
            cf.alphas(i, j) = upwind(cf.Fs(i, j));
        }
    }

    for (int i = 0; i <= imax; ++i) {
        for (int j = 0; j < jmax; ++j) {
            cf.dF(i, j) = cf.Fe(i, j) - cf.Fw(i, j) + cf.Fn(i, j) - cf.Fs(i, j);
            cf.A0_p(i, j) = controlVar.flow_steady ? 0.0 : domain.dV_v(i, j) / variables.dt;
        }
    }

    return cf;
}

ConvFluxCoeffs computeConvFluxPhi(const StateVar &stateVar, const ControlVar &controlVar,
                                   const Domain &domain, const Variables &variables,
                                   const BC & /*bc*/) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const Field2D &U = stateVar.U;
    const Field2D &V = stateVar.V;

    ConvFluxCoeffs cf;
    cf.Fe = Field2D(imax + 1, jmax + 1);
    cf.Fw = Field2D(imax + 1, jmax + 1);
    cf.Fn = Field2D(imax + 1, jmax + 1);
    cf.Fs = Field2D(imax + 1, jmax + 1);
    cf.alphae = Field2D(imax + 1, jmax + 1);
    cf.alphaw = Field2D(imax + 1, jmax + 1);
    cf.alphan = Field2D(imax + 1, jmax + 1);
    cf.alphas = Field2D(imax + 1, jmax + 1);
    cf.dF = Field2D(imax + 1, jmax + 1);
    cf.A0_p = Field2D(imax + 1, jmax + 1);

    // Interior: i in [1, imax-1], j in [1, jmax-1] -- same range as
    // computeDiffFluxPhi's interior (all four faces share one raw block
    // here, unlike momentum's split ranges). No CoEW/CoNS interpolation
    // at all -- U/V already live exactly at a P-cell's own faces, so
    // each face flux reads the matching U/V value directly, with no
    // BC-dependent special case (ConvFlux.m's flag==0 branch has none).
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            cf.Fe(i, j) = U(i, j) * domain.dyv[j - 1];
            cf.Fw(i, j) = U(i - 1, j) * domain.dyv[j - 1];
            cf.Fn(i, j) = V(i, j) * domain.dxu[i - 1];
            cf.Fs(i, j) = V(i, j - 1) * domain.dxu[i - 1];

            cf.alphae(i, j) = upwind(cf.Fe(i, j));
            cf.alphaw(i, j) = upwind(cf.Fw(i, j));
            cf.alphan(i, j) = upwind(cf.Fn(i, j));
            cf.alphas(i, j) = upwind(cf.Fs(i, j));
        }
    }

    for (int i = 0; i <= imax; ++i) {
        for (int j = 0; j <= jmax; ++j) {
            cf.dF(i, j) = cf.Fe(i, j) - cf.Fw(i, j) + cf.Fn(i, j) - cf.Fs(i, j);
            cf.A0_p(i, j) = controlVar.transport_steady ? 0.0 : domain.dV_p(i, j) / variables.dt;
        }
    }

    return cf;
}
