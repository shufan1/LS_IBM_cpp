#include "RhsPhiADRE.h"

void rhsPhiADRE(const StateVar &stateVar, const BC &bc, const Flux &flux,
                 const ConvFluxCoeffs &convFluxPhi, const Domain &domain, const IBMCoeff &ibmCoeffPhi,
                 const IBM &ibm, const Variables &variables, const LS &ls, const Field2D &ap_p,
                 const std::vector<double> &A1_g, int speciesIndex, Vec RHS) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const double alpha_q = variables.alpha_q;

    const Field2D &A0_p = convFluxPhi.A0_p;
    const Field2D &Fe = convFluxPhi.Fe;
    const Field2D &Fw = convFluxPhi.Fw;
    const Field2D &Fn = convFluxPhi.Fn;
    const Field2D &Fs = convFluxPhi.Fs;
    const Field2D &alphae = convFluxPhi.alphae;
    const Field2D &alphaw = convFluxPhi.alphaw;
    const Field2D &alphan = convFluxPhi.alphan;
    const Field2D &alphas = convFluxPhi.alphas;
    const std::vector<double> &D2_a = flux.Diffu_Phi.D2_a;
    const Field2D &phi = stateVar.phi[speciesIndex];
    const Field2D &phi_prev = stateVar.phi_prev[speciesIndex];
    const std::vector<double> &phi_a = bc.phi_a[speciesIndex];

    Field2D S0(imax + 1, jmax + 1);
    Field2D S_e(imax + 1, jmax + 1);
    Field2D S_w(imax + 1, jmax + 1);
    Field2D S_n(imax + 1, jmax + 1);
    Field2D S_s(imax + 1, jmax + 1);
    Field2D S_ur(imax + 1, jmax + 1);

    // 1.compute S0: accumulation from last time (n-1)
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            S0(i, j) = A0_p(i, j) * phi_prev(i, j);
        }
    }

    // 2.compute S_ur, under-relaxation
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            S_ur(i, j) = ap_p(i, j) * (1.0 - alpha_q) * phi(i, j);
        }
    }

    // 3.compute Se,Sw,Sn,Ss, quick iteration, use phi from previous Quick iteration (k)
    // Each face's positive-/negative-flow branches (alphae/alphaw/alphan/
    // alphas switch between them, same upwind switches ConvFlux.cpp
    // already computed) use domain.g1c_*/g2c_* (Coordinates.cpp) as the
    // 3-point QUICK interpolation weights.
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double g1p = domain.g1c_e_p[i], g2p = domain.g2c_e_p[i];
            double g1n = domain.g1c_e_n[i], g2n = domain.g2c_e_n[i];
            S_e(i, j) = Fe(i, j) * (alphae(i, j) * (phi(i, j) * (g1p - g2p) + phi(i - 1, j) * g2p -
                                                      phi(i + 1, j) * g1p) +
                                     (1.0 - alphae(i, j)) * (phi(i + 1, j) * (g1n - g2n) +
                                                              phi(i + 2, j) * g2n - phi(i, j) * g1n));
        }
    }

    for (int i = 2; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double g1p = domain.g1c_w_p[i], g2p = domain.g2c_w_p[i];
            double g1n = domain.g1c_w_n[i], g2n = domain.g2c_w_n[i];
            S_w(i, j) = Fw(i, j) * (alphaw(i, j) * (phi(i - 1, j) * (g2p - g1p) + phi(i, j) * g1p -
                                                      phi(i - 2, j) * g2p) +
                                     (1.0 - alphaw(i, j)) * (phi(i, j) * (g2n - g1n) +
                                                              phi(i - 1, j) * g1n - phi(i + 1, j) * g2n));
        }
    }

    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            double g1p = domain.g1c_n_p[j], g2p = domain.g2c_n_p[j];
            double g1n = domain.g1c_n_n[j], g2n = domain.g2c_n_n[j];
            S_n(i, j) = Fn(i, j) * (alphan(i, j) * (phi(i, j) * (g1p - g2p) + phi(i, j - 1) * g2p -
                                                      phi(i, j + 1) * g1p) +
                                     (1.0 - alphan(i, j)) * (phi(i, j + 1) * (g1n - g2n) +
                                                              phi(i, j + 2) * g2n - phi(i, j) * g1n));
        }
    }

    // i in [2,imax-2], not the full [1,imax-1] -- narrower than S_n's
    // i-range even though nothing in the formula itself depends on i;
    // kept faithful to RHSPHIADRE.m's own S_s(3:imax-1,3:jmax,...) range
    // rather than "fixed", same spirit as COEFFU.cpp's south-boundary
    // quirk.
    for (int i = 2; i <= imax - 2; ++i) {
        for (int j = 2; j <= jmax - 1; ++j) {
            double g1p = domain.g1c_s_p[j], g2p = domain.g2c_s_p[j];
            double g1n = domain.g1c_s_n[j], g2n = domain.g2c_s_n[j];
            S_s(i, j) = Fs(i, j) * (alphas(i, j) * (phi(i, j - 1) * (g2p - g1p) + phi(i, j) * g1p -
                                                      phi(i, j - 2) * g2p) +
                                     (1.0 - alphas(i, j)) * (phi(i, j) * (g2n - g1n) +
                                                              phi(i, j - 1) * g1n - phi(i, j + 1) * g2n));
        }
    }

    // West/inlet Dirichlet row's S_w override (RHSPHIADRE.m:114) --
    // reuses the same Fw+D2_a coefficient as coeffPhiADRE's own S_p(1,j)
    // (see CoeffPhiADRE.cpp), now multiplied by the actual known
    // boundary value instead of folded into ap. Overwrites row i=1,
    // which the general S_w loop above never touches (it starts at i=2).
    for (int j = 1; j <= jmax - 1; ++j) {
        S_w(1, j) = (Fw(1, j) + D2_a[j]) * phi_a[j] * (ls.psi(0, j) > 0.0);
    }

    // 4.overwrite solid cell values
    for (size_t idx = 0; idx < ibmCoeffPhi.I_solid.size(); ++idx) {
        int is = ibmCoeffPhi.I_solid[idx];
        int js = ibmCoeffPhi.J_solid[idx];
        S_e(is, js) = 0.0;
        S_w(is, js) = 0.0;
        S_n(is, js) = 0.0;
        S_s(is, js) = 0.0;
        S0(is, js) = ibm.phi_inside_psi;
        S_ur(is, js) = 0.0;
    }

    // 5.updsate ghost cell rhs values
    for (size_t idx = 0; idx < ibmCoeffPhi.I_g.size(); ++idx) {
        int ig = ibmCoeffPhi.I_g[idx];
        int jg = ibmCoeffPhi.J_g[idx];
        S_e(ig, jg) = 0.0;
        S_w(ig, jg) = 0.0;
        S_n(ig, jg) = 0.0;
        S_s(ig, jg) = 0.0;
        S0(ig, jg) = A1_g[idx];
        S_ur(ig, jg) = 0.0;
    }

    // ---- Fill RHS (caller-allocated, same k(i,j) bijection as
    // coeffPhiADRE's CM).
    auto k = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double S = S_e(i, j) + S_w(i, j) + S_n(i, j) + S_s(i, j) + S0(i, j) + S_ur(i, j);
            VecSetValue(RHS, k(i, j), S, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);
}
