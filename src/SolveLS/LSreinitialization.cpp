#include "LSreinitialization.h"
#include "WENODerivative.h"
#include "godunovGradientNorm.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// solve  ∂ψ/∂τ = sign(ψ) · (1 − ‖∇ψ‖)
Field2D LSreinitialization(const Field2D &psi_in, const Field2D &psi_prev, const Domain &domain,
                           const Variables &variables) {

    // Four sweeps restore ‖∇ψ‖ = 1 out to about two cells from the interface
    // n_iter × dtau  =  4 × 0.5·dx  =  2 cells
    // solve until covergence () / we say up till 4 iterations
        // 1 estimate ∇ψ by weno
        // 2 estiamte ‖∇ψ‖
        // 3. smoothed form of sign function: on orignal psi
        // 4. Runge–Kutta to solve get corrected psi
    // flatten out anything beyond, far from psi = 0 we dont care

    const int Nx = psi_in.nx();
    const int Ny = psi_in.ny();
    const size_t N = psi_in.data().size();

    const double dtau = variables.dtau;
    const double h = variables.LSgamma;
    const int n_iter = variables.n_iter_ReLS;
    const std::string &scheme = variables.TimeSchemeRLS;
    const double epsSign = *std::min_element(domain.dxp.begin(), domain.dxp.end());

    Field2D psi = psi_in;
    const Field2D psi_n = psi_in;
    const std::vector<double> &pn = psi_n.data();

    // 3. smoothed form of sign function: on orignal psi
    std::vector<double> signPsi(N);
    for (size_t k = 0; k < N; ++k)
        signPsi[k] = pn[k] / std::sqrt(pn[k] * pn[k] + epsSign * epsSign);

    for (int it = 0; it < n_iter; ++it) {
        const std::vector<double> &p0 = psi.data();

        // 1 estimate ∇ψ by weno
        const LSGradient g1 = wenoDerivative(psi, psi_prev, domain, h, LSEquation::ReinitializationEqn);
        // 2 estiamte ‖∇ψ‖
        const Field2D G1 = godunovGradientNorm(g1, psi_n);
        const std::vector<double> &gg1 = G1.data();

        // 4. Runge–Kutta to solve get corrected psi
        Field2D psi_m1(Nx, Ny);
        std::vector<double> &pm1 = psi_m1.data();
        for (size_t n = 0; n < N; ++n)
            pm1[n] = p0[n] + dtau * signPsi[n] * (1.0 - gg1[n]) * (gg1[n] != 0.0);

        if (scheme == "RK1") { psi = psi_m1; continue; }

        const LSGradient g2 = wenoDerivative(psi_m1, psi_prev, domain, h, LSEquation::ReinitializationEqn);
        const Field2D G2 = godunovGradientNorm(g2, psi_n);
        const std::vector<double> &gg2 = G2.data();

        Field2D psi_m2(Nx, Ny);
        std::vector<double> &pm2 = psi_m2.data();
        for (size_t n = 0; n < N; ++n)
            pm2[n] = pm1[n] + dtau * signPsi[n] * (1.0 - gg2[n]) * (gg2[n] != 0.0);

        if (scheme == "RK2") {
            Field2D next(Nx, Ny);
            std::vector<double> &r = next.data();
            for (size_t n = 0; n < N; ++n) r[n] = 0.5 * (p0[n] + pm2[n]);
            psi = next;
            continue;
        }

        if (scheme == "RK3") {
            Field2D psi_m12(Nx, Ny);
            std::vector<double> &pm12 = psi_m12.data();
            for (size_t n = 0; n < N; ++n) pm12[n] = 0.75 * p0[n] + 0.25 * pm2[n];

            const LSGradient g3 = wenoDerivative(psi_m12, psi_prev, domain, h, LSEquation::ReinitializationEqn);
            const Field2D G3 = godunovGradientNorm(g3, psi_n);
            const std::vector<double> &gg3 = G3.data();

            Field2D psi_m32(Nx, Ny);
            std::vector<double> &pm32 = psi_m32.data();
            for (size_t n = 0; n < N; ++n)
                pm32[n] = pm12[n] + dtau * signPsi[n] * (1.0 - gg3[n]) * (gg3[n] != 0.0);

            Field2D next(Nx, Ny);
            std::vector<double> &r = next.data();
            for (size_t n = 0; n < N; ++n) r[n] = (p0[n] + 2.0 * pm32[n]) / 3.0;
            psi = next;
            continue;
        }

        throw std::runtime_error("LSreinitialization: unknown TimeSchemeRLS \"" + scheme +
                                  "\" -- expected \"RK1\", \"RK2\" or \"RK3\"");
    }

    // flatten out anything beyond, far from psi = 0 we dont care
    std::vector<double> &p = psi.data();
    for (size_t n = 0; n < N; ++n) p[n] = std::max(-h, std::min(h, p[n]));

    return psi;
}
