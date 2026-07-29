#include "LSnormals.h"
#include <algorithm>
#include <cmath>

LSNormals computeLSNormals(const Field2D &psi, const Domain &domain) {
    const int Nx = psi.nx();
    const int Ny = psi.ny();

    LSNormals result;
    result.nx = Field2D(Nx, Ny);
    result.ny = Field2D(Nx, Ny);

    const double dx = domain.dxp.empty() ? 0.0 : *std::min_element(domain.dxp.begin(), domain.dxp.end());
    const double dy = domain.dyp.empty() ? 0.0 : *std::min_element(domain.dyp.begin(), domain.dyp.end());
    const double epsi = 20.0 * dx;

    for (int i = 1; i <= Nx - 2; ++i) {
        for (int j = 1; j <= Ny - 2; ++j) {
            if (std::abs(psi(i, j)) >= epsi) continue;

            double nxv = (psi(i + 1, j) - psi(i - 1, j)) / (2.0 * dx);
            double nyv = (psi(i, j + 1) - psi(i, j - 1)) / (2.0 * dy);
            const double grad = std::sqrt(nxv * nxv + nyv * nyv);
            if (grad > 1e-15) {
                nxv /= grad;
                nyv /= grad;
            }
            result.nx(i, j) = nxv;
            result.ny(i, j) = nyv;
        }
    }
    return result;
}
