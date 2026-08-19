#include "CalculateQG.h"
#include "../Utilities/Interp.h"
#include <algorithm>
#include <cmath>

std::vector<std::vector<double>> calculateQG(const Variables &variables,
                                              const std::vector<Field2D> &phi,
                                              const Domain &domain, const LS &ls,
                                              const std::vector<int> &I_g,
                                              const std::vector<int> &J_g) {
    const int Np = variables.Np;
    int num_g = static_cast<int>(I_g.size());
    std::vector<std::vector<double>> q_G(Np, std::vector<double>(num_g, 0.0));

    double dx = *std::min_element(domain.dxp.begin(), domain.dxp.end());
    double Delta = std::sqrt(2.0) * dx;

    for (int k = 0; k < num_g; k++) {
        int i_k = I_g[k], j_k = J_g[k];
        if (std::abs(ls.psi(i_k, j_k)) >= variables.LSgamma) continue;

        // get boundary point coordinates
        double x_g = domain.xp[i_k] + ls.nx(i_k, j_k) * std::abs(ls.psi(i_k, j_k));
        double y_g = domain.yp[j_k] + ls.ny(i_k, j_k) * std::abs(ls.psi(i_k, j_k));

        // get the two probe points in the fluid domain
        double x_d = x_g + Delta * ls.nx(i_k, j_k);
        double y_d = y_g + Delta * ls.ny(i_k, j_k);
        double x_2d = x_d + Delta * ls.nx(i_k, j_k);
        double y_2d = y_d + Delta * ls.ny(i_k, j_k);

        // use bilinear inerpolation get their values
        std::vector<double> phi_d(Np);
        std::vector<double> phi_2d(Np);
        for (int i_s = 0; i_s < Np; i_s++) {
            phi_d[i_s] = bilinearInterp(domain.xp, domain.yp, phi[i_s], x_d, y_d);
            phi_2d[i_s] = bilinearInterp(domain.xp, domain.yp, phi[i_s], x_2d, y_2d);
        }

        // for each species, get estimated phi at interface:
        // second-order one-sided (forward) finite difference
        // dC/dn = (-3C_gamma+4C1-C2)/2d = A C_gamma -> solve C_gamma
        // C_gamma = (3I-2dA)^-1 (4C1-C2)
        std::vector<double> phi_gamma(Np);
        for (int i_s = 0; i_s < Np; i_s++) {
            // get interpolated c at interface per species
            double phi_gamma_s = 0.0;
            for (int j_s = 0; j_s < Np; j_s++) {
                // variables.inv_A already defined as 3I-2dA^-1
                phi_gamma_s += variables.inv_A[i_s * Np + j_s] * (4 * phi_d[j_s] - phi_2d[j_s]);
            }
            phi_gamma[i_s] = phi_gamma_s;
        }

        // compute off diagonal contribution as q_k
        for (int i_s = 0; i_s < Np; i_s++) {
            double q_k = 0.0;
            for (int j_s = 0; j_s < Np; j_s++) {
                if (i_s != j_s) {
                    q_k += variables.A[i_s * Np + j_s] * phi_gamma[j_s];
                }
            }
            q_G[i_s][k] = -q_k;
        }
    }

    return q_G;
}
