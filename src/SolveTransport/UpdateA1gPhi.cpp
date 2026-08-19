#include "UpdateA1gPhi.h"


std::vector<double> updateA1gPhi(const IBM &ibm, const IBMCoeff &ibmCoeffPhi,
                                   const std::vector<double> &q_G_species) {
    const int numg = static_cast<int>(ibmCoeffPhi.I_g.size());
    std::vector<double> A1_g(numg);
    const double alpha = ibm.alpha_phi;
    const double d = ibmCoeffPhi.Delta;

    for (int k = 0; k < numg; ++k) {
        double beta_use = ibmCoeffPhi.betaG[k];
        double r_g = ibmCoeffPhi.r_g[k];
        double q_use = q_G_species[k];
        A1_g[k] = q_use * (d + 2.0 * r_g + r_g * r_g / d) / (2.0 * alpha - beta_use * d);
    }

    return A1_g;
}
