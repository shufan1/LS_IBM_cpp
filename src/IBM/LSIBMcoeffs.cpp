#include "LSIBMcoeffs.h"
#include "LSPointIdent.h"

std::vector<IBMCoeff> LSIBMcoeffsUV(const IBM &ibm, const Domain &domain, const LS &ls) {
    std::vector<IBMCoeff> result;
    result.reserve(2);
    result.push_back(  // U
        LSPointIdent(domain, ibm.alpha, ibm.beta, ibm.q, ibm.BQu, ls, -1, ibm.treshold));
    result.push_back(  // V
        LSPointIdent(domain, ibm.alpha, ibm.beta, ibm.q, ibm.BQv, ls, 0, ibm.treshold));
    return result;
}

std::vector<IBMCoeff> LSIBMcoeffsPhi(const IBM &ibm, const Domain &domain, const LS &ls, int Np) {
    // Deliberately NOT a literal mirror of LSIBMcoeffs.m, which seeds
    // this with a dead single-scalar beta_phi and relies on
    // SolveTransportADRE.m's update_A1g to overwrite every species'
    // landa_g_k/A1_g with the real beta_G before first use. Reading the
    // real per-species value straight from ibm.beta_phi[i_s] (populated
    // by main.cpp from -diag(variables.A) once A is loaded -- see
    // IBM::beta_phi's comment) instead skips that redundant round-trip
    // and produces already-correct geometry for every species directly
    // -- safe even for a future nonlinear-BC model, since that path
    // re-derives beta every QUICK iteration regardless of what this seed
    // used (see SolveTransportADRE.cpp's nonlinearReactionBC branch).
    std::vector<IBMCoeff> result;
    result.reserve(Np);
    for (int i_s = 0; i_s < Np; ++i_s) {
        // computeA1g=false: based on q_G which needs I_g/J_g
        // updated when q_G is computed per Quick iteration in solving transport
        result.push_back(LSPointIdent(domain, ibm.alpha_phi, ibm.beta_phi[i_s], ibm.q_phi[i_s], ibm.BQp,
                                       ls, 1, ibm.treshold, /*computeA1g=*/false));
    }
    return result;
}
