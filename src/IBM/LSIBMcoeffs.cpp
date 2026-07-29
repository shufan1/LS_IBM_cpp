#include "LSIBMcoeffs.h"
#include "LSPointIdent.h"

std::vector<IBMCoeff> LSIBMcoeffs(const IBM &ibm, const Domain &domain, const LS &ls,
                                   const std::vector<Field2D> &phi) {
    std::vector<IBMCoeff> result;
    result.reserve(2);
    result.push_back(  // U
        LSPointIdent(domain, ibm.alpha, ibm.beta, ibm.q, ibm.BQu, ls, -1, phi, ibm.treshold));
    result.push_back(  // V
        LSPointIdent(domain, ibm.alpha, ibm.beta, ibm.q, ibm.BQv, ls, 0, phi, ibm.treshold));
    return result;
}
