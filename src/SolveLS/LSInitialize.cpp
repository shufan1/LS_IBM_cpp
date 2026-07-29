#include "LSInitialize.h"
#include <cmath>
#include <stdexcept>

Field2D computeLSInitialize(const Domain &domain, const LSCase &lsCase) {
    
    Field2D psi(domain.imax + 1, domain.jmax + 1);
    // Case 1 (grain/circle): signed distance to a circle centered at
    // (xc,yc) with radius diamcyl/2.
    if (lsCase.caseId == 1) {
        for (int i = 0; i <= domain.imax; ++i) {
            for (int j = 0; j <= domain.jmax; ++j) {
                double dx = domain.xp[i] - lsCase.xc;
                double dy = domain.yp[j] - lsCase.yc;
                psi(i, j) = std::sqrt(dx * dx + dy * dy) - lsCase.diamcyl / 2.0;
            }
        }
    } else{
      
        throw std::runtime_error(
            "computeLSInitialize: only LSCase.caseId==1 (grain) is implemented; got caseId=" +
            std::to_string(lsCase.caseId));
    }
    return psi;
}
