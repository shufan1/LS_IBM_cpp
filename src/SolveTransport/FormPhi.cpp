#include "FormPhi.h"

Field2D formPhi(const Domain &domain, Vec /*phiVec*/, const BC & /*bc*/, const Field2D & /*phiA*/,
                const Field2D & /*flag*/, double /*phiInside*/) {
    return Field2D(domain.imax + 1, domain.jmax + 1);
}
