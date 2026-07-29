#include "RHSP_PISO.h"
#include "SolveUVP.h"  // pinnedPressureCell()

namespace {
// corrU/corrV: the neighbor-coupling correction term, computed over
// sys.d_piso's own valid interior range and left at 0 elsewhere (Field2D
// default) -- see RHSP_PISO.h for why that zero-padding is exactly what
// reproduces RHSP_PISO.m's own boundary behavior.
Field2D correctionTerm(const MomentumCoeffs &sys, const Field2D &delta, int nx, int ny,
                        int iBegin, int iEnd, int jBegin, int jEnd) {
    Field2D corr(nx, ny);
    for (int i = iBegin; i <= iEnd; ++i) {
        for (int j = jBegin; j <= jEnd; ++j) {
            double neighborSum = sys.aw(i, j) * delta(i - 1, j) + sys.ae(i, j) * delta(i + 1, j) +
                                  sys.as(i, j) * delta(i, j - 1) + sys.an(i, j) * delta(i, j + 1);
            corr(i, j) = -sys.d_piso(i, j) * neighborSum;
        }
    }
    return corr;
}
}  // namespace

void rhsPPiso(const MomentumCoeffs &sysU, const MomentumCoeffs &sysV,
              const Field2D &deltaU, const Field2D &deltaV, const Domain &domain,
              const ControlVar &controlVar, Vec RHS_P2) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // Same interior range as coeffU()/coeffV()'s own ap/d_SIMPLC/d_piso.
    Field2D corrU = correctionTerm(sysU, deltaU, imax, jmax + 1, 1, imax - 2, 1, jmax - 1);
    Field2D corrV = correctionTerm(sysV, deltaV, imax + 1, jmax, 1, imax - 1, 1, jmax - 2);

    // Same reference-cell exclusion as coeffP()/rhsP() -- keeps this
    // consistently indexed against CM_p.
    PressurePin pin = pinnedPressureCell(imax, jmax, controlVar.imposePresBC);
    auto kFull = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    const int kPin = kFull(pin.i0, pin.j0);
    auto kReduced = [kPin](int kf) { return kf < kPin ? kf : kf - 1; };

    // Same face convention as rhsP(): west/east U-faces are
    // corrU(i-1,j)/corrU(i,j); south/north V-faces are
    // corrV(i,j-1)/corrV(i,j).
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            if (i == pin.i0 && j == pin.j0) continue;
            // add west, east, south and north face contribution
            // west face contribution: from u-momentum, involves four of u neighbors
            double S0 = (corrU(i - 1, j) - corrU(i, j)) + (corrV(i, j - 1) - corrV(i, j));
            VecSetValue(RHS_P2, kReduced(kFull(i, j)), S0, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(RHS_P2);
    VecAssemblyEnd(RHS_P2);
}
