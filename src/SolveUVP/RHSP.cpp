#include "RHSP.h"
#include "SolveUVP.h"  // pinnedPressureCell()

void rhsP(const Field2D &U_star, const Field2D &V_star, const Domain &domain,
          const ControlVar &controlVar, const IBM &ibm, Vec RHS_P2) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // Same reference-cell exclusion as coeffP()'s CM2 -- see COEFFP.cpp
    // for why this keeps CM_p and RHS_P2 consistently indexed.
    PressurePin pin = pinnedPressureCell(imax, jmax, controlVar.imposePresBC);
    auto kFull = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    const int kPin = kFull(pin.i0, pin.j0);
    auto kReduced = [kPin](int kf) { return kf < kPin ? kf : kf - 1; };

    // Ghost-cell placeholder substitution (RHSP.h): reads ibm.flag_u/
    // flag_v instead of mutating U_star/V_star -- see RHSP.h for why
    // that matches MATLAB's own local-copy scoping without needing an
    // actual copy. phi_inside is hardcoded to 0 here, matching RHSP.m's
    // own local `phi_inside=0` (not ibm.u_inside_psi -- a separate,
    // coincidentally-equal constant in the active config).
    const double phi_inside = 0.0;
    auto uMasked = [&](int i, int j) { return ibm.flag_u(i, j) == 1.0 ? phi_inside : U_star(i, j); };
    auto vMasked = [&](int i, int j) { return ibm.flag_v(i, j) == 1.0 ? 0.0 : V_star(i, j); };

    // P-cell (i,j)'s west/east U-faces are U(i-1,j)/U(i,j); south/north
    // V-faces are V(i,j-1)/V(i,j) -- same face convention already
    // established for S_Pres in coeffU()/coeffV() (P(i,j)/P(i+1,j) there
    // is this same relationship viewed from the U-cell's side).
    // every index rewrite value
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            if (i == pin.i0 && j == pin.j0) continue;
            double S0 = (uMasked(i - 1, j) - uMasked(i, j)) * domain.dyv[j - 1] +
                        (vMasked(i, j - 1) - vMasked(i, j)) * domain.dxu[i - 1];
            VecSetValue(RHS_P2, kReduced(kFull(i, j)), S0, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(RHS_P2);
    VecAssemblyEnd(RHS_P2);
}
