#include "COEFFP.h"
#include "SolveUVP.h"  // pinnedPressureCell()

void coeffP(const Field2D &d_u, const Field2D &d_v, const Domain &domain,
            const ControlVar &controlVar, Mat CM2) {

    // coeff d_u = A_iJ / a_iJ, a_iJ coeff in u-momentum
    // coeff d_v = A_Ij / a_Ij, a_Ij coeff in v-momentum

    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // ---- Scratch coefficients (not returned -- see COEFFP.h) ----
    Field2D a_e(imax + 1, jmax + 1), a_w(imax + 1, jmax + 1);
    Field2D a_n(imax + 1, jmax + 1), a_s(imax + 1, jmax + 1);
    Field2D a_p(imax + 1, jmax + 1);

    // a_e/a_w only cover i=1..imax-2 / i=2..imax-1 respectively (one shy
    // of the full i=1..imax-1 interior on each side); a_n/a_s similarly
    // cover j=1..jmax-2 / j=2..jmax-1 -- same "one edge is missing, the
    // other equation's own boundary treatment covers it" pattern as
    // coeffU()/coeffV(), except here there's no boundary correction term
    // at all (S_p is provably always zero -- see COEFFP.h) that
    // formPCor()'s unconditional zero-gradient BCs handles instead.
    for (int i = 1; i <= imax - 2; ++i)
        for (int j = 1; j <= jmax - 1; ++j) a_e(i, j) = -d_u(i, j) * domain.dyv[j - 1];
    for (int i = 2; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax - 1; ++j) a_w(i, j) = -d_u(i - 1, j) * domain.dyv[j - 1];
    for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax - 2; ++j) a_n(i, j) = -d_v(i, j) * domain.dxu[i - 1];
    for (int i = 1; i <= imax - 1; ++i)
        for (int j = 2; j <= jmax - 1; ++j) a_s(i, j) = -d_v(i, j - 1) * domain.dxu[i - 1];

    // a_p = -(sum of neighbors); any exact zero is forced to 1 (COEFFP.m:
    // `AP(AP==0)=1`) as a safeguard against a singular diagonal entry.
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double ap = -(a_w(i, j) + a_e(i, j) + a_s(i, j) + a_n(i, j));
            a_p(i, j) = (ap == 0.0) ? 1.0 : ap;
        }
    }

    // ---- Fill the reduced system CM2 (caller-allocated), excluding the
    // pinned cell ----
    PressurePin pin = pinnedPressureCell(imax, jmax, controlVar.imposePresBC);

    auto kFull = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    const int kPin = kFull(pin.i0, pin.j0);
    // Every full-system index after the pinned one shifts down by 1 to
    // keep the reduced system's row numbering dense (0..L-2).
    auto kReduced = [kPin](int kf) { return kf < kPin ? kf : kf - 1; };

    MatZeroEntries(CM2);

    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            if (i == pin.i0 && j == pin.j0) continue;  // pinned cell has no row of its own
            int row = kReduced(kFull(i, j));
            MatSetValue(CM2, row, row, a_p(i, j), INSERT_VALUES);

            auto link = [&](int ni, int nj, double coeff) {
                if (ni == pin.i0 && nj == pin.j0) return;  // pinned column doesn't exist either
                MatSetValue(CM2, row, kReduced(kFull(ni, nj)), coeff, INSERT_VALUES);
            };
            if (i > 1) link(i - 1, j, a_w(i, j));
            if (i < imax - 1) link(i + 1, j, a_e(i, j));
            if (j > 1) link(i, j - 1, a_s(i, j));
            if (j < jmax - 1) link(i, j + 1, a_n(i, j));
        }
    }
    MatAssemblyBegin(CM2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CM2, MAT_FINAL_ASSEMBLY);
}
