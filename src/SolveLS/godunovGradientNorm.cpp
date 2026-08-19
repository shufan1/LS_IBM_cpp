#include "godunovGradientNorm.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// Not sqrt(psi_x^2 + psi_y^2): after the reinitialization WENO pass there
// are TWO candidates per direction -- a = psi_xn (left-biased) and
// b = psi_xp (right-biased) -- and averaging or picking arbitrarily
// converges to the wrong (non-viscosity) solution. The Godunov flux
// selects between them.
//
// The selection is written branchlessly as max of squared one-sided
// parts. For psi_n >= 0, in x, with aP = max(a,0) and bM = min(b,0):
//
//    a (backward)  b (forward)   aP    bM    max(aP^2,bM^2)   meaning
//   ------------- ------------- ----- ----- ---------------- ---------------
//        > 0           > 0        a     0         a^2         psi rising;
//                                                             upwind is the
//                                                             BACKWARD diff
//        < 0           < 0        0     b         b^2         psi falling;
//                                                             upwind is the
//                                                             FORWARD diff
//        > 0           < 0        a     b     max(a^2,b^2)    local MAXIMUM
//                                                             (kink) -- take
//                                                             the steeper side
//        < 0           > 0        0     0          0          local MINIMUM
//                                                             -- gradient is
//                                                             genuinely zero
//
// Gn is the mirror image, using aM = min(a,0) and bP = max(b,0), because
// for psi_n < 0 the reinitialization correction travels the other way.
//
// sign(psi_n) is what picks the branch. The Hamiltonian
// sign(psi)*(|grad psi| - 1) is nonlinear, so there is no velocity whose
// sign gives the upwind direction -- but the correction always propagates
// OUTWARD from the zero contour, so the sign of psi does tell you which
// way. It replaces sign(u) from the advection branch.
//
// The two components then combine as an ordinary Euclidean norm, so G
// really is ||grad psi||, just with each component chosen upwind rather
// than centred.
Field2D godunovGradientNorm(const LSGradient &grad, const Field2D &psi_n) {
    if (grad.psi_xn.nx() == 0 || grad.psi_xp.nx() == 0 || grad.psi_yn.nx() == 0 ||
        grad.psi_yp.nx() == 0) {
        throw std::runtime_error(
            "godunovGradientNorm: needs the four one-sided derivatives -- call "
            "wenoDerivative with LSEquation::ReinitializationEqn, not LevelSetEqn");
    }

    const int nx = psi_n.nx();
    const int ny = psi_n.ny();

    Field2D G(nx, ny);
    std::vector<double> &g = G.data();

    const std::vector<double> &pn = psi_n.data();
    const std::vector<double> &pxn = grad.psi_xn.data();
    const std::vector<double> &pxp = grad.psi_xp.data();
    const std::vector<double> &pyn = grad.psi_yn.data();
    const std::vector<double> &pyp = grad.psi_yp.data();

    const size_t N = g.size();

    for (size_t n = 0; n < N; ++n) {
        const double a = pxn[n], b = pxp[n], c = pyn[n], d = pyp[n];

        double gx, gy;
        if (pn[n] >= 0.0) {
            const double aP = std::max(a, 0.0), bM = std::min(b, 0.0);
            const double cP = std::max(c, 0.0), dM = std::min(d, 0.0);
            gx = std::max(aP * aP, bM * bM);
            gy = std::max(cP * cP, dM * dM);
        } else {
            const double aM = std::min(a, 0.0), bP = std::max(b, 0.0);
            const double cM = std::min(c, 0.0), dP = std::max(d, 0.0);
            gx = std::max(aM * aM, bP * bP);
            gy = std::max(cM * cM, dP * dP);
        }

        g[n] = std::sqrt(gx + gy);
    }

    return G;
}
