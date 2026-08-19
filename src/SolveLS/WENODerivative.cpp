#include "WENODerivative.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// Both branches of all three functions are implemented.

WenoDifferences wenoDifferences(const Field2D &psi, const Field2D &psi_prev, const Domain &domain,
                                 double h, LSEquation equation, const Field2D *u,
                                 const Field2D *v) {
    const int nx = psi.nx();
    const int ny = psi.ny();

    // LSdirDerivates.m:13-14. Both from the P-grid spacing, and both the
    // minimum over it -- the WENO reconstruction assumes a uniform mesh
    // (the file's own header comment says so: "so far this is limited to
    // uniform grids").
    const double dx = *std::min_element(domain.dxp.begin(), domain.dxp.end());
    const double dy = *std::min_element(domain.dyp.begin(), domain.dyp.end());

    // The tube is fixed by psi_prev, the pre-advection level set, so it
    // cannot drift between RK stages (LSdirDerivates.m:31).
    const Field2D &psi_tube = psi_prev;

    WenoDifferences dir;

    if (equation == LSEquation::LevelSetEqn) {
        if (!u || !v)
            throw std::runtime_error(
                "wenoDifferences: LevelSetEqn needs u and v -- the stencil side is chosen from "
                "their sign");

        dir.vx1 = Field2D(nx, ny);
        dir.vx2 = Field2D(nx, ny);
        dir.vx3 = Field2D(nx, ny);
        dir.vx4 = Field2D(nx, ny);
        dir.vx5 = Field2D(nx, ny);
        dir.vy1 = Field2D(nx, ny);
        dir.vy2 = Field2D(nx, ny);
        dir.vy3 = Field2D(nx, ny);
        dir.vy4 = Field2D(nx, ny);
        dir.vy5 = Field2D(nx, ny);

        // ---- x-direction (lines 32-52) ----
        // i in [3, nx-4] 0-based == MATLAB 4:nx-3. The margin is the
        // 6-point stencil's reach: the u>=0 branch reaches i-3, the u<0
        // branch i+3.
        for (int i = 3; i <= nx - 4; ++i) {
            for (int j = 1; j <= ny - 2; ++j) {
                if (std::abs(psi_tube(i, j)) >= h) continue;

                if ((*u)(i, j) >= 0.0) {
                    // Upwind: information arrives from the left, lean left.
                    dir.vx1(i, j) = (psi(i - 2, j) - psi(i - 3, j)) / dx;
                    dir.vx2(i, j) = (psi(i - 1, j) - psi(i - 2, j)) / dx;
                    dir.vx3(i, j) = (psi(i, j) - psi(i - 1, j)) / dx;
                    dir.vx4(i, j) = (psi(i + 1, j) - psi(i, j)) / dx;
                    dir.vx5(i, j) = (psi(i + 2, j) - psi(i + 1, j)) / dx;
                } else {
                    // Mirror image. Everything downstream -- smoothness
                    // indicators, weights, blend -- is identical;
                    // flipping the stencil IS the upwinding.
                    dir.vx1(i, j) = (psi(i + 3, j) - psi(i + 2, j)) / dx;
                    dir.vx2(i, j) = (psi(i + 2, j) - psi(i + 1, j)) / dx;
                    dir.vx3(i, j) = (psi(i + 1, j) - psi(i, j)) / dx;
                    dir.vx4(i, j) = (psi(i, j) - psi(i - 1, j)) / dx;
                    dir.vx5(i, j) = (psi(i - 1, j) - psi(i - 2, j)) / dx;
                }
            }
        }

        // ---- y-direction (lines 53-75) ----
        // Its own loop with the margins swapped. Not merged with the x
        // loop -- the ranges genuinely differ, and a node near a boundary
        // can have vx* set while vy* stays zero.
        for (int i = 1; i <= nx - 2; ++i) {
            for (int j = 3; j <= ny - 4; ++j) {
                if (std::abs(psi_tube(i, j)) >= h) continue;

                if ((*v)(i, j) >= 0.0) {
                    dir.vy1(i, j) = (psi(i, j - 2) - psi(i, j - 3)) / dy;
                    dir.vy2(i, j) = (psi(i, j - 1) - psi(i, j - 2)) / dy;
                    dir.vy3(i, j) = (psi(i, j) - psi(i, j - 1)) / dy;
                    dir.vy4(i, j) = (psi(i, j + 1) - psi(i, j)) / dy;
                    dir.vy5(i, j) = (psi(i, j + 2) - psi(i, j + 1)) / dy;
                } else {
                    dir.vy1(i, j) = (psi(i, j + 3) - psi(i, j + 2)) / dy;
                    dir.vy2(i, j) = (psi(i, j + 2) - psi(i, j + 1)) / dy;
                    dir.vy3(i, j) = (psi(i, j + 1) - psi(i, j)) / dy;
                    dir.vy4(i, j) = (psi(i, j) - psi(i, j - 1)) / dy;
                    dir.vy5(i, j) = (psi(i, j - 1) - psi(i, j - 2)) / dy;
                }
            }
        }

    } else if (equation == LSEquation::ReinitializationEqn) {
        // ---- ReinitializationEqn (LSdirDerivates.m:87-172) ----
        // Same five differences as above, but computed BOTH ways at every
        // node instead of picking a side. There is no u to test: the
        // characteristic speed sign(psi)*psi_x/|grad psi| is not known
        // until psi_x is, so the choice is deferred to
        // godunovGradientNorm's Godunov flux. vxn* is the left-biased
        // run psi(i-3)..psi(i+2) -- identical to the u>=0 stencil above --
        // and vxp* the right-biased psi(i-2)..psi(i+3), identical to u<0.
        dir.vxn1 = Field2D(nx, ny);
        dir.vxn2 = Field2D(nx, ny);
        dir.vxn3 = Field2D(nx, ny);
        dir.vxn4 = Field2D(nx, ny);
        dir.vxn5 = Field2D(nx, ny);
        dir.vxp1 = Field2D(nx, ny);
        dir.vxp2 = Field2D(nx, ny);
        dir.vxp3 = Field2D(nx, ny);
        dir.vxp4 = Field2D(nx, ny);
        dir.vxp5 = Field2D(nx, ny);
        dir.vyn1 = Field2D(nx, ny);
        dir.vyn2 = Field2D(nx, ny);
        dir.vyn3 = Field2D(nx, ny);
        dir.vyn4 = Field2D(nx, ny);
        dir.vyn5 = Field2D(nx, ny);
        dir.vyp1 = Field2D(nx, ny);
        dir.vyp2 = Field2D(nx, ny);
        dir.vyp3 = Field2D(nx, ny);
        dir.vyp4 = Field2D(nx, ny);
        dir.vyp5 = Field2D(nx, ny);

        // ---- x-direction (lines 89-115) ----
        for (int i = 3; i <= nx - 4; ++i) {
            for (int j = 1; j <= ny - 2; ++j) {
                // Dilated tube: the minimum of |psi_prev| over the node
                // and its 8 neighbours, so the active set is one cell
                // wider than LevelSetEqn's pointwise test. The reinit
                // correction propagates outward, and a node at the band
                // edge needs computed neighbours to lean on.
                double N0 = std::abs(psi_tube(i, j));
                for (int di = -1; di <= 1; ++di)
                    for (int dj = -1; dj <= 1; ++dj)
                        N0 = std::min(N0, std::abs(psi_tube(i + di, j + dj)));
                if (N0 >= h) continue;

                dir.vxn1(i, j) = (psi(i - 2, j) - psi(i - 3, j)) / dx;
                dir.vxn2(i, j) = (psi(i - 1, j) - psi(i - 2, j)) / dx;
                dir.vxn3(i, j) = (psi(i, j) - psi(i - 1, j)) / dx;
                dir.vxn4(i, j) = (psi(i + 1, j) - psi(i, j)) / dx;
                dir.vxn5(i, j) = (psi(i + 2, j) - psi(i + 1, j)) / dx;

                dir.vxp1(i, j) = (psi(i + 3, j) - psi(i + 2, j)) / dx;
                dir.vxp2(i, j) = (psi(i + 2, j) - psi(i + 1, j)) / dx;
                dir.vxp3(i, j) = (psi(i + 1, j) - psi(i, j)) / dx;
                dir.vxp4(i, j) = (psi(i, j) - psi(i - 1, j)) / dx;
                dir.vxp5(i, j) = (psi(i - 1, j) - psi(i - 2, j)) / dx;
            }
        }

        // ---- y-direction (lines 117-141) ----
        for (int i = 1; i <= nx - 2; ++i) {
            for (int j = 3; j <= ny - 4; ++j) {
                double N0 = std::abs(psi_tube(i, j));
                for (int di = -1; di <= 1; ++di)
                    for (int dj = -1; dj <= 1; ++dj)
                        N0 = std::min(N0, std::abs(psi_tube(i + di, j + dj)));
                if (N0 >= h) continue;

                dir.vyn1(i, j) = (psi(i, j - 2) - psi(i, j - 3)) / dy;
                dir.vyn2(i, j) = (psi(i, j - 1) - psi(i, j - 2)) / dy;
                dir.vyn3(i, j) = (psi(i, j) - psi(i, j - 1)) / dy;
                dir.vyn4(i, j) = (psi(i, j + 1) - psi(i, j)) / dy;
                dir.vyn5(i, j) = (psi(i, j + 2) - psi(i, j + 1)) / dy;

                dir.vyp1(i, j) = (psi(i, j + 3) - psi(i, j + 2)) / dy;
                dir.vyp2(i, j) = (psi(i, j + 2) - psi(i, j + 1)) / dy;
                dir.vyp3(i, j) = (psi(i, j + 1) - psi(i, j)) / dy;
                dir.vyp4(i, j) = (psi(i, j) - psi(i, j - 1)) / dy;
                dir.vyp5(i, j) = (psi(i, j - 1) - psi(i, j - 2)) / dy;
            }
        }

    } else {
        throw std::runtime_error("wenoDifferences: unknown LSEquation");
    }

    return dir;
}

WenoWeights wenoWeights(const WenoDifferences &dir, double eps, LSEquation equation) {
    // Ideal (linear) weights. At exactly these values the three
    // 3rd-order candidates' errors cancel and the blend is 5th order --
    // the nonlinear weights below only drift off them near a kink.
    const double d1 = 1.0 / 10.0;
    const double d2 = 6.0 / 10.0;
    const double d3 = 3.0 / 10.0;

    const double C = 13.0 / 12.0;

    WenoWeights w;

    if (equation == LSEquation::LevelSetEqn) {
        const int nx = dir.vx1.nx();
        const int ny = dir.vx1.ny();

        w.wx1 = Field2D(nx, ny);
        w.wx2 = Field2D(nx, ny);
        w.wx3 = Field2D(nx, ny);
        w.wy1 = Field2D(nx, ny);
        w.wy2 = Field2D(nx, ny);
        w.wy3 = Field2D(nx, ny);

        const std::vector<double> &vx1 = dir.vx1.data(), &vx2 = dir.vx2.data();
        const std::vector<double> &vx3 = dir.vx3.data(), &vx4 = dir.vx4.data();
        const std::vector<double> &vx5 = dir.vx5.data();
        const std::vector<double> &vy1 = dir.vy1.data(), &vy2 = dir.vy2.data();
        const std::vector<double> &vy3 = dir.vy3.data(), &vy4 = dir.vy4.data();
        const std::vector<double> &vy5 = dir.vy5.data();

        std::vector<double> &wx1 = w.wx1.data(), &wx2 = w.wx2.data(), &wx3 = w.wx3.data();
        std::vector<double> &wy1 = w.wy1.data(), &wy2 = w.wy2.data(), &wy3 = w.wy3.data();

        const size_t N = wx1.size();
        for (size_t n = 0; n < N; ++n) {
            // ---- x (LSWeights.m:29-31, 38-40, 46-48) ----
            const double a1 = vx1[n], a2 = vx2[n], a3 = vx3[n], a4 = vx4[n], a5 = vx5[n];
            const double px1 = a1 - 2.0 * a2 + a3, qx1 = a1 - 4.0 * a2 + 3.0 * a3;
            const double px2 = a2 - 2.0 * a3 + a4, qx2 = a2 - a4;
            const double px3 = a3 - 2.0 * a4 + a5, qx3 = 3.0 * a3 - 4.0 * a4 + a5;
            const double tx1 = eps + C * px1 * px1 + 0.25 * qx1 * qx1;
            const double tx2 = eps + C * px2 * px2 + 0.25 * qx2 * qx2;
            const double tx3 = eps + C * px3 * px3 + 0.25 * qx3 * qx3;
            const double ax1 = d1 / (tx1 * tx1), ax2 = d2 / (tx2 * tx2), ax3 = d3 / (tx3 * tx3);
            const double axSum = ax1 + ax2 + ax3;
            wx1[n] = ax1 / axSum;
            wx2[n] = ax2 / axSum;
            wx3[n] = ax3 / axSum;

            // ---- y (LSWeights.m:33-35, 42-44, 50-52) ----
            const double b1 = vy1[n], b2 = vy2[n], b3 = vy3[n], b4 = vy4[n], b5 = vy5[n];
            const double py1 = b1 - 2.0 * b2 + b3, qy1 = b1 - 4.0 * b2 + 3.0 * b3;
            const double py2 = b2 - 2.0 * b3 + b4, qy2 = b2 - b4;
            const double py3 = b3 - 2.0 * b4 + b5, qy3 = 3.0 * b3 - 4.0 * b4 + b5;
            const double ty1 = eps + C * py1 * py1 + 0.25 * qy1 * qy1;
            const double ty2 = eps + C * py2 * py2 + 0.25 * qy2 * qy2;
            const double ty3 = eps + C * py3 * py3 + 0.25 * qy3 * qy3;
            const double ay1 = d1 / (ty1 * ty1), ay2 = d2 / (ty2 * ty2), ay3 = d3 / (ty3 * ty3);
            const double aySum = ay1 + ay2 + ay3;
            wy1[n] = ay1 / aySum;
            wy2[n] = ay2 / aySum;
            wy3[n] = ay3 / aySum;
        }

    } else if (equation == LSEquation::ReinitializationEqn) {
        // ---- LSWeights.m:63-155 ----
        // Identical arithmetic to the LevelSetEqn branch, run on four
        // groups of five differences instead of two: the left- and
        // right-biased runs in each direction. No group is preferred
        // here -- all four sets of weights survive, and the choice
        // between the resulting derivatives happens later, in
        // godunovGradientNorm's Godunov flux.
        const int nx = dir.vxn1.nx();
        const int ny = dir.vxn1.ny();

        w.wxn1 = Field2D(nx, ny);
        w.wxn2 = Field2D(nx, ny);
        w.wxn3 = Field2D(nx, ny);
        w.wxp1 = Field2D(nx, ny);
        w.wxp2 = Field2D(nx, ny);
        w.wxp3 = Field2D(nx, ny);
        w.wyn1 = Field2D(nx, ny);
        w.wyn2 = Field2D(nx, ny);
        w.wyn3 = Field2D(nx, ny);
        w.wyp1 = Field2D(nx, ny);
        w.wyp2 = Field2D(nx, ny);
        w.wyp3 = Field2D(nx, ny);

        const std::vector<double> &vxn1 = dir.vxn1.data(), &vxn2 = dir.vxn2.data();
        const std::vector<double> &vxn3 = dir.vxn3.data(), &vxn4 = dir.vxn4.data();
        const std::vector<double> &vxn5 = dir.vxn5.data();
        const std::vector<double> &vxp1 = dir.vxp1.data(), &vxp2 = dir.vxp2.data();
        const std::vector<double> &vxp3 = dir.vxp3.data(), &vxp4 = dir.vxp4.data();
        const std::vector<double> &vxp5 = dir.vxp5.data();
        const std::vector<double> &vyn1 = dir.vyn1.data(), &vyn2 = dir.vyn2.data();
        const std::vector<double> &vyn3 = dir.vyn3.data(), &vyn4 = dir.vyn4.data();
        const std::vector<double> &vyn5 = dir.vyn5.data();
        const std::vector<double> &vyp1 = dir.vyp1.data(), &vyp2 = dir.vyp2.data();
        const std::vector<double> &vyp3 = dir.vyp3.data(), &vyp4 = dir.vyp4.data();
        const std::vector<double> &vyp5 = dir.vyp5.data();

        std::vector<double> &wxn1 = w.wxn1.data(), &wxn2 = w.wxn2.data(), &wxn3 = w.wxn3.data();
        std::vector<double> &wxp1 = w.wxp1.data(), &wxp2 = w.wxp2.data(), &wxp3 = w.wxp3.data();
        std::vector<double> &wyn1 = w.wyn1.data(), &wyn2 = w.wyn2.data(), &wyn3 = w.wyn3.data();
        std::vector<double> &wyp1 = w.wyp1.data(), &wyp2 = w.wyp2.data(), &wyp3 = w.wyp3.data();

        const size_t N = wxn1.size();
        for (size_t n = 0; n < N; ++n) {
            // ---- x, left-biased (LSWeights.m:90-92, 107-109, 123-125) ----
            const double a1 = vxn1[n], a2 = vxn2[n], a3 = vxn3[n], a4 = vxn4[n], a5 = vxn5[n];
            const double pn1 = a1 - 2.0 * a2 + a3, qn1 = a1 - 4.0 * a2 + 3.0 * a3;
            const double pn2 = a2 - 2.0 * a3 + a4, qn2 = a2 - a4;
            const double pn3 = a3 - 2.0 * a4 + a5, qn3 = 3.0 * a3 - 4.0 * a4 + a5;
            const double txn1 = eps + C * pn1 * pn1 + 0.25 * qn1 * qn1;
            const double txn2 = eps + C * pn2 * pn2 + 0.25 * qn2 * qn2;
            const double txn3 = eps + C * pn3 * pn3 + 0.25 * qn3 * qn3;
            const double axn1 = d1 / (txn1 * txn1), axn2 = d2 / (txn2 * txn2),
                         axn3 = d3 / (txn3 * txn3);
            const double axnSum = axn1 + axn2 + axn3;
            wxn1[n] = axn1 / axnSum;
            wxn2[n] = axn2 / axnSum;
            wxn3[n] = axn3 / axnSum;

            // ---- x, right-biased (LSWeights.m:94-96, 111-113, 127-129) ----
            const double b1 = vxp1[n], b2 = vxp2[n], b3 = vxp3[n], b4 = vxp4[n], b5 = vxp5[n];
            const double pp1 = b1 - 2.0 * b2 + b3, qp1 = b1 - 4.0 * b2 + 3.0 * b3;
            const double pp2 = b2 - 2.0 * b3 + b4, qp2 = b2 - b4;
            const double pp3 = b3 - 2.0 * b4 + b5, qp3 = 3.0 * b3 - 4.0 * b4 + b5;
            const double txp1 = eps + C * pp1 * pp1 + 0.25 * qp1 * qp1;
            const double txp2 = eps + C * pp2 * pp2 + 0.25 * qp2 * qp2;
            const double txp3 = eps + C * pp3 * pp3 + 0.25 * qp3 * qp3;
            const double axp1 = d1 / (txp1 * txp1), axp2 = d2 / (txp2 * txp2),
                         axp3 = d3 / (txp3 * txp3);
            const double axpSum = axp1 + axp2 + axp3;
            wxp1[n] = axp1 / axpSum;
            wxp2[n] = axp2 / axpSum;
            wxp3[n] = axp3 / axpSum;

            // ---- y, left-biased (LSWeights.m:98-100, 115-117, 131-133) ----
            const double c1 = vyn1[n], c2 = vyn2[n], c3 = vyn3[n], c4 = vyn4[n], c5 = vyn5[n];
            const double rn1 = c1 - 2.0 * c2 + c3, sn1 = c1 - 4.0 * c2 + 3.0 * c3;
            const double rn2 = c2 - 2.0 * c3 + c4, sn2 = c2 - c4;
            const double rn3 = c3 - 2.0 * c4 + c5, sn3 = 3.0 * c3 - 4.0 * c4 + c5;
            const double tyn1 = eps + C * rn1 * rn1 + 0.25 * sn1 * sn1;
            const double tyn2 = eps + C * rn2 * rn2 + 0.25 * sn2 * sn2;
            const double tyn3 = eps + C * rn3 * rn3 + 0.25 * sn3 * sn3;
            const double ayn1 = d1 / (tyn1 * tyn1), ayn2 = d2 / (tyn2 * tyn2),
                         ayn3 = d3 / (tyn3 * tyn3);
            const double aynSum = ayn1 + ayn2 + ayn3;
            wyn1[n] = ayn1 / aynSum;
            wyn2[n] = ayn2 / aynSum;
            wyn3[n] = ayn3 / aynSum;

            // ---- y, right-biased (LSWeights.m:102-104, 119-121, 135-137) ----
            const double e1 = vyp1[n], e2 = vyp2[n], e3 = vyp3[n], e4 = vyp4[n], e5 = vyp5[n];
            const double rp1 = e1 - 2.0 * e2 + e3, sp1 = e1 - 4.0 * e2 + 3.0 * e3;
            const double rp2 = e2 - 2.0 * e3 + e4, sp2 = e2 - e4;
            const double rp3 = e3 - 2.0 * e4 + e5, sp3 = 3.0 * e3 - 4.0 * e4 + e5;
            const double typ1 = eps + C * rp1 * rp1 + 0.25 * sp1 * sp1;
            const double typ2 = eps + C * rp2 * rp2 + 0.25 * sp2 * sp2;
            const double typ3 = eps + C * rp3 * rp3 + 0.25 * sp3 * sp3;
            const double ayp1 = d1 / (typ1 * typ1), ayp2 = d2 / (typ2 * typ2),
                         ayp3 = d3 / (typ3 * typ3);
            const double aypSum = ayp1 + ayp2 + ayp3;
            wyp1[n] = ayp1 / aypSum;
            wyp2[n] = ayp2 / aypSum;
            wyp3[n] = ayp3 / aypSum;
        }

    } else {
        throw std::runtime_error("wenoWeights: unknown LSEquation");
    }

    return w;
}


LSGradient wenoDerivative(const Field2D &psi, const Field2D &psi_prev, const Domain &domain, double h,
                           LSEquation equation, const Field2D *u, const Field2D *v) {
    // Same three-step shape as LSFindDerivative.m:16-48 -- differences,
    // weights, blend. This function owns none of the numerics itself.
    const WenoDifferences dir = wenoDifferences(psi, psi_prev, domain, h, equation, u, v);

    // LSFindDerivative.m:17. Guards the 1/(eps+S)^2 division when a
    // sub-stencil is perfectly smooth (S == 0), which is the normal case
    // outside the tube where every difference is zero.
    const double eps = 1e-6;
    const WenoWeights w = wenoWeights(dir, eps, equation);

    LSGradient grad;
    const int nx = psi.nx();
    const int ny = psi.ny();

    if (equation == LSEquation::LevelSetEqn) {
        // ---- LSFindDerivative.m:44-48 ----
        // The three bracketed terms are the three candidate 3rd-order
        // estimates of the derivative, each from its own 3-point
        // sub-stencil.
        grad.psi_x = Field2D(nx, ny);
        grad.psi_y = Field2D(nx, ny);

        std::vector<double> &px = grad.psi_x.data();
        std::vector<double> &py = grad.psi_y.data();

        const std::vector<double> &vx1 = dir.vx1.data(), &vx2 = dir.vx2.data();
        const std::vector<double> &vx3 = dir.vx3.data(), &vx4 = dir.vx4.data();
        const std::vector<double> &vx5 = dir.vx5.data();
        const std::vector<double> &wx1 = w.wx1.data(), &wx2 = w.wx2.data();
        const std::vector<double> &wx3 = w.wx3.data();

        const std::vector<double> &vy1 = dir.vy1.data(), &vy2 = dir.vy2.data();
        const std::vector<double> &vy3 = dir.vy3.data(), &vy4 = dir.vy4.data();
        const std::vector<double> &vy5 = dir.vy5.data();
        const std::vector<double> &wy1 = w.wy1.data(), &wy2 = w.wy2.data();
        const std::vector<double> &wy3 = w.wy3.data();

        const size_t N = px.size();
        for (size_t n = 0; n < N; ++n) {
            px[n] = wx1[n] * (vx1[n] / 3.0 - 7.0 / 6.0 * vx2[n] + 11.0 / 6.0 * vx3[n]) +
                    wx2[n] * (-vx2[n] / 6.0 + 5.0 / 6.0 * vx3[n] + vx4[n] / 3.0) +
                    wx3[n] * (vx3[n] / 3.0 + 5.0 / 6.0 * vx4[n] - vx5[n] / 6.0);

            py[n] = wy1[n] * (vy1[n] / 3.0 - 7.0 / 6.0 * vy2[n] + 11.0 / 6.0 * vy3[n]) +
                    wy2[n] * (-vy2[n] / 6.0 + 5.0 / 6.0 * vy3[n] + vy4[n] / 3.0) +
                    wy3[n] * (vy3[n] / 3.0 + 5.0 / 6.0 * vy4[n] - vy5[n] / 6.0);
        }

    } else if (equation == LSEquation::ReinitializationEqn) {
        // ---- LSFindDerivative.m:95-110 ----
        // The same blend as the LevelSetEqn branch above, run four times
        // instead of twice: once per side, per direction. Nothing is
        // selected here -- all four survive, and godunovGradientNorm
        // picks between them pointwise.
        grad.psi_xn = Field2D(nx, ny);
        grad.psi_xp = Field2D(nx, ny);
        grad.psi_yn = Field2D(nx, ny);
        grad.psi_yp = Field2D(nx, ny);

        std::vector<double> &pxn = grad.psi_xn.data();
        std::vector<double> &pxp = grad.psi_xp.data();
        std::vector<double> &pyn = grad.psi_yn.data();
        std::vector<double> &pyp = grad.psi_yp.data();

        const std::vector<double> &vxn1 = dir.vxn1.data(), &vxn2 = dir.vxn2.data();
        const std::vector<double> &vxn3 = dir.vxn3.data(), &vxn4 = dir.vxn4.data();
        const std::vector<double> &vxn5 = dir.vxn5.data();
        const std::vector<double> &wxn1 = w.wxn1.data(), &wxn2 = w.wxn2.data();
        const std::vector<double> &wxn3 = w.wxn3.data();

        const std::vector<double> &vxp1 = dir.vxp1.data(), &vxp2 = dir.vxp2.data();
        const std::vector<double> &vxp3 = dir.vxp3.data(), &vxp4 = dir.vxp4.data();
        const std::vector<double> &vxp5 = dir.vxp5.data();
        const std::vector<double> &wxp1 = w.wxp1.data(), &wxp2 = w.wxp2.data();
        const std::vector<double> &wxp3 = w.wxp3.data();

        const std::vector<double> &vyn1 = dir.vyn1.data(), &vyn2 = dir.vyn2.data();
        const std::vector<double> &vyn3 = dir.vyn3.data(), &vyn4 = dir.vyn4.data();
        const std::vector<double> &vyn5 = dir.vyn5.data();
        const std::vector<double> &wyn1 = w.wyn1.data(), &wyn2 = w.wyn2.data();
        const std::vector<double> &wyn3 = w.wyn3.data();

        const std::vector<double> &vyp1 = dir.vyp1.data(), &vyp2 = dir.vyp2.data();
        const std::vector<double> &vyp3 = dir.vyp3.data(), &vyp4 = dir.vyp4.data();
        const std::vector<double> &vyp5 = dir.vyp5.data();
        const std::vector<double> &wyp1 = w.wyp1.data(), &wyp2 = w.wyp2.data();
        const std::vector<double> &wyp3 = w.wyp3.data();

        const size_t Nr = pxn.size();
        for (size_t n = 0; n < Nr; ++n) {
            pxn[n] = wxn1[n] * (vxn1[n] / 3.0 - 7.0 / 6.0 * vxn2[n] + 11.0 / 6.0 * vxn3[n]) +
                     wxn2[n] * (-vxn2[n] / 6.0 + 5.0 / 6.0 * vxn3[n] + vxn4[n] / 3.0) +
                     wxn3[n] * (vxn3[n] / 3.0 + 5.0 / 6.0 * vxn4[n] - vxn5[n] / 6.0);

            pxp[n] = wxp1[n] * (vxp1[n] / 3.0 - 7.0 / 6.0 * vxp2[n] + 11.0 / 6.0 * vxp3[n]) +
                     wxp2[n] * (-vxp2[n] / 6.0 + 5.0 / 6.0 * vxp3[n] + vxp4[n] / 3.0) +
                     wxp3[n] * (vxp3[n] / 3.0 + 5.0 / 6.0 * vxp4[n] - vxp5[n] / 6.0);

            pyn[n] = wyn1[n] * (vyn1[n] / 3.0 - 7.0 / 6.0 * vyn2[n] + 11.0 / 6.0 * vyn3[n]) +
                     wyn2[n] * (-vyn2[n] / 6.0 + 5.0 / 6.0 * vyn3[n] + vyn4[n] / 3.0) +
                     wyn3[n] * (vyn3[n] / 3.0 + 5.0 / 6.0 * vyn4[n] - vyn5[n] / 6.0);

            pyp[n] = wyp1[n] * (vyp1[n] / 3.0 - 7.0 / 6.0 * vyp2[n] + 11.0 / 6.0 * vyp3[n]) +
                     wyp2[n] * (-vyp2[n] / 6.0 + 5.0 / 6.0 * vyp3[n] + vyp4[n] / 3.0) +
                     wyp3[n] * (vyp3[n] / 3.0 + 5.0 / 6.0 * vyp4[n] - vyp5[n] / 6.0);
        }

    } else {
        throw std::runtime_error("wenoDerivative: unknown LSEquation");
    }

    return grad;
}
