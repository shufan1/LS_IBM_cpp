#include "Coordinates.h"
#include <json-c/json.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace {

// MATLAB colon operator start:step:stop, restricted to the monotonic
// increasing ranges Coordinates.m actually uses. The +eps guard absorbs
// floating-point round-off in (stop-start)/step so an exact multiple
// (the common case here -- dx/dy are chosen so ranges land exactly on
// their endpoints) doesn't drop or add a spurious last point.
std::vector<double> colonRange(double start, double step, double stop) {
    double span = (stop - start) / step;
    long n = static_cast<long>(std::floor(span + 1e-9 * std::max(1.0, std::fabs(span))));
    std::vector<double> v(n + 1);
    for (long i = 0; i <= n; ++i) v[i] = start + i * step;
    return v;
}

std::vector<double> diff(const std::vector<double> &x) {
    std::vector<double> d(x.size() - 1);
    for (size_t i = 0; i + 1 < x.size(); ++i) d[i] = x[i + 1] - x[i];
    return d;
}

std::vector<double> concat(std::vector<double> a, const std::vector<double> &b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

std::vector<double> addScalar(const std::vector<double> &x, double s) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) y[i] = x[i] + s;
    return y;
}

// Mirrors Utilities/QuickInterp.m: the two 2nd-order QUICK interpolation
// weights for a face at x, given the nearest upstream node xU, the next
// upstream node xUU, and the downstream node xD.
void quickInterp(double x, double xU, double xUU, double xD, double &g1, double &g2) {
    g1 = ((x - xU) * (x - xUU)) / ((xD - xU) * (xD - xUU));
    g2 = ((x - xU) * (xD - x)) / ((xU - xUU) * (xD - xUU));
}

}  // namespace

Domain computeCoordinates(double lx, double ly, double diamcyl, const Grid &grid) {
    const double r = grid.r;
    const bool expon = grid.expon;
    const double A = grid.A;
    const double dvdxdy = grid.dvdxdy;
    const bool uniform = grid.uniform;
    const double lengthUnit = grid.lengthUnit;
    const double zoomedAreax = grid.zoomedAreax;
    const double zoomedAreay = grid.zoomedAreay;
    const double Lx_l = grid.Lx_l;
    const double Ly_b = grid.Ly_b;

    const int imax_fine = static_cast<int>(std::ceil(dvdxdy * lx / lengthUnit)) + 1;
    const int jmax_fine = static_cast<int>(std::ceil(dvdxdy * ly / lengthUnit)) + 1;
    const double dx = lx / (imax_fine - 1);
    const double dy = ly / (jmax_fine - 1);

    std::vector<double> xu, xv, yu, yv;

    if (!uniform) {
        const int n = static_cast<int>(std::ceil(std::log(A) / std::log(r)));
        std::vector<double> indx(n);
        for (int i = 0; i < n; ++i) indx[i] = i + 1;  // 1:n

        const double dx_r = dx / A;
        const double dy_r = dy / A;
        std::vector<double> xuRef = addScalar(colonRange(dx_r, dx_r, zoomedAreax * diamcyl), Lx_l);
        std::vector<double> yvRef = addScalar(colonRange(dy_r, dy_r, zoomedAreay * diamcyl), Ly_b);

        std::vector<double> x_trans, y_trans;
        double resX, resY;
        if (expon) {
            x_trans.resize(n);
            y_trans.resize(n);
            for (int i = 0; i < n; ++i) {
                double rp = std::pow(r, indx[i]);
                x_trans[i] = ((1.0 - rp) / (1.0 - r)) * dx_r;
                y_trans[i] = ((1.0 - rp) / (1.0 - r)) * dy_r;
            }
            resX = Lx_l - x_trans.back();
            resY = Ly_b - y_trans.back();
            if (resY < dx || resX < dx) {
                throw std::runtime_error(
                    "Coordinates: the transition grid exceeds the coarse grid area length. "
                    "Please reduce the coarse grid cell or A.");
            }
        } else {
            resX = Lx_l;
            resY = Ly_b;
        }

        // ---- x-axis (xu primary, xv derived) ----
        {
            int m = static_cast<int>(std::floor(resX / dx));
            double dx_mod = resX / m;
            std::vector<double> x_coarse_L(m);
            for (int i = 0; i < m; ++i) x_coarse_L[i] = i * dx_mod;

            std::vector<double> flip0 = x_trans;
            std::reverse(flip0.begin(), flip0.end());
            flip0.push_back(0.0);
            std::vector<double> x_transL(flip0.size());
            for (size_t i = 0; i < flip0.size(); ++i) x_transL[i] = Lx_l - flip0[i];

            xu = concat(concat(x_coarse_L, x_transL), xuRef);
            xu = concat(xu, addScalar(x_trans, xu.back()));

            double res = lx - xu.back();
            m = static_cast<int>(std::floor(res / dx));
            dx_mod = res / m;
            std::vector<double> x_coarse_R(m);
            for (int i = 1; i <= m; ++i) x_coarse_R[i - 1] = xu.back() + i * dx_mod;
            xu = concat(xu, x_coarse_R);

            std::vector<double> mid(xu.size() - 1);
            for (size_t i = 0; i + 1 < xu.size(); ++i) mid[i] = (xu[i] + xu[i + 1]) / 2.0;
            xv.push_back(-dx_mod / 2.0);
            xv.insert(xv.end(), mid.begin(), mid.end());
            xv.push_back(mid.back() + dx_mod);
        }

        // ---- y-axis (yv primary, yu derived) ----
        {
            int m = static_cast<int>(std::floor(resY / dy));
            double dy_mod = resY / m;
            std::vector<double> y_coarse_B(m);
            for (int i = 0; i < m; ++i) y_coarse_B[i] = i * dy_mod;

            std::vector<double> flip0 = y_trans;
            std::reverse(flip0.begin(), flip0.end());
            flip0.push_back(0.0);
            std::vector<double> y_transB(flip0.size());
            for (size_t i = 0; i < flip0.size(); ++i) y_transB[i] = Ly_b - flip0[i];

            yv = concat(concat(y_coarse_B, y_transB), yvRef);
            yv = concat(yv, addScalar(y_trans, yv.back()));

            double res = ly - yv.back();
            m = static_cast<int>(std::floor(res / dy));
            dy_mod = res / m;
            std::vector<double> y_coarse_U(m);
            for (int i = 1; i <= m; ++i) y_coarse_U[i - 1] = yv.back() + i * dy_mod;
            yv = concat(yv, y_coarse_U);

            std::vector<double> mid(yv.size() - 1);
            for (size_t i = 0; i + 1 < yv.size(); ++i) mid[i] = (yv[i] + yv[i + 1]) / 2.0;
            yu.push_back(-dy_mod / 2.0);
            yu.insert(yu.end(), mid.begin(), mid.end());
            yu.push_back(mid.back() + dy_mod);
        }
    } else {
        xu = colonRange(0.0, dx, lx);
        xv = colonRange(-dx / 2.0, dx, lx + dx / 2.0);
        yu = colonRange(-dy / 2.0, dy, ly + dy / 2.0);
        yv = colonRange(0.0, dy, ly);
    }

    Domain domain;
    domain.lx = lx;
    domain.ly = ly;
    domain.xu = xu;
    domain.yu = yu;
    domain.xv = xv;
    domain.yv = yv;
    domain.xp = xv;  // scalar nodes (p, psi, phi) coincide with xv/yu
    domain.yp = yu;

    domain.imax = static_cast<int>(domain.xp.size()) - 1;
    domain.jmax = static_cast<int>(domain.yp.size()) - 1;

    domain.dxu = diff(domain.xu);
    domain.dxv = diff(domain.xv);
    domain.dxp = diff(domain.xp);
    domain.dyu = diff(domain.yu);
    domain.dyv = diff(domain.yv);
    domain.dyp = diff(domain.yp);

    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // Central-difference interpolation weights:
    //   CoEWu=(xp(2:end-1)-xu(1:end-1))./dxu;   CoNSu=(xu-xv(1:end-1))./dxv;
    //   CoEWv=(yv(1:end)-yp(1:end-1))./dyu;     CoNSv=(yp(2:end-1)-yv(1:end-1))./dyv;
    //   CoEWp=(xu(1:end)-xp(1:end-1))./dxp;     CoNSp=(yv(1:end)-yp(1:end-1))./dyp;
    domain.CoEWu.resize(imax - 1);
    for (int i = 0; i <= imax - 2; ++i) domain.CoEWu[i] = (domain.xp[i + 1] - domain.xu[i]) / domain.dxu[i];

    domain.CoNSu.resize(imax);
    for (int i = 0; i <= imax - 1; ++i) domain.CoNSu[i] = (domain.xu[i] - domain.xv[i]) / domain.dxv[i];

    domain.CoEWv.resize(jmax);
    for (int j = 0; j <= jmax - 1; ++j) domain.CoEWv[j] = (domain.yv[j] - domain.yp[j]) / domain.dyu[j];

    domain.CoNSv.resize(jmax - 1);
    for (int j = 0; j <= jmax - 2; ++j) domain.CoNSv[j] = (domain.yp[j + 1] - domain.yv[j]) / domain.dyv[j];

    domain.CoEWp.resize(imax);
    for (int i = 0; i <= imax - 1; ++i) domain.CoEWp[i] = (domain.xu[i] - domain.xp[i]) / domain.dxp[i];

    domain.CoNSp.resize(jmax);
    for (int j = 0; j <= jmax - 1; ++j) domain.CoNSp[j] = (domain.yv[j] - domain.yp[j]) / domain.dyp[j];

    // Per-cell volumes -- each an outer product of two zero-padded 1D
    // vectors:
    //   dV_u=[0 dxv(2:end)]'*[0 dyv 0];
    //   dV_v=[0 dxu 0]'*[0 dyu(2:end-1) 0];
    //   dV_p=[0 dxu 0]'*[0 dyv 0];
    {
        std::vector<double> A(imax, 0.0);      // [0, dxv(2:end)]
        for (int i = 1; i <= imax - 1; ++i) A[i] = domain.dxv[i];

        std::vector<double> B(jmax + 1, 0.0);  // [0, dyv, 0]
        for (int j = 1; j <= jmax - 1; ++j) B[j] = domain.dyv[j - 1];

        domain.dV_u = Field2D(imax, jmax + 1);
        for (int i = 0; i < imax; ++i)
            for (int j = 0; j <= jmax; ++j)
                domain.dV_u(i, j) = A[i] * B[j];
    }
    {
        std::vector<double> Av(imax + 1, 0.0);  // [0, dxu, 0]
        for (int i = 1; i <= imax - 1; ++i) Av[i] = domain.dxu[i - 1];

        std::vector<double> Bv(jmax, 0.0);      // [0, dyu(2:end-1), 0]
        for (int j = 1; j <= jmax - 2; ++j) Bv[j] = domain.dyu[j];

        domain.dV_v = Field2D(imax + 1, jmax);
        for (int i = 0; i <= imax; ++i)
            for (int j = 0; j < jmax; ++j)
                domain.dV_v(i, j) = Av[i] * Bv[j];
    }
    {
        std::vector<double> Ap(imax + 1, 0.0);  // [0, dxu, 0] -- same as Av above
        for (int i = 1; i <= imax - 1; ++i) Ap[i] = domain.dxu[i - 1];

        std::vector<double> Bp(jmax + 1, 0.0);  // [0, dyv, 0] -- same as B above
        for (int j = 1; j <= jmax - 1; ++j) Bp[j] = domain.dyv[j - 1];

        domain.dV_p = Field2D(imax + 1, jmax + 1);
        for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax; ++j)
                domain.dV_p(i, j) = Ap[i] * Bp[j];
    }

    // Phi/scalar-transport QUICK weights (Coordinates.m's "Transport"
    // section, g1c_*/g2c_*) -- each face's 3-point upstream-biased
    // stencil, positive- and negative-flow variants. East/west vary
    // along i (a P-cell's east face sits at xu[i], west face at
    // xu[i-1]); north/south vary along j the same way with yv/yp.
    // West/south's formula is exactly east/north's evaluated one index
    // lower (a cell's west face IS its west neighbor's east face) --
    // kept as separate loops rather than derived from each other, to
    // stay a direct, checkable translation of Coordinates.m.
    domain.g1c_e_p.assign(imax + 1, 0.0);
    domain.g2c_e_p.assign(imax + 1, 0.0);
    domain.g1c_e_n.assign(imax + 1, 0.0);
    domain.g2c_e_n.assign(imax + 1, 0.0);
    for (int i = 1; i <= imax - 2; ++i) {
        quickInterp(domain.xu[i], domain.xp[i], domain.xp[i - 1], domain.xp[i + 1], domain.g1c_e_p[i],
                    domain.g2c_e_p[i]);
        quickInterp(domain.xu[i], domain.xp[i + 1], domain.xp[i + 2], domain.xp[i], domain.g1c_e_n[i],
                    domain.g2c_e_n[i]);
    }

    domain.g1c_w_p.assign(imax + 1, 0.0);
    domain.g2c_w_p.assign(imax + 1, 0.0);
    domain.g1c_w_n.assign(imax + 1, 0.0);
    domain.g2c_w_n.assign(imax + 1, 0.0);
    for (int i = 2; i <= imax - 1; ++i) {
        quickInterp(domain.xu[i - 1], domain.xp[i - 1], domain.xp[i - 2], domain.xp[i], domain.g1c_w_p[i],
                    domain.g2c_w_p[i]);
        quickInterp(domain.xu[i - 1], domain.xp[i], domain.xp[i + 1], domain.xp[i - 1], domain.g1c_w_n[i],
                    domain.g2c_w_n[i]);
    }

    domain.g1c_n_p.assign(jmax + 1, 0.0);
    domain.g2c_n_p.assign(jmax + 1, 0.0);
    domain.g1c_n_n.assign(jmax + 1, 0.0);
    domain.g2c_n_n.assign(jmax + 1, 0.0);
    for (int j = 1; j <= jmax - 2; ++j) {
        quickInterp(domain.yv[j], domain.yp[j], domain.yp[j - 1], domain.yp[j + 1], domain.g1c_n_p[j],
                    domain.g2c_n_p[j]);
        quickInterp(domain.yv[j], domain.yp[j + 1], domain.yp[j + 2], domain.yp[j], domain.g1c_n_n[j],
                    domain.g2c_n_n[j]);
    }

    domain.g1c_s_p.assign(jmax + 1, 0.0);
    domain.g2c_s_p.assign(jmax + 1, 0.0);
    domain.g1c_s_n.assign(jmax + 1, 0.0);
    domain.g2c_s_n.assign(jmax + 1, 0.0);
    for (int j = 2; j <= jmax - 1; ++j) {
        quickInterp(domain.yv[j - 1], domain.yp[j - 1], domain.yp[j - 2], domain.yp[j], domain.g1c_s_p[j],
                    domain.g2c_s_p[j]);
        quickInterp(domain.yv[j - 1], domain.yp[j], domain.yp[j + 1], domain.yp[j - 1], domain.g1c_s_n[j],
                    domain.g2c_s_n[j]);
    }

    return domain;
}

namespace {
// json-c's JSON_C_TO_STRING_PRETTY puts one array element per line, which
// is unreadable for coordinate axes with hundreds of points -- there's no
// flag to keep the top-level object pretty while collapsing just arrays,
// so this renders one array as a single compact ("[1, 2, 3]") line, using
// json-c purely for its double-formatting (consistent with everything
// else this project writes/reads via json-c).
std::string arrayLine(const std::vector<double> &v) {
    json_object *arr = json_object_new_array();
    for (double x : v) json_object_array_add(arr, json_object_new_double(x));
    std::string s = json_object_to_json_string_ext(arr, JSON_C_TO_STRING_SPACED);
    json_object_put(arr);
    return s;
}
}  // namespace

void saveCoordinatesJson(const Domain &domain, const std::string &path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("saveCoordinatesJson: failed to write " + path);
    }

    out << "{\n";
    out << "  \"imax\": " << domain.imax << ",\n";
    out << "  \"jmax\": " << domain.jmax << ",\n";
    out << "  \"lx\": " << domain.lx << ",\n";
    out << "  \"ly\": " << domain.ly << ",\n";
    out << "  \"xu\": " << arrayLine(domain.xu) << ",\n";
    out << "  \"yu\": " << arrayLine(domain.yu) << ",\n";
    out << "  \"xv\": " << arrayLine(domain.xv) << ",\n";
    out << "  \"yv\": " << arrayLine(domain.yv) << ",\n";
    out << "  \"xp\": " << arrayLine(domain.xp) << ",\n";
    out << "  \"yp\": " << arrayLine(domain.yp) << "\n";
    out << "}\n";

    if (!out) {
        throw std::runtime_error("saveCoordinatesJson: failed to write " + path);
    }
}
