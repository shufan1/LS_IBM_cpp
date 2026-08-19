#include "LSmirPointsBQ.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// LAPACK LU factorization + inverse-from-factorization -- the same pair
// of routines MATLAB's own `A = A^(-1)` calls internally, used here
// instead of a hand-written elimination so this port's ghost-cell mirror-
// point interpolation weights (N=4 bilinear / N=6 biquadratic) round the
// same way MATLAB's do, U ghostnot just to the same mathematical answer.
extern "C" {
void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
}

namespace {
// In-place inversion of a small dense N x N matrix (N=4 or 6 here).
// LAPACK is column-major and this project's matrices are row-major
// std::vector<std::vector<double>>, but no explicit transpose is needed:
// row-major storage of A is byte-for-byte identical to column-major
// storage of A^T, so handing LAPACK a row-major flatten of A makes it
// compute inv(A^T) = inv(A)^T in column-major -- which is, by the same
// row/column-major equivalence, byte-for-byte the row-major storage of
// inv(A) itself. Reading the output buffer back as row-major therefore
// gives the correct inverse of A directly.
void invertMatrix(std::vector<std::vector<double>> &A) {
    int n = static_cast<int>(A.size());
    std::vector<double> flat(static_cast<size_t>(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) flat[i * n + j] = A[i][j];

    std::vector<int> ipiv(n);
    int info = 0;
    dgetrf_(&n, &n, flat.data(), &n, ipiv.data(), &info);
    if (info != 0) {
        throw std::runtime_error(
            "LSmirPointsBQ: dgetrf_ failed (singular or near-singular mirror-point interpolation matrix)");
    }

    int lwork = n * n;
    std::vector<double> work(static_cast<size_t>(lwork));
    dgetri_(&n, flat.data(), &n, ipiv.data(), work.data(), &lwork, &info);
    if (info != 0) {
        throw std::runtime_error("LSmirPointsBQ: dgetri_ failed (singular mirror-point interpolation matrix)");
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i][j] = flat[i * n + j];
}
}  // namespace

void LSmirPointsBQ(const std::vector<double> &x, const std::vector<double> &y, double alpha,
                    double beta, double q, const std::vector<double> &X_g,
                    const std::vector<double> &Y_g, int BQ, double dx, const Field2D &psi,
                    const Field2D &nx, const Field2D &ny, IBMCoeff &ibm_coeff, bool computeA1g,
                    bool applyEdgeSafeguard, double x_0, double lx) {
    const int numg = static_cast<int>(X_g.size());
    const double Delta = std::sqrt(2.0) * dx;

    ibm_coeff.numg = numg;
    ibm_coeff.Delta = Delta;
    ibm_coeff.betaG.resize(numg);
    ibm_coeff.r_g.resize(numg);
    ibm_coeff.I_m.resize(numg);
    ibm_coeff.J_m.resize(numg);
    ibm_coeff.I1.resize(numg);
    ibm_coeff.J1.resize(numg);
    ibm_coeff.I2.resize(numg);
    ibm_coeff.J2.resize(numg);
    ibm_coeff.I3.resize(numg);
    ibm_coeff.J3.resize(numg);
    ibm_coeff.I4.resize(numg);
    ibm_coeff.J4.resize(numg);
    ibm_coeff.A1_g.resize(numg);
    ibm_coeff.lambda_g_1.resize(numg);
    ibm_coeff.lambda_g_2.resize(numg);
    ibm_coeff.lambda_g_3.resize(numg);
    ibm_coeff.lambda_g_4.resize(numg);
    if (BQ == 1) {
        ibm_coeff.I5.resize(numg);
        ibm_coeff.J5.resize(numg);
        ibm_coeff.I6.resize(numg);
        ibm_coeff.J6.resize(numg);
        ibm_coeff.lambda_g_5.resize(numg);
        ibm_coeff.lambda_g_6.resize(numg);
    }

    for (int k = 0; k < numg; ++k) {
        // -------------------------------------------------
        // 1. Project to the mirror point, store it to (Im, Jm): d = sqrt(2) * dx
        int i = ibm_coeff.I_g[k];
        int j = ibm_coeff.J_g[k];
        double nx_g = nx(i, j);
        double ny_g = ny(i, j);

        double X_ib_g = X_g[k] + nx_g * Delta;
        double X_m_g = X_ib_g + nx_g * Delta;
        double Y_ib_g = Y_g[k] + ny_g * Delta;
        double Y_m_g = Y_ib_g + ny_g * Delta;

        // find the last cell with center still less than the mirror point (x is
        int I = static_cast<int>(std::lower_bound(x.begin(), x.end(), X_m_g) - x.begin()) - 1;
        int J = static_cast<int>(std::lower_bound(y.begin(), y.end(), Y_m_g) - y.begin()) - 1;
        ibm_coeff.I_m[k] = I;
        ibm_coeff.J_m[k] = J;
    
        // -------------------------------------------------
        // 2. Find that mirror point's stencil, store it to (I1,J1), (I2,J2), ....
        // if bilinear: the 4 corners of the grid cell (I,J)-(I+1,J+1) 
        //      does not depend on which way the normal point ( as a simplication, fine for now)
        // if biquadratic
            // I1,I2: two nearest cell horiztonally away from boundary
            // I3,I4: two nearest cell vertically away from boundary
            // I5: one unit away from boundary along digital direction
            // I6: nearest to the mirror point between A) one horiztonal & two veritical OR B) two horiztonal & one veritical
        if (BQ == 0) {
            // bilinear: the 4 corners of the grid cell (I,J)-(I+1,J+1) --
            ibm_coeff.I1[k] = I;
            ibm_coeff.J1[k] = J;
            ibm_coeff.I2[k] = I + 1;
            ibm_coeff.J2[k] = J;
            ibm_coeff.I3[k] = I;
            ibm_coeff.J3[k] = J + 1;
            ibm_coeff.I4[k] = I + 1;
            ibm_coeff.J4[k] = J + 1;
        } else {  // BQ == 1, biquadratic
            int I0 = (nx_g > 0) ? I : I + 1;
            int I1v = (nx_g > 0) ? I0 + 1 : I0 - 1;
            int I2v = (nx_g > 0) ? I0 + 2 : I0 - 2;

            int J0 = (ny_g > 0) ? J : J + 1;
            int J3v = (ny_g > 0) ? J0 + 1 : J0 - 1;
            int J4v = (ny_g > 0) ? J0 + 2 : J0 - 2;

            ibm_coeff.I1[k] = I1v;
            ibm_coeff.J1[k] = J0;
            ibm_coeff.I2[k] = I2v;
            ibm_coeff.J2[k] = J0;
            ibm_coeff.I3[k] = I0;
            ibm_coeff.J3[k] = J3v;
            ibm_coeff.I4[k] = I0;
            ibm_coeff.J4[k] = J4v;
            ibm_coeff.I5[k] = I1v;
            ibm_coeff.J5[k] = J3v;

            double df = std::hypot(X_m_g - x[I1v], Y_m_g - y[J4v]);
            double ds = std::hypot(X_m_g - x[I2v], Y_m_g - y[J3v]);
            if (df < ds) {
                ibm_coeff.I6[k] = I1v;
                ibm_coeff.J6[k] = J4v;
            } else {
                ibm_coeff.I6[k] = I2v;
                ibm_coeff.J6[k] = J3v;
            }
        }

        // -------------------------------------------------
        // // 3. from bi-quadratic interpolation get vec b_k and e_k
            //  bi-quadratic interpolation from neighbors' local coordinates
            //  rows like [1, x, y, xy] (bilinear) or [1, x, y, xy, x², y²] (biquadratic)
            //  phi = c0 + c1*x + c2*y + ...
            //  neighbor phi_k = mat_A vec_c ==> c = inv(A) phi_k
            //  Phi_m = c00 + c01xm + … = Σ(inv(A)_1k phi_k) + Σ(inv(A)_2k phi_k)  xm+ Σ(inv(A)_3k phi_k) ym + ..
            //                          = Σ(inv(A)_1k+inv(A)_2k*xm + … ) phi_k
            //                          = Σb_k phi_k
            // similar for d_phi/dr_n at mirror point: dphi/dr_n = nx*r_gphi/dx + ny* dphi/dy
            //  dphi/dn_m = Σe_k * phi_k
        // -------------------------------------------------
        std::vector<double> b, e;
        if (BQ == 0) {
            // get relative coordiate around I1,J1 .
            int Ia = ibm_coeff.I1[k], Ja = ibm_coeff.J1[k];
            double x_m = X_m_g - x[Ia];
            double y_m = Y_m_g - y[Ja];

            double x1 = x[ibm_coeff.I1[k]] - x[Ia], y1 = y[ibm_coeff.J1[k]] - y[Ja];
            double x2 = x[ibm_coeff.I2[k]] - x[Ia], y2 = y[ibm_coeff.J2[k]] - y[Ja];
            double x3 = x[ibm_coeff.I3[k]] - x[Ia], y3 = y[ibm_coeff.J3[k]] - y[Ja];
            double x4 = x[ibm_coeff.I4[k]] - x[Ia], y4 = y[ibm_coeff.J4[k]] - y[Ja];

            std::vector<std::vector<double>> A = {
                {1.0, x1, y1, x1 * y1},
                {1.0, x2, y2, x2 * y2},
                {1.0, x3, y3, x3 * y3},
                {1.0, x4, y4, x4 * y4},
            };
            invertMatrix(A);  // A is now inv(A)

            b.resize(4);
            e.resize(4);
            for (int c = 0; c < 4; ++c) {
                b[c] = A[0][c] + A[1][c] * x_m + A[2][c] * y_m + A[3][c] * x_m * y_m;
                e[c] = nx_g * (A[1][c] + A[3][c] * y_m) + ny_g * (A[2][c] + A[3][c] * x_m);
            }
        } else {  // BQ == 1
            // get relative coordiate around I0,J0 .
            int I0 = (nx_g > 0) ? ibm_coeff.I_m[k] : ibm_coeff.I_m[k] + 1;
            int J0 = (ny_g > 0) ? ibm_coeff.J_m[k] : ibm_coeff.J_m[k] + 1;
            double x_m = X_m_g - x[I0];
            double y_m = Y_m_g - y[J0];

            double x1 = x[ibm_coeff.I1[k]] - x[I0], y1 = y[ibm_coeff.J1[k]] - y[J0];
            double x2 = x[ibm_coeff.I2[k]] - x[I0], y2 = y[ibm_coeff.J2[k]] - y[J0];
            double x3 = x[ibm_coeff.I3[k]] - x[I0], y3 = y[ibm_coeff.J3[k]] - y[J0];
            double x4 = x[ibm_coeff.I4[k]] - x[I0], y4 = y[ibm_coeff.J4[k]] - y[J0];
            double x5 = x[ibm_coeff.I5[k]] - x[I0], y5 = y[ibm_coeff.J5[k]] - y[J0];
            double x6 = x[ibm_coeff.I6[k]] - x[I0], y6 = y[ibm_coeff.J6[k]] - y[J0];

            std::vector<std::vector<double>> A = {
                {1.0, x1, y1, x1 * y1, x1 * x1, y1 * y1},
                {1.0, x2, y2, x2 * y2, x2 * x2, y2 * y2},
                {1.0, x3, y3, x3 * y3, x3 * x3, y3 * y3},
                {1.0, x4, y4, x4 * y4, x4 * x4, y4 * y4},
                {1.0, x5, y5, x5 * y5, x5 * x5, y5 * y5},
                {1.0, x6, y6, x6 * y6, x6 * x6, y6 * y6},
            };
            invertMatrix(A);  // A is now inv(A)

            b.resize(6);
            e.resize(6);
            for (int c = 0; c < 6; ++c) {
                b[c] = A[0][c] + A[1][c] * x_m + A[2][c] * y_m + A[3][c] * x_m * y_m +
                       A[4][c] * x_m * x_m + A[5][c] * y_m * y_m;
                e[c] = nx_g * (A[1][c] + A[3][c] * y_m + 2 * A[4][c] * x_m) +
                       ny_g * (A[2][c] + A[3][c] * x_m + 2 * A[5][c] * y_m);
            }
        }
        
        // -------------------------------------------------
        // 4. Combines with the ghost cell polynomial to produce
        //       the ghost cell equation coefficients and rhs: lambda_g_k1, lambda_g_k2,.. and A1_g
            //  polynomial around ghost cell: P(r) = a0+a1r+a2r^2+O(r^3)
            //  using robin boundary condition: -alpha*P'(r=0) = beta*P(r=0) + q
            //                                   P'(r=0) = a1, P(r=0) = a0 ->
            //                                     alpha*a1- beta*a0 = q
            //  mirror point value, mirror point derivative
            //                                     P(d) = a0 + a1*d + a2*d^2 = phi_m = Σ b_k phi_k
            //                                     P'(d) = a0 + a1*d + a2*d^2 = dphi/dn_m = Σ e_k phi_k
            //  a0 + a1*d + a2*d^2 = phi_m
            //  a0 + a1*d + a2*d^2 = dphi/dn_m
            //  alpha*a1- beta*a0 = q
            //  3 by 3 matrix solve for a0, a1, a2
            //                  -> ai = C1i'_i*phi_m + C2i'_i*r_gphi/dn_m + C3_i,
            //                      i = 0, 1, 2 each coefficient, C1, C2, C3 = some coefficients
            // at ghost cell: P(r_g) = a0 - a1*r_g + a2*r_g^2, C
            //                 phi_g = B*phi_m + E*r_gphi/dn_m + A_1 (B,E,A are scalar computed from C1i,C2i,C3i)
            // plug in intepolation on phi_m and dphi/dn_m:
            //              phi_g = B*(Σb_k phi_k) + E*(Σe_k phi_k) + A_1
            //              phi_g = Σ lambda_g_k *phi_k + A_1
        double r_g = std::abs(psi(i, j));  // distance from the ghost cell itself to the true interface
        double d = Delta;
        ibm_coeff.r_g[k] = r_g;


        // if the mirror poitn is too close to the one of the four outer edges of the computational domain
        // use 0 reactivity 
        double betaUse = beta;
        if (applyEdgeSafeguard) {
            int yLen = static_cast<int>(y.size());
            if (J >= yLen - 5 || J <= 4 || x[I] <= x_0 || x[I] >= lx - x_0) betaUse = 0.0;
        }
        ibm_coeff.betaG[k] = betaUse;

        double denom = 2.0 * alpha - betaUse * d;

        double B = (2.0 * alpha + 2.0 * betaUse * r_g + betaUse * r_g * r_g / d) / denom;
        double E = (-alpha * d - betaUse * d * r_g + (alpha / d - betaUse) * r_g * r_g) / denom;

        // MATLAB: A1_g's own formula genuinely differs between BQ==1 and BQ==0.
        // Skipped entirely when computeA1g is false -- see LSmirPointsBQ.h's
        // comment for why callers that only need lambda_g_k/geometry (not a
        // real, current q) ask for that.
        if (computeA1g) {
            ibm_coeff.A1_g[k] = (BQ == 1) ? 4.0 * Delta * q / (-betaUse * Delta + 2.0 * alpha)
                                           : q * (d + 2.0 * r_g + r_g * r_g / d) / denom;
        }

        std::vector<double> lambda(b.size());
        for (size_t c = 0; c < b.size(); ++c) lambda[c] = B * b[c] + E * e[c];

        ibm_coeff.lambda_g_1[k] = lambda[0];
        ibm_coeff.lambda_g_2[k] = lambda[1];
        ibm_coeff.lambda_g_3[k] = lambda[2];
        ibm_coeff.lambda_g_4[k] = lambda[3];
        if (BQ == 1) {
            ibm_coeff.lambda_g_5[k] = lambda[4];
            ibm_coeff.lambda_g_6[k] = lambda[5];
        }
        // -------------------------------------------------
        // ( 5. when building coefficient matrix: [ ..., -lambda_g_1,..., -lambda_g_2, .. 1 ], rhs = A_1)
    }
}
