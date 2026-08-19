#include "ConvergenceResiduals.h"
#include "Debug.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace {
// max|CM*x-RHS| / max|diag(CM).*x|, with x built from phi's interior
// range via the same bijective k(i,j)=(i-iBegin)+(j-jBegin)*stride
// coeffU()/coeffV()/solveMomentumSystem() use.
double momentumResidual(Mat CM, Vec RHS, const Field2D &phi, int iBegin, int iEnd,
                         int jBegin, int jEnd, int stride, const char *label = nullptr) {
    auto k = [iBegin, jBegin, stride](int i, int j) { return (i - iBegin) + (j - jBegin) * stride; };

    Vec x;
    VecDuplicate(RHS, &x);
    for (int i = iBegin; i <= iEnd; ++i)
        for (int j = jBegin; j <= jEnd; ++j)
            VecSetValue(x, k(i, j), phi(i, j), INSERT_VALUES);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    Vec r;
    VecDuplicate(RHS, &r);
    MatMult(CM, x, r);
    VecAXPY(r, -1.0, RHS);  // r = CM*x - RHS
    VecAbs(r);
    double numerator;
    VecMax(r, nullptr, &numerator);

    Vec diag, prod;
    VecDuplicate(RHS, &diag);
    VecDuplicate(RHS, &prod);
    MatGetDiagonal(CM, diag);
    VecPointwiseMult(prod, diag, x);
    VecAbs(prod);
    double denominator;
    VecMax(prod, nullptr, &denominator);

    VecDestroy(&x);
    VecDestroy(&r);
    VecDestroy(&diag);
    VecDestroy(&prod);

    if (debug::residual && label) {
        PetscPrintf(PETSC_COMM_WORLD, "      [%s] numerator=%e denominator=%e\n", label, numerator,
                    denominator);
    }
    return numerator / denominator;
}
}  // namespace

ConvergenceResult convergenceResiduals(const Field2D &U, const Field2D &V,
                                        Vec RHS_P2, int ii, const Field2D &PCOR,
                                        Mat CM_u, Mat CM_v, Vec RHS_U, Vec RHS_V,
                                        const Domain &domain, int disc_scheme,
                                        double tol, int PISO) {
    // Unlike coeffU()/coeffV(), disc_scheme doesn't gate the residual
    // computation itself -- ConvergenceResiduals.m computes resi1-4
    // identically regardless of disc_scheme; only the message text
    // below depends on it.
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    ConvergenceResult r;
    r.resi1 = momentumResidual(CM_u, RHS_U, U, 1, imax - 2, 1, jmax - 1, imax - 2, "U");
    r.resi2 = momentumResidual(CM_v, RHS_V, V, 1, imax - 1, 1, jmax - 2, imax - 1, "V");

    VecMax(RHS_P2, nullptr, &r.resi3);

    r.resi4 = 0.0;
    for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax - 1; ++j) r.resi4 = std::max(r.resi4, std::abs(PCOR(i, j)));

    r.resi_max = std::max({r.resi1, r.resi2, r.resi3});

    if (disc_scheme != 1 && disc_scheme != 2) {
        throw std::runtime_error(
            "convergenceResiduals: disc_scheme==3 (QUICK) message formatting is not implemented -- "
            "only disc_scheme==1/2's shared text is (needs iter_qq_u/v, err_q_u/v, which nothing "
            "in this port produces)");
    }

    // Two independent checks (not if/else), matching ConvergenceResiduals.m:
    // both can in principle set `message`, though in practice they're
    // mutually exclusive except right at ii==1. MATLAB shares this exact
    // text between disc_scheme==1 and ==2 -- only ==3 differs.
    char buf[256];
    if (r.resi_max > tol || ii < 2) {
        if (PISO) {
            // MATLAB's own text calls this "converged" unconditionally on
            // PISO's one and only call, regardless of resi vs tol -- kept
            // faithful rather than "fixed", since PISO always breaks right
            // after this regardless of the actual residual.
            std::snprintf(buf, sizeof(buf),
                           "The PISO method converged, Residual U = %8.8f, V = %8.8f and mass = %8.8f ",
                           r.resi1, r.resi2, r.resi3);
        } else {
            std::snprintf(buf, sizeof(buf),
                           "At SIMPLE iter : %3d , Residual for U = %8.8f, V = %8.8f and mass = %8.8f",
                           ii, r.resi1, r.resi2, r.resi3);
        }
        r.message = buf;
    }
    if (r.resi_max < tol && ii > 1) {
        std::snprintf(
            buf, sizeof(buf),
            "The SIMPLE method converged with %3d iterations, residual for U = %8.8f, V = %8.8f and mass = %8.8f .",
            ii, r.resi1, r.resi2, r.resi3);
        r.message = buf;
    }

    return r;
}
