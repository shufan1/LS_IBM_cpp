#include "solveHJEq.h"
#include "WENODerivative.h"
#include <stdexcept>
#include <vector>

Field2D solveHJEq(const Field2D &psi, const Field2D &u, const Field2D &v, const Field2D &psi_prev,
                  const Domain &domain, const Variables &variables) {
    // solveHJEq.m:3-5. MATLAB overwrites VARIABLES.dtLS with VARIABLES.dt
    // on entry, discarding whatever defineLSvariables() had put there
    // (the nLSupdate-scaled alternative is present but commented out).
    // At the point LSeqSolve is reached, VARIABLES.dt is still the BIG dt
    // -- the /nLSupdate division happens later, around the flow/transport
    // call, and is undone afterwards.
    const double dt = variables.dt;
    const double h = variables.LSgamma;
    const std::string &scheme = variables.TimeSchemeLS;

    // Every field here has the identical (imax+1)x(jmax+1) shape, so the
    // stage updates run straight over the flat buffer -- the 2D indices
    // are never needed, and Field2D has no arithmetic operators.
    const int Nx = psi.nx();
    const int Ny = psi.ny();
    const std::vector<double> &psi_in = psi.data();
    const std::vector<double> &uu = u.data();
    const std::vector<double> &vv = v.data();
    const size_t N = psi_in.size();

    // ---- Stage 1 (lines 16-19), common to all three schemes ----
    //   psi_n1 = psi_n - dt*( u*psi_x + v*psi_y )
    LSGradient grad1 = wenoDerivative(psi, psi_prev, domain, h, LSEquation::LevelSetEqn, &u, &v);
    const std::vector<double> &gx1 = grad1.psi_x.data();
    const std::vector<double> &gy1 = grad1.psi_y.data();

    Field2D psi_n1(Nx, Ny);
    std::vector<double> &p1 = psi_n1.data();
    for (size_t n = 0; n < N; ++n) p1[n] = psi_in[n] - dt * (uu[n] * gx1[n] + vv[n] * gy1[n]);

    if (scheme == "RK1") return psi_n1;  // line 23

    // ---- Stage 2 (lines 27-31 for RK2, 37-39 for RK3) ----
    // Identical in both schemes, so hoisted above the split. RK1 returns
    // before reaching it, so nothing extra is computed there.
    //   psi_n2 = psi_n1 - dt*( u*psi_x + v*psi_y )   evaluated at psi_n1
    LSGradient grad2 = wenoDerivative(psi_n1, psi_prev, domain, h, LSEquation::LevelSetEqn, &u, &v);
    const std::vector<double> &gx2 = grad2.psi_x.data();
    const std::vector<double> &gy2 = grad2.psi_y.data();

    Field2D psi_n2(Nx, Ny);
    std::vector<double> &p2 = psi_n2.data();
    for (size_t n = 0; n < N; ++n) p2[n] = p1[n] - dt * (uu[n] * gx2[n] + vv[n] * gy2[n]);

    if (scheme == "RK2") {
        //   psi = 0.5*( psi_n + psi_n2 )                          line 33
        Field2D result(Nx, Ny);
        std::vector<double> &r = result.data();
        for (size_t n = 0; n < N; ++n) r[n] = 0.5 * (psi_in[n] + p2[n]);
        return result;
    }

    if (scheme == "RK3") {
        //   psi_n12 = 0.75*psi_n + 0.25*psi_n2                    line 41
        Field2D psi_n12(Nx, Ny);
        std::vector<double> &p12 = psi_n12.data();
        for (size_t n = 0; n < N; ++n) p12[n] = 0.75 * psi_in[n] + 0.25 * p2[n];

        // BUG-FOR-BUG -- see the block comment in solveHJEq.h.
        //
        // MATLAB line 43 evaluates the third stage's derivative at
        // psi_n1, an exact duplicate of the line-37 call above, when
        // SSP-RK3 requires L(psi_n12). Reproduced on purpose so the port
        // matches the reference bit-for-bit during validation.
        //
        // grad2 already holds exactly that duplicate call's result --
        // same u, v, psi_n1, psi_prev, domain, h -- so gx2/gy2 are reused
        // below instead of recomputed. That skips a full WENO sweep for
        // provably identical values; it is not a behavioural divergence.
        //
        // THE FIX, once validation passes -- evaluate at psi_n12 here,
        // and fix solveHJEq.m:43 in the same commit:
        //     LSGradient grad3 = wenoDerivative(psi_n12, psi_prev, domain, h, LSEquation::LevelSetEqn, &u, &v);
        //   then use grad3.psi_x/psi_y in the loop below instead of gx2/gy2.
        //
        //   psi_n32 = psi_n12 - dt*( u*psi_x + v*psi_y )          line 45
        Field2D psi_n32(Nx, Ny);
        std::vector<double> &p32 = psi_n32.data();
        for (size_t n = 0; n < N; ++n) p32[n] = p12[n] - dt * (uu[n] * gx2[n] + vv[n] * gy2[n]);

        //   psi = ( psi_n + 2*psi_n32 )/3                         line 47
        Field2D result(Nx, Ny);
        std::vector<double> &r = result.data();
        for (size_t n = 0; n < N; ++n) r[n] = (psi_in[n] + 2.0 * p32[n]) / 3.0;
        return result;
    }

    // MATLAB falls through its if/elseif chain and returns whatever `psi`
    // happened to be -- silently, and in a way that looks like a
    // converged answer. Fail loudly instead.
    throw std::runtime_error("solveHJEq: unknown TimeSchemeLS \"" + scheme +
                              "\" -- expected \"RK1\", \"RK2\" or \"RK3\"");
}
