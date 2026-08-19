#include "COEFFV.h"
#include <cstdio>
#include <stdexcept>

MomentumCoeffs coeffV(const StateVar &stateVar, const Field2D &V_star_old,
                       const Domain &domain, const DiffFluxCoeffs &diffuFluxV,
                       const ConvFluxCoeffs &convFluxV, const IBM &ibm,
                       const IBMCoeff &ibmCoeffV, const Variables &variables, const BC &bc,
                       int disc_scheme, Mat CM, Vec RHS) {
    if (disc_scheme != 2) {
        throw std::runtime_error(
            "coeffV: disc_scheme==1 or 3 (upwind/QUICK) is not implemented -- "
            "only disc_scheme==2 (central difference) is");
    }

    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const Field2D &De = diffuFluxV.De;
    const Field2D &Dw = diffuFluxV.Dw;
    const Field2D &Dn = diffuFluxV.Dn;
    const Field2D &Ds = diffuFluxV.Ds;
    const Field2D &Fe = convFluxV.Fe;
    const Field2D &Fw = convFluxV.Fw;
    const Field2D &Fn = convFluxV.Fn;
    const Field2D &Fs = convFluxV.Fs;
    const Field2D &A0_p = convFluxV.A0_p;
    const Field2D &dF = convFluxV.dF;
    const Field2D &P = stateVar.P;
    const Field2D &V_prev = stateVar.V_prev;
    const double alpha_v = variables.alpha_v;
    const double coef = 1.0 - alpha_v;

    MomentumCoeffs c;
    c.aw = Field2D(imax + 1, jmax);
    c.ae = Field2D(imax + 1, jmax);
    c.as = Field2D(imax + 1, jmax);
    c.an = Field2D(imax + 1, jmax);
    c.ap = Field2D(imax + 1, jmax);
    c.S = Field2D(imax + 1, jmax);
    c.d_SIMPLC = Field2D(imax + 1, jmax);
    c.d_piso = Field2D(imax + 1, jmax);

    // ---- Interior coefficients (central difference) ----
    // ae/aw only cover i=1..imax-2 / i=2..imax-1 respectively (one shy
    // of the full i=1..imax-1 interior on each side) -- the missing row
    // on each is exactly where the boundary correction below takes
    // over. Same story for an/as with j=1..jmax-3 / j=2..jmax-2.
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            c.ae(i, j) = -(De(i, j) - domain.CoNSu[i] * Fe(i, j));
        }
    }
    for (int i = 2; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            c.aw(i, j) = -(Dw(i, j) + (1.0 - domain.CoNSu[i]) * Fw(i, j));
        }
    }
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 3; ++j) {
            c.an(i, j) = -(Dn(i, j) - domain.CoNSv[j] * Fn(i, j));
        }
    }
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 2; j <= jmax - 2; ++j) {
            c.as(i, j) = -(Ds(i, j) + (1.0 - domain.CoNSv[j]) * Fs(i, j));
        }
    }

    // ---- Boundary corrections ----
    // West/east carry a factor of 2 (V borrows the ghost-padded U-axis
    // there, same reason coeffU's south/north do); south/north don't
    // (V's own real boundary rows, same reason coeffU's west/east
    // don't). Ported exactly as COEFFV.m has it, not symmetrized.
    Field2D S_p(imax + 1, jmax);
    std::vector<double> S_p_w(jmax, 0.0), S_p_e(jmax, 0.0);
    std::vector<double> S_p_s(imax + 1, 0.0), S_p_n(imax + 1, 0.0);

    for (int j = 1; j <= jmax - 2; ++j) {
        S_p_w[j] = 2.0 * (Dw(1, j) + (1.0 - domain.CoNSu[1]) * Fw(1, j)) * (bc.BC_w_v == 1 ? 1.0 : 0.0);
        S_p(1, j) += S_p_w[j];

        S_p_e[j] = 2.0 * (De(imax - 1, j) - domain.CoNSu[imax - 1] * Fe(imax - 1, j)) * (bc.BC_e_v == 1 ? 1.0 : 0.0);
        S_p(imax - 1, j) += S_p_e[j];
    }
    for (int i = 1; i <= imax - 1; ++i) {
        S_p_s[i] = (Ds(i, 1) + (1.0 - domain.CoNSv[0]) * Fs(i, 1)) * (bc.BC_s_v == 1 ? 1.0 : 0.0);
        S_p(i, 1) += S_p_s[i];

        S_p_n[i] = (Dn(i, jmax - 2) - domain.CoNSv[jmax - 2] * Fn(i, jmax - 2)) * (bc.BC_n_v == 1 ? 1.0 : 0.0);
        S_p(i, jmax - 2) += S_p_n[i];
    }

    // ---- ap / d_SIMPLC / d_piso ----
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            double sumA = c.aw(i, j) + c.ae(i, j) + c.as(i, j) + c.an(i, j);
            c.ap(i, j) = (-sumA + A0_p(i, j) + dF(i, j) + S_p(i, j)) / alpha_v;
            c.d_SIMPLC(i, j) = domain.dxu[i - 1] / (c.ap(i, j) + sumA);
            c.d_piso(i, j) = domain.dxu[i - 1] / c.ap(i, j);
        }
    }

    // ---- Source term: S = S0 + S_Pres + S_ur + S_u ----
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            double S0 = A0_p(i, j) * V_prev(i, j);
            double S_Pres = (P(i, j) - P(i, j + 1)) * domain.dxu[i - 1];
            double S_ur = c.ap(i, j) * coef * V_star_old(i, j);
            c.S(i, j) = S0 + S_Pres + S_ur;
        }
    }
    for (int j = 1; j <= jmax - 2; ++j) {
        c.S(1, j) += bc.V_a[j] * S_p_w[j];
        c.S(imax - 1, j) += bc.V_b[j] * S_p_e[j];
    }
    for (int i = 1; i <= imax - 1; ++i) {
        c.S(i, 1) += bc.V_c[i] * S_p_s[i];
        c.S(i, jmax - 2) += bc.V_d[i] * S_p_n[i];
    }

    // ---- Immersed boundary: solid cells (COEFFV.m's "Immersed Boundary
    // Treating" section, the I_solid loop). Pin ap=1, zero the
    // off-diagonal coefficients and both velocity-correction
    // sensitivities, and set the source to a fixed constant -- MATLAB
    // sets S0=phi_inside (== ibm.u_inside_psi -- COEFFV.m reuses this
    // same constant, not a separate V-specific one) with
    // S_Pres=S_ur=S_u=0 there; c.S already holds S0+S_Pres+S_ur+S_u
    // summed directly (not 4 separate arrays like MATLAB), so this
    // collapses to just c.S = u_inside_psi.
    for (size_t idx = 0; idx < ibmCoeffV.I_solid.size(); ++idx) {
        int is = ibmCoeffV.I_solid[idx];
        int js = ibmCoeffV.J_solid[idx];
        c.ae(is, js) = 0.0;
        c.aw(is, js) = 0.0;
        c.an(is, js) = 0.0;
        c.as(is, js) = 0.0;
        c.ap(is, js) = 1.0;
        c.d_SIMPLC(is, js) = 0.0;
        c.d_piso(is, js) = 0.0;
        c.S(is, js) = ibm.u_inside_psi;
    }

    // ---- Immersed boundary: ghost cells (COEFFV.m's I_g loop, BC_e_p!=1
    // branch only -- BC_e_p==1 uses a differently-sized linear system
    // this port's own k(i,j) bijection doesn't support, matching every
    // other BC_e_p==1 omission in this project). Same collapsing as the
    // solid loop: MATLAB sets S0=RHS_V_g (=ibmCoeffV.A1_g), S_Pres=
    // S_ur=S_u=0, so c.S becomes A1_g directly. The mirror point's
    // neighbor coupling (-lambda_g_k at each of I1..I6/J1..J6) can't be
    // added here yet -- it needs k(i,j), which isn't defined until
    // below, and MatZeroEntries (right below) would wipe anything set
    // here first anyway. See the "Ghost-cell mirror-point coupling"
    // block after the fill loop for the actual insertion -- both that
    // loop and the main fill loop above use ADD_VALUES (CM is freshly
    // zeroed, so a single write is unaffected) rather than
    // INSERT_VALUES, matching MATLAB's own A_g_sparse construction
    // (literally a sum of one sparse(...) matrix per corner): any
    // mirror-point neighbor that happens to coincide with a standard
    // band column, or with another of the same ghost cell's own
    // corners, gets correctly summed instead of one contribution
    // silently overwriting the other.
    for (size_t idx = 0; idx < ibmCoeffV.I_g.size(); ++idx) {
        int ig = ibmCoeffV.I_g[idx];
        int jg = ibmCoeffV.J_g[idx];
        c.ae(ig, jg) = 0.0;
        c.aw(ig, jg) = 0.0;
        c.an(ig, jg) = 0.0;
        c.as(ig, jg) = 0.0;
        c.ap(ig, jg) = 1.0;
        c.d_SIMPLC(ig, jg) = 0.0;
        c.d_piso(ig, jg) = 0.0;
        c.S(ig, jg) = ibmCoeffV.A1_g[idx];
    }

    // ---- Fill CM/RHS (caller-allocated -- see COEFFV.h) ----
    // Same bijective-row-per-cell scheme as coeffU(), just over V's own
    // interior range (i=1..imax-1, j=1..jmax-2) and shape (imax-1 rows).
    auto k = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };

    MatZeroEntries(CM);
    // every index in vec RHS rewrite value
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            int row = k(i, j);
            MatSetValue(CM, row, row, c.ap(i, j), ADD_VALUES);
            if (i > 1) MatSetValue(CM, row, k(i - 1, j), c.aw(i, j), ADD_VALUES);
            if (i < imax - 1) MatSetValue(CM, row, k(i + 1, j), c.ae(i, j), ADD_VALUES);
            if (j > 1) MatSetValue(CM, row, k(i, j - 1), c.as(i, j), ADD_VALUES);
            if (j < jmax - 2) MatSetValue(CM, row, k(i, j + 1), c.an(i, j), ADD_VALUES);
            VecSetValue(RHS, row, c.S(i, j), INSERT_VALUES);
        }
    }

    // Ghost-cell mirror-point coupling (MATLAB's A_g_sparse, summed into
    // CM_v): each ghost row's standard band is already all-zero except
    // the diagonal (set above), so these extra -lambda_g_k entries at
    // the stencil neighbors' own rows are the only nonzero off-diagonal
    // terms that row ends up with -- matches
    // phi_ghost = Sum(lambda_g_k*phi_k) + A1_g.
    for (size_t idx = 0; idx < ibmCoeffV.I_g.size(); ++idx) {
        int row = k(ibmCoeffV.I_g[idx], ibmCoeffV.J_g[idx]);
        // DIAGNOSTIC (temporary): verify each mirror-point neighbor falls
        // within k(i,j)'s valid domain [1,imax-1]x[1,jmax-2] before using
        // it as a column index -- an out-of-range neighbor would silently
        // corrupt unrelated heap memory via MatSetValue on a
        // --with-debugging=0 PETSc build instead of erroring cleanly.
        auto checkAndSet = [&](int ni, int nj, double val) {
            if (ni < 1 || ni > imax - 1 || nj < 1 || nj > jmax - 2) {
                fprintf(stderr,
                        "coeffV: ghost cell (%d,%d) mirror-point neighbor (%d,%d) OUT OF RANGE "
                        "[1,%d]x[1,%d]\n",
                        ibmCoeffV.I_g[idx], ibmCoeffV.J_g[idx], ni, nj, imax - 1, jmax - 2);
                return;
            }
            MatSetValue(CM, row, k(ni, nj), val, ADD_VALUES);
        };
        checkAndSet(ibmCoeffV.I1[idx], ibmCoeffV.J1[idx], -ibmCoeffV.lambda_g_1[idx]);
        checkAndSet(ibmCoeffV.I2[idx], ibmCoeffV.J2[idx], -ibmCoeffV.lambda_g_2[idx]);
        checkAndSet(ibmCoeffV.I3[idx], ibmCoeffV.J3[idx], -ibmCoeffV.lambda_g_3[idx]);
        checkAndSet(ibmCoeffV.I4[idx], ibmCoeffV.J4[idx], -ibmCoeffV.lambda_g_4[idx]);
        if (ibm.BQv == 1) {
            checkAndSet(ibmCoeffV.I5[idx], ibmCoeffV.J5[idx], -ibmCoeffV.lambda_g_5[idx]);
            checkAndSet(ibmCoeffV.I6[idx], ibmCoeffV.J6[idx], -ibmCoeffV.lambda_g_6[idx]);
        }
    }

    MatAssemblyBegin(CM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CM, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);

    return c;
}
