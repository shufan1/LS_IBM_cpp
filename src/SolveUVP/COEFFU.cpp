#include "COEFFU.h"
#include "../Utilities/Debug.h"
#include <cstdio>
#include <stdexcept>

MomentumCoeffs coeffU(const StateVar &stateVar, const Field2D &U_star_old,
                       const Domain &domain, const DiffFluxCoeffs &diffuFluxU,
                       const ConvFluxCoeffs &convFluxU, const IBM &ibm,
                       const IBMCoeff &ibmCoeffU, const Variables &variables, const BC &bc,
                       int disc_scheme, Mat CM, Vec RHS) {
    if (disc_scheme != 2) {
        throw std::runtime_error(
            "coeffU: disc_scheme==1 or 3 (upwind/QUICK) is not implemented -- "
            "only disc_scheme==2 (central difference) is");
    }

    const int imax = domain.imax;
    const int jmax = domain.jmax;
    const Field2D &De = diffuFluxU.De;
    const Field2D &Dw = diffuFluxU.Dw;
    const Field2D &Dn = diffuFluxU.Dn;
    const Field2D &Ds = diffuFluxU.Ds;
    const Field2D &Fe = convFluxU.Fe;
    const Field2D &Fw = convFluxU.Fw;
    const Field2D &Fn = convFluxU.Fn;
    const Field2D &Fs = convFluxU.Fs;
    const Field2D &A0_p = convFluxU.A0_p;
    const Field2D &dF = convFluxU.dF;
    const Field2D &P = stateVar.P;
    const Field2D &U_prev = stateVar.U_prev;
    const double alpha_u = variables.alpha_u;


    MomentumCoeffs c;
    c.aw = Field2D(imax, jmax + 1);
    c.ae = Field2D(imax, jmax + 1);
    c.as = Field2D(imax, jmax + 1);
    c.an = Field2D(imax, jmax + 1);
    c.ap = Field2D(imax, jmax + 1);
    c.S = Field2D(imax, jmax + 1);
    c.d_SIMPLC = Field2D(imax, jmax + 1);
    c.d_piso = Field2D(imax, jmax + 1);

    // ---- Interior coefficients (central difference) ----
    // ae/aw only cover i=1..imax-3 / i=2..imax-2 respectively (one shy of
    // the full i=1..imax-2 interior on each side) -- the missing row on
    // each is exactly where the boundary correction below (S_p_e/S_p_w)
    // takes over instead of an interior link. Same story for an/as with
    // j=jmax-1/j=1 (north/south).
    for (int i = 1; i <= imax - 3; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            c.ae(i, j) = -(De(i, j) - domain.CoEWu[i] * Fe(i, j));
        }
    }
    for (int i = 2; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            c.aw(i, j) = -(Dw(i, j) + (1.0 - domain.CoEWu[i - 1]) * Fw(i, j));
        }
    }
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            c.an(i, j) = -(Dn(i, j) - domain.CoEWv[j] * Fn(i, j));
        }
    }
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 2; j <= jmax - 1; ++j) {
            c.as(i, j) = -(Ds(i, j) + (1.0 - domain.CoEWv[j - 1]) * Fs(i, j));
        }
    }

    // ---- Boundary corrections: move the missing interior link into the
    // source term S_p, gated by whether that edge is actually Dirichlet.
    // West/east act on all j; south/north act on all i.
    Field2D S_p(imax, jmax + 1);
    std::vector<double> S_p_w(jmax + 1, 0.0), S_p_e(jmax + 1, 0.0);
    std::vector<double> S_p_s(imax, 0.0), S_p_n(imax, 0.0);

    for (int j = 1; j <= jmax - 1; ++j) {
        S_p_w[j] = (Dw(1, j) + (1.0 - domain.CoEWu[0]) * Fw(1, j)) * (bc.BC_w_u == 1 ? 1.0 : 0.0);
        S_p(1, j) += S_p_w[j];

        S_p_e[j] = (De(imax - 2, j) - domain.CoEWu[imax - 2] * Fe(imax - 2, j)) * (bc.BC_e_u == 1 ? 1.0 : 0.0);
        S_p(imax - 2, j) += S_p_e[j];
    }
    // South's "1 - CoEWv*Fs" (rather than "(1-CoEWv)*Fs", like the
    // interior `as` and the north term below use) is ported exactly as
    // written in COEFFU.m -- MATLAB's `.*` binds tighter than binary `-`,
    // so `(1- CoEWv(...) .*Fs(...))` really does parse as
    // `1 - (CoEWv.*Fs)`, not `(1-CoEWv).*Fs`. Kept faithful rather than
    // "fixed" since this is meant to validate bit-for-bit against the
    // MATLAB reference.
    for (int i = 1; i <= imax - 2; ++i) {
        S_p_s[i] = 2.0 * (Ds(i, 1) + (1.0 - domain.CoEWv[1] * Fs(i, 1))) * (bc.BC_s_u == 1 ? 1.0 : 0.0);
        S_p(i, 1) += S_p_s[i];

        S_p_n[i] = 2.0 * (Dn(i, jmax - 1) - domain.CoEWv[jmax - 1] * Fn(i, jmax - 1)) * (bc.BC_n_u == 1 ? 1.0 : 0.0);
        S_p(i, jmax - 1) += S_p_n[i];
    }

    // ---- ap / d_SIMPLC / d_piso ----
    //  diagonal index, and the coefficient used in pressure correction term
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double sumA = c.aw(i, j) + c.ae(i, j) + c.as(i, j) + c.an(i, j);
            c.ap(i, j) = (-sumA + A0_p(i, j) + dF(i, j) + S_p(i, j)) / alpha_u;
            c.d_SIMPLC(i, j) = domain.dyv[j - 1] / (c.ap(i, j) + sumA); //SIMPLEC
            c.d_piso(i, j) = domain.dyv[j - 1] / c.ap(i, j);
        }
    }
    // DIAGNOSTIC (temporary): watch the last interior column's own ap and
    // dF for a sign flip -- Fe(imax-2,j) reads U(imax-1,j) (the outlet,
    // just overwritten by formUV's M_in/M_out rescale), so if that rescale
    // ever swings U(imax-1,j) far enough, dF(imax-2,j) can go large and
    // negative, driving ap(imax-2,j) non-positive and breaking diagonal
    // dominance right at this row.
    // {
    //     int i = imax - 2;
    //     int jmid = jmax / 2;
    //     fprintf(stderr, "  [coeffU] i=%d(outlet-adjacent) j=%d: dF=%e ap=%e Fe=%e U_outlet=%e\n", i,
    //             jmid, dF(i, jmid), c.ap(i, jmid), Fe(i, jmid), stateVar.U(imax - 1, jmid));
    // }

    // ---- Source term: S = S0 (old timestep) + S_Pres (pressure gradient) + S_ur (under-relaxation) + S_u (boundary values) ----
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double S0 = A0_p(i, j) * U_prev(i, j);
            double S_Pres = (P(i, j) - P(i + 1, j)) * domain.dyv[j - 1];
            double S_ur = c.ap(i, j) * (1-alpha_u) * U_star_old(i, j);
            c.S(i, j) = S0 + S_Pres + S_ur;
        }
    }
    for (int j = 1; j <= jmax - 1; ++j) {
        c.S(1, j) += bc.U_a[j] * S_p_w[j];
        c.S(imax - 2, j) += bc.U_b[j] * S_p_e[j];
    }
    for (int i = 1; i <= imax - 2; ++i) {
        c.S(i, 1) += bc.U_c[i] * S_p_s[i];
        c.S(i, jmax - 1) += bc.U_d[i] * S_p_n[i];
    }

    // DEBUG (temporary): full coefficient breakdown at the specific U
    // cells being traced (biggest coeffU discrepancies from
    // compare_matrices.m) -- ae/aw/an/as/ap plus the S0/S_Pres/S_ur pieces
    // (recomputed here since the loop above only keeps their sum) that add
    // up to the final source term c.S, BEFORE the immersed-boundary
    // solid/ghost overwrite below. Fires every SIMPLE iteration.
    // {
    //     // (i,j) here are C++ 0-based -- MATLAB-native equivalent is (i+1,j+1).
    //     int wanted[][2] = {{214, 64}, {186, 137}};
    //     for (auto &w : wanted) {
    //         int i = w[0], j = w[1];
    //         double S0 = A0_p(i, j) * U_prev(i, j);
    //         double S_Pres = (P(i, j) - P(i + 1, j)) * domain.dyv[j - 1];
    //         double S_ur = c.ap(i, j) * (1 - alpha_u) * U_star_old(i, j);
    //         printf("  [COEFFU] (%d,%d) [0-based]: ae=%.8e aw=%.8e an=%.8e as=%.8e ap=%.8e | S0=%.8e "
    //                "S_Pres=%.8e S_ur=%.8e S(final)=%.8e\n",
    //                i, j, c.ae(i, j), c.aw(i, j), c.an(i, j), c.as(i, j), c.ap(i, j), S0, S_Pres, S_ur,
    //                c.S(i, j));
    //     }
    // }

    // ---- Immersed boundary: solid cells (COEFFU.m's "Immersed Boundary
    // Treating" section, the I_solid loop). Pin ap=1, zero the
    // off-diagonal coefficients and both velocity-correction sensitivities, and set the source to a fixed constant 
    // sets S0=u_inside_psi with S_Pres=S_ur=S_u=0 there; 
    //  u_inside_psi.
    //

    // solid cell row [ 0,0,... 1,0,0...], rhs = value_inside
    // ghost cell row [ 0, -lambda_g1, 0, ..., -lambfa_g2, 0, .... 1 (at ghost cell diagonal), ...] = A1_g from IBMcoeff
    
    // update solid cell
    for (size_t idx = 0; idx < ibmCoeffU.I_solid.size(); ++idx) {
        int is = ibmCoeffU.I_solid[idx];
        int js = ibmCoeffU.J_solid[idx];
        c.ae(is, js) = 0.0;
        c.aw(is, js) = 0.0;
        c.an(is, js) = 0.0;
        c.as(is, js) = 0.0;
        c.ap(is, js) = 1.0;
        c.d_SIMPLC(is, js) = 0.0;
        c.d_piso(is, js) = 0.0;
        c.S(is, js) = ibm.u_inside_psi;
    }

    // ---- Immersed boundary: ghost cells (COEFFU.m:551-657's I_g loop,
    // update ghost cell row
    // zero out all neighbor, set diagonal = 1, 
    // set S0=RHS_U_g (=ibmCoeffU.A1_g)
    for (size_t idx = 0; idx < ibmCoeffU.I_g.size(); ++idx) {
        int ig = ibmCoeffU.I_g[idx];
        int jg = ibmCoeffU.J_g[idx];
        c.ae(ig, jg) = 0.0;
        c.aw(ig, jg) = 0.0;
        c.an(ig, jg) = 0.0;
        c.as(ig, jg) = 0.0;
        c.ap(ig, jg) = 1.0;
        c.d_SIMPLC(ig, jg) = 0.0;
        c.d_piso(ig, jg) = 0.0;
        c.S(ig, jg) = ibmCoeffU.A1_g[idx];
    }
    // The mirror point's neighbor coupling (-lambda_g_k at each of
    // I1..I6/J1..J6) update directly on PESTc CM object later

    // ---- Fill CM/RHS (caller-allocated -- see COEFFU.h) ----
    // k(i,j) is a bijection over the interior range, so each loop
    // iteration owns exactly one row of the linear system: no two
    // iterations ever write the same (row,col) entry, even though a
    // given column (e.g. k(i,j) itself) legitimately appears in several
    // different rows -- once for each neighbor whose equation
    // references it.
    auto k = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 2); };

    MatZeroEntries(CM);
    // every index in vec RHS rewrite value
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            int row = k(i, j);
            MatSetValue(CM, row, row, c.ap(i, j), ADD_VALUES);
            if (i > 1) MatSetValue(CM, row, k(i - 1, j), c.aw(i, j), ADD_VALUES);
            if (i < imax - 2) MatSetValue(CM, row, k(i + 1, j), c.ae(i, j), ADD_VALUES);
            if (j > 1) MatSetValue(CM, row, k(i, j - 1), c.as(i, j), ADD_VALUES);
            if (j < jmax - 1) MatSetValue(CM, row, k(i, j + 1), c.an(i, j), ADD_VALUES);
            VecSetValue(RHS, row, c.S(i, j), INSERT_VALUES);
        }
    }

    // DEBUG (temporary): print the (row,col,val) tuples just written to CM
    // for the wanted cells, straight from c.ap/ae/aw/an/as -- i.e. exactly
    // what the main fill loop above inserted, before the ghost-cell
    // mirror-point coupling below adds anything more. Reading this back
    // via a live MatGetRow would need an extra intermediate assembly this
    // function doesn't otherwise do; printing the already-known values
    // avoids that risk entirely.
    // {
    //     int wanted[][2] = {{214, 64}, {186, 137}};
    //     for (auto &w : wanted) {
    //         int i = w[0], j = w[1];
    //         int row = k(i, j);
    //         printf("  [CM main-loop] (%d,%d) row=%d (0-based): (col=%d,val=%.8e)", i, j, row, row,
    //                c.ap(i, j));
    //         if (i > 1) printf(" (col=%d,val=%.8e)", k(i - 1, j), c.aw(i, j));
    //         if (i < imax - 2) printf(" (col=%d,val=%.8e)", k(i + 1, j), c.ae(i, j));
    //         if (j > 1) printf(" (col=%d,val=%.8e)", k(i, j - 1), c.as(i, j));
    //         if (j < jmax - 1) printf(" (col=%d,val=%.8e)", k(i, j + 1), c.an(i, j));
    //         printf("\n");
    //     }
    // }

    // Ghost-cell mirror-point coupling (MATLAB's A_g_sparse, summed into
    // CM_u): each ghost row's standard band is already all-zero except
    // the diagonal (set above), so these extra -lambda_g_k entries at
    // the stencil neighbors' own rows are the only nonzero off-diagonal
    // terms that row ends up with -- matches
    // phi_ghost = Sum(lambda_g_k*phi_k) + A1_g.
    //
    // All MatSetValue calls on CM in this function use ADD_VALUES (CM is
    // freshly MatZeroEntries()'d above, so a single write behaves
    // identically to INSERT_VALUES) rather than INSERT_VALUES, to match
    // MATLAB's own construction: COEFFU.m's A_g_sparse is literally the
    // SUM of one sparse(...) matrix per corner (A_g_sparse =
    // sparse(A1_I_g,A1_J_g,A1_g,L,L) + sparse(A2_I_g,...) + ...), so any
    // two corners (of the same ghost cell, or of the ghost row vs the
    // standard 5-point band) that ever land on the same (row,col) get
    // correctly summed there. INSERT_VALUES would instead silently keep
    // only the last write and drop the other's contribution.
    //
    // DIAGNOSTIC (temporary): still checks whether this geometry actually
    // has any such coincidence within a single ghost cell's own corners
    // -- confirms whether the ADD_VALUES fix above changes anything for
    // this run, or is just defensive-for-other-geometries.
    if (debug::ibm_stencil) {
        int nDup = 0;
        for (size_t idx = 0; idx < ibmCoeffU.I_g.size(); ++idx) {
            int cI[6] = {ibmCoeffU.I1[idx], ibmCoeffU.I2[idx], ibmCoeffU.I3[idx], ibmCoeffU.I4[idx],
                         ibm.BQu == 1 ? ibmCoeffU.I5[idx] : -999, ibm.BQu == 1 ? ibmCoeffU.I6[idx] : -999};
            int cJ[6] = {ibmCoeffU.J1[idx], ibmCoeffU.J2[idx], ibmCoeffU.J3[idx], ibmCoeffU.J4[idx],
                         ibm.BQu == 1 ? ibmCoeffU.J5[idx] : -999, ibm.BQu == 1 ? ibmCoeffU.J6[idx] : -999};
            int n = ibm.BQu == 1 ? 6 : 4;
            for (int a = 0; a < n; ++a)
                for (int b = a + 1; b < n; ++b)
                    if (cI[a] == cI[b] && cJ[a] == cJ[b]) {
                        printf("  [coeffU DUP] ghost(%d,%d) corners %d and %d both land on (%d,%d)\n",
                               ibmCoeffU.I_g[idx], ibmCoeffU.J_g[idx], a + 1, b + 1, cI[a], cJ[a]);
                        ++nDup;
                    }
        }
        printf("  U ghost-cell self-duplicate-corner count: %d\n", nDup);
    }

    for (size_t idx = 0; idx < ibmCoeffU.I_g.size(); ++idx) {
        int row = k(ibmCoeffU.I_g[idx], ibmCoeffU.J_g[idx]);
        // DIAGNOSTIC (temporary): verify each mirror-point neighbor falls
        // within k(i,j)'s valid domain [1,imax-2]x[1,jmax-1] before using
        // it as a column index -- an out-of-range neighbor would silently
        // corrupt unrelated heap memory via MatSetValue on a
        // --with-debugging=0 PETSc build instead of erroring cleanly.
        auto checkAndSet = [&](const char *label, int ni, int nj, double val) {
            if (ni < 1 || ni > imax - 2 || nj < 1 || nj > jmax - 1) {
                fprintf(stderr,
                        "coeffU: ghost cell (%d,%d) mirror-point neighbor (%d,%d) OUT OF RANGE "
                        "[1,%d]x[1,%d]\n",
                        ibmCoeffU.I_g[idx], ibmCoeffU.J_g[idx], ni, nj, imax - 2, jmax - 1);
                return;
            }
            // DEBUG (temporary): trace exactly what gets written into
            // C++ 0-based ghost cell (186,64) / stencil corner (184,62) --
            // i.e. which lambda_g_k contribution (if any) targets this
            // specific entry. row/col printed here are PETSc's own raw
            // 0-based indices (not shifted to match MATLAB's 1-based
            // convention).
            // if (ibmCoeffU.I_g[idx] == 186 && ibmCoeffU.J_g[idx] == 64 && ni == 184 && nj == 62) {
            //     printf("  [CM contribution] ghost(%d,%d) %s -> (row=%d,col=%d, 0-based) val=%.8e\n",
            //            ibmCoeffU.I_g[idx], ibmCoeffU.J_g[idx], label, row, k(ni, nj), val);
            // }
            MatSetValue(CM, row, k(ni, nj), val, ADD_VALUES);
        };
        checkAndSet("lambda_g_1", ibmCoeffU.I1[idx], ibmCoeffU.J1[idx], -ibmCoeffU.lambda_g_1[idx]);
        checkAndSet("lambda_g_2", ibmCoeffU.I2[idx], ibmCoeffU.J2[idx], -ibmCoeffU.lambda_g_2[idx]);
        checkAndSet("lambda_g_3", ibmCoeffU.I3[idx], ibmCoeffU.J3[idx], -ibmCoeffU.lambda_g_3[idx]);
        checkAndSet("lambda_g_4", ibmCoeffU.I4[idx], ibmCoeffU.J4[idx], -ibmCoeffU.lambda_g_4[idx]);
        if (ibm.BQu == 1) {
            checkAndSet("lambda_g_5", ibmCoeffU.I5[idx], ibmCoeffU.J5[idx], -ibmCoeffU.lambda_g_5[idx]);
            checkAndSet("lambda_g_6", ibmCoeffU.I6[idx], ibmCoeffU.J6[idx], -ibmCoeffU.lambda_g_6[idx]);
        }
    }

    MatAssemblyBegin(CM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CM, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);

    // DEBUG (temporary): same wanted-cell row dump, now after the
    // ghost-cell mirror-point coupling above and final assembly -- diff
    // against the "[CM main-loop]" print above to see exactly what the
    // IBM step added/changed for each row.
    // {
    //     int wanted[][2] = {{214, 64}, {186, 137}};
    //     for (auto &w : wanted) {
    //         int row = k(w[0], w[1]);
    //         PetscInt ncols;
    //         const PetscInt *cols;
    //         const PetscScalar *vals;
    //         MatGetRow(CM, row, &ncols, &cols, &vals);
    //         printf("  [CM after-IBM] (%d,%d) row=%d (0-based):", w[0], w[1], row);
    //         for (PetscInt c2 = 0; c2 < ncols; ++c2)
    //             printf(" (col=%d,val=%.8e)", (int)cols[c2], (double)vals[c2]);
    //         printf("\n");
    //         MatRestoreRow(CM, row, &ncols, &cols, &vals);
    //     }
    // }

    // DEBUG (temporary): the single entry at C++ 0-based row=k(186,64),
    // col=k(184,62) -- MATLAB native equivalent is CM_u(25323,24523).
    // {
    //     PetscInt rowIdx = k(186, 64), colIdx = k(184, 62);
    //     double val;
    //     MatGetValues(CM, 1, &rowIdx, 1, &colIdx, &val);
    //     printf("  [CM value] CM_u(row=%d,col=%d, 0-based) = %.8e\n", (int)rowIdx, (int)colIdx, val);
    // }

    return c;
}
