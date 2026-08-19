#include "CoeffPhiADRE.h"

TransportCoeffs coeffPhiADRE(const Domain &domain, const Flux &flux, const ConvFluxCoeffs &convFluxPhi,
                              const IBMCoeff &ibmCoeffPhi, const IBM &ibm, const Variables &variables,
                              Mat CM) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    const Field2D &De = flux.Diffu_Phi.De;
    const Field2D &Dw = flux.Diffu_Phi.Dw;
    const Field2D &Dn = flux.Diffu_Phi.Dn;
    const Field2D &Ds = flux.Diffu_Phi.Ds;
    const std::vector<double> &D1_a = flux.Diffu_Phi.D1_a;
    const std::vector<double> &D2_a = flux.Diffu_Phi.D2_a;
    const Field2D &Fe = convFluxPhi.Fe;
    const Field2D &Fw = convFluxPhi.Fw;
    const Field2D &Fn = convFluxPhi.Fn;
    const Field2D &Fs = convFluxPhi.Fs;
    const Field2D &A0_p = convFluxPhi.A0_p;
    const Field2D &dF = convFluxPhi.dF;
    const Field2D &alphae = convFluxPhi.alphae;
    const Field2D &alphaw = convFluxPhi.alphaw;
    const Field2D &alphan = convFluxPhi.alphan;
    const Field2D &alphas = convFluxPhi.alphas;
    const double alpha_q = variables.alpha_q;

    Field2D ae(imax + 1, jmax + 1);
    Field2D aw(imax + 1, jmax + 1);
    Field2D an(imax + 1, jmax + 1);
    Field2D as(imax + 1, jmax + 1);
    Field2D ap(imax + 1, jmax + 1);

    // ---- Interior coefficients (COEFFPHIADRE.m's "Inner Points" block) ----
    // interior block is i∈[1,imax-1], j∈[1,jmax-1]
    // ae/aw: i in [2, imax-2], j in [1, jmax-1] (0-indexed) -- one shy of
    // the full [1, imax-1] interior on each side; the missing west
    // (i=1) and east (i=imax-1) columns are the boundary-corrected rows
    // (D1_a/D2_a west correction, AW_e east correction) filled below.
    for (int i = 2; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            ae(i, j) = -(De(i, j) - (1.0 - alphae(i, j)) * Fe(i, j));
            // will add ae on the right most interior later:  i = 2
            aw(i, j) = -(Dw(i, j) + alphaw(i, j) * Fw(i, j));
            // will add aw on the left most interior later: imax-1
        }
    }
    // an/as: i in [1, imax-1], j in [2, jmax-2] -- likewise one shy on
    // the south (j=1) and north (j=jmax-1) rows (AN_s/AS_n corrections,
    // filled below).
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 2; j <= jmax - 2; ++j) {
            an(i, j) = -(Dn(i, j) - (1.0 - alphan(i, j)) * Fn(i, j));
            // will add an on the bottom most interior later:  j = 1
            as(i, j) = -(Ds(i, j) + alphas(i, j) * Fs(i, j));
            // will add an on the top most interior later:  j = jmax-1
        }
    }

    // ---- West boundary (i=1): Dirichlet inlet (BC=1). No west neighbor
    // exists (aw(1,j) stays 0), so its diffusive contribution is instead
    // folded into ae(1,j) via the one-sided D1_a correction, and the part
    // that multiplies the known boundary value (D2_a, plus the ordinary
    // convective Fw) is captured in S_p for the ap assembly below --
    // RHSPHIADRE.m/rhsPhiADRE() reuses the same Fw+D2_a combination,
    // multiplied by the actual boundary phi value, for the RHS.
    Field2D S_p(imax + 1, jmax + 1);
    for (int j = 1; j <= jmax - 1; ++j) {
        ae(1, j) = -(De(1, j) - (1.0 - alphae(1, j)) * Fe(1, j) - D1_a[j]);
        S_p(1, j) = Fw(1, j) + D2_a[j]; // a boundary-only correction field
    }

    // ---- East boundary (i=imax-1): zero-diffusive outlet (BC=3). No
    // east neighbor (ae(imax-1,j) stays 0, no correction needed); aw at
    // this row just extends the ordinary interior formula one column
    // further, since its own west neighbor (i=imax-2) is a real interior
    // node.
    for (int j = 1; j <= jmax - 1; ++j) {
        aw(imax - 1, j) = -(Dw(imax - 1, j) + alphaw(imax - 1, j) * Fw(imax - 1, j));
    }

    // ---- South boundary (j=1) / north boundary (j=jmax-1): both just
    // extend the ordinary interior an/as formula to the row the interior
    // loop above skipped -- no one-sided correction or S_p contribution
    // needed here (as(i,1) and an(i,jmax-1) stay 0: no south/north
    // neighbor exists at either row).
    for (int i = 1; i <= imax - 1; ++i) {
        an(i, 1) = -(Dn(i, 1) - (1.0 - alphan(i, 1)) * Fn(i, 1));
        as(i, jmax - 1) = -(Ds(i, jmax - 1) + alphas(i, jmax - 1) * Fs(i, jmax - 1));
    }

    // ---- ap assembly ----
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            double sumA = aw(i, j) + ae(i, j) + as(i, j) + an(i, j);
            ap(i, j) = (-sumA + A0_p(i, j) + dF(i, j) + S_p(i, j)) / alpha_q;
        }
    }

    // ---- Immersed boundary: solid cells. Pin ap=1, zero every
    // off-diagonal coefficient -- the RHS value (phi_inside_psi) is set
    // separately by rhsPhiADRE(), which is called every QUICK iteration
    // (this function only builds CM, never a RHS).
    for (size_t idx = 0; idx < ibmCoeffPhi.I_solid.size(); ++idx) {
        int is = ibmCoeffPhi.I_solid[idx];
        int js = ibmCoeffPhi.J_solid[idx];
        ae(is, js) = 0.0;
        aw(is, js) = 0.0;
        an(is, js) = 0.0;
        as(is, js) = 0.0;
        ap(is, js) = 1.0;
    }

    // ---- Immersed boundary: ghost cells. Same pinning; the -lambda_g_k
    // mirror-point coupling itself is added below, right before the CM
    // fill (updateCoeffGhost() in UpdateGhostReactionPhi.h is what
    // re-pushes a freshly-rederived version of this for the nonlinear-
    // reaction-BC case, every QUICK iteration -- this function is not
    // rerun mid-solve).
    for (size_t idx = 0; idx < ibmCoeffPhi.I_g.size(); ++idx) {
        int ig = ibmCoeffPhi.I_g[idx];
        int jg = ibmCoeffPhi.J_g[idx];
        ae(ig, jg) = 0.0;
        aw(ig, jg) = 0.0;
        an(ig, jg) = 0.0;
        as(ig, jg) = 0.0;
        ap(ig, jg) = 1.0;
    }

    // ---- Fill CM (caller-allocated -- see CoeffPhiADRE.h). No RHS here:
    // the transient/pressure-like source terms change every QUICK
    // iteration, so rhsPhiADRE() builds that Vec separately, reusing
    // this ap.
    // k(i,j) is a bijection over the interior range, so each loop
    // iteration owns exactly one row of the linear system: no two
    // iterations ever write the same (row,col) entry, even though a
    // given column (e.g. k(i,j) itself) legitimately appears in several
    // different rows -- once for each neighbor whose equation
    // references it.
    auto k = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };

    MatZeroEntries(CM);
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            int row = k(i, j);
            MatSetValue(CM, row, row, ap(i, j), ADD_VALUES);
            if (i > 1) MatSetValue(CM, row, k(i - 1, j), aw(i, j), ADD_VALUES);
            if (i < imax - 1) MatSetValue(CM, row, k(i + 1, j), ae(i, j), ADD_VALUES);
            if (j > 1) MatSetValue(CM, row, k(i, j - 1), as(i, j), ADD_VALUES);
            if (j < jmax - 1) MatSetValue(CM, row, k(i, j + 1), an(i, j), ADD_VALUES);
        }
    }
    

    // add ghostcell coefficients
    for (size_t idx = 0; idx < ibmCoeffPhi.I_g.size(); ++idx) {
        int row = k(ibmCoeffPhi.I_g[idx], ibmCoeffPhi.J_g[idx]);

        MatSetValue(CM, row, k(ibmCoeffPhi.I1[idx], ibmCoeffPhi.J1[idx]), -ibmCoeffPhi.lambda_g_1[idx], ADD_VALUES);
        MatSetValue(CM, row, k(ibmCoeffPhi.I2[idx], ibmCoeffPhi.J2[idx]), -ibmCoeffPhi.lambda_g_2[idx], ADD_VALUES);
        MatSetValue(CM, row, k(ibmCoeffPhi.I3[idx], ibmCoeffPhi.J3[idx]), -ibmCoeffPhi.lambda_g_3[idx], ADD_VALUES);
        MatSetValue(CM, row, k(ibmCoeffPhi.I4[idx], ibmCoeffPhi.J4[idx]), -ibmCoeffPhi.lambda_g_4[idx], ADD_VALUES);

        if (ibm.BQp == 1) {
            MatSetValue(CM, row, k(ibmCoeffPhi.I5[idx], ibmCoeffPhi.J5[idx]), -ibmCoeffPhi.lambda_g_5[idx], ADD_VALUES);
            MatSetValue(CM, row, k(ibmCoeffPhi.I6[idx], ibmCoeffPhi.J6[idx]), -ibmCoeffPhi.lambda_g_6[idx], ADD_VALUES);
        }
    }

    MatAssemblyBegin(CM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CM, MAT_FINAL_ASSEMBLY);

    TransportCoeffs c;
    c.ap = ap;
    return c;
}
