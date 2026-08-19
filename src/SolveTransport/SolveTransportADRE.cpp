#include "SolveTransportADRE.h"
#include "CoeffPhiADRE.h"
#include "RhsPhiADRE.h"
#include "CalculateQG.h"
#include "UpdateA1gPhi.h"
#include "UpdateGhostReactionPhi.h"
#include "../Utilities/ConvFlux.h"
#include "../Utilities/Timer.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace {
// DEBUG (temporary): dumps a Mat/Vec pair as a MATLAB-loadable .m script
// (spconvert-style ASCII, via PETSc's own MATLAB viewer format) --
// duplicated from SolveUVP.cpp's own dumpForMatlabComparison() (internal
// linkage there, not reusable across files) for the phi coefficient
// matrix/RHS, so compare_phi_matrices.m can diff CM_phi/RHS_phi against
// SolveTransportADRE.m's own Soln.CM_phi{i_s}/RHS_phi(:,i_s) the same
// way compare_matrices.m already does for U/V/P.
void dumpForMatlabComparison(Mat M, Vec b, const std::string &matName, const std::string &vecName,
                              const std::string &path) {
    PetscObjectSetName((PetscObject)M, matName.c_str());
    PetscObjectSetName((PetscObject)b, vecName.c_str());
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_SELF, path.c_str(), &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(M, viewer);
    VecView(b, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
}

// Unpacks one species' solved KSP vector into a full (imax+1, jmax+1)
// phi field + boundary values, then overwrites solid cells -- the
// phi-transport counterpart of FormPhi.h's own (still-unimplemented)
// formPhi(), kept local here rather than there since it sources the
// west/inlet Dirichlet value from bc.phi_a[i_s] (the caller's per-
// species inlet vector, see VariableNonDim.cpp's getBC()), masked here
// by the level set (ls.psi(0,j)>0) the same way rhsPhiADRE.cpp's own
// S_w(1,j) inlet term is.
Field2D formPhiADRE(const Domain &domain, Vec phiVec, const LS &ls, const Field2D &flag, double phiInside,
                     const std::vector<double> &phi_inlet) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;
    Field2D phi(imax + 1, jmax + 1);

    // ---- unpack the solved interior -- same k(i,j) bijection as
    // coeffPhiADRE()/rhsPhiADRE().
    auto k = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    const PetscScalar *arr;
    VecGetArrayRead(phiVec, &arr);
    for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax - 1; ++j) phi(i, j) = arr[k(i, j)];
    VecRestoreArrayRead(phiVec, &arr);

    // ---- boundary conditions (FORMPHI.m) ----
    for (int i = 1; i <= imax - 1; ++i) {
        phi(i, 0) = phi(i, 1);            // south: zero-gradient mirror
        phi(i, jmax) = phi(i, jmax - 1);  // north: zero-gradient mirror
    }
    for (int j = 0; j <= jmax; ++j) {
        phi(0, j) = phi_inlet[j] * (ls.psi(0, j) > 0.0);  // west: Dirichlet inlet
        phi(imax, j) = phi(imax - 1, j);                  // east: zero-gradient mirror
    }

    // ---- immersed boundary: solid cells (flag==2) ----
    for (int i = 0; i <= imax; ++i)
        for (int j = 0; j <= jmax; ++j)
            if (flag(i, j) == 2.0) phi(i, j) = phiInside;

    return phi;
}

// DEBUG (temporary): dumps every species' per-ghost-cell q_G/A1_g (this
// QUICK pass) and lambda_g_1..4 (geometry, frozen for the whole run) to
// JSON, keyed by (I_g,J_g) index order -- see SolveTransportADRE.h's own
// comment. Hand-rolled JSON (no json-c object graph) to match
// VariableNonDim.cpp's saveCurrentStateJson() style.
void dumpGhostDebug(const std::string &path, double time, const std::vector<int> &I_g,
                     const std::vector<int> &J_g, const std::vector<std::vector<double>> &q_G,
                     const std::vector<std::vector<double>> &A1_g_all,
                     const std::vector<IBMCoeff> &ibmCoeffPhi) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("dumpGhostDebug: failed to write " + path);

    auto darr = [](const std::vector<double> &v) {
        std::string s = "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s += ",";
            s += std::to_string(v[i]);
        }
        return s + "]";
    };
    auto iarr = [](const std::vector<int> &v) {
        std::string s = "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s += ",";
            s += std::to_string(v[i]);
        }
        return s + "]";
    };

    out << "{\n";
    out << "  \"time\": " << time << ",\n";
    out << "  \"I_g\": " << iarr(I_g) << ",\n";
    out << "  \"J_g\": " << iarr(J_g) << ",\n";
    out << "  \"species\": {\n";
    for (size_t i_s = 0; i_s < q_G.size(); ++i_s) {
        out << "    \"phi_" << i_s << "\": {\n";
        out << "      \"q_G\": " << darr(q_G[i_s]) << ",\n";
        out << "      \"A1_g\": " << darr(A1_g_all[i_s]) << ",\n";
        out << "      \"lambda_g_1\": " << darr(ibmCoeffPhi[i_s].lambda_g_1) << ",\n";
        out << "      \"lambda_g_2\": " << darr(ibmCoeffPhi[i_s].lambda_g_2) << ",\n";
        out << "      \"lambda_g_3\": " << darr(ibmCoeffPhi[i_s].lambda_g_3) << ",\n";
        out << "      \"lambda_g_4\": " << darr(ibmCoeffPhi[i_s].lambda_g_4) << "\n";
        out << "    }" << (i_s + 1 < q_G.size() ? "," : "") << "\n";
    }
    out << "  }\n";
    out << "}\n";

    if (!out) throw std::runtime_error("dumpGhostDebug: failed to write " + path);
}
}  // namespace

// Mirrors SolveTransport/SolveTransportADRE.m: one unsteady ADRE step
// for phi, all Np species, via the Hayase QUICK loop below. The only
// remaining placeholder is the nonlinear-reaction-BC path (step 1.2) --
// updateGhostibmLambda()/updateCoeffGhost() are still empty stubs, so
// controlVar.nonlinearReactionBC must stay false until those are filled
// in (see UpdateGhostReactionPhi.h).
void solveTransportADRE(ControlVar &controlVar, const Domain &domain,
                        const Variables &variables, StateVar &stateVar,
                        const IBM &ibm, std::vector<IBMCoeff> &ibmCoeffPhi,
                        const BC &bc, const LS &ls, const Flux &diffu_flux,
                        int subIdx, const std::string &outputAdreDir) {
    const int Np = variables.Np;

    // get convective flux coeffcient, uA -- once per call, not per QUICK iteration
    ConvFluxCoeffs convFluxPhi = computeConvFluxPhi(stateVar, controlVar, domain, variables, bc);

    // CM_phi/RHS_phi/x_phi: allocated ONCE per species here (same
    // discipline as SolveUVP.cpp's CM_u/RHS_u/x_u -- reused every QUICK
    // iteration below rather than malloc/free per pass for the same-
    // shaped memory). nz=11 per row for the same reason SolveUVP.cpp's
    // momentum matrices use it: a ghost row's standard 5-point band still
    // consumes 5 structural slots even though coeffPhiADRE()'s IBM
    // treatment zeroed their *values*, plus up to 6 more genuinely
    // distinct columns for the -lambda_g_k mirror-point corners.
    const int L_phi = (domain.imax - 1) * (domain.jmax - 1);
    std::vector<Mat> CM_phi(Np);
    std::vector<Vec> RHS_phi(Np), x_phi(Np);
    for (int i_s = 0; i_s < Np; ++i_s) {
        MatCreateSeqAIJ(PETSC_COMM_SELF, L_phi, L_phi, 11, nullptr, &CM_phi[i_s]);
        VecCreateSeq(PETSC_COMM_SELF, L_phi, &RHS_phi[i_s]);
        VecDuplicate(RHS_phi[i_s], &x_phi[i_s]);
    }

    // CM_phi: built ONCE per species here, not inside the QUICK loop --
    // beta_G is fixed under the current linear-reaction-matrix model
    // (ibmCoeffPhi[i_s]'s geometry/lambda_g_k already reflect that
    // species' real beta_G, from LSIBMcoeffsPhi() -- no need to re-derive
    // it here too), see step 1.2 below for the nonlinear-BC exception.
    // coeffPhiADRE() fills the whole band itself, including the ghost
    // rows' pinning + -lambda_g_k coupling; only the nonlinear-BC case
    // (step 1.2) needs to re-push a fresh lambda_g_k into CM afterward.
    std::vector<TransportCoeffs> sysPhi(Np);
    {
        // "assembly_transport" -- this one-time CM_phi build is C++'s
        // counterpart to SolveTransportADRE.m's own per-QUICK-iteration
        // COEFFPHIADRE loop (see the QUICK loop's own assembly_transport
        // scope below for why call frequency differs between languages,
        // even though the category name matches).
        ScopedTimer t("assembly_transport");
        for (int i_s = 0; i_s < Np; ++i_s) {
            sysPhi[i_s] = coeffPhiADRE(domain, diffu_flux, convFluxPhi, ibmCoeffPhi[i_s], ibm, variables, CM_phi[i_s]);
        }
    }

    // best-iterate checkpoint (mirrors SolveTransportADRE.m's phi_conv/
    // err_q tracking) -- in case iter_qq hits 250 without ever converging
    // below tol_q, return whichever pass had the smallest max iterate-to-
    // iterate change, not just whatever the last pass happened to produce.
    std::vector<Field2D> phi_conv = stateVar.phi;
    double bestErrSoFar = 100.0;

    // Reset every call (ControlVar persists across timesteps) -- without
    // this, err_q/iter_qq stay at whatever value the PREVIOUS call left
    // them at (err_q<=tol_q, since that's why that call's loop exited),
    // so the while condition below would be false immediately and the
    // QUICK loop would never run again after the first successful call.
    controlVar.err_q = 1.0;
    controlVar.iter_qq = 0;

    // ===========================  The Hayase's Quick Loop ====================
    // Hayse's loop, keeps 5 stencil structure by using lagged values for the second neighbor
    // so they end up on the RHS

    while ((controlVar.err_q > controlVar.tol_q && controlVar.iter_qq < 250) ) {
        // 0. get q_G from current phi_k (cross-species via reaction matrix A, not per-species)
        std::vector<std::vector<double>> q_G =
            calculateQG(variables, stateVar.phi, domain, ls, ibmCoeffPhi[0].I_g, ibmCoeffPhi[0].J_g);

        double maxDelta = 0.0;
        // A1_g per species, only actually populated/used when this pass's
        // ghost-cell values are about to be dumped (see below) --
        // collected here instead of inside dumpGhostDebug's call site
        // since A1_g itself is scoped to one species' iteration below.
        std::vector<std::vector<double>> A1_g_all(Np);
        // for each phi_i
        for (int i_s = 0; i_s < Np; i_s ++){
            std::vector<double> A1_g;
            {
                // "assembly_transport": mirrors SolveTransportADRE.m's
                // own t_timer span around its update_A1g/COEFFPHIADRE/
                // RHSPHIADRE block -- per-species here (vs. MATLAB's one
                // sample covering all Np species per iteration), so
                // "calls" won't match, but summed "total(s)" is directly
                // comparable.
                ScopedTimer t("assembly_transport");

                // 1. compute ghost cell rhs a1_g from q_G(:,i_s)
                A1_g = updateA1gPhi(ibm, ibmCoeffPhi[i_s], q_G[i_s]);
                A1_g_all[i_s] = A1_g;

                // 1.2 if nonlinear reaction BC: lambda_g_k depends on phi too,
                //     not just A1_g -- re-derive it and update ONLY the ghost
                //     rows' coupling entries in CM_phi[i_s] (the diffusive/convective band is still frozen).
                if (controlVar.nonlinearReactionBC) {
                    // TODO: betaSpecies would need to be re-derived from the
                    // current phi/geochemistry here, not read from
                    // ibm.beta_phi -- the whole point of this branch is that
                    // beta stops being fixed. Placeholder reuses the fixed,
                    // -diag(A)-derived value for now.
                    updateGhostibmLambda(ibmCoeffPhi[i_s], stateVar.phi);
                    updateCoeffGhost(domain, ibm, ls, ibm.beta_phi[i_s], ibmCoeffPhi[i_s], CM_phi[i_s]);
                }

                // 2. compute rhs
                rhsPhiADRE(stateVar, bc, diffu_flux, convFluxPhi, domain, ibmCoeffPhi[i_s], ibm, variables, ls,
                           sysPhi[i_s].ap, A1_g, i_s, RHS_phi[i_s]);
            }

            // DEBUG: dump this species' CM_phi/RHS_phi, first QUICK pass
            // only (controlVar.iter_qq is still 0 here) -- see
            // compare_phi_matrices.m.
            if (!outputAdreDir.empty() && controlVar.iter_qq == 0) {
                std::string base = outputAdreDir + "/cpp_phi_" + std::to_string(i_s) + "_" +
                                    std::to_string(controlVar.iTime) + "_" + std::to_string(subIdx) + ".m";
                dumpForMatlabComparison(CM_phi[i_s], RHS_phi[i_s], "CM_phi_" + std::to_string(i_s),
                                         "RHS_phi_" + std::to_string(i_s), base);
                printf("  [ADRE debug] saved %s\n", base.c_str());
            }

            Field2D phi_i_new;
            {
                // "solve_transport": mirrors SolveTransportADRE.m's own
                // t_timer span around its ilu/bicgstab/FORMPHI block.
                ScopedTimer t("solve_transport");

                // 3. solve each phi_i at k+1 iteration -- mirrors
                // SolveTransportADRE.m's ilu(CM_phi{i_s},'nofill','off') +
                // bicgstab(CM_phi{i_s},RHS_phi(:,i_s),tolbicg_c,maxit_c,L,U):
                // PCILU defaults to zero fill-in ('nofill') and no modified-
                // ILU diagonal compensation ('milu off'), so no extra PC
                // calls are needed beyond selecting "ilu" -- KSPSetType/
                // PCSetType pick the same bicgstab+ILU combination, and
                // KSPSetTolerances mirrors the (tol,maxit) pair MATLAB passes
                // positionally into bicgstab(). A fresh KSP per call (like
                // SolveUVP.cpp's momentum solves) since CM_phi[i_s]'s ILU
                // factors would otherwise go stale the moment the ghost rows'
                // lambda_g_k get re-pushed under the nonlinear-BC case.
                KSP ksp_phi;
                KSPCreate(PETSC_COMM_SELF, &ksp_phi);
                KSPSetOperators(ksp_phi, CM_phi[i_s], CM_phi[i_s]);
                // set up "bcgs" solver
                KSPSetType(ksp_phi, controlVar.ksp_type_scalar.c_str());
                KSPSetTolerances(ksp_phi, controlVar.tolbicg_c, PETSC_DEFAULT, PETSC_DEFAULT, controlVar.maxit_c);
                // ilu preconditioner
                PC pc_phi;
                KSPGetPC(ksp_phi, &pc_phi);
                PCSetType(pc_phi, controlVar.pc_type_scalar.c_str());
                // solve
                KSPSolve(ksp_phi, RHS_phi[i_s], x_phi[i_s]);
                KSPDestroy(&ksp_phi);

                // add boundary values i = 0, imax, j = 0, jmax, dirichlet and 0 gradient BC
                phi_i_new = formPhiADRE(domain, x_phi[i_s], ls, ibmCoeffPhi[i_s].flag, ibm.phi_inside_psi,
                                         bc.phi_a[i_s]);
            }

            // 4. delta_phi = phi_i - phi_i_k
            // TEMP (validation only, revert after use): SolveTransportADRE.m
            // line ~126 assigns ControlVar.err_q = max(max(abs(...))) INSIDE
            // this same per-species loop, overwriting rather than
            // max-accumulating -- so MATLAB's QUICK-loop exit test only ever
            // reflects the LAST species (i_s=Np-1)'s delta, silently ignoring
            // species 0..Np-2. Reproducing that overwrite (instead of the
            // correct std::max accumulation) here to test whether it explains
            // the persistent 1e-2/1e-3-level phi disagreement against MATLAB.
            double delta = 0.0;
            for (size_t k = 0; k < phi_i_new.data().size(); ++k) {
                delta = std::max(delta, std::fabs(phi_i_new.data()[k] - stateVar.phi[i_s].data()[k]));
            }
            maxDelta = delta;

            stateVar.phi[i_s] = phi_i_new;
        }

        // DEBUG: dump this call's ghost-cell q_G/A1_g/lambda_g_k, first
        // QUICK pass only (controlVar.iter_qq is still 0 here -- it's
        // incremented below) -- see SolveTransportADRE.h's comment.
        if (!outputAdreDir.empty() && controlVar.iter_qq == 0) {
            std::string path = outputAdreDir + "/dataRDE" + std::to_string(controlVar.iTime) + "_" +
                                std::to_string(subIdx) + "_ghost.json";
            dumpGhostDebug(path, controlVar.time, ibmCoeffPhi[0].I_g, ibmCoeffPhi[0].J_g, q_G, A1_g_all,
                           ibmCoeffPhi);
            printf("  [ADRE debug] saved %s\n", path.c_str());
        }

        controlVar.err_q = maxDelta;
        // 5. update phi_conv if maxDelta is the smallest so far
        if (maxDelta < bestErrSoFar) {
            bestErrSoFar = maxDelta;
            phi_conv = stateVar.phi;
        }
        controlVar.iter_qq++;
        printf("    [QUICK] iter_qq=%d err_q=%.6e (tol_q=%.3e)%s\n", controlVar.iter_qq, maxDelta,
               controlVar.tol_q, (maxDelta == bestErrSoFar) ? "  <- best" : "");
    }

    if (controlVar.err_q <= controlVar.tol_q) {
        printf("  [QUICK] converged after %d iterations, err_q=%.6e\n", controlVar.iter_qq, controlVar.err_q);
    } else {
        printf("  [QUICK] did NOT converge in %d iterations (err_q=%.6e > tol_q=%.3e) -- using best iterate, "
               "bestErr=%.6e\n",
               controlVar.iter_qq, controlVar.err_q, controlVar.tol_q, bestErrSoFar);
    }

    // Post-loop ghost+solid masking (SolveTransportADRE.m:116-118):
    // `phi(:,:,i_s) = (~IBM_coeffP(i_s).flag_p) .* phi(:,:,i_s)` --
    // zeroes every cell where flag_p~=0 (ghost=1 or solid=2), for every
    // species, on phi_conv (the best/final QUICK iterate) before it's
    // written back to stateVar.phi. Same pattern as SolveUVP.cpp's own
    // post-SIMPLE-loop fU/fV masking (SolveUVP.m:328-332).
    for (int i_s = 0; i_s < Np; ++i_s) {
        Field2D &ph = phi_conv[i_s];
        const Field2D &flag = ibmCoeffPhi[i_s].flag;
        for (int i = 0; i < ph.nx(); ++i)
            for (int j = 0; j < ph.ny(); ++j)
                if (flag(i, j) != 0.0) ph(i, j) = 0.0;
    }

    stateVar.phi = phi_conv;
    stateVar.phi_prev = stateVar.phi;

    for (int i_s = 0; i_s < Np; ++i_s) {
        MatDestroy(&CM_phi[i_s]);
        VecDestroy(&RHS_phi[i_s]);
        VecDestroy(&x_phi[i_s]);
    }
}
