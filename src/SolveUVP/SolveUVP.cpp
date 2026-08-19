#include "SolveUVP.h"
#include "COEFFU.h"
#include "COEFFV.h"
#include "COEFFP.h"
#include "RHSP.h"
#include "RHSP_PISO.h"
#include "../Utilities/ConvFlux.h"
#include "../Utilities/ConvergenceResiduals.h"
#include "../Utilities/Debug.h"
#include "../Utilities/Timer.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <json-c/json.h>

// ============================================================================

namespace {
// Sum of U_star over j=[1, jmax-1] (0-indexed) at fixed streamwise index
// `i` -- i.e. the total flow rate through the vertical cross-section at
// that x-location, used to compare inflow (i=0) against outflow
// (i=imax-1) for the mass-flux rebalancing below.
double sumCrossSectionFlow(const Field2D &U_star, int i) {
    double total = 0.0;
    for (int j = 1; j <= U_star.ny() - 2; ++j) total += U_star(i, j);
    return total;
}

// Solves CM * x = RHS with a fresh KSP, warm-started from `sol`'s own
// current interior values (so each SIMPLE iteration resumes from the
// last one's solution instead of from zero), then writes the result
// back into `sol` in place.
//
// CM/RHS/x are all caller-allocated and outlive this call -- solveUVP()
// allocates them once before its SIMPLE loop and reuses them every
// iteration (their sparsity pattern/size never changes, only their
// values do), rather than paying malloc/free every pass for the same-
// shaped memory. Only the KSP itself is still created/destroyed per
// call.
//
// k(i,j) = (i-iBegin) + (j-jBegin)*stride must be the exact same
// bijection coeffU()/coeffV() used to fill CM/RHS.
void solveMomentumSystem(Mat CM, Vec RHS, Vec x, Field2D &sol, int iBegin, int iEnd,
                          int jBegin, int jEnd, int stride, const ControlVar &controlVar,
                          const char *label = nullptr) {
    auto k = [iBegin, jBegin, stride](int i, int j) { return (i - iBegin) + (j - jBegin) * stride; };

    KSP ksp;
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, CM, CM);
    KSPSetType(ksp, controlVar.ksp_type_momentum.c_str());
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, controlVar.pc_type_momentum.c_str());
    KSPSetTolerances(ksp, controlVar.tolbicg, PETSC_DEFAULT, PETSC_DEFAULT, controlVar.maxit);

    // Tested reverting this to match MATLAB's always-x0=0 bicgstab() call
    // (suspecting the warm start was the source of the ~5e-4 U_star
    // mismatch near the grain) -- ruled out: the mismatch was essentially
    // unchanged, and cold-starting every SIMPLE iteration at this loose
    // rtol=1e-5 measurably destabilized the outer M_in/M_out feedback
    // loop (reintroduced the runaway-mass-imbalance blowup). Restored.
    for (int i = iBegin; i <= iEnd; ++i)
        for (int j = jBegin; j <= jEnd; ++j)
            VecSetValue(x, k(i, j), sol(i, j), INSERT_VALUES);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

    KSPSolve(ksp, RHS, x);

    if (debug::ksp && label) {
        KSPConvergedReason reason;
        PetscInt its;
        PetscReal rnorm;
        KSPGetConvergedReason(ksp, &reason);
        KSPGetIterationNumber(ksp, &its);
        KSPGetResidualNorm(ksp, &rnorm);
        PetscPrintf(PETSC_COMM_WORLD,
                    "      [%s KSP] reason=%d its=%d rnorm=%e\n", label, (int)reason, (int)its,
                    (double)rnorm);
    }

    const double *vals;
    VecGetArrayRead(x, &vals);
    for (int i = iBegin; i <= iEnd; ++i)
        for (int j = jBegin; j <= jEnd; ++j)
            sol(i, j) = vals[k(i, j)];
    VecRestoreArrayRead(x, &vals);

    KSPDestroy(&ksp);
}

// Solves CM_p * x_p = RHS_p and unpacks the result into PCOR's non-
// pinned interior cells. Shared by TODO 3's first-corrector solve and
// PISO's second-corrector solve below -- both need the identical KSP
// setup and the identical reduced-index unpacking, just with a
// different RHS. Leaves PCOR's pinned cell and boundary edges alone --
// formPCor() fills those separately, right after this call. No warm
// start here (unlike solveMomentumSystem()) -- not asked for, and
// MATLAB's own pcg() call doesn't pass one either.
void solvePressureSystem(Mat CM_p, Vec RHS_p, Vec x_p, Field2D &PCOR,
                          const Domain &domain, const ControlVar &controlVar) {
    KSP ksp;
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, CM_p, CM_p);
    KSPSetType(ksp, controlVar.ksp_type_pressure.c_str());
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, controlVar.pc_type_pressure.c_str());
    KSPSetTolerances(ksp, controlVar.tol_q, PETSC_DEFAULT, PETSC_DEFAULT, controlVar.maxit_c);

    KSPSolve(ksp, RHS_p, x_p);

    const int imax = domain.imax;
    const int jmax = domain.jmax;
    PressurePin pin = pinnedPressureCell(imax, jmax, controlVar.imposePresBC);
    auto kFull = [imax](int i, int j) { return (i - 1) + (j - 1) * (imax - 1); };
    const int kPin = kFull(pin.i0, pin.j0);
    auto kReduced = [kPin](int kf) { return kf < kPin ? kf : kf - 1; };

    const double *vals;
    VecGetArrayRead(x_p, &vals);
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            if (i == pin.i0 && j == pin.j0) continue;
            PCOR(i, j) = vals[kReduced(kFull(i, j))];
        }
    }
    VecRestoreArrayRead(x_p, &vals);

    KSPDestroy(&ksp);
}

// DEBUG (temporary): dump a Mat/Vec pair to a MATLAB-loadable .m script --
// `run(path)` in MATLAB reconstructs `matName` (via spconvert) and `vecName`
// as plain variables, for direct numeric comparison against SolveUVP.m's
// own Soln.CM_u/RHS_U etc. at the same (first) SIMPLE iteration.
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

// DEBUG (temporary): dump U_star/V_star -- the just-solved,
// pre-pressure-correction velocity fields rhsP() consumes -- to a JSON
// file, same flat row-major convention as VariableNonDim::
// saveCurrentStateJson(), so debug_compare/compare_state.m can reshape
// and diff them against SolveUVP.m's own pre-RHSP snapshot. Isolates
// whether a downstream RHS_P discrepancy comes from U_star/V_star
// themselves (upstream, already-solved momentum fields) or from rhsP()'s
// own formula.
void dumpPrePressureState(const Field2D &U_star, const Field2D &V_star, const std::string &path) {
    json_object *root = json_object_new_object();
    auto addField = [&](const char *name, const Field2D &f) {
        json_object *arr = json_object_new_array();
        for (double v : f.data()) json_object_array_add(arr, json_object_new_double(v));
        json_object_object_add(root, name, arr);
    };
    addField("U_star", U_star);
    addField("V_star", V_star);
    FILE *fp = fopen(path.c_str(), "w");
    if (fp) {
        fputs(json_object_to_json_string_ext(root, JSON_C_TO_STRING_SPACED), fp);
        fclose(fp);
    }
    json_object_put(root);
}

// DEBUG (temporary): appends one SIMPLE iteration's worth of traced-cell
// U values to a growing JSON array, rewriting the whole file each call
// -- tracks whether the ghost-cell U discrepancy compare_output.m found
// at the final, post-SIMPLE-loop state (MATLAB exactly 0, C++ ~-0.5 at
// these cells) is already present at iteration 1, or compounds over the
// SIMPLE loop's ~24 passes for this timestep. Only ever called while
// controlVar.iTime==1 (see call sites below), so the static array below
// is safe to accumulate across this function's whole lifetime -- it's
// only ever fed by that one solveUVP() call.
// wanted[][2] here are C++ 0-based, matching main.cpp's ghost-cell
// diagnostic; MATLAB-native equivalent is (i+1,j+1): (194,63) (208,63)
// (194,140) (208,140) (190,64).
void dumpIterTraceEntry(int ii, bool piso, const Field2D &U, const std::string &path) {
    static const int wanted[5][2] = {{193, 62}, {207, 62}, {193, 139}, {207, 139}, {189, 63}};
    static json_object *arr = json_object_new_array();

    printf("  [iterTrace] ii=%d piso=%d U(193,62)=%.8e U(207,62)=%.8e U(193,139)=%.8e "
           "U(207,139)=%.8e U(189,63)=%.8e\n",
           ii, (int)piso, U(193, 62), U(207, 62), U(193, 139), U(207, 139), U(189, 63));

    json_object *entry = json_object_new_object();
    json_object_object_add(entry, "ii", json_object_new_int(ii));
    json_object_object_add(entry, "piso", json_object_new_boolean(piso));
    json_object *vals = json_object_new_array();
    for (auto &w : wanted) json_object_array_add(vals, json_object_new_double(U(w[0], w[1])));
    json_object_object_add(entry, "U", vals);
    json_object_array_add(arr, entry);

    FILE *fp = fopen(path.c_str(), "w");
    if (fp) {
        fputs(json_object_to_json_string_ext(arr, JSON_C_TO_STRING_SPACED), fp);
        fclose(fp);
    }
}
}  // namespace

// ============================================================================

void formUV(Field2D &U_star, Field2D &V_star, const BC &bc, int iter) {
    const int imax = U_star.nx();
    const int jmax = V_star.ny();

    // Type codes: 1 = Dirichlet (prescribed value), 3 = Neumann
    // (zero-gradient).

    // ---- U ----

    // MATLAB's FORMUV.m rebuilds U as a fresh zeros(imax,jmax+1) array on
    // every call, so any cell the reshape/boundary lines below never touch
    // implicitly stays zero. U_star here is mutated in place instead, so
    // the one cell that MATLAB's coverage genuinely misses -- (imax-1,
    // jmax-1), the outlet's last row just below the north wall, excluded
    // by both the East loop (stops at jmax-2) and the corner copy below
    // (which only writes (imax-1,jmax), not (imax-1,jmax-1)) -- would
    // otherwise keep whatever stale/initial value it already held. Zero it
    // explicitly so the corner copy below (which reads this cell) matches
    // MATLAB's FORMUV.m exactly.
    U_star(imax - 1, jmax - 1) = 0.0;

    // West/inlet: set to the prescribed inflow value at every row.
    for (int j = 1; j <= jmax - 1; ++j) {
        U_star(0, j) = bc.U_a[j];
    }

    // East/outlet.
    if (bc.BC_e_p != 1) {
        // Copy the last interior column across (zero streamwise gradient).
        for (int j = 1; j <= jmax - 2; ++j) {
            U_star(imax - 1, j) = U_star(imax - 2, j);
        }
        // From the second pass onward, rescale the whole outlet row so its
        // total flow rate matches the inlet's total flow rate exactly.
        if (iter > 1) {
            double M_in = sumCrossSectionFlow(U_star, 0);
            double M_out = sumCrossSectionFlow(U_star, imax - 1);
            PetscPrintf(PETSC_COMM_WORLD, "      [formUV] iter=%d M_in=%e M_out=%e ratio=%e\n", iter,
                        M_in, M_out, M_in / M_out);
            for (int j = 1; j <= jmax - 2; ++j) {
                U_star(imax - 1, j) = (M_in / M_out) * U_star(imax - 2, j);
            }
        }
    } else if (iter > 1) {
        // Outlet row is left as whatever it already holds; only rescale
        // it (in place) to match the inlet's total flow rate.
        double M_in = sumCrossSectionFlow(U_star, 0);
        double M_out = sumCrossSectionFlow(U_star, imax - 1);
        for (int j = 1; j <= jmax - 2; ++j) {
            U_star(imax - 1, j) = (M_in / M_out) * U_star(imax - 1, j);
        }
    }

    // North/south walls: either pin the wall value exactly (antisymmetric
    // mirror through the ghost point) or copy the adjacent interior value
    // (zero-gradient), depending on the BC type at each wall.
    for (int i = 0; i < imax; ++i) {
        if (bc.BC_n_u == 1) {
            U_star(i, jmax) = 2 * bc.U_d[i] - U_star(i, jmax - 1);
        } else if (bc.BC_n_u == 3) {
            U_star(i, jmax) = U_star(i, jmax - 1);
        }
        if (bc.BC_s_u == 1) {
            U_star(i, 0) = 2 * bc.U_c[i] - U_star(i, 1);
        } else if (bc.BC_s_u == 3) {
            U_star(i, 0) = U_star(i, 1);
        }
    }

    // Corner points, set last so they override whatever the wall formulas
    // above just wrote at these two locations.
    U_star(imax - 1, 0) = U_star(imax - 1, 1);
    U_star(imax - 1, jmax) = U_star(imax - 1, jmax - 1);

    // ---- V ----

    // Same reasoning as U's corner above: MATLAB's V = zeros(imax+1,jmax)
    // reset means these four corners -- never covered by the N/S wall loop
    // (i=1..imax-1) or the W/E edge loop (j=1..jmax-2) in either language
    // -- are always exactly zero in MATLAB. Match that explicitly.
    V_star(0, 0) = 0.0;
    V_star(0, jmax - 1) = 0.0;
    V_star(imax, 0) = 0.0;
    V_star(imax, jmax - 1) = 0.0;

    // South/north walls: either the prescribed wall value directly, or a
    // zero-gradient copy, depending on the BC type at each wall.
    for (int i = 1; i <= imax - 1; ++i) {
        if (bc.BC_n_v == 1) {
            V_star(i, jmax - 1) = bc.V_d[i];
        } else if (bc.BC_n_v == 3) {
            V_star(i, jmax - 1) = V_star(i, jmax - 2);
        }
        if (bc.BC_s_v == 1) {
            V_star(i, 0) = bc.V_c[i];
        } else if (bc.BC_s_v == 3) {
            V_star(i, 0) = V_star(i, 1);
        }
    }

    // West/east edges: always an antisymmetric mirror (pinning V=0 at the
    // inlet) and a zero-gradient copy (outlet) -- no BC-type check on
    // either edge. Corner points are left untouched.
    for (int j = 1; j <= jmax - 2; ++j) {
        V_star(0, j) = -V_star(1, j);
        V_star(imax, j) = V_star(imax - 1, j);
    }
}

// ============================================================================

PressurePin pinnedPressureCell(int imax, int jmax, const std::string &imposePresBC) {
    int j0;
    if (imposePresBC == "middle") {
        int J = static_cast<int>(std::ceil((jmax - 1) / 2.0));  // 1-indexed row
        j0 = J - 1;
    } else if (imposePresBC == "top") {
        j0 = jmax - 1;
    } else {  // "bottom"
        j0 = 1;
    }
    return {imax - 1, j0};
}

void formPCor(Field2D &PCOR, const Domain &domain, const ControlVar &controlVar) {
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    // Pin the reference point to zero (see pinnedPressureCell()'s comment
    // in SolveUVP.h for how it's chosen).
    PressurePin pin = pinnedPressureCell(imax, jmax, controlVar.imposePresBC);
    PCOR(pin.i0, pin.j0) = 0.0;

    // Boundary conditions: a zero-gradient copy on all four edges,
    // unconditionally (no BC-type check on any edge). Corner points are
    // left untouched.
    for (int i = 1; i <= imax - 1; ++i) {
        PCOR(i, 0) = PCOR(i, 1);
        PCOR(i, jmax) = PCOR(i, jmax - 1);
    }
    for (int j = 1; j <= jmax - 1; ++j) {
        PCOR(0, j) = PCOR(1, j);
        PCOR(imax, j) = PCOR(imax - 1, j);
    }
}

// ============================================================================
// overwrite U_star, V_star and P_star in place
void newUVP(Field2D &U_star, Field2D &V_star, Field2D &P_star,
            const Field2D &PCOR, const Field2D &d_u, const Field2D &d_v,
            double alpha_p, const BC &bc, int iter) {
    const int imax = U_star.nx();
    const int jmax = V_star.ny();

    // ---- U: interior correction, then boundary (NEWUVP.m:63-153) ----
    // U(i,j) = U_star(i,j) + d_u(i,j) * (PCOR(i,j) - PCOR(i+1,j)) -- nudges
    // U toward divergence-free using the pressure-correction gradient
    // across the cell it sits between.
    for (int i = 1; i <= imax - 2; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            U_star(i, j) = U_star(i, j) + d_u(i, j) * (PCOR(i, j) - PCOR(i + 1, j));
        }
    }

    // West/inlet: the full column range (0..jmax, including both
    // corners) -- wider than formUV()'s own west treatment (1..jmax-1),
    // which leaves the corners to the north/south wall formulas instead.
    // (NEWUVP.m also applies an IBM mask here, U_a.*(flag_u~=2) -- not
    // implemented, defaults to always-true.)
    for (int j = 0; j <= jmax; ++j) {
        U_star(0, j) = bc.U_a[j];
    }

    // East/outlet: if not P boundary not Dirichlet, manufacture a value on the right so gradient = 
    //              if dirchilet, rescale flux so flux in = flux out
    if (bc.BC_e_p != 1) {
        for (int j = 1; j <= jmax - 2; ++j) {
            U_star(imax - 1, j) = U_star(imax - 2, j);
        }
        if (iter > 1) {
            double M_in = sumCrossSectionFlow(U_star, 0);
            double M_out = sumCrossSectionFlow(U_star, imax - 1);
            for (int j = 1; j <= jmax - 2; ++j) {
                U_star(imax - 1, j) = (M_in / M_out) * U_star(imax - 2, j);
            }
        }
    } else if (iter > 1) {
        double M_in = sumCrossSectionFlow(U_star, 0);
        double M_out = sumCrossSectionFlow(U_star, imax - 1);
        for (int j = 1; j <= jmax - 2; ++j) {
            U_star(imax - 1, j) = (M_in / M_out) * U_star(imax - 1, j);
        }
    }

    // North & South: no slip conditon, manufactor U, so interplolation between boundary pacth and intertior patch give u = 0 at the wall
    if (bc.BC_n_u == 1) {
        for (int i = 0; i < imax; ++i) U_star(i, jmax) = 2 * bc.U_d[i] - U_star(i, jmax - 1);
    } else if (bc.BC_n_u == 3) {
        for (int i = 1; i <= imax - 2; ++i) U_star(i, jmax) = U_star(i, jmax - 1);
    }

    // South: full row range for both branches
    if (bc.BC_s_u == 1) {
        for (int i = 0; i < imax; ++i) U_star(i, 0) = 2 * bc.U_c[i] - U_star(i, 1);
    } else if (bc.BC_s_u == 3) {
        for (int i = 0; i < imax; ++i) U_star(i, 0) = U_star(i, 1);
    }

    // Corner points, set last -- same as formUV()'s own corner formula.
    U_star(imax - 1, 0) = U_star(imax - 1, 1);
    U_star(imax - 1, jmax) = U_star(imax - 1, jmax - 1);

    // ---- V: interior correction, then boundary (NEWUVP.m:155-195) ----
    // Same correction as U, using the north-south PCOR gradient instead
    // of east-west, matching V's own orientation.
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 2; ++j) {
            V_star(i, j) = V_star(i, j) + d_v(i, j) * (PCOR(i, j) - PCOR(i, j + 1));
        }
    }

    // North/south walls and west/east edges: 1) dirchilet BC condition 3) symmetry plane or free-slip wall
    for (int i = 1; i <= imax - 1; ++i) {
        if (bc.BC_n_v == 1) {
            V_star(i, jmax - 1) = bc.V_d[i];
        } else if (bc.BC_n_v == 3) {
            V_star(i, jmax - 1) = V_star(i, jmax - 2);
        }
        if (bc.BC_s_v == 1) {
            V_star(i, 0) = bc.V_c[i];
        } else if (bc.BC_s_v == 3) {
            V_star(i, 0) = V_star(i, 1);
        }
    }

    // west/inlet boundary of V: the assumed inflow is simple and straight, only u component so v = 0
    // east/outlet boudnary: nothing changes velocity so v = v on the left
    for (int j = 1; j <= jmax - 2; ++j) {
        V_star(0, j) = -V_star(1, j);
        V_star(imax, j) = V_star(imax - 1, j);
    }

    // ---- P: interior update, then boundary (NEWUVP.m:202-213) ----
    // P(i,j) = P_star(i,j) + alpha_p * PCOR(i,j) -- interior only
    // (i=1..imax-1, j=1..jmax-1, matching P's own established range,
    // same as coeffP()'s -- NEWUVP.m's own loop is `for j=2:jmax, for
    // i=2:imax`, not the full array).
    for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
            P_star(i, j) = P_star(i, j) + alpha_p * PCOR(i, j);
        }
    }

    // Boundary: a zero-gradient copy of the just-computed interior on
    // south/north/west, but a sign flip (not a copy) on east. Not the
    // same formula as the interior update, and not derivable from it.
    for (int i = 1; i <= imax - 1; ++i) {
        P_star(i, 0) = P_star(i, 1);            // south, dP/dy=0
        P_star(i, jmax) = P_star(i, jmax - 1);  // north, dP/dy=0
    }
    for (int j = 1; j <= jmax - 1; ++j) {
        P_star(0, j) = P_star(1, j);             // west, dP/dx=0
        P_star(imax, j) = -P_star(imax - 1, j);  // east, p=0
    }
}

// ============================================================================

void newUVPPiso(Field2D &U_star, Field2D &V_star, Field2D &P_star, const Field2D &PCOR,
                const MomentumCoeffs &sysU, const MomentumCoeffs &sysV,
                const Field2D &deltaU, const Field2D &deltaV, double alpha_p,
                const BC &bc, int iter) {
    
                    // Same PCOR-gradient correction for U/V, boundary conditions for
    // U/V/P, and alpha_p*PCOR update for P that newUVP() always applies
    // -- PISO doesn't change any of this part.
    newUVP(U_star, V_star, P_star, PCOR, sysU.d_SIMPLC, sysV.d_SIMPLC, alpha_p, bc, iter);

    // PISO's own extra term: 
    // account the neighbor-coupling contribution  (from this iteration's first-correction velocity change) divided by the raw diagonal ap 
    const int imaxU = U_star.nx();
    const int jmaxU = U_star.ny();
    for (int i = 1; i <= imaxU - 2; ++i) {
        for (int j = 1; j <= jmaxU - 1; ++j) {
            double neighborSum = sysU.aw(i, j) * deltaU(i - 1, j) + sysU.ae(i, j) * deltaU(i + 1, j) +
                                  sysU.as(i, j) * deltaU(i, j - 1) + sysU.an(i, j) * deltaU(i, j + 1);
            U_star(i, j) -= neighborSum / sysU.ap(i, j);
        }
    }
    // same PISO extra term for V
    const int imaxV = V_star.nx();
    const int jmaxV = V_star.ny();
    for (int i = 1; i <= imaxV - 1; ++i) {
        for (int j = 1; j <= jmaxV - 2; ++j) {
            double neighborSum = sysV.aw(i, j) * deltaV(i - 1, j) + sysV.ae(i, j) * deltaV(i + 1, j) +
                                  sysV.as(i, j) * deltaV(i, j - 1) + sysV.an(i, j) * deltaV(i, j + 1);
            V_star(i, j) -= neighborSum / sysV.ap(i, j);
        }
    }
}

// ============================================================================

void solveUVP(ControlVar &controlVar, const Domain &domain,
              const Variables &variables, StateVar &stateVar, const IBM &ibm,
              const IBMCoeff &ibmCoeffU, const IBMCoeff &ibmCoeffV, const BC &bc,
              const Flux &diffu_flux) {
    // U_star/V_star: the momentum-only (pre-pressure-correction) velocity
    // guess for the current iteration -- updated by solveMomentumSystem(),
    // then corrected in place by newUVP()/newUVPPiso().
    //
    // U_star/V_star: current iteration value -> potentially be the solution
    // U_star_old/V_star_old: value from previous iteration
    // U_star_uncorrected/V_star_uncorrected: this iteration's OWN pre-correction guess
    //
    // U_star/V_star/U_star_old/V_star_old all start at the current stateVar.U/V 


    Field2D U_star = stateVar.U;
    Field2D V_star = stateVar.V;
    Field2D U_star_old = U_star;
    Field2D V_star_old = V_star;
    Field2D PCOR(domain.imax + 1, domain.jmax + 1);

    // U/V-momentum CM/RHS/x: allocated once here, update value each iteration
    // Fe/Fw/Fn/Fs (hence ae/aw/an/as/ap/S) are relinearized around the current velocity field every pass
    const int Lu = (domain.imax - 2) * (domain.jmax - 1);
    const int Lv = (domain.imax - 1) * (domain.jmax - 2);
    Mat CM_u, CM_v;
    Vec RHS_u, x_u, RHS_v, x_v;
    // nz=11, not 5: a plain interior row only ever needs 5 (diagonal +
    // the 4 orthogonal neighbors), but a ghost row's own standard-band
    // MatSetValue calls still consume 5 structural slots even though
    // coeffU()/coeffV()'s IBM treatment zeroed their *values* (PETSc
    // counts (row,col) pairs that were ever written, not which ones are
    // nonzero, unless MAT_IGNORE_ZERO_ENTRIES is set) -- then the extra
    // -lambda_g_k mirror-point couplings add up to 6 more, genuinely
    // distinct columns. Worst case is 5+6=11, not 5 or even 7. With too
    // small an nz, PETSc's default MAT_NEW_NONZERO_ALLOCATION_ERR
    // rejects the extra ghost-row entries outright ("Argument out of
    // range: New nonzero ... caused a malloc") instead of silently
    // growing.
    MatCreateSeqAIJ(PETSC_COMM_SELF, Lu, Lu, 11, nullptr, &CM_u);
    VecCreateSeq(PETSC_COMM_SELF, Lu, &RHS_u);
    VecDuplicate(RHS_u, &x_u);
    MatCreateSeqAIJ(PETSC_COMM_SELF, Lv, Lv, 11, nullptr, &CM_v);
    VecCreateSeq(PETSC_COMM_SELF, Lv, &RHS_v);
    VecDuplicate(RHS_v, &x_v);

    // Pressure-correction CM_p2: allocate once
    // its entries depend on sysU.d_SIMPLC/sysV.d_SIMPLC, which change every SIMPLE iteration
    const int Lp = (domain.imax - 1) * (domain.jmax - 1) - 1;
    Mat CM_p;
    Vec RHS_p, x_p;
    MatCreateSeqAIJ(PETSC_COMM_SELF, Lp, Lp, 5, nullptr, &CM_p);
    VecCreateSeq(PETSC_COMM_SELF, Lp, &RHS_p);
    VecDuplicate(RHS_p, &x_p);

    // ii<2 forces at least 2 passes through the loop before the residual
    // is even checked. The ii>100 cap below stops it from running forever.
    controlVar.ii = 0;
    controlVar.resi = 1e10;

    while ((controlVar.resi > controlVar.tol || controlVar.ii < 2) && controlVar.f == 0) {
        controlVar.ii += 1;
        

        // Timing instrumentation below mirrors SolveUVP.m's own tic/toc
        // blocks one-for-one, same category names (assembly_u/_v/_p,
        // solve_u/_v/_p, solve_p_piso, reconstruction_uvp), so the C++
        // and MATLAB TIMING SUMMARY tables line up row by row. See
        // Utilities/Timer.h. Each ScopedTimer records on scope exit, so
        // the declarations that outlive a timed block (sysU/sysV, needed
        // later by coeffP and the PISO branch) are hoisted above it.

        // ---- U-momentum ----
        MomentumCoeffs sysU;
        {
            // ConvFlux + COEFFU, matching SolveUVP.m's assembly_u block.
            ScopedTimer t("assembly_u");
            ConvFluxCoeffs convFluxU = computeConvFluxU(stateVar, controlVar, domain, variables, bc);
            sysU = coeffU(stateVar, U_star_old, domain, diffu_flux.Diffu_U, convFluxU,
                           ibm, ibmCoeffU, variables, bc, controlVar.disc_scheme_vel,
                           CM_u, RHS_u);
        }
        {
            // x_u is the Vec object solvedm U_star is updated with x_u values
            // Covers KSP setup + solve, i.e. MATLAB's ilu()+bicgstab() pair.
            // Also covers the Vec pack/unpack around the solve, which has no
            // MATLAB counterpart (there the solution vector is reshaped later,
            // inside FORMUV) -- negligible next to the solve itself.
            ScopedTimer t("solve_u");
            solveMomentumSystem(CM_u, RHS_u, x_u, U_star, 1, domain.imax - 2, 1, domain.jmax - 1,
                                 domain.imax - 2, controlVar, "U");
        }

        // DEBUG (temporary, inline): traced-cell U_star values right out
        // of the KSP solve, before formUV()/the pressure correction touch
        // anything -- same 5 cells as dumpIterTraceEntry's post-correction
        // print below, so the two together show whether the ghost-cell
        // discrepancy is already present pre-formUV, or only appears
        // after formUV/newUVP run.
        if (controlVar.iTime == 1) {
            printf("  [iterTrace postsolve] ii=%d U(193,62)=%.8e U(207,62)=%.8e U(193,139)=%.8e "
                   "U(207,139)=%.8e U(189,63)=%.8e\n",
                   controlVar.ii, U_star(193, 62), U_star(207, 62), U_star(193, 139),
                   U_star(207, 139), U_star(189, 63));
        }

        // ---- V-momentum ----
        MomentumCoeffs sysV;
        {
            ScopedTimer t("assembly_v");
            ConvFluxCoeffs convFluxV = computeConvFluxV(stateVar, controlVar, domain, variables, bc);
            sysV = coeffV(stateVar, V_star_old, domain, diffu_flux.Diffu_V, convFluxV,
                           ibm, ibmCoeffV, variables, bc, controlVar.disc_scheme_vel,
                           CM_v, RHS_v);
        }
        {
            ScopedTimer t("solve_v");
            solveMomentumSystem(CM_v, RHS_v, x_v, V_star, 1, domain.imax - 1, 1, domain.jmax - 2,
                                 domain.imax - 1, controlVar, "V");
        }

        // DEBUG (temporary): dump U_star/V_star right after both KSP
        // solves, before formUV() touches any boundary values -- isolates
        // whether the interior discrepancy already exists right out of the
        // linear solve, or gets introduced/amplified by formUV().
        if (controlVar.ii == 1) {
            dumpPrePressureState(
                U_star, V_star,
                "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/cpp_postsolve_state.json");
        }

        {
            // FORMUV + RHSP -- SolveUVP.m books both under
            // reconstruction_uvp (its first of two samples per iteration).
            ScopedTimer t("reconstruction_uvp");

            // Fill in U_star/V_star's boundary conditions (their interior
            // values were just solved for above).
            formUV(U_star, V_star, bc, controlVar.ii);

            // DEBUG (temporary, inline): same 5 traced cells, right after
            // formUV() -- isolates whether formUV's boundary/outlet-
            // rescale logic changes them (it shouldn't; none of these 5
            // cells are on a formUV-touched row/column).
            if (controlVar.iTime == 1) {
                printf("  [iterTrace preRHSP] ii=%d U(193,62)=%.8e U(207,62)=%.8e U(193,139)=%.8e "
                       "U(207,139)=%.8e U(189,63)=%.8e\n",
                       controlVar.ii, U_star(193, 62), U_star(207, 62), U_star(193, 139),
                       U_star(207, 139), U_star(189, 63));
            }

            // DEBUG (temporary): dump U_star/V_star right before rhsP()
            // consumes them -- see debug_compare/compare_state.m.
            if (controlVar.ii == 1) {
                dumpPrePressureState(
                    U_star, V_star,
                    "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/cpp_prepressure_state.json");
            }

            // ---- Pressure correction ----

            // fill RHS_p fills allocated-once/reset-every-iteration.
            rhsP(U_star, V_star, domain, controlVar, ibmCoeffU, ibmCoeffV, RHS_p);
        }

        {
            //  CM_p fills in place, allocated-once/reuse-every-iteration pattern as CM_u/CM_v.
            //  Reused unmodified by PISO's second corrector below.
            ScopedTimer t("assembly_p");
            coeffP(sysU.d_SIMPLC, sysV.d_SIMPLC, domain, controlVar, CM_p);
        }

        // DEBUG (temporary): dump the first SIMPLE iteration's CM/RHS for
        // U/V/P so they can be loaded into MATLAB and compared numerically
        // against SolveUVP.m's own Soln.CM_u/RHS_U, Soln.CM_v/RHS_V,
        // Soln.CM_p2/RHS_P2 at the same iteration.
        if (controlVar.ii == 1) {
            const std::string dir = "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/";
            dumpForMatlabComparison(CM_u, RHS_u, "CM_u", "RHS_U", dir + "cpp_u.m");
            dumpForMatlabComparison(CM_v, RHS_v, "CM_v", "RHS_V", dir + "cpp_v.m");
            dumpForMatlabComparison(CM_p, RHS_p, "CM_p2", "RHS_P2", dir + "cpp_p.m");
        }

        {
            //  solves CM_p*x_p=RHS_p and unpacks the result into PCOR's non-pinned interior cells.
            ScopedTimer t("solve_p");
            solvePressureSystem(CM_p, RHS_p, x_p, PCOR, domain, controlVar);
        }

        // ---- Velcoity correction ------
        // Snapshot the pre-correction guess for PISO's own delta term below
        // (its own quantity, distinct from U_star_old/V_star_old: those now
        // hold the *post*-correction value once newUVP() runs, for
        // coeffU()/coeffV()'s under-relaxation term next SIMPLE iteration --
        // PISO's delta needs this iteration's correction step in isolation,
        // i.e. after minus before, so it needs its own separate "before"
        // copy). Hoisted above the timed block below like sysU/sysV, since
        // PISO (further down, outside that block) is what reads it.
        Field2D U_star_uncorrected = U_star;
        Field2D V_star_uncorrected = V_star;
        {
            // FORMPCOR + NEWUVP -- SolveUVP.m's second reconstruction_uvp
            // sample for this iteration.
            ScopedTimer t("reconstruction_uvp");

            // Fills PCOR's pinned cell (=0) and boundary edges -- the part
            // solvePressureSystem() deliberately leaves alone.
            formPCor(PCOR, domain, controlVar);

            // Correct U_star/V_star into stateVar.U/V using PCOR's gradient,
            // and update U_star, V_star and stateVar.P
            newUVP(U_star, V_star, stateVar.P, PCOR, sysU.d_SIMPLC, sysV.d_SIMPLC, variables.alpha_p, bc, controlVar.ii);
            U_star_old = U_star;
            V_star_old = V_star;
            stateVar.U = U_star;
            stateVar.V = V_star;

            if (controlVar.iTime == 1) {
                dumpIterTraceEntry(
                    controlVar.ii, false, stateVar.U,
                    "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/cpp_iter_trace.json");
            }
        }

        // ---- Convergence check ----
        // compute residuals 1) CM_u * U - RHS_U, 2) CM_v * V - RHS_V,  3) divergnece check, resuse RHS_p from pressure equation 
        // the divergence check could be done by recomputed RHS_p on the corrected U and V. Matlab does not do that, current version matches the orignal matlab
        {
            ConvergenceResult conv =
                convergenceResiduals(stateVar.U, stateVar.V, RHS_p, controlVar.ii, PCOR, CM_u, CM_v,
                                      RHS_u, RHS_v, domain, controlVar.disc_scheme_vel, controlVar.tol,
                                      /*PISO=*/0);
            controlVar.resi = conv.resi_max;
            controlVar.messageFlow = conv.message;
            PetscPrintf(PETSC_COMM_WORLD,
                        "    SIMPLE ii=%d resi=%e (U=%e, V=%e, mass=%e)\n", controlVar.ii,
                        controlVar.resi, conv.resi1, conv.resi2, conv.resi3);
        }

        // ----- PISO: second pressure-correction pass ---------
        // PISO: reusing the momentum coefficients from above without rebuilding them, 
        // PISO mode always does exactly 2 corrector passes per call, never iterates SIMPLE further.
        if (controlVar.PISO == 1) {
            // deltaU/deltaV: the velocity change from this iteration's first pressure correction alone
            Field2D deltaU(stateVar.U.nx(), stateVar.U.ny());
            for (int i = 0; i < deltaU.nx(); ++i)
                for (int j = 0; j < deltaU.ny(); ++j) deltaU(i, j) = U_star(i, j) - U_star_uncorrected(i, j);
            Field2D deltaV(stateVar.V.nx(), stateVar.V.ny());
            for (int i = 0; i < deltaV.nx(); ++i)
                for (int j = 0; j < deltaV.ny(); ++j) deltaV(i, j) = V_star(i, j) - V_star_uncorrected(i, j);
    
            {
                // RHSP_PISO + the second-corrector pressure solve, matching
                // SolveUVP.m's solve_p_piso block. Note MATLAB's own
                // agmg() call there is commented-out dead code (no
                // Linux-compatible AGMG build exists on Sherlock, per
                // SolveUVP.m's comment above it) -- PISO's P_cor_vec is
                // never actually re-solved in MATLAB, only PISO=0's path
                // is exercised there. This port does perform the real
                // solve here (via the same ilu+cg-equivalent
                // solvePressureSystem() as the first corrector), matching
                // MATLAB's intent rather than its broken/dead execution.
                // deltaU/deltaV above stay outside the window, same as
                // MATLAB's dU_star/dV_star.
                ScopedTimer t("solve_p_piso");

                // Second corrector's RHS -- the neighbor-coupling term the first (diagonal-only) correction ignored,
                // mass imbalance already corrected in SMPLE
                rhsPPiso(sysU, sysV, deltaU, deltaV, domain, controlVar, RHS_p);

                // CM_p is reused unmodified from the first corrector
                solvePressureSystem(CM_p, RHS_p, x_p, PCOR, domain, controlVar);
            }

            {
                // FORMPCOR + NEWUVP(PISO) -- SolveUVP.m's third
                // reconstruction_uvp sample, PISO-only.
                ScopedTimer t("reconstruction_uvp");
                formPCor(PCOR, domain, controlVar);

                // Second correction layers on top of the first, plus PISO's own neighbor-coupling term.
                newUVPPiso(U_star, V_star, stateVar.P, PCOR, sysU, sysV, deltaU, deltaV, variables.alpha_p,
                           bc, controlVar.ii);
                U_star_old = U_star;
                V_star_old = V_star;
                stateVar.U = U_star;
                stateVar.V = V_star;

                if (controlVar.iTime == 1) {
                    dumpIterTraceEntry(
                        controlVar.ii, true, stateVar.U,
                        "/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/debug_compare/cpp_iter_trace.json");
                }
            }

            // Deviates from ConvergenceResiduals.m here: MATLAB reuses the
            // stale, pre-either-corrector Soln.RHS_P2 for this call's mass residual. 
            // Recompute it fresh from the field PISO returned, so resi3 (and hence resi) reflects the real mpost-correction divergence
            // Left untimed, like the convergence check itself: SolveUVP.m has
            // no timer around either, so folding this extra (C++-only) rhsP
            // into a category would make the two summaries non-comparable.
            rhsP(stateVar.U, stateVar.V, domain, controlVar, ibmCoeffU, ibmCoeffV, RHS_p);
            // compute residuals 1) CM_u * U - RHS_U, 2) CM_v * V - RHS_V,  3) divergnece check, resuse RHS_p from pressure equation
            ConvergenceResult convPiso =
                convergenceResiduals(stateVar.U, stateVar.V, RHS_p, controlVar.ii, PCOR, CM_u,
                                      CM_v, RHS_u, RHS_v, domain, controlVar.disc_scheme_vel,
                                      controlVar.tol, /*PISO=*/1);
            controlVar.resi = convPiso.resi_max;
            controlVar.messageFlow = convPiso.message;
            PetscPrintf(PETSC_COMM_WORLD, "    SIMPLE ii=%d resi=%e (PISO)\n", controlVar.ii, controlVar.resi);

            break; // PISO break after one SIMPLE iteration

        }

        if (controlVar.ii > 50) {
            PetscPrintf(PETSC_COMM_WORLD, "    SIMPLE/PISO did not converge after %d iterations (resi=%e)\n",
                        controlVar.ii, controlVar.resi);
            break;
        }
    }

    // Match SolveUVP.m:328-332 (`fU = double(IBM_coeffU.flag_u == 0);
    // StateVar.U = StateVar.U.*fU;`, same for V/flag_v): zero every non-
    // fluid (ghost flag==1 or solid flag==2) cell's U/V before returning.
    // The ghost cell's internally-reconstructed value is genuinely needed
    // during the SIMPLE iterations -- it couples correctly into
    // neighboring fluid cells' own equations, and tracks MATLAB's own
    // internal value almost exactly throughout the loop (confirmed via
    // the iterTrace diagnostic above) -- but MATLAB doesn't consider it a
    // physically meaningful velocity, so it discards it here before
    // reporting/saving the final state. Without this, C++ was instead
    // saving the raw, unmasked ghost value.
    for (int i = 0; i < stateVar.U.nx(); ++i)
        for (int j = 0; j < stateVar.U.ny(); ++j)
            if (ibmCoeffU.flag(i, j) != 0.0) stateVar.U(i, j) = 0.0;
    for (int i = 0; i < stateVar.V.nx(); ++i)
        for (int j = 0; j < stateVar.V.ny(); ++j)
            if (ibmCoeffV.flag(i, j) != 0.0) stateVar.V(i, j) = 0.0;

    // Freeze this call's converged velocity as the next call's "old
    // timestep" reference for coeffU()/coeffV()'s S0 term (SolveUVP.m:267,
    // `StateVar.U_old = StateVar.U;`, run once per solveUVP() call, right
    // after its own SIMPLE loop -- not once per SIMPLE iteration, and not
    // just once at simulation start).
    stateVar.U_prev = stateVar.U;
    stateVar.V_prev = stateVar.V;

    MatDestroy(&CM_u);
    VecDestroy(&RHS_u);
    VecDestroy(&x_u);
    MatDestroy(&CM_v);
    VecDestroy(&RHS_v);
    VecDestroy(&x_v);
    MatDestroy(&CM_p);
    VecDestroy(&RHS_p);
    VecDestroy(&x_p);
}
