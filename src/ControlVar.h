#pragma once
#include <string>

// Mirrors beta_AD_correctu/setUpControlVar.m. Fields below are listed in
// the same order as that file, with the same defaults.
//
// Note: setUpControlVar.m also takes a DOMAIN argument, but only uses it
// to compute a StateVar field that the function never actually returns
// (dead code in the original) -- so DOMAIN isn't needed here.
//
// This only covers the fields setUpControlVar.m sets directly (the
// config/-loadable ones). Several more ControlVar fields get added
// dynamically elsewhere in the MATLAB code as the solve progresses
// (e.g. SolveUVP.m sets .resi/.ii/.messageFlow, SolveTransportADRE.m sets
// .err_q/.iter_qq, modelSimulation.m sets .PISO/.noTime/.endTime) --
// those aren't added yet since the modules that use them aren't ported.
class ControlVar {
public:
    std::string output_folder = "Output_fracture_uni_mineral_rough/";

    double time = 0.0;

    // The simulation timestep -- MATLAB's dt_man (setUpVariablesNonDim.m
    // :259), moved here because it is a run-control knob, not a physical
    // property: it sits alongside noLStime/savedat/tolerances rather than
    // with Re/Pe/geometry. main.cpp copies it into Variables::dt, which
    // is what the solver modules actually read.
    //
    // It is the geometry (level-set) step. The transport sub-step is
    // dt_man/nLSupdate -- see main.cpp's model-4 block. Run length is
    // noLStime * dt_man in simulated time.
    //
    // (MATLAB's dt_diff/dt_cour, the diffusion and Courant limits, are
    // computed there and discarded: `dt = dt_man` unconditionally. Not
    // ported for the same reason.)
    double dt_man = 0.05;

    int savedat = 1;
    // MATLAB: floor(0.01 / VARIABLES.dt) -- depends on VARIABLES, which
    // isn't ported yet. Left at a plain default here; set it explicitly
    // once VARIABLES.dt exists, rather than baking that dependency into
    // this class's constructor.
    int rat = 0;

    double tol = 5e-3;
    int f = 0;  // checked as a flag (`ControlVar.f==0`) in SolveUVP.m's outer loop
    double tolbicg = 1e-5;  // max tol for bicgstab
    int maxit = 2000;       // max iter for bicgstab

    int disc_scheme_vel = 2;
    // Which of MATLAB's 4 numbered models (1-4, each a different
    // steady/unsteady combination for flow+transport+LS) is active,
    // expressed directly as these two flags instead of a model number.
    // flow_steady=true, transport_steady=false is "model 4": steady flow,
    // unsteady transport and LS -- the only one this project implements
    // (see main.cpp).
    bool flow_steady = true;
    std::string imposePresBC = "top";  // "bottom", "top", or "middle"

    bool transport_steady = false;

    int noLStime = 5;
    int iStart = 0;
    double endTime = 0;

    // Outer timestep counter (main.cpp's own `iTime` loop variable,
    // mirrored here so SolveUVP.cpp's per-SIMPLE-iteration debug dumps
    // can gate on "first timestep only" the same way other debug dumps
    // already gate on controlVar.ii==1 for "first SIMPLE iteration only".
    int iTime = 0;

    double tol_q = 5e-4;
    double tolbicg_c = 5e-5;
    int maxit_c = 2000;

    // SolveTransportADRE's own Hayase QUICK-loop bookkeeping (mirrors
    // SolveUVP's ii/resi pattern above). err_q is the max iterate-to-
    // iterate change in phi across all species this pass (NOT a true
    // matrix residual -- see SolveTransportADRE.cpp); iter_qq is the pass
    // counter, capped at 250.
    double err_q = 1.0;
    int iter_qq = 0;

    // True only if the reactive ghost-cell boundary condition's beta
    // depends on phi/geochemistry itself (e.g. real pH-dependent
    // kinetics), rather than the fixed linear reaction-rate matrix this
    // project currently uses (variables.A, loaded once). Gates whether
    // SolveTransportADRE's QUICK loop needs to re-derive lambda_g_k via
    // updateGhostibmLambda()/updateCoeffGhost() every iteration --
    // always false for now.
    bool nonlinearReactionBC = false;

    bool ADRE = true;  // advection-diffusion-reaction vs. diffusion-reaction transport

    // Progress/diagnostic printing. true (default): the one-time setup
    // parameter dump, the per-geometry-step banner, and the per-step [LS]
    // band/CFL/solid-cell lines. false: only warnings, errors and the
    // final timing summary -- for batch runs where the log is noise.
    //
    // Distinct from `debug` above, which controls what is written to DISK
    // (the per-ADRE-substep dumps), not what is written to stdout. Also
    // distinct from the compile-time switches in Utilities/Debug.h, which
    // gate the heavyweight per-ghost-cell/per-KSP traces.
    //
    // main.cpp copies this into Variables::verbose, because the solver
    // modules (LSeqSolve and below) take Variables but not ControlVar.
    bool verbose = true;

    // true: skip LSeqSolve() entirely every outer iTime -- ls.psi never
    // advances, matching modelSimulation_no_LS.m's own frozen-geometry
    // behavior (Milestone 1's validation scope: flow+transport alone,
    // against a fixed immersed boundary). false (default): call
    // LSeqSolve() every outer iTime, the real moving-geometry model
    // (modelSimulation.m). MATLAB has this split as two different
    // top-level scripts (RunADRE_no_LS.m vs RunADRE.m calling two
    // different modelSimulation* functions); this project has one
    // binary/main.cpp, so it's a config flag instead.
    bool freezeLS = false;

    // DEBUG ONLY, controls how much gets saved to disk for a run --
    // false (default): only the normal per-controlVar.savedat
    // dataRDE<iTime>dt.json snapshot, to output_folder as configured.
    // true: ALSO the per-ADRE-substep dumps (dumpGhostDebug/
    // dumpForMatlabComparison in SolveTransportADRE.cpp, plus the
    // per-substep state dump in flowSStransportUS.cpp) to
    // <projectDir>/output_ADRE, AND redirects the normal per-savedat
    // snapshot to output_folder with "_debug" appended, so a long
    // diagnostic run's saves never mix with/overwrite a normal run's
    // output_folder. See main.cpp's own outputDir/outputAdreDir setup.
    bool debug = false;

    // Set dynamically by modelSimulation/flowSStransportUS in MATLAB
    // (flow_steady/transport_steady already declared above).
    int PISO = 0;
    int noTime = 10;

    // Set by SolveUVP once the flow solve is ported; printed by
    // flowSStransportUS. Empty until then.
    std::string messageFlow;

    // SolveUVP's own outer-loop bookkeeping (SolveUVP.m dynamically sets
    // ControlVar.ii/.resi -- these mirror that). ii is the iteration
    // counter (loop requires ii>=2 before resi is even checked); resi is
    // the worst-of-{U-momentum, V-momentum, continuity} residual, driving
    // the `resi>tol` convergence check.
    int ii = 0; // SIMPLE/PISO iteration
    double resi = 1e10;

    // Linear solver choice per equation, as plain PETSc KSPType/PCType
    // strings ("bcgs", "cg", "gmres", ... / "ilu", "gamg", "hypre",
    // "jacobi", ...) -- passed straight to KSPSetType()/PCSetType() once
    // the real solves are wired up (SolveUVP.cpp's TODOs). Not hardcoded
    // to MATLAB's bicgstab/pcg choice: momentum's matrix is non-symmetric
    // (convection + IBM ghost coupling), so it needs a general Krylov
    // method (bicgstab -- PETSc's registered name is "bcgs", not
    // "bicgstab" -- is a reasonable default: fixed, low memory cost per
    // iteration, unlike GMRES's growing Krylov basis). Pressure's matrix
    // is SPD (a Poisson-type operator), so CG is the right *method* --
    // but plain/ILU-preconditioned CG's iteration count scales with the
    // condition number (~h^-2 for a Poisson operator), a bottleneck
    // already identified in this project's MATLAB profiling work.
    // Defaulting pressure's preconditioner to "gamg" (algebraic
    // multigrid) rather than "ilu" reflects that -- override via
    // config.json if you want to compare against ILU/BoomerAMG/etc.

    // Krylov Subpace (method) PETSc's name for:  iterative linear solver itself
    std::string ksp_type_momentum = "bcgs";
    std::string pc_type_momentum = "ilu";
    std::string ksp_type_pressure = "cg";
    std::string pc_type_pressure = "gamg";

    // Scalar transport (phi): MATLAB's own bicgstab()+ilu(...,'nofill')
    // choice, unmodified -- unlike pressure, there's no known conditioning
    // problem here motivating a deviation from the MATLAB reference.
    std::string ksp_type_scalar = "bcgs";
    std::string pc_type_scalar = "ilu";

    // All fields at their MATLAB defaults (matches setUpControlVar.m as-is).
    ControlVar() = default;

    // Same defaults, then overridden field-by-field by whatever keys are
    // present under config.json's top-level "ControlVar" section
    // (missing section or keys keep their default). Throws
    // std::runtime_error if the file can't be read.
    explicit ControlVar(const std::string &configPath);
};
