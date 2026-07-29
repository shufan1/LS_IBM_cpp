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
    double timedt = 50000.0;
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

    double tol_q = 5e-4;
    double tolbicg_c = 5e-5;
    int maxit_c = 2000;

    bool ADRE = true;  // advection-diffusion-reaction vs. diffusion-reaction transport

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

    // All fields at their MATLAB defaults (matches setUpControlVar.m as-is).
    ControlVar() = default;

    // Same defaults, then overridden field-by-field by whatever keys are
    // present under config.json's top-level "ControlVar" section
    // (missing section or keys keep their default). Throws
    // std::runtime_error if the file can't be read.
    explicit ControlVar(const std::string &configPath);
};
