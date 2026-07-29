// Mirrors beta_AD_correctu/Models/modelSimulation.m (and, for Milestone 1,
// modelSimulation_no_LS.m). Called once from main.cpp as runModel(...);
// owns the outer time loop and all solver detail -- main.cpp stays thin.
//
// Roadmap only below -- nothing implemented yet. Order matches the staged
// plan: Milestone 1 (flow, frozen geometry) first, validated against the
// MATLAB _no_LS reference trace, before Milestone 2 (+transport) and
// Milestone 3 (+level-set/moving geometry). See
// beta_AD_correctu/PARALLELIZATION_PLAN.md for the full rationale.

// ================================================================
// Milestone 1: flow only, on a FROZEN geometry (no LS/IBM porting
// yet). Goal: match the MATLAB _no_LS reference trace on 1 rank,
// then confirm rank-count-independent results + speedup on N ranks.
// ================================================================

// 1. Domain decomposition / grid setup
//    - mirrors Utilities/Coordinates.m
//    - PETSc DMStag, 2D block decomposition (per the plan: comm
//      volume ~sqrt(P), 4 neighbors/rank)
//    - staggered MAC layout: u/v on faces, p on centers
//    - ghost width 1 for flow/pressure (level-set band needs width 3,
//      not relevant until Milestone 3)

// 2. Load fixed geometry + IBM coefficients (Milestone 1/2 only)
//    - read IBM_coeffU/V/P (or the LS.psi they came from) exported
//      from the MATLAB _no_LS harness (RunADRE_no_LS.m)
//    - no LSInitialize/LSIBMcoeffs ported yet -- this is static input

// 3. Load boundary conditions + fluid/flow parameters
//    - mirrors BC struct + relevant VARIABLES fields from
//      setUpVariablesNonDim.m (Re, uinflow, BC_*_u/v/p, etc.)

// 4. Assemble momentum coefficients (per outer Picard/SIMPLE iterate)
//    - mirrors SolveUVP/COEFFU.m, COEFFV.m -> PETSc Mat via
//      DMStag + MatSetValuesStencil
//    - halo exchange (non-blocking, overlapped with interior compute)

// 5. Assemble + solve pressure Poisson
//    - mirrors SolveUVP/COEFFP.m, FORMPCOR.m
//    - PETSc KSP: CG + GAMG/BoomerAMG, constant null space (pure
//      Neumann), reuse AMG hierarchy across outer iterations
//    - this is the plan's make-or-break scaling test

// 6. Solve momentum (u, v)
//    - mirrors SolveUVP/NEWUVP.m, FORMUV.m
//    - PETSc KSP: BiCGStab/GMRES + block-Jacobi/ILU(0) or ASM

// 7. Outer SIMPLE/PISO Picard loop
//    - mirrors SolveUVP/SolveUVP.m's `while (resi>tol || ii<2)`
//    - re-assemble coefficients each pass with the latest u/v/p
//      iterate (lagged-coefficient nonlinearity handling, same as
//      MATLAB -- no Jacobian needed)
//    - global residual norms via PETSc VecNorm/MPI_Allreduce,
//      batched to minimize synchronization points

// 8. Checkpoint + compare
//    - write U/V/P via PETSc binary Vec output
//    - diff against the MATLAB _no_LS reference trace (1 rank first,
//      then confirm rank-count independence)
//    - benchmark strong scaling once correctness is confirmed

// ================================================================
// Milestone 2 (later): + scalar transport, still frozen geometry
//    - mirrors SolveTransport/{COEFFPHIADRE,RHSPHIADRE,
//      SolveTransportADRE}.m -> PETSc Mat/KSP (BiCGStab + ILU/ASM)
//    - same QUICK-loop Picard structure as MATLAB
// ================================================================

// ================================================================
// Milestone 3 (later): + level-set update + IBM recomputation
//    - mirrors SolveLS/LSeqSolve.m + IBM/LSIBMcoeffs.m
//    - ghost width 3 (WENO5 reinit stencil), separate from the
//      width-1 flow/pressure halo
//    - trickiest in parallel (IBM interpolation across rank
//      boundaries, load imbalance near the moving front) -- saved
//      for last per the plan
// ================================================================

// void runModel(...) { }
