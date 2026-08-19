#pragma once

// Compile-time switches for the diagnostics added while validating this port
// against the MATLAB reference.
//
// These all fire inside the time/SIMPLE loop, so leaving them on buries the
// solver's real output. They are switched off rather than deleted because
// they are the instrumentation the MATLAB-vs-C++ comparison depends on, and
// re-deriving the exact probe cells and index conventions is the expensive
// part.
//
// constexpr, not #define, so the guarded code is still parsed and
// type-checked on every build and cannot silently rot while switched off --
// the optimiser drops it when the flag is false.
namespace debug {

// Per-ghost-cell mirror-point stencil dumps: which cells are classified as
// ghosts, which four corners each interpolates from, their lambda_g_k
// weights, and whether any corner is non-fluid or duplicated.
// Sites: main.cpp, SolveUVP/COEFFU.cpp
inline constexpr bool ibm_stencil = false;

// Per-solve KSP convergence: converged reason, iteration count, residual
// norm. Two lines per SIMPLE iteration (U and V).
// Sites: SolveUVP/SolveUVP.cpp
inline constexpr bool ksp = false;

// The numerator and denominator that go into each equation's scaled
// residual, for when the residual itself looks wrong and the question is
// which half moved.
// Sites: Utilities/ConvergenceResiduals.cpp
inline constexpr bool residual = false;

}  // namespace debug
