#pragma once
#include <vector>
#include <petsc.h>
#include "../VariableNonDim.h"
#include "../IBM/IBMCoeff.h"

// Nonlinear-reaction-BC counterpart to UpdateA1gPhi.h's updateA1gPhi():
// used only when controlVar.nonlinearReactionBC is true, i.e. beta
// itself is a function of the current phi field, not the fixed
// -diag(variables.A) value the linear model uses. In that case
// lambda_g_k (not just A1_g) depends on phi, so both need re-deriving
// every QUICK iteration -- see SolveTransportADRE.cpp step 1.2.
//
// Both functions REUSE this species' already-found ghost-cell geometry
// (ibmCoeffPhi.I_m/J_m, I1..I4 or I1..I6 when ibm.BQp==1, all populated
// once by LSPointIdent()/LSmirPointsBQ() -- see IBMCoeff.h) rather than
// redoing the mirror-point search or stencil-corner lookup. TODO once
// implemented: LSmirPointsBQ.cpp's b/e interpolation vectors (the
// geometry-only half of lambda_g_k = B*b + E*e -- see its own comment)
// aren't currently cached anywhere, only the already-combined
// lambda_g_k is. Since b/e depend only on stencil-corner geometry, never
// on beta, updateGhostibmLambda() only strictly needs to redo the
// (cheap) B/E part if b/e get cached too; without that, it has to redo
// the interpolation-matrix build + inversion (the expensive part) as
// well, same cost as the one-time LSmirPointsBQ() call itself.

// Recomputes this species' ghost cells' beta from the current phi field
// (nonlinear reaction rate law -- not yet defined/ported), then
// lambda_g_1..6 (and, since it depends on the same beta, A1_g) from that
// new beta plus the reused stencil geometry. Mutates ibmCoeffPhi in
// place. phi is the *whole* stateVar.phi (every species), not just this
// one, since a real reaction-rate law may need every species'
// concentration, not only this one's (mirrors calculateQG()'s own
// cross-species dependence via variables.A).
void updateGhostibmLambda(IBMCoeff &ibmCoeffPhi, const std::vector<Field2D> &phi);

// Pushes the just-updated lambda_g_1..6/A1_g (from updateGhostibmLambda()
// above) into CM's ghost rows -- same ap=1 pinning + -lambda_g_k
// mirror-point coupling coeffPhiADRE() does once for the linear case,
// but callable every QUICK iteration instead of once. domain/ls/
// betaSpecies are here for whatever additional
// re-derivation this ends up needing beyond a plain re-push (e.g. if the
// nonlinear beta also has to be re-evaluated at this point rather than
// solely inside updateGhostibmLambda) -- finalize once the nonlinear
// reaction rate law itself is defined.
void updateCoeffGhost(const Domain &domain, const IBM &ibm, const LS &ls, double betaSpecies,
                       IBMCoeff &ibmCoeffPhi, Mat CM);
