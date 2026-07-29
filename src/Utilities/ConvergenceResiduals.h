#pragma once
#include <petsc.h>
#include <string>
#include "../VariableNonDim.h"

// resi1/resi2: normalized max residual of the U-/V-momentum equations
// (max|CM*x-RHS| / max|diag(CM).*x|), evaluated at the post-correction
// velocity field. resi3: max of the mass-imbalance RHS (continuity
// residual) passed in as RHS_P2 -- see that parameter's own comment
// below for what it should hold. resi4: max|PCOR| over the interior --
// logged but not part of the convergence test. resi_max = max(resi1,
// resi2,resi3) is what drives solveUVP()'s loop.
struct ConvergenceResult {
    double resi1, resi2, resi3, resi4;
    double resi_max;
    std::string message;
};

// Ported from Utilities/ConvergenceResiduals.m. Unlike coeffU()/coeffV(),
// disc_scheme doesn't gate the residual computation itself -- resi1-4/
// resi_max are computed unconditionally, exactly as MATLAB does. Only
// the message text depends on it: disc_scheme==1 and ==2 share the same
// wording in the source (both implemented here); ==3 (QUICK) has its
// own distinct wording that needs iter_qq_u/v, err_q_u/v -- nothing in
// this port produces those, so that case throws std::runtime_error.
// BC.BC_e_p==1's alternate sizing isn't implemented either, matching
// coeffU()/coeffV(), which never build that system in the first place.
//
// U/V: stateVar.U/V *after* newUVP()/newUVPPiso() has applied this
// iteration's correction. RHS_P2: the mass-imbalance RHS to report as
// resi3 -- solveUVP()'s plain-SIMPLE call passes rhsP()'s pre-correction
// output, matching ConvergenceResiduals.m's own Soln.RHS_P2 reuse there;
// its PISO call instead passes a *fresh* rhsP(stateVar.U, stateVar.V,...)
// call on the just-corrected field, a deliberate deviation from MATLAB
// (which reuses that same stale pre-correction value there too) so PISO's
// resi3 reflects the actual post-correction divergence rather than a
// number from before either corrector ran. CM_u/RHS_U, CM_v/RHS_V: this
// iteration's already-assembled momentum systems (coeffU()/coeffV()'s
// own CM/RHS, unmodified since that call). PCOR: this iteration's
// pressure correction, after formPCor().
ConvergenceResult convergenceResiduals(const Field2D &U, const Field2D &V,
                                        Vec RHS_P2, int ii, const Field2D &PCOR,
                                        Mat CM_u, Mat CM_v, Vec RHS_U, Vec RHS_V,
                                        const Domain &domain, int disc_scheme,
                                        double tol, int PISO);
