#pragma once
#include "../VariableNonDim.h"
#include "../Utilities/Field2D.h"

// WENO5 spatial derivatives of the level set.
//
// This file does NOT map 1:1 onto a MATLAB file. It consolidates three:
//
//     SolveLS/LSdirDerivates.m   -> wenoDifferences()
//     SolveLS/LSWeights.m        -> wenoWeights()
//     SolveLS/LSFindDerivative.m -> wenoDerivative()
//
// Each dispatches on `equation`, exactly as MATLAB's own `equation`
// string argument does, because the two PDEs in SolveLS/ need different
// derivatives:
//
//   LevelSetEqn         d(psi)/dt + u*psi_x + v*psi_y = 0
//                       H(grad psi) = u*psi_x + v*psi_y, LINEAR.
//                       dH/d(psi_x) = u -- the characteristic speed is a
//                       field you already have. So the upwind side can be
//                       chosen up front from sign(u), and ONE derivative
//                       per direction is enough: psi_x, psi_y.
//
//   ReinitializationEqn d(psi)/d(tau) + sign(psi)*(|grad psi| - 1) = 0
//                       H(grad psi) = sign(psi)*(sqrt(psi_x^2+psi_y^2)-1),
//                       NONLINEAR.
//                       dH/d(psi_x) = sign(psi)*psi_x/|grad psi| -- the
//                       characteristic speed depends on psi_x ITSELF, the
//                       very thing being computed. There is no sign to
//                       branch on before the fact, so BOTH directions are
//                       computed -- psi_xn/psi_xp, psi_yn/psi_yp -- and
//                       the choice is deferred to godunovGradientNorm,
//                       which resolves it pointwise with a Godunov flux.
//
// Everything is restricted to a narrow tube around the interface, keyed
// on `psi_prev` -- the level set at the START of the geometry step, not the
// current RK stage value, so the tube cannot drift between stages. The
// two branches test it differently: LevelSetEqn pointwise, |psi_prev| < h;
// ReinitializationEqn dilated over the 9-point neighbourhood
// (LSdirDerivates.m:91-94), one cell wider in every direction, because
// its correction propagates outward and a node at the band edge needs
// computed neighbours to lean on.
enum class LSEquation { LevelSetEqn, ReinitializationEqn };

// The undivided first differences the WENO reconstruction is built from
// -- LSdirDerivates.m's `dir` struct. Only the fields belonging to the
// requested branch are filled; the rest stay default-constructed (empty,
// no allocation), same as MATLAB returning a struct with different
// fieldnames per branch.
//
//   LevelSetEqn          vx1..vx5, vy1..vy5              (10, upwind side)
//   ReinitializationEqn  vxn/vxp/vyn/vyp 1..5            (20, both sides)
struct WenoDifferences {
    Field2D vx1, vx2, vx3, vx4, vx5;
    Field2D vy1, vy2, vy3, vy4, vy5;

    Field2D vxn1, vxn2, vxn3, vxn4, vxn5;
    Field2D vxp1, vxp2, vxp3, vxp4, vxp5;
    Field2D vyn1, vyn2, vyn3, vyn4, vyn5;
    Field2D vyp1, vyp2, vyp3, vyp4, vyp5;
};

// in : psi (current stage value), psi_prev (tube reference), domain, h,
//      equation; u/v only for LevelSetEqn, which needs their sign to pick
//      the stencil. Pass nothing for ReinitializationEqn -- MATLAB takes
//      u/v there too but never reads them.
// out: the differences above, zero outside the tube and outside the
//      6-point stencil's own reach at the array edges.
WenoDifferences wenoDifferences(const Field2D &psi, const Field2D &psi_prev, const Domain &domain,
                                 double h, LSEquation equation, const Field2D *u = nullptr,
                                 const Field2D *v = nullptr);

// The normalized WENO weights -- LSWeights.m's `Weight` struct. Same
// per-branch split: 6 for LevelSetEqn, 12 for ReinitializationEqn. They
// sum to 1 across each group of three.
struct WenoWeights {
    Field2D wx1, wx2, wx3;
    Field2D wy1, wy2, wy3;

    Field2D wxn1, wxn2, wxn3;
    Field2D wxp1, wxp2, wxp3;
    Field2D wyn1, wyn2, wyn3;
    Field2D wyp1, wyp2, wyp3;
};

// Each group of five differences contains three overlapping 3-point
// sub-stencils, each a 3rd-order estimate of the same derivative. Every
// one gets a Jiang-Shu smoothness indicator -- roughly its squared
// variation of psi:
//     S1 = 13/12*(v1-2*v2+v3)^2 + 1/4*(v1-4*v2+3*v3)^2
//     S2 = 13/12*(v2-2*v3+v4)^2 + 1/4*(v2-v4)^2
//     S3 = 13/12*(v3-2*v4+v5)^2 + 1/4*(3*v3-4*v4+v5)^2
// then the ideal weights 1/10, 6/10, 3/10 are divided by (eps+S)^2 and
// renormalized.
//
// That division is the point of WENO. Where psi is smooth all three S are
// comparably small, the weights sit near their ideal values and the blend
// is 5th order. Where a kink falls inside one sub-stencil its S blows up,
// 1/(eps+S)^2 drives its weight to ~0, and the scheme falls back to a
// 3rd-order estimate built only from smooth data. Level sets grow kinks
// constantly -- corners, merging fronts, a grain pinching off.
//
// Nodes outside the tube have all differences zero, hence S = 0 and a
// finite weight d_k/eps^2, which renormalizes back to the ideal d_k and
// multiplies zeros. No division by zero; the derivative comes out 0.
//
// `eps` is the caller's, 1e-6 at LSFindDerivative.m:17. In MATLAB that
// shadows the builtin `eps` (2.2e-16) -- it is the WENO regularizer, six
// orders of magnitude larger, not machine epsilon.
//
// ReinitializationEqn runs the identical arithmetic on four groups of
// five differences instead of two -- the left- and right-biased runs in
// each direction. No group is preferred here; all four sets of weights
// survive, and the choice between the resulting derivatives is made
// later by godunovGradientNorm.
WenoWeights wenoWeights(const WenoDifferences &dir, double eps, LSEquation equation);

// The finished reconstruction. Again only the requested branch's fields
// are filled -- MATLAB's `psinp`.
//
//   LevelSetEqn          psi_x, psi_y
//   ReinitializationEqn  psi_xn, psi_xp, psi_yn, psi_yp
struct LSGradient {
    Field2D psi_x, psi_y;
    Field2D psi_xn, psi_xp, psi_yn, psi_yp;
};

// Orchestrator, as MATLAB has it: wenoDifferences, then wenoWeights, then
// blend the three candidates. Owns no numerics itself.
LSGradient wenoDerivative(const Field2D &psi, const Field2D &psi_prev, const Domain &domain, double h,
                           LSEquation equation, const Field2D *u = nullptr,
                           const Field2D *v = nullptr);
