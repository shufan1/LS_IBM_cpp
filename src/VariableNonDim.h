#pragma once
#include <array>
#include <string>
#include <vector>
#include "Utilities/Field2D.h"
// Mirrors the LSCase struct built in setUpVariablesNonDim.m. `case` is a
// C++ keyword, so LSCase.case -> caseId here. These fields are DERIVED
// (computed from geometry + VariableNonDim's config), not JSON-configurable
// directly.
//
// Only case 1 (grain/circle) is populated right now -- xc/yc/diamcyl are
// exactly what LSInitialize.m's case==1 branch needs. Other cases (2:
// flat fracture, 3: uni-mineral rough fracture, 4: bi-mineral, 5/6:
// square/square_rot) need more fields (h, BoundaryCurve, etc.) once
// those geometry branches are ported -- not added yet since they'd be
// unused dead fields until then.
struct LSCase {
    int caseId = 0;  // 1=circle/grain, 2=flat fracture, 3=uni-mineral rough,
                      // 4=bi-mineral, 5=square, 6=square_rot

    // Grain-only fields (case 1).
    double xc = 0.0;
    double yc = 0.0;
    double diamcyl = 0.0;
};

// Mirrors the Grid struct built in setUpVariablesNonDim.m, right before
// the Coordinates(lx, ly, diamcyl, Grid) call. Not geometry-dependent --
// same fields regardless of which case is active. All of these ARE
// directly JSON-configurable (see VariableNonDim below).
struct Grid {
    double r = 1.1;             // stretching ratio for the non-uniform mesh
    bool expon = true;          // exponential vs. linear transition zone
    double A = 2.0;             // stretch factor
    double dvdxdy = 40.0;       // grid points per lengthUnit
    bool uniform = true;        // uniform vs. graded/stretched mesh
    double zoomedAreax = 1.5;
    double zoomedAreay = 1.5;
    double Lx_l = 1.5 / 50.0;
    double Ly_b = 0.5 / 50.0;
    double lengthUnit = 0.0;    // set from VariableNonDim::lengthUnit, not
                                 // an independent JSON key (avoids the two
                                 // going out of sync)
};

// Mirrors the DOMAIN struct built by Coordinates.m -- grid-point
// coordinates and spacings are real (see Utilities/Coordinates.h/.cpp);
// everything array-sized downstream (BC boundary arrays, StateVar.U/V/P,
// LS.psi) can now be sized correctly from imax/jmax.
//
// xu/yu/xv/yv are the staggered MAC grid axis vectors (u on vertical
// faces, v on horizontal faces); xp/yp (== xv/yu) are the scalar
// (pressure/level-set/species) nodes -- psi and phi live here too, so
// they don't need their own coordinate arrays.
//
// NOT ported: the QUICK advection-interpolation coefficients -- dead
// code for momentum (COEFFU.m/COEFFV.m never read them) and only ever
// consumed as a scalar-transport RHS correction elsewhere, so not needed
// here. The full 2D meshgrids (Xu/Yu/...) aren't stored either -- they're
// just the tensor product of these 1D axes.
struct Domain {
    int imax = 0;
    int jmax = 0;
    double lx = 0.0;
    double ly = 0.0;

    std::vector<double> xu, yu;  // length imax,   jmax+1
    std::vector<double> xv, yv;  // length imax+1, jmax
    std::vector<double> xp, yp;  // length imax+1, jmax+1 (== xv, yu)

    std::vector<double> dxu, dxv, dxp;  // diff(xu)/diff(xv)/diff(xp)
    std::vector<double> dyu, dyv, dyp;  // diff(yu)/diff(yv)/diff(yp)

    // Central-difference interpolation weights -- used both to
    // interpolate velocity onto a face (needed by ConvFlux regardless of
    // discretization scheme) and directly by COEFFU/COEFFV/COEFFP's
    // active central-difference branch. Each is 1D (varies along one
    // axis only), same storage convention as dxu/dyv/etc.
    std::vector<double> CoEWu, CoNSu;  // length imax-1, imax
    std::vector<double> CoEWv, CoNSv;  // length jmax,   jmax-1
    std::vector<double> CoEWp, CoNSp;  // length imax,   jmax

    // Per-cell volumes (U's/V's/P's own control volume), needed for each
    // equation's transient coefficient A0_p=dV/dt. Genuinely 2D (built as
    // outer products in Coordinates.m), unlike everything else above.
    Field2D dV_u;  // shape (imax,   jmax+1) -- matches U
    Field2D dV_v;  // shape (imax+1, jmax)   -- matches V
    Field2D dV_p;  // shape (imax+1, jmax+1) -- matches P

    // Phi/scalar-transport QUICK-scheme interpolation weights (Coordinates.m's
    // g1c_*/g2c_* -- computed once via QuickInterp(), consumed by
    // RhsPhiADRE.cpp's deferred-correction S_e/S_w/S_n/S_s terms). Each
    // pair only varies along the axis its face is normal to (e/w vary with
    // i, n/s vary with j) -- MATLAB stores them as full 2D arrays via an
    // outer product with an all-ones vector, but since the other axis is a
    // pure broadcast, this port keeps them 1D, same convention as
    // CoEWu/CoNSu/etc above. "_p"/"_n" suffix is the flow direction
    // (positive/negative), not P-grid/node -- matches MATLAB's own naming.
    // Zero outside the interior range the QUICK stencil is valid for
    // (needs 2 real neighbors on the upstream side).
    std::vector<double> g1c_e_p, g2c_e_p, g1c_e_n, g2c_e_n;  // length imax+1, nonzero for i in [1,imax-2]
    std::vector<double> g1c_w_p, g2c_w_p, g1c_w_n, g2c_w_n;  // length imax+1, nonzero for i in [2,imax-1]
    std::vector<double> g1c_n_p, g2c_n_p, g1c_n_n, g2c_n_n;  // length jmax+1, nonzero for j in [1,jmax-2]
    std::vector<double> g1c_s_p, g2c_s_p, g1c_s_n, g2c_s_n;  // length jmax+1, nonzero for j in [2,jmax-1]
};

// Mirrors the BC struct.
struct BC {
    // Boundary-value arrays (need DOMAIN for sizing, StateVar/LS for
    // values -- empty until then).
    std::vector<double> U_a, U_b, U_c, U_d;
    std::vector<double> V_a, V_b, V_c, V_d;

    // phi_a/b/c/d, matching MATLAB's (1, jmax+1[, Np]) shapes -- one
    // vector<double> (length jmax+1) per species, same pattern for all
    // four (phi_b/c/d only carry an Np dimension in MATLAB for phi_a,
    // but keeping the same container type here for consistency, even
    // though b/c/d stay all-zero: confirmed by tracing every consumer
    // that MATLAB never assigns them past their initial zeros()
    // allocation -- BC_e_phi/BC_n_phi/BC_s_phi are all Neumann,
    // implemented by zeroing the outward-face coefficient during
    // assembly, not by reading a stored boundary value).
    std::vector<std::vector<double>> phi_a, phi_b, phi_c, phi_d;

    // Domain-edge type codes (1=Dirichlet, 3=Neumann) -
    int BC_e_u = 3, BC_w_u = 1, BC_n_u = 1, BC_s_u = 1;
    int BC_e_v = 3, BC_w_v = 1, BC_n_v = 1, BC_s_v = 1;
    int BC_e_phi = 3, BC_w_phi = 1, BC_n_phi = 3, BC_s_phi = 3;
    int BC_n_p = 3, BC_s_p = 3, BC_w_p = 3, BC_e_p = 3;
    double P0_e = 0.0;
};

// Mirrors the IBM struct.
struct IBM {
    // Robin-BC for velocity: -alpha * dphi/dn - beta * phi = q, no-slip, no-penetration condition
    double q = 0.0, alpha = 0.0, beta = 1.0;
    // Robin-BC for scalar: -alpha_phi * dphi/dn - beta_phi(i_s) * phi =
    // q_phi(i_s), same shape as velocity's. In MATLAB, alpha_phi and
    // q_phi are one constant shared by every species (see
    // SolveTransportADRE.m:47's update_A1g call, which passes IBM.q_phi
    // with no species index) -- this port doesn't keep that restriction:
    // q_phi is per-species here since nothing about the Robin BC itself
    // requires species to share a value, only MATLAB's config happened
    // to. alpha_phi stays a shared scalar (no config knob to vary it per
    // species yet; add one the same way if that's ever needed).
    //
    // beta_phi is per-species too, populated from -diag(variables.A)
    // right after Variables::defineReactivity() runs (see main.cpp) --
    // NOT by getIBM(), which runs before A is loaded. MATLAB's own
    // beta_phi is a dead single scalar (LSIBMcoeffs.m seeds an initial
    // landa_g_k/A1_g with it that update_A1g immediately overwrites
    // before first use, and its other consumer, Da/interfaceVelocityCoeff
    // in setUpVariablesNonDim.m, is assigned and never read anywhere) --
    // this field intentionally does NOT mirror that; it stores the real
    // value instead.
    double alpha_phi = -1.0;
    std::vector<double> beta_phi;  // size Np, beta_phi[i_s] = -A[i_s*Np+i_s]
    std::vector<double> q_phi;     // size Np, independent per species
    double xc = 0.0, yc = 0.0, diamcyl = 0.0;
    int nrgrainx = 1, nrgrainy = 1;

    // values of velocity and scalar inside the solids 
    double u_inside_psi = 0.0, phi_inside_psi = 0.0;
    int BQu = 0, BQv = 0, BQp = 0; // how to interpolate mirror point value: 0 = bilinear; 1=bilinear-quadratic
    double treshold = 100.0 * 2.220446049250313e-16;  // 100*eps

    // Per-cell IBM classification for U/V, matching LSPointIdent.m's
    // convention: 0=fluid, 1=ghost, 2=solid Built by LSIBMcoeffs.m 
    Field2D flag_u;  // shape (imax,   jmax+1) -- matches U
    Field2D flag_v;  // shape (imax+1, jmax)   -- matches V
};


// Mirrors the StateVar struct. All array fields need DOMAIN for sizing
// -- empty until Coordinates() is ported.
struct StateVar {
    // U, V, P: staggered MAC grid, flattened Field2D (index = i*ny + j).
    Field2D U, V, P;

    // phi, phi_prev: MATLAB's (imax+1, jmax+1, Np) -- one Field2D
    // (imax+1 x jmax+1) per species, same "vector-of-per-species"
    // pattern as BC.phi_a, rather than trying to bolt a 3rd dimension
    // onto Field2D itself.
    std::vector<Field2D> phi, phi_prev;

    std::vector<double> P_cor_vec;

    // Old-timestep values (feed coeffU()/coeffV()'s transient S0 term,
    // A0_p*U_prev) -- distinct from coeffU()/coeffV()'s own U_star_old
    // parameter, which is the previous SIMPLE-iteration's momentum
    // guess, not the previous timestep's converged state.
    Field2D U_prev, V_prev, P_prev;
};

// Mirrors the LS struct (uni-mineral case only -- LS1/LS2/LS_s for
// bi-mineral not added, see VariableNonDim's uniMineral discussion).
struct LS {
    Field2D psi;  // signed distance to the immersed boundary, on the (xp,yp)
                   // scalar grid -- real for case 1 (see LSInitialize.cpp)

    // psi interpolated onto the U-grid (xu,yu) and V-grid (xv,yv) via
    // bilinear interpolation (Utilities/Interp.h), mirroring the interp2
    // calls in setUpVariablesNonDim.m. Used as the fluid-region mask for
    // U/V: psiU>0 / psiV>0 marks fluid nodes.
    Field2D psiU, psiV;

    // Level-set surface normals on the (xp,yp) scalar grid (nx = dpsi/dx,
    // ny = dpsi/dy of the normalized gradient), from computeLSNormals()
    // -- see VariableNonDim.cpp's getLS().
    Field2D nx, ny;

    // ---- Filled by LSeqSolve() (SolveLS/LSeqSolve.h), stage 1 ----
    // Interface extension velocity on the (xp,yp) grid: the local
    // recession speed times the surface normal, so u = speed*nx,
    // v = speed*ny. Nonzero only inside the narrow band |psi| < LSgamma;
    // zero until LSeqSolve() runs, which in the frozen-geometry milestone
    // is never. Consumed by solveHJEq() to advect psi.
    Field2D u, v;

    // MATLAB's LS.q_out and LS.beta_out are deliberately NOT ported.
    // LSVelocityExtrapolation.m writes them, and LSPointIdentnew.m:236-240
    // reads them back into q_G/beta_G to seed the scalar IBM coefficients
    // -- but SolveTransportADRE.m throws that seed away before it is ever
    // used: it computes its own q_G via calculate_qG (line 35, then again
    // every QUICK iteration at line 131), builds beta_G inline from
    // -diag(A) (lines 36-41), and update_A1g() (line 47) overwrites every
    // IBM_coeffP field LSPointIdentnew had just produced on the first
    // QUICK pass. beta_out is additionally just -diag(A), which already
    // exists here as IBM::beta_phi (see main.cpp). Carrying them would
    // only reproduce MATLAB's ordering coupling, where skipping LSeqSolve
    // makes LSPointIdentnew fail with "Unrecognized field name q_out".

    int caseId = 0;
};

// Mirrors the VARIABLES struct. Np/phi_inlet/phi_init/Re/Pe/D/density/
// dimensional/dissolution are real now. inv_A/A come from
// defineReactivity() (loads the reaction-rate matrix, see VariableNonDim
// .cpp); dtau still needs DOMAIN and is TODO.
//
struct Variables {
    double D = 0.0, Re = 0.0, density = 0.0;
    bool dimensional = false;
    double alpha_u = 0.7, alpha_v = 0.7, alpha_p = 1.0, alpha_q = 1.0;
    double dt = 0.05;  // dt_man
    double Pe = 0.0;
    int Np = 3;
    std::array<double, 3> phi_inlet{};
    std::array<double, 3> phi_init{};
    bool dissolution = true;

    // Copied from ControlVar::verbose by main.cpp -- the solver modules
    // take Variables but not ControlVar. Gates progress printing only.
    bool verbose = true;
    int n_iter_ReLS = 4;
    std::string TimeSchemeLS = "RK3";
    std::string TimeSchemeRLS = "RK3";

    double dtau = 0.0; //used for the level-set reinitialization procedure.

    // Dimensional reference scales, from setUpVariablesNonDim.m:405,457-458.
    // Used by LSVelocityExtrapolation to turn a nondimensional molar flux
    // into a physical interface recession speed.
    double u_real = 0.12e-2;   // m/s
    double c_real = 10.0;      // mol/m^3
    double L_real = 0.5e-2;    // m
    double molarVol = 36.9;    // cm^3/mol (calcite). NOTE the unit --
                                // LSVelocityExtrapolation.m:139 writes the
                                // literal 36.9e-6, i.e. this value
                                // converted to m^3/mol.

    // Pd: per-species diffusivity-like normalization, hardcoded in MATLAB
    // too (setUpVariablesNonDim.m:463). defineReactivity() divides row i
    // of A by Pd[i], so every use of A carries that scaling.
    //
    // Exposed here because LSVelocityExtrapolation has to undo it. The
    // scaled A is right for the transport BC -- each species' equation is
    // nondimensionalized by its own diffusivity -- but the interface
    // recession speed needs the TRUE molar flux, which is
    // diffusivity-independent. So that one term multiplies its species'
    // Pd back in. MATLAB writes the undo as a bare `/1.261` literal
    // (= 1/0.793 rounded); reading it from here instead keeps the two
    // ends of the cancellation provably consistent.
    std::array<double, 3> Pd{9.3, 0.793, 1.91};

    // A: the (rescaled-by-Pd) Np x Np reaction-rate matrix,
    // LSVelocityExtrapolation's biquadratic extrapolation scheme. Both
    // flattened row-major (index = i*Np + j) -- filled by defineReactivity()
    std::vector<double> A;
    std::vector<double> inv_A;

    int nLSupdate = 10;
    double dtLS =  dt;
    double LSgamma = 0.025;
    double LSbeta = 0.015;
    double Pe_vel = Pe;

    // Grid-spacing multipliers behind dtau and LSgamma. Set from
    // config.json by getVariables(); defineLSvariables() needs them
    // because it is a Variables method and never sees VariableNonDim.
    double dtau_fac = 0.5;
    double LSgamma_fac = 5.0;


    // Fills in A/inv_A . loaded from A_0.2.json. inv_A = inv(3*I - 2*l*A), used by
    void defineReactivity(const Domain &domain, const std::string &aMatJsonPath);

    // sets the Variables- side level-set parameters (nLSupdate, dtLS, LSgamma,
    // LSbeta, Pe_vel, also dtau, dimesionless dt) 
    void defineLSvariables(const Domain &domain);
};

// Mirrors setUpVariablesNonDim.m in full: the config fields (translated
// up through the Coordinates() call, JSON-configurable, same two-
// constructor pattern as ControlVar) plus one method per remaining
// MATLAB return value (LSCase, DOMAIN, VARIABLES, BC, IBM, StateVar, LS).
// See VariableNonDim.cpp for exactly what's real vs. TODO in each method
// -- everything downstream of DOMAIN stays empty until Coordinates()
// (grid generation) is actually ported; the fill logic is otherwise
// already written as if DOMAIN were real, so no further changes should
// be needed once it is.
//
// Only geometry == "grain" is implemented; "fracture"/"square"/
// "square_rot" leave lsCase/lx/ly/nrgrainx/nrgrainy at their defaults
// (see the TODO branches in VariableNonDim.cpp).
class VariableNonDim {
public:
    // ---- Domain size / geometry selection (raw config) ----
    std::string geometry = "grain";  // "grain" | "fracture" | "square" | "square_rot"
    double lengthUnit = 0.2;
    double buffer_dist = 0.0;

    // ---- Additional variables ----
    double Re = 0.6;           // kinematic viscosity (nondimensional)
    double uinflow = 1.0;
    double Pe = 600.0;         // Peclet number
    int Np = 3;                // number of chemical species tracked
    std::array<double, 3> phi_inlet{1.0, 0.0, 0.0};
    std::array<double, 3> phi_init{1.0, 0.0, 0.0};
    double density = 1.0;
    bool dimensional = false;

    // ---- BCs: domain edges (u/v/p) ----
    // Simple type codes for the four edges of the rectangular domain:
    // 1 = Dirichlet, 3 = Neumann. Different from q/alpha/beta below,
    // which are the immersed-boundary (solid surface) condition instead.
    int BQu = 0;
    int BQv = 0;
    int BC_e_u = 3, BC_w_u = 1, BC_n_u = 1, BC_s_u = 1;
    int BC_n_v = 1, BC_s_v = 1, BC_e_v = 3, BC_w_v = 1;
    int BC_n_p = 3, BC_s_p = 3, BC_e_p = 3, BC_w_p = 3;
    double P0_e = 0.0;

    // ---- Immersed boundary: Robin condition at the solid/mineral
    // surface, -alpha*(dphi/dn) - beta*phi = q, where phi is whichever
    // field (u, v, or the scalar) is being constrained. This is a
    // *different* boundary than the domain-edge BCs above -- it's
    // applied at the object's surface via the IBM ghost-cell method.
    //
    // For velocity: alpha=0, beta=1, q=0 collapses to phi=0 -- i.e. a
    // plain no-slip condition at the solid wall.
    double q = 0.0;
    double alpha = 0.0;
    double beta = 1.0;

    // For the scalar: q_phi=0, alpha_phi=-1 -- the per-species reactive-
    // flux coefficient the transport solve's ghost-cell BC actually uses
    // is beta_G(:,i_s) = -diag(variables.A)(i_s), from the reaction-rate
    // matrix (see LSIBMcoeffsPhi()). MATLAB also has a beta_phi=178 here,
    // but it's dead in the active model: LSIBMcoeffs.m only uses it to
    // seed an initial landa_g_k/A1_g that update_A1g immediately
    // overwrites before first use, and its other consumer (Da/
    // interfaceVelocityCoeff in setUpVariablesNonDim.m) is assigned and
    // never read anywhere. Confirmed by tracing every consumer; not
    // ported here.
    bool uniMineral = true;   // hardcoded true in the MATLAB source too
                               // (the bi-mineral else-branch is unreachable
                               // dead code there); not translated here.
    int BQp = 0;
    // Per-species, unlike MATLAB's own shared scalar (see IBM::q_phi's
    // comment) -- this project isn't constrained to a literal MATLAB
    // mirror here, and ibm.q_phi is already per-species internally, so
    // config.json can set each species independently instead of only
    // ever broadcasting one shared value.
    std::array<double, 3> q_phi{};
    double alpha_phi = -1.0;
    bool dissolution = true;
    int BC_e_phi = 3, BC_w_phi = 1, BC_n_phi = 3, BC_s_phi = 3;  // domain-edge BCs for phi

    // ---- Objects ----
    double u_inside_psi = 0.0;
    double phi_inside_psi = 0.0;

    // ---- Grid (feeds getDomain()) ----
    Grid grid;

    // ---- Level-set solver knobs (feed getVariables() and
    // Variables::defineLSvariables()) ----
    // All hardcoded in setUpVariablesNonDim.m:422-429,527-531. Exposed
    // here because they are genuine tuning parameters, not physics.
    //
    // CAUTION -- n_iter_ReLS and dtau_fac are NOT independent. The
    // reinitialization equation propagates outward from the interface at
    // unit pseudo-time speed, so what matters is their product:
    //
    //     repair radius (in cells) = n_iter_ReLS * dtau_fac = 4 * 0.5 = 2
    //
    // Two cells is what the defaults buy, and it is chosen to cover the
    // only things that read the distance property -- r_g = abs(psi) at
    // ghost cells (1 cell) and computeLSNormals' central difference
    // (1 cell) -- with 2x margin, then stop. Halving dtau_fac without
    // doubling n_iter_ReLS silently halves the repair radius; raising
    // n_iter_ReLS buys reach at the cost of MORE interface drift, since
    // the discrete scheme does not exactly preserve the zero contour.
    // Change them as a pair, and watch the product.
    int n_iter_ReLS = 4;                  // reinitialization sweeps per LS step
    std::string TimeSchemeLS = "RK3";     // advection      (RK1 | RK2 | RK3)
    std::string TimeSchemeRLS = "RK3";    // reinitialization
    double dtau_fac = 0.5;                // dtau   = dtau_fac   * min(dxp)
    double LSgamma_fac = 5.0;             // LSgamma = LSgamma_fac * min(dxp)
    int nLSupdate = 10;                   // transport sub-steps per LS update

    // dt_man lives on ControlVar now -- it is a run-control knob, not a
    // physical property. main.cpp copies ControlVar::dt_man into
    // Variables::dt, which is what the solver modules read.

    // ---- Derived (computed by computeDerivedGeometry(), not from JSON) ----
    LSCase lsCase;
    double lx = 0.0;
    double ly = 0.0;
    int nrgrainx = 1;
    int nrgrainy = 1;
    double D = 1.0 / 600.0;  // 1/Pe, recomputed after Pe is finalized

    // All fields at their MATLAB defaults, then derived geometry computed.
    VariableNonDim();

    // Same defaults, overridden field-by-field by whatever keys are
    // present under config.json's "geometry"/"physicalParameters"/
    // "boundaryConditions" sections (plus top-level "dt_man") -- see
    // VariableNonDim.cpp's header comment for exactly which section each
    // field lives under. Derived geometry is recomputed afterward from
    // the (possibly overridden) values. Throws std::runtime_error if the
    // file can't be read.
    explicit VariableNonDim(const std::string &configPath);

    // ---- One method per setUpVariablesNonDim.m return value besides
    // LSCase's own fields (already public above). Call in this order --
    // each one after getDomain() takes the Domain it depends on
    // explicitly, rather than this class caching it as hidden state. ----

    // Thin accessor -- lsCase is already computed by the constructor.
    LSCase getLSCase() const { return lsCase; }

    // Real now -- delegates to Utilities/Coordinates.h's
    // computeCoordinates(), mirroring Coordinates.m.
    Domain getDomain() const;

    // Real for the scalar fields, including dtau (needs DOMAIN.dxp).
    // A/inv_A are filled separately by defineReactivity(), not this
    // method.
    Variables getVariables(const Domain &domain) const;

    // Real for q/alpha/beta/q_phi/alpha_phi/phi_inside_*/xc/yc/diamcyl/
    // nrgrainx/nrgrainy/BQ*/treshold. Np sizes the copy from this class's
    // fixed-3 q_phi array into ibm.q_phi's per-species vector (q_phi
    // doesn't depend on A, unlike beta_phi, so this doesn't need to wait
    // for Variables::defineReactivity() the way beta_phi does -- see
    // IBM::beta_phi's comment).
    IBM getIBM(const Domain &domain, int Np) const;

    // Real now: sizes AND fills U/V/P/phi with their masked initial
    // values (uniform inflow / phi_init masked to the fluid region via
    // ls.psiU/psiV/psi). Must be called AFTER getLS() (needs its masks)
    // and BEFORE getBC() (which reads this).
    StateVar getStateVar(const Domain &domain, const LS &ls) const;

    // psi is real for case 1/grain (see SolveLS/LSInitialize.cpp); nx/ny
    // from computeLSNormals() (SolveLS/LSnormals.cpp). Takes lsCase
    // explicitly (even though it's also a member) for consistency with
    // every other getX() method here taking its dependencies as
    // parameters rather than reading `this->` state.
    LS getLS(const Domain &domain, const LSCase &lsCase) const;

    // Real for the BC_*/P0_e type codes. U_a/V_a/phi_a/... boundary-value
    // arrays mirror MATLAB's `U_a(:) = U(1,:)`-style fill, so they need
    // StateVar's (masked, psi-dependent) U/V/phi -- TODO until getLS()
    // is real too, since StateVar itself isn't fully real without psi.
    BC getBC(const Domain &domain, const StateVar &stateVar, const LS &ls) const;

    // Dumps the current flow/transport/level-set state (U/V/P/phi, psi)
    // plus the simulation time to a JSON file, mirroring
    // modelSimulation.m's periodic CurrentStateVar save (dataRDE<iTime>
    // dt.mat there, .json here). Call from rank 0 only. Throws
    // std::runtime_error if the file can't be written.
    void saveCurrentStateJson(const StateVar &stateVar, const LS &ls, double time,
                               const std::string &path) const;

private:
    // Fills in lsCase/lx/ly/nrgrainx/nrgrainy/D from geometry + the
    // fields above. Mirrors the "Domain size" and "Objects" (grain
    // centers) sections of setUpVariablesNonDim.m; only geometry=="grain"
    // is implemented (see VariableNonDim.cpp).
    void computeDerivedGeometry();
};
