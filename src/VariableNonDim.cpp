// Mirrors beta_AD_correctu/setUpVariablesNonDim.m (LS_IBM_matlab's copy,
// which currently has geometry = "grain" active) in full: the
// constructor/computeDerivedGeometry() handle everything up through the
// Coordinates() call (JSON-configurable, see the header comment below
// for the config.json schema); the getX() methods mirror the 7 MATLAB
// return values, one function each, called in dependency order from
// main.cpp.
//
// config.json layout for the constructor's fields (see ControlVar.cpp
// for the separate "ControlVar" section):
//   "geometry":            geometry, lengthUnit, buffer_dist, grid.*
//   "physicalParameters":  Re, uinflow, Pe, Np, phi_inlet, phi_init,
//                          density, dimensional
//   "boundaryConditions":  u/v/p/phi domain-edge BC codes, plus
//                          "immersedBoundary" (the Robin-BC coefficients
//                          at the solid/mineral surface -- q/alpha/beta
//                          for velocity, q_phi/alpha_phi/beta_phi for
//                          the scalar, and u_inside_psi/phi_inside_psi)
//   "dt_man":              top-level (not under any section above)
#include "VariableNonDim.h"
#include "Utilities/Coordinates.h"
#include "Utilities/Field2D.h"
#include "Utilities/Interp.h"
#include "SolveLS/LSInitialize.h"
#include <json-c/json.h>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace {
// Small helper: get object[key] as an int/double/bool/string if present,
// leave `field` untouched otherwise. Cuts down on repetition below.
void getIfPresent(json_object *obj, const char *key, double &field) {
    json_object *v;
    if (obj && json_object_object_get_ex(obj, key, &v)) field = json_object_get_double(v);
}
void getIfPresent(json_object *obj, const char *key, int &field) {
    json_object *v;
    if (obj && json_object_object_get_ex(obj, key, &v)) field = json_object_get_int(v);
}
void getIfPresent(json_object *obj, const char *key, bool &field) {
    json_object *v;
    if (obj && json_object_object_get_ex(obj, key, &v)) field = json_object_get_boolean(v);
}
void getIfPresent(json_object *obj, const char *key, std::string &field) {
    json_object *v;
    if (obj && json_object_object_get_ex(obj, key, &v)) field = json_object_get_string(v);
}
void getIfPresent(json_object *obj, const char *key, std::array<double, 3> &field) {
    json_object *v;
    if (obj && json_object_object_get_ex(obj, key, &v) && json_object_is_type(v, json_type_array)) {
        for (int k = 0; k < 3 && k < (int)json_object_array_length(v); ++k)
            field[k] = json_object_get_double(json_object_array_get_idx(v, k));
    }
}
// Returns obj[key] as a json_object*, or nullptr if absent -- for
// stepping into nested sections without crashing on a missing one.
json_object *getSection(json_object *obj, const char *key) {
    json_object *v = nullptr;
    if (obj) json_object_object_get_ex(obj, key, &v);
    return v;
}
}  // namespace

VariableNonDim::VariableNonDim() {
    computeDerivedGeometry();
}

VariableNonDim::VariableNonDim(const std::string &configPath) : VariableNonDim() {
    json_object *root = json_object_from_file(configPath.c_str());
    if (!root) {
        throw std::runtime_error("VariableNonDim: failed to read config file: " + configPath);
    }

    json_object *geomSec = getSection(root, "geometry");
    getIfPresent(geomSec, "geometry", geometry);
    getIfPresent(geomSec, "lengthUnit", lengthUnit);
    getIfPresent(geomSec, "buffer_dist", buffer_dist);
    json_object *gridSec = getSection(geomSec, "grid");
    getIfPresent(gridSec, "r", grid.r);
    getIfPresent(gridSec, "expon", grid.expon);
    getIfPresent(gridSec, "A", grid.A);
    getIfPresent(gridSec, "dvdxdy", grid.dvdxdy);
    getIfPresent(gridSec, "uniform", grid.uniform);
    getIfPresent(gridSec, "zoomedAreax", grid.zoomedAreax);
    getIfPresent(gridSec, "zoomedAreay", grid.zoomedAreay);
    getIfPresent(gridSec, "Lx_l", grid.Lx_l);
    getIfPresent(gridSec, "Ly_b", grid.Ly_b);

    json_object *physSec = getSection(root, "physicalParameters");
    getIfPresent(physSec, "Re", Re);
    getIfPresent(physSec, "uinflow", uinflow);
    getIfPresent(physSec, "Pe", Pe);
    getIfPresent(physSec, "Np", Np);
    getIfPresent(physSec, "phi_inlet", phi_inlet);
    getIfPresent(physSec, "phi_init", phi_init);
    getIfPresent(physSec, "density", density);
    getIfPresent(physSec, "dimensional", dimensional);

    json_object *bcSec = getSection(root, "boundaryConditions");
    json_object *bcU = getSection(bcSec, "u");
    getIfPresent(bcU, "BC_e", BC_e_u);
    getIfPresent(bcU, "BC_w", BC_w_u);
    getIfPresent(bcU, "BC_n", BC_n_u);
    getIfPresent(bcU, "BC_s", BC_s_u);

    json_object *bcV = getSection(bcSec, "v");
    getIfPresent(bcV, "BC_n", BC_n_v);
    getIfPresent(bcV, "BC_s", BC_s_v);
    getIfPresent(bcV, "BC_e", BC_e_v);
    getIfPresent(bcV, "BC_w", BC_w_v);

    json_object *bcP = getSection(bcSec, "p");
    getIfPresent(bcP, "BC_n", BC_n_p);
    getIfPresent(bcP, "BC_s", BC_s_p);
    getIfPresent(bcP, "BC_e", BC_e_p);
    getIfPresent(bcP, "BC_w", BC_w_p);
    getIfPresent(bcP, "P0_e", P0_e);

    json_object *bcPhi = getSection(bcSec, "phi");
    getIfPresent(bcPhi, "BC_e", BC_e_phi);
    getIfPresent(bcPhi, "BC_w", BC_w_phi);
    getIfPresent(bcPhi, "BC_n", BC_n_phi);
    getIfPresent(bcPhi, "BC_s", BC_s_phi);

    json_object *ibSec = getSection(bcSec, "immersedBoundary");
    json_object *ibVel = getSection(ibSec, "velocity");
    getIfPresent(ibVel, "q", q);
    getIfPresent(ibVel, "alpha", alpha);
    getIfPresent(ibVel, "beta", beta);
    getIfPresent(ibVel, "BQu", BQu);
    getIfPresent(ibVel, "BQv", BQv);

    json_object *ibScalar = getSection(ibSec, "scalar");
    getIfPresent(ibScalar, "q_phi", q_phi);
    getIfPresent(ibScalar, "alpha_phi", alpha_phi);
    getIfPresent(ibScalar, "beta_phi", beta_phi);
    getIfPresent(ibScalar, "BQp", BQp);
    getIfPresent(ibScalar, "uniMineral", uniMineral);
    getIfPresent(ibScalar, "dissolution", dissolution);

    getIfPresent(ibSec, "u_inside_psi", u_inside_psi);
    getIfPresent(ibSec, "phi_inside_psi", phi_inside_psi);

    getIfPresent(root, "dt_man", dt_man);

    json_object_put(root);

    // geometry/lengthUnit/Pe (and anything else derived depends on) may
    // have just been overridden -- recompute.
    computeDerivedGeometry();
}

void VariableNonDim::computeDerivedGeometry() {
    D = 1.0 / Pe;
    grid.lengthUnit = lengthUnit;  // MATLAB: Grid.lengthUnit = lengthUnit
                                    // (not an independent JSON key)

    // ================== Domain size ==================
    double diamcyl = 0.0;

    if (geometry == "grain") {
        // ======= circular grain =======
        lsCase.caseId = 1;  // 1 ==> circles

        diamcyl = 2 * lengthUnit;  // Diameter of cylinders

        nrgrainx = 1;  // Number of cylinders in one row
        nrgrainy = 1;

        double S = 2 * diamcyl;  // Space between centers of objects

        double freeEast = 2 * diamcyl;   // exit length after cylinders
        double freeWest = 2 * diamcyl;   // entrance length after
        double freeNorth = 0.75 * diamcyl;
        double freeSouth = 0.75 * diamcyl;

        lx = freeWest + nrgrainx * diamcyl + (nrgrainx - 1) * S + freeEast;
        ly = freeSouth + nrgrainy * diamcyl + (nrgrainy - 1) * S + freeNorth;

        // ---- Objects: define grain centers ----
        // MATLAB computes xcent/ycent arrays here too, but for a single
        // grain (nrgrainx=nrgrainy=1) the code immediately overrides xc/yc
        // with plain domain-center scalars -- that's the only path taken,
        // so xcent/ycent themselves are never actually used.
        lsCase.xc = lx / 2.0;
        lsCase.yc = ly / 2.0;
        lsCase.diamcyl = diamcyl;

    } else if (geometry == "fracture") {
        // TODO: rough/flat fracture geometry -- not ported yet.
        // MATLAB: LSCase.case=2, b/lambda/a/lx/x_0/ly derived from b;
        // diamcyl/xc/yc/nrgrainx/nrgrainy left empty.
    } else if (geometry == "square") {
        // TODO: square grain geometry + grain-center placement -- not
        // ported yet.
    } else if (geometry == "square_rot") {
        // TODO: rotated square grain geometry + grain-center placement --
        // not ported yet.
    }

}

Domain VariableNonDim::getDomain() const {
    return computeCoordinates(lx, ly, lsCase.diamcyl, grid);
}

Variables VariableNonDim::getVariables(const Domain &domain) const {
    Variables v;
    v.D = D;
    v.Re = Re;
    v.density = density;
    v.dimensional = dimensional;
    // intVelCoeff = beta_phi * molarVol / D (MATLAB: molarVol = 36.9,
    // hardcoded there too for the non-bi-mineral case).
    v.intVelCoeff = beta_phi * 36.9 / D;
    v.intVelMethod = 2;
    v.alpha_u = 0.7;
    v.alpha_v = 0.7;  // MATLAB: alpha_v = alpha_u
    v.alpha_p = 1.0;
    v.alpha_q = 1.0;
    v.dt = dt_man;
    v.Pe = Pe;
    v.Np = Np;
    v.phi_inlet = phi_inlet;
    v.phi_init = phi_init;
    v.dissolution = dissolution;
    v.n_iter_ReLS = 4;
    v.TimeSchemeLS = "RK3";
    v.TimeSchemeRLS = "RK3";
    // v.A/v.inv_A: filled by defineReactivity(), not here.
    return v;
}



BC VariableNonDim::getBC(const Domain &domain, const StateVar &stateVar, const LS &ls) const {
    BC bc;
    bc.BC_e_u = BC_e_u; bc.BC_w_u = BC_w_u; bc.BC_n_u = BC_n_u; bc.BC_s_u = BC_s_u;
    bc.BC_e_v = BC_e_v; bc.BC_w_v = BC_w_v; bc.BC_n_v = BC_n_v; bc.BC_s_v = BC_s_v;
    bc.BC_e_phi = BC_e_phi; bc.BC_w_phi = BC_w_phi; bc.BC_n_phi = BC_n_phi; bc.BC_s_phi = BC_s_phi;
    bc.BC_n_p = BC_n_p; bc.BC_s_p = BC_s_p; bc.BC_w_p = BC_w_p; bc.BC_e_p = BC_e_p;
    bc.P0_e = P0_e;

    // U boundary values (setUpVariablesNonDim.m: `U_a(:) = U(1,:);
    // U_c(:) = 0; U_d(:) = 0;`) -- a=west, b=east, c=south, d=north.
    // U_a  hold the actual uinflow value once getStateVar() fills 
    // U_b is never assigned (BC_e_u is Neumann, no stored value needed); 
    // U_c/U_d are explicitly zero (no-slip south/north walls).
    bc.U_a.resize(domain.jmax + 1);
    for (int j = 0; j < domain.jmax + 1; ++j) bc.U_a[j] = stateVar.U(0, j);
    bc.U_b.assign(domain.jmax + 1, 0.0);
    bc.U_c.assign(domain.imax, 0.0);
    bc.U_d.assign(domain.imax, 0.0);

    // V boundary values -- MATLAB never reassigns these past their
    // initial zeros() allocation (V's BCs -- no-slip north/south, no
    // inflow at west/east -- are already satisfied by V==0 everywhere),
    // so this is purely correct sizing, not a stand-in for missing logic.
    bc.V_a.assign(domain.jmax, 0.0);
    bc.V_b.assign(domain.jmax, 0.0);
    bc.V_c.assign(domain.imax + 1, 0.0);
    bc.V_d.assign(domain.imax + 1, 0.0);

    // phi_a: (1, jmax+1, Np) in MATLAB -- one zero-filled vector<double>
    // (length jmax+1) per species. Real per-species inlet values need
    // LS.psi (masks by psi(1,:)>0) -- still TODO, needs LSInitialize().
    bc.phi_a.assign(Np, std::vector<double>(domain.jmax + 1, 0.0));

    // phi_b/c/d: MATLAB never assigns these past their initial zeros()
    // allocation (BC_e_phi/BC_n_phi/BC_s_phi are Neumann, applied by
    // zeroing a coefficient during assembly, not by reading a stored
    // value) 
    bc.phi_b.assign(Np, std::vector<double>(domain.jmax + 1, 0.0));
    bc.phi_c.assign(Np, std::vector<double>(domain.jmax + 1, 0.0));
    bc.phi_d.assign(Np, std::vector<double>(domain.jmax + 1, 0.0));
    return bc;
}

IBM VariableNonDim::getIBM(const Domain &domain) const {
    IBM ibm;
    ibm.q = q; ibm.alpha = alpha; ibm.beta = beta;
    ibm.q_phi = q_phi; ibm.alpha_phi = alpha_phi; ibm.beta_phi = beta_phi;
    ibm.u_inside_psi = u_inside_psi;
    ibm.phi_inside_psi = phi_inside_psi;
    ibm.xc = lsCase.xc; ibm.yc = lsCase.yc; ibm.diamcyl = lsCase.diamcyl;
    ibm.nrgrainx = nrgrainx; ibm.nrgrainy = nrgrainy;
    ibm.BQu = BQu; ibm.BQv = BQv; ibm.BQp = BQp;
    ibm.treshold = 100.0 * 2.220446049250313e-16;

    // flag_u/flag_v: sized correctly, all-zero (every cell fluid) --
    // Field2D's fill default already gives this, LSIBMcoeffs.m isn't
    // ported so there's nothing real to fill them with yet.
    ibm.flag_u = Field2D(domain.imax, domain.jmax + 1);
    ibm.flag_v = Field2D(domain.imax + 1, domain.jmax);

    // // beta_phi (MATLAB: beta_rand = ones(imax,jmax)*beta_phi) -- real now
    // // that domain.imax/jmax come from computeCoordinates().
    // ibm.beta_phi.assign(static_cast<size_t>(domain.imax) * static_cast<size_t>(domain.jmax), beta_phi);
    return ibm;
}



LS VariableNonDim::getLS(const Domain &domain, const LSCase &lsCase) const {
    LS ls;
    ls.caseId = lsCase.caseId;

    // psi = LSInitialize(DOMAIN, LSCase): real now for case 1 (grain),
    // see SolveLS/LSInitialize.cpp.
    ls.psi = computeLSInitialize(domain, lsCase);

    // psiU/psiV = interp2(psi onto the U-/V-grid), mirroring
    // setUpVariablesNonDim.m:367-371. psi lives on the (xp,yp) scalar
    // grid; interpolate it to each U-node (xu[i], yu[j]) and V-node
    // (xv[i], yv[j]). Used as the fluid-region mask for U/V downstream.
    ls.psiU = Field2D(domain.imax, domain.jmax + 1);
    for (int i = 0; i < domain.imax; ++i)
        for (int j = 0; j < domain.jmax + 1; ++j)
            ls.psiU(i, j) = bilinearInterp(domain.xp, domain.yp, ls.psi, domain.xu[i], domain.yu[j]);

    ls.psiV = Field2D(domain.imax + 1, domain.jmax);
    for (int i = 0; i < domain.imax + 1; ++i)
        for (int j = 0; j < domain.jmax; ++j)
            ls.psiV(i, j) = bilinearInterp(domain.xp, domain.yp, ls.psi, domain.xv[i], domain.yv[j]);

    // TODO: nx/ny need LSnormals() -- not ported yet.
    return ls;
}

StateVar VariableNonDim::getStateVar(const Domain &domain, const LS &ls) const {
    // Mirrors setUpVariablesNonDim.m's "Define matrix variables /
    // storage" + "initialize U, V, phi" sections. U/V/P/phi are sized AND
    // given their masked initial values, using ls.psiU/psiV/psi as the
    // fluid-region masks.
    const int imax = domain.imax;
    const int jmax = domain.jmax;

    StateVar state_var;
    state_var.U = Field2D(imax, jmax + 1);        // (imax, jmax+1)
    state_var.V = Field2D(imax + 1, jmax);        // (imax+1, jmax), stays 0
    state_var.P = Field2D(imax + 1, jmax + 1);    // (imax+1, jmax+1), stays 0

    // ---- U:  U starts as uniform inflow everywhere in the fluid, zero inside the grain.
    for (int i = 0; i < imax; ++i)
        for (int j = 1; j < jmax; ++j)           // interior cols == MATLAB 2:end-1
            state_var.U(i, j) = uinflow;
    for (int i = 0; i < imax; ++i)
        for (int j = 0; j <= jmax; ++j)
            if (!(ls.psiU(i, j) > 0.0)) state_var.U(i, j) = 0.0;  // solid -> 0


    // ---- V: start at all 0 . ----

    // ---- phi: per species, phi_init masked to the fluid region (psi>=0), 
    // then the west-inlet row set to phi_inlet where psi>0.
    state_var.phi.assign(Np, Field2D(imax + 1, jmax + 1));
    for (int s = 0; s < Np; ++s) {
        Field2D &ph = state_var.phi[s];
        for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax; ++j)
                ph(i, j) = (ls.psi(i, j) >= 0.0) ? phi_init[s] : 0.0;
        for (int j = 0; j <= jmax; ++j)          // west inlet row (i==0)
            ph(0, j) = (ls.psi(0, j) > 0.0) ? phi_inlet[s] : 0.0;
    }
    
    state_var.phi_old = state_var.phi;           // phi_old = phi
    state_var.U_prev = state_var.U;
    state_var.V_prev = state_var.V;
    state_var.P_prev = state_var.P;
    // P_cor_vec: pressure-correction solution vector from the SIMPLE/PISO pressure-velocity coupling scheme
    // boundary/ghost pressure nodes are excluded because they're set by BCs, not solved for
    state_var.P_cor_vec.assign(static_cast<size_t>(imax - 1) * (jmax - 1), 0.0);

    return state_var;
}





namespace {
// Same compact-array trick as Utilities/Coordinates.cpp's
// saveCoordinatesJson: json-c's pretty-printer puts one element per
// line, unreadable for U/V/P/phi/psi, so render each array as a single
// line via json-c (for consistent double formatting) and write the
// surrounding object structure by hand.
std::string arrayLine(const std::vector<double> &v) {
    json_object *arr = json_object_new_array();
    for (double x : v) json_object_array_add(arr, json_object_new_double(x));
    std::string s = json_object_to_json_string_ext(arr, JSON_C_TO_STRING_SPACED);
    json_object_put(arr);
    return s;
}
}  // namespace

void VariableNonDim::saveCurrentStateJson(const StateVar &stateVar, const LS &ls, double time,
                                           const std::string &path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("saveCurrentStateJson: failed to write " + path);
    }

    out << "{\n";
    out << "  \"time\": " << time << ",\n";
    out << "  \"U\": " << arrayLine(stateVar.U.data()) << ",\n";
    out << "  \"V\": " << arrayLine(stateVar.V.data()) << ",\n";
    out << "  \"P\": " << arrayLine(stateVar.P.data()) << ",\n";

    // phi: one species per Field2D (see StateVar's declaration comment)
    // -- written as a sub-object { "phi_0": [...], "phi_1": [...], ... },
    // one compact line per species keyed by species index.
    out << "  \"phi\": {\n";
    for (size_t s = 0; s < stateVar.phi.size(); ++s) {
        out << "    \"phi_" << s << "\": " << arrayLine(stateVar.phi[s].data());
        if (s + 1 < stateVar.phi.size()) out << ",";
        out << "\n";
    }
    out << "  },\n";

    out << "  \"psi\": " << arrayLine(ls.psi.data()) << "\n";
    out << "}\n";

    if (!out) {
        throw std::runtime_error("saveCurrentStateJson: failed to write " + path);
    }
}


namespace {
// Inverts a flat, row-major 3x3 matrix via the analytic adjugate
// formula. Hardcoded to 3x3 (not a general NxN inverse) since Pd below
// and A_0.2.json's matrix are both fixed at 3 species.
std::vector<double> invert3x3(const std::vector<double> &m) {
    double a = m[0], b = m[1], c = m[2];
    double d = m[3], e = m[4], f = m[5];
    double g = m[6], h = m[7], i = m[8];

    double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if (det == 0.0) {
        throw std::runtime_error("invert3x3: singular matrix");
    }
    double invDet = 1.0 / det;

    return {
        (e * i - f * h) * invDet, (c * h - b * i) * invDet, (b * f - c * e) * invDet,
        (f * g - d * i) * invDet, (a * i - c * g) * invDet, (c * d - a * f) * invDet,
        (d * h - e * g) * invDet, (b * g - a * h) * invDet, (a * e - b * d) * invDet,
    };
}
}  // namespace

void Variables::defineReactivity(const Domain &domain, const std::string &aMatJsonPath) {
    // Mirrors setUpVariablesNonDim.m:460-465. A_0.2.mat isn't readable
    // directly from C++, so this loads the one-time JSON conversion
    // instead (see scripts/convert_A_mat_to_json.m) -- same {"A": [[...
    // ]]} shape jsonencode(struct('A', A)) produces in MATLAB.
    json_object *root = json_object_from_file(aMatJsonPath.c_str());
    if (!root) {
        throw std::runtime_error("defineReactivity: failed to read " + aMatJsonPath);
    }
    json_object *aArr = nullptr;
    if (!json_object_object_get_ex(root, "A", &aArr) || !json_object_is_type(aArr, json_type_array) ||
        (int)json_object_array_length(aArr) != Np) {
        json_object_put(root);
        throw std::runtime_error("defineReactivity: " + aMatJsonPath + " missing a " +
                                  std::to_string(Np) + "x" + std::to_string(Np) + " \"A\" array");
    }
    std::vector<double> aRaw(static_cast<size_t>(Np) * Np, 0.0);
    for (int i = 0; i < Np; ++i) {
        json_object *row = json_object_array_get_idx(aArr, i);
        for (int j = 0; j < Np; ++j) {
            aRaw[static_cast<size_t>(i) * Np + j] =
                json_object_get_double(json_object_array_get_idx(row, j));
        }
    }
    json_object_put(root);

    // Pd: per-species diffusivity-like normalization, hardcoded in
    // MATLAB too (setUpVariablesNonDim.m:463) -- fixed at 3 species.
    static const double Pd[3] = {9.3, 0.793, 1.91};

    // A = inv(Pd)*A: Pd is diagonal, so this just divides row i by Pd[i].
    std::vector<double> a(aRaw.size());
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < Np; ++j) a[i * Np + j] = aRaw[i * Np + j] / Pd[i];

    // l = dxp(1,1)*sqrt(2) -- domain.dxp[0] is the same value MATLAB's
    // dxp(1,1) picks out (dxp is uniform across y even where MATLAB
    // broadcasts it to a 2D array).
    double l = domain.dxp.empty() ? 0.0 : domain.dxp[0] * std::sqrt(2.0);

    // inv_A = inv(3*eye(3) - 2*l*A)
    std::vector<double> m(a.size());
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < Np; ++j) m[i * Np + j] = (i == j ? 3.0 : 0.0) - 2.0 * l * a[i * Np + j];

    A = a;
    inv_A = invert3x3(m);
}

void Variables::defineLSvariables(const Domain &domain) {
   // smallest spacing
    double min_dxp = domain.dxp.empty()
                         ? 0.0
                         : *std::min_element(domain.dxp.begin(), domain.dxp.end());
    LSband = 10.0 * min_dxp;
    nLSupdate = 10;
    dtLS = dt;
    LSgamma = 5.0 * min_dxp;
    LSbeta = 3.0 * min_dxp;
    Pe_vel = Pe;
    dtau = 0.5 * min_dxp;  // fac = 0.5
}