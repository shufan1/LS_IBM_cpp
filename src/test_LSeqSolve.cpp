// Isolated unit test for LSeqSolve -- the C++ counterpart of
// debug_compare/test_LSeqSolve.m. Given the SAME (psi, phi) input (read
// from a JSON file test_LSeqSolve.m already produced -- see its own
// comment for the flat row-major convention), calls LSeqSolve() ONCE
// and dumps the result, so debug_compare/compare_LSeqSolve.m can diff
// it directly against MATLAB's own output. No long simulation run
// involved -- this tests the level-set step itself, in isolation from
// any accumulated drift.
//
// Usage: test_ls_eq_solve <project_directory> <input_json_path> <output_json_path>
//   <project_directory>: same convention as main.cpp -- must contain
//     config.json (for Domain/Variables -- see the assumption this
//     matches Output_Siqin's own setup, documented in
//     test_LSeqSolve.m's own header comment).
//   <input_json_path>: the LSeqSolve_input_<tstep>.json file
//     test_LSeqSolve.m wrote.
//   <output_json_path>: where to write this run's LS output.
#include <petsc.h>
#include <mpi.h>
#include <json-c/json.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "VariableNonDim.h"
#include "SolveLS/LSnormals.h"
#include "SolveLS/LSeqSolve.h"

namespace {

// Reads {time, psi: [...], phi: {phi_0: [...], ...}} -- test_LSeqSolve.m's
// own JSON schema. psi/phi arrays are flat, Field2D::data()-convention
// (index = i*ny+j), same as saveCurrentStateJson()'s own arrays.
struct LoadedState {
    double time = 0.0;
    std::vector<double> psi;
    std::vector<std::vector<double>> phi;  // [species][flat index]
};

std::vector<double> readDoubleArray(json_object *arr) {
    std::vector<double> out;
    int n = json_object_array_length(arr);
    out.reserve(n);
    for (int i = 0; i < n; ++i) out.push_back(json_object_get_double(json_object_array_get_idx(arr, i)));
    return out;
}

LoadedState loadStateJson(const std::string &path) {
    json_object *root = json_object_from_file(path.c_str());
    if (!root) throw std::runtime_error("test_LSeqSolve: failed to read " + path);

    LoadedState s;
    json_object *v;
    if (json_object_object_get_ex(root, "time", &v)) s.time = json_object_get_double(v);
    if (json_object_object_get_ex(root, "psi", &v)) s.psi = readDoubleArray(v);

    json_object *phiObj;
    if (json_object_object_get_ex(root, "phi", &phiObj)) {
        // Species keys are phi_0, phi_1, ... -- read in that order, not
        // json-c's own (insertion) iteration order, so species index is
        // never accidentally scrambled.
        for (int i_s = 0;; ++i_s) {
            json_object *sp;
            std::string key = "phi_" + std::to_string(i_s);
            if (!json_object_object_get_ex(phiObj, key.c_str(), &sp)) break;
            s.phi.push_back(readDoubleArray(sp));
        }
    }
    json_object_put(root);
    return s;
}

// Unflattens a Field2D::data()-convention array (index = i*ny+j) into a
// Field2D of the given shape -- the read-side counterpart of Field2D's
// own flat storage, not a generic reshape.
Field2D toField2D(const std::vector<double> &flat, int nx, int ny) {
    Field2D f(nx, ny);
    if (static_cast<int>(flat.size()) != nx * ny) {
        throw std::runtime_error("test_LSeqSolve: flat array size " + std::to_string(flat.size()) +
                                  " != nx*ny (" + std::to_string(nx) + "*" + std::to_string(ny) + ")");
    }
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) f(i, j) = flat[static_cast<size_t>(i) * ny + j];
    return f;
}

void saveOutputJson(const std::string &path, const LS &ls) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("test_LSeqSolve: failed to write " + path);

    auto arrayLine = [](const std::vector<double> &v) {
        std::ostringstream s;
        s << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s << ",";
            s << v[i];
        }
        s << "]";
        return s.str();
    };

    out << "{\n";
    out << "  \"psi\": " << arrayLine(ls.psi.data()) << ",\n";
    out << "  \"nx\": " << arrayLine(ls.nx.data()) << ",\n";
    out << "  \"ny\": " << arrayLine(ls.ny.data()) << ",\n";
    out << "  \"u\": " << arrayLine(ls.u.data()) << ",\n";
    out << "  \"v\": " << arrayLine(ls.v.data()) << "\n";
    out << "}\n";
}

}  // namespace

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 4) {
        if (rank == 0)
            fprintf(stderr, "usage: %s <project_directory> <input_json_path> <output_json_path>\n", argv[0]);
        PetscFinalize();
        return 1;
    }
    std::string projectDir = argv[1];
    std::string inputPath = argv[2];
    std::string outputPath = argv[3];
    std::string configPath = projectDir + "/config.json";

    // Domain/Variables from this project's own config -- NOT loaded from
    // Output_Siqin's coordinate.mat. Assumes config.json describes the
    // same "grain" test case Output_Siqin was generated from; see
    // test_LSeqSolve.m's own sanity check (imax/jmax/xp/yp comparison)
    // for the MATLAB-side half of that verification. No equivalent
    // assertion here yet -- if this ever gets used against a differently
    // -configured project directory, the mismatch would show up as
    // nonsense output, not a clear error.
    VariableNonDim varSetup(configPath);
    Domain domain = varSetup.getDomain();
    Variables variables = varSetup.getVariables(domain);
    variables.defineReactivity(domain, "A_0.2.json");
    variables.defineLSvariables(domain);
    variables.verbose = false;

    LoadedState s = loadStateJson(inputPath);
    if (rank == 0) printf("loaded %s, time=%.6f, Np=%zu\n", inputPath.c_str(), s.time, s.phi.size());

    const int nx = domain.imax + 1;
    const int ny = domain.jmax + 1;

    LS ls;
    ls.caseId = 1;  // 1 = circle/grain, matches test_LSeqSolve.m's LS.case=1
    ls.psi = toField2D(s.psi, nx, ny);
    LSNormals normals = computeLSNormals(ls.psi, domain);
    ls.nx = normals.nx;
    ls.ny = normals.ny;

    StateVar stateVar;
    stateVar.phi.resize(s.phi.size());
    for (size_t i_s = 0; i_s < s.phi.size(); ++i_s) stateVar.phi[i_s] = toField2D(s.phi[i_s], nx, ny);

    LSeqSolve(ls, stateVar, domain, variables);

    if (rank == 0) {
        saveOutputJson(outputPath, ls);
        printf("saved %s\n", outputPath.c_str());
    }

    PetscFinalize();
    return 0;
}
