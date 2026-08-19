#include "ControlVar.h"
#include <json-c/json.h>
#include <stdexcept>

ControlVar::ControlVar(const std::string &configPath) : ControlVar() {
    json_object *root = json_object_from_file(configPath.c_str());
    if (!root) {
        throw std::runtime_error("ControlVar: failed to read config file: " + configPath);
    }

    // All of this class's fields live under the top-level "ControlVar"
    // section of config.json (see the "geometry"/"physicalParameters"/
    // "boundaryConditions" sections for VariableNonDim's fields instead).
    // Missing section -> keep every default, same as a missing key.
    json_object *section;
    if (json_object_object_get_ex(root, "ControlVar", &section)) {
        json_object *v;
        if (json_object_object_get_ex(section, "output_folder", &v)) output_folder = json_object_get_string(v);
        if (json_object_object_get_ex(section, "time", &v)) time = json_object_get_double(v);
        if (json_object_object_get_ex(section, "timedt", &v)) timedt = json_object_get_double(v);
        if (json_object_object_get_ex(section, "savedat", &v)) savedat = json_object_get_int(v);
        if (json_object_object_get_ex(section, "rat", &v)) rat = json_object_get_int(v);
        if (json_object_object_get_ex(section, "tol", &v)) tol = json_object_get_double(v);
        if (json_object_object_get_ex(section, "f", &v)) f = json_object_get_int(v);
        if (json_object_object_get_ex(section, "tolbicg", &v)) tolbicg = json_object_get_double(v);
        if (json_object_object_get_ex(section, "maxit", &v)) maxit = json_object_get_int(v);
        if (json_object_object_get_ex(section, "disc_scheme_vel", &v)) disc_scheme_vel = json_object_get_int(v);
        if (json_object_object_get_ex(section, "flow_steady", &v)) flow_steady = json_object_get_boolean(v);
        if (json_object_object_get_ex(section, "imposePresBC", &v)) imposePresBC = json_object_get_string(v);
        if (json_object_object_get_ex(section, "transport_steady", &v)) transport_steady = json_object_get_boolean(v);
        if (json_object_object_get_ex(section, "noLStime", &v)) noLStime = json_object_get_int(v);
        if (json_object_object_get_ex(section, "iStart", &v)) iStart = json_object_get_int(v);
        if (json_object_object_get_ex(section, "tol_q", &v)) tol_q = json_object_get_double(v);
        if (json_object_object_get_ex(section, "tolbicg_c", &v)) tolbicg_c = json_object_get_double(v);
        if (json_object_object_get_ex(section, "maxit_c", &v)) maxit_c = json_object_get_int(v);
        if (json_object_object_get_ex(section, "ADRE", &v)) ADRE = json_object_get_boolean(v);

        // Linear solver choice, plain PETSc KSPType/PCType strings (see
        // ControlVar.h's field comments for why these aren't hardcoded).
        if (json_object_object_get_ex(section, "ksp_type_momentum", &v)) ksp_type_momentum = json_object_get_string(v);
        if (json_object_object_get_ex(section, "pc_type_momentum", &v)) pc_type_momentum = json_object_get_string(v);
        if (json_object_object_get_ex(section, "ksp_type_pressure", &v)) ksp_type_pressure = json_object_get_string(v);
        if (json_object_object_get_ex(section, "pc_type_pressure", &v)) pc_type_pressure = json_object_get_string(v);
    }

    json_object_put(root);
}
