cmake_minimum_required(VERSION 3.18)
project(LS_IBM_cpp LANGUAGES CXX) # set PROJECT_SOURCE_DIR

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

find_package(MPI REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)

find_package(json-c CONFIG REQUIRED)

# LAPACK: used by LSmirPointsBQ.cpp's dgetrf_/dgetri_-based matrix inverse
# (LU decomposition with partial pivoting) -- the same algorithm MATLAB's
# own inv()/A^-1 uses internally, so this is what keeps the ghost-cell
# mirror-point interpolation weights numerically close to the MATLAB
# reference. PETSc on this cluster is already linked against OpenBLAS
# (which bundles LAPACK), so this doesn't add a new runtime dependency --
# `ml math openblas/0.3.10` (same version PETSc itself uses) just needs to
# be loaded so CMake can find it at configure time.
find_package(LAPACK REQUIRED)

# Environment sanity check (MPI + PETSc + multi-node) -- keep this around to
# re-verify the toolchain later; it is not the solver itself.
add_executable(ls_ibm_smoketest src/smoketest.cpp)
target_include_directories(ls_ibm_smoketest PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(ls_ibm_smoketest PRIVATE MPI::MPI_CXX PkgConfig::PETSC)

# Isolated repro for the "Unable to find requested KSP type bicgstab" crash --
# tests only PetscInitialize + KSPSetType(bicgstab) on a tiny hand-built
# matrix, independent of any of this project's own matrix-assembly code.
add_executable(ksp_smoketest src/ksp_smoketest.cpp)
target_include_directories(ksp_smoketest PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(ksp_smoketest PRIVATE MPI::MPI_CXX PkgConfig::PETSC)

# Isolated unit test for LSeqSolve -- given the same (psi, phi) input as
# debug_compare/test_LSeqSolve.m (a JSON file, not a long simulation
# run), does this produce the same output? See test_LSeqSolve.cpp's own
# header comment. Needs everything LSeqSolve.cpp's own call chain
# touches, same list as ls_ibm below minus the flow/transport/IBM
# modules this test never calls (VariableNonDim.cpp still references
# LSInitialize.cpp's computeLSInitialize() inside getLS(), even though
# this test doesn't call getLS() -- the symbol still needs linking).
add_executable(test_ls_eq_solve
    src/test_LSeqSolve.cpp
    src/VariableNonDim.cpp
    src/Utilities/Coordinates.cpp
    src/Utilities/Interp.cpp
    src/Utilities/Timer.cpp
    src/SolveLS/LSInitialize.cpp
    src/SolveLS/LSnormals.cpp
    src/SolveLS/LSVelocityExtrapolation.cpp
    src/SolveLS/WENODerivative.cpp
    src/SolveLS/solveHJEq.cpp
    src/SolveLS/godunovGradientNorm.cpp
    src/SolveLS/LSreinitialization.cpp
    src/SolveLS/LSeqSolve.cpp
)
target_include_directories(test_ls_eq_solve PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(test_ls_eq_solve PRIVATE MPI::MPI_CXX PkgConfig::PETSC json-c::json-c)


# The actual solver, built as modules get ported from MATLAB (mirrors the
# beta_AD_correctu/<Module>/<File>.m layout). Add each new .cpp here by hand.
add_executable(ls_ibm
    src/main.cpp
    src/ControlVar.cpp
    src/VariableNonDim.cpp
    src/Utilities/Coordinates.cpp
    src/Utilities/DiffuFlux.cpp
    src/Utilities/ConvFlux.cpp
    src/Utilities/Timer.cpp
    src/Utilities/Interp.cpp
    src/Utilities/ConvergenceResiduals.cpp
    src/SolveLS/LSInitialize.cpp
    src/SolveLS/LSnormals.cpp
    src/IBM/LSPointIdent.cpp
    src/IBM/LSmirPointsBQ.cpp
    src/IBM/LSIBMcoeffs.cpp
    src/SolveUVP/SolveUVP.cpp
    src/SolveUVP/COEFFU.cpp
    src/SolveUVP/COEFFV.cpp
    src/SolveUVP/COEFFP.cpp
    src/SolveUVP/RHSP.cpp
    src/SolveUVP/RHSP_PISO.cpp
    src/SolveTransport/CalculateQG.cpp
    src/SolveTransport/CoeffPhiADRE.cpp
    src/SolveTransport/RhsPhiADRE.cpp
    src/SolveTransport/UpdateA1gPhi.cpp
    src/SolveTransport/UpdateGhostReactionPhi.cpp
    src/SolveTransport/SolveTransportADRE.cpp
    src/SolveLS/LSVelocityExtrapolation.cpp
    src/SolveLS/WENODerivative.cpp
    src/SolveLS/solveHJEq.cpp
    src/SolveLS/godunovGradientNorm.cpp
    src/SolveLS/LSreinitialization.cpp
    src/SolveLS/LSeqSolve.cpp
    src/Models/flowSStransportUS.cpp
)
target_include_directories(ls_ibm PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(ls_ibm PRIVATE MPI::MPI_CXX PkgConfig::PETSC json-c::json-c LAPACK::LAPACK)

# Example of how the source list grows as modules get ported -- uncomment
# and extend as real files land in these folders:
# add_executable(ls_ibm
#     src/main.cpp
#     src/SolveTransport/SolveTransportADRE.cpp
#     src/IBM/LSIBMcoeffs.cpp
# )
# if two .cpp files in the same add_executable() both define main(), you get a hard linker error — multiple  definition of 'main' 

# 1. target_include_directories(transport_solver PRIVATE ${MFEM_SOURCE_DIR})
# "When compiling this target's .cpp files, where should the compiler look for headers?" This
#  becomes a -I${MFEM_SOURCE_DIR} flag on every compile command for transport_solver. It's what lets
#  #include "mfem.hpp" actually resolve to a real file instead of erroring "file not found."


#  2. target_link_directories(transport_solver PRIVATE ${MFEM_SOURCE_DIR})
#  Answers: "When linking the final executable, what directories should the linker search for library
#  files?" This becomes a -L${MFEM_SOURCE_DIR} flag. On its own it doesn't link anything — it just adds a
#  directory to the search path.

#  3. target_link_libraries(transport_solver PRIVATE libmfem.a)
#  Answers: "Which actual library file(s) should get linked in?" Since line 2 already told the linker where
#  to search, line 3 only needs the bare filename — the linker checks every directory from
#  target_link_directories (plus system defaults) until it finds a libmfem.a.