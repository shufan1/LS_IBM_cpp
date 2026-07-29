#!/bin/bash
# Run this from an interactive shell (sh_dev / salloc), not the login node.
# Configure is cheap/idempotent; the build step only recompiles changed files.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

ml purge
ml math petsc/3.18.5 cmake/3.24.2 openblas/0.3.10
ml devel json-c/0.18

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=1
cmake --build build -j "${SLURM_CPUS_ON_NODE:-$(nproc)}"
