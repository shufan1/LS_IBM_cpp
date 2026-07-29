#!/bin/bash
# Runs ls_ibm directly as a single process -- no srun/sbatch, no MPI
# launcher, just the plain binary. Only good for quick single-rank
# testing (main.cpp doesn't do any real solving yet, so this is cheap for
# now) -- run it from an sh_dev/salloc shell, never the login node, since
# it'll stop being cheap once the solver's actual time loop is ported.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

source ./config.sh

ml purge
ml math petsc/3.18.5 openblas/0.3.10
ml devel json-c/0.18

"$LS_IBM_CPP_DIR/bin/ls_ibm" "$PROJECT_DIR"
