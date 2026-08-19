#!/bin/bash
#SBATCH -p serc
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J ls_ibm_smoketest
#SBATCH -o logs/smoketest.%j.out
#SBATCH -e logs/smoketest.%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# Build already done via ./build.sh -- this job only runs the binary.
ml math petsc/3.18.5
ml cmake/3.24.2

srun --mem-per-cpu=2GB ./bin/ls_ibm_smoketest
