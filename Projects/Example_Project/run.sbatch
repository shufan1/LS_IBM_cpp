#!/bin/bash
#SBATCH -p serc
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH -J ls_ibm_example_project
#SBATCH -o sherlock_log/run.%j.out
#SBATCH -e sherlock_log/run.%j.err

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Guard against stale Slurm env vars leaking in from an outer interactive
# allocation (e.g. sh_dev) that don't match this job's own --cpus-per-task
# -- srun fatals if SLURM_CPUS_PER_TASK disagrees with SLURM_TRES_PER_TASK.
unset SLURM_CPUS_PER_TASK

source ./config.sh

ml purge
ml math petsc/3.18.5 openblas/0.3.10
ml devel json-c/0.18

# main.cpp takes the project directory as its one positional argument
# and reads $PROJECT_DIR/config.json from it.
srun "$LS_IBM_CPP_DIR/bin/ls_ibm" "$PROJECT_DIR"
