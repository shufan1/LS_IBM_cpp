#!/bin/bash
#SBATCH -p normal
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH -J ls_ibm_example_project
#SBATCH -o sherlock_log/run.%j.out
#SBATCH -e sherlock_log/run.%j.err

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p sherlock_log

# Guard against stale Slurm env vars leaking in from an outer interactive
# allocation (e.g. sh_dev) that don't match this job's own --cpus-per-task
# -- srun fatals if SLURM_CPUS_PER_TASK disagrees with SLURM_TRES_PER_TASK.
unset SLURM_CPUS_PER_TASK

source ./config.sh

ml purge
ml math petsc/3.18.5
ml devel json-c/0.18

# Temporary diagnostic for the libjson-c.so.5 "cannot open shared object
# file" error -- remove once resolved. Confirms whether LD_LIBRARY_PATH is
# actually set on the node(s) srun launches on, not just in this script's
# own shell.
echo "=== diagnostic: batch-script node ($(hostname)) ==="
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
ldd "$LS_IBM_CPP_DIR/bin/ls_ibm" | grep -i json
echo "=== diagnostic: on every allocated node, via srun ==="
srun --ntasks-per-node=1 --nodes="$SLURM_JOB_NUM_NODES" bash -c \
    'echo "$(hostname): LD_LIBRARY_PATH=$LD_LIBRARY_PATH"'

# main.cpp takes the project directory as its one positional argument
# and reads $PROJECT_DIR/config.json from it.
srun "$LS_IBM_CPP_DIR/bin/ls_ibm" "$PROJECT_DIR"
