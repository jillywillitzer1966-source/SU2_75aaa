#!/bin/bash
#SBATCH -J deto_test
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -t 02:00:00
#SBATCH -p cpu
#SBATCH -o SU2_output_%j.log
#SBATCH -e SU2_error_%j.log

set -euo pipefail

CFG_PATH="${CFG_PATH:-${1:-}}"
if [[ -z "$CFG_PATH" ]]; then
  echo "Usage: sbatch sub_deto.sh /absolute/or/relative/path/to/case.cfg" >&2
  exit 2
fi

if [[ ! -f "$CFG_PATH" ]]; then
  echo "Config file not found: $CFG_PATH" >&2
  exit 2
fi

CFG_ABS=$(readlink -f "$CFG_PATH")
CASE_DIR=$(dirname "$CFG_ABS")
CFG_NAME=$(basename "$CFG_ABS")

cd "$CASE_DIR"
source /home/jmyang/detonationFoam/env_su2_deto.sh

NTASKS=${SLURM_NTASKS:-2}

echo "[sub_deto] job_id=${SLURM_JOB_ID:-unknown} host=$(hostname) pwd=$PWD"
echo "[sub_deto] cfg=$CFG_ABS ntasks=$NTASKS"
echo "[sub_deto] mpirun=$(which mpirun)"
echo "[sub_deto] SU2_CFD=$(which SU2_CFD)"

srun -n "$NTASKS" SU2_CFD "$CFG_NAME"
