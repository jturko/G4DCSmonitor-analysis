#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00               # now sized to ONE macro, not 5
#SBATCH --job-name=g4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-user=j.turko@hzdr.de
#SBATCH --mail-type=FAIL                 # was END,FAIL -> that's now 224 emails!

set -eu
MANIFEST="${1:?usage: sbatch submit_macro.sh <manifest>}"
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

MAC="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")"
echo "[array ${SLURM_ARRAY_TASK_ID}] $MAC on $(hostname) at $(date)"
G4DCSmonitor "$MAC"

