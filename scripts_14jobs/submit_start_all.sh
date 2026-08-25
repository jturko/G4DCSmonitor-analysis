#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00              # PER TASK now — tune to one base fuel's runtime
#SBATCH --array=0-13                   # 14 base fuels -> 14 independent tasks
#SBATCH --job-name=start_all           # fallback; overridden by run.sh
#SBATCH --output=logs/%x_%A_%a.out     # %A=array job id, %a=task index
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-user=j.turko@hzdr.de
#SBATCH --mail-type=END,FAIL           # NB: fires per-task -> up to 14 emails

set -eu
MACRO_DIR="${1:?usage: sbatch submit_start_all.sh macros/<dir>}"

cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Same 14-fuel list as start_all.sh, indexed by the array task id
BASE_FUELS=(0 1 2 3 5 6 7 8 13 14 15 22 23 32)
B="${BASE_FUELS[$SLURM_ARRAY_TASK_ID]}"

echo "[array ${SLURM_ARRAY_TASK_ID}] base fuel ${B} on $(hostname) at $(date)"
./start_one_base_fuel.sh "$B" "$MACRO_DIR"



###  ## version 1 
###  #!/bin/bash
###  #SBATCH --nodes=1
###  #SBATCH --ntasks=1
###  #SBATCH --cpus-per-task=128
###  #SBATCH --mem=128G
###  #SBATCH --time=5-00:00:00
###  #SBATCH --job-name=start_all          # fallback; overridden by run.sh
###  #SBATCH --output=logs/%x_%j.out       # %x=job name (has the dir), %j=jobid
###  #SBATCH --error=logs/%x_%j.err
###  #SBATCH --mail-user=j.turko@hzdr.de
###  #SBATCH --mail-type=BEGIN,END,FAIL
###  
###  set -eu
###  MACRO_DIR="${1:?usage: sbatch submit_start_all.sh macros/<dir>}"
###  
###  cd "$SLURM_SUBMIT_DIR"
###  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
###  ./start_all.sh "$MACRO_DIR"
