#!/bin/bash
#SBATCH --job-name=rtF
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --ntasks=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=rtF-%A_%a.out
set -uo pipefail
export DEPOSIT_ROOT="$HOME/co2-basalt"
source "$DEPOSIT_ROOT/revision/slurm/env.sh"
CASE=p32_150_s1481
TAGS=(30d 10d 1d)
D="$DEPOSIT_ROOT/revision/runs/F_${TAGS[$SLURM_ARRAY_TASK_ID]}__$CASE"
cd "$D" || exit 1
echo "$D"; grep -A8 "RATE LIST" pflotran_co2.in
rm -f pflotran_co2.h5 pflotran_co2-mas.dat pflotran.log
srun --mpi=pmix -n "${SLURM_NTASKS:-4}" pflotran -input_prefix pflotran_co2 \
     > pflotran.log 2>&1
echo "exit=$?"; tail -5 pflotran.log
