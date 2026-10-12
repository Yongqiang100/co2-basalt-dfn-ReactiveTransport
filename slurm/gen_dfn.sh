#!/bin/bash
# =============================================================================
# DFN generation for the revision blocks. Cheap (~70 s/realisation at ~94k
# elements) and independent of the reactive-transport budget, so start this
# immediately regardless of how the benchmark comes out.
#
#   sbatch --array=0-49 --export=ALL,BLOCK=A 02_gen_dfn.sh   # expansion
#   sbatch --array=0-19 --export=ALL,BLOCK=B 02_gen_dfn.sh   # held-out
#   sbatch --array=0-15 --export=ALL,BLOCK=D 02_gen_dfn.sh   # domain size
#
# Block D needs 01_matrix_ext.patch applied first. Serialise the array with
# %4 (e.g. --array=0-49%4) if LaGriT contends for memory.
# =============================================================================
#SBATCH --job-name=gendfn
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
# Without an explicit --mem, Slurm applies DefMemPerNode = the WHOLE node
# (122 GB), so only one job can ever run and the array serialises. LaGriT is
# the memory consumer here; 6G covers the 250k-element x2.00 networks.
# Override:  sbatch --mem=<N>G ...
#SBATCH --time=04:00:00
#SBATCH --output=gendfn-%A_%a.out

set -euo pipefail
source /etc/profile.d/z99-local-modules.sh 2>/dev/null || true
# dfnWorks needs its own module; env.sh handles PFLOTRAN + python
module load dfnworks/2.7 || { echo 'FATAL: dfnworks module'; exit 1; }
source /opt/miniforge3/etc/profile.d/conda.sh
conda activate /opt/sw/conda/envs/dfnworks
export PMIX_MCA_psec=native

# ---- resolve paths, whether sbatch ran from revision/ or the deposit root ----
if   [[ -n "${DEPOSIT_ROOT:-}" ]]; then :
elif [[ -f "$SLURM_SUBMIT_DIR/prepare_dfn.py" ]];    then DEPOSIT_ROOT="$SLURM_SUBMIT_DIR"
elif [[ -f "$SLURM_SUBMIT_DIR/../prepare_dfn.py" ]]; then DEPOSIT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
else echo "FATAL: prepare_dfn.py not found from $SLURM_SUBMIT_DIR; set DEPOSIT_ROOT"; exit 1; fi
export DEPOSIT_ROOT
REVDIR="$DEPOSIT_ROOT/revision"
cd "$DEPOSIT_ROOT"
echo "DEPOSIT_ROOT=$DEPOSIT_ROOT"
# NOT sourcing env.sh here: it begins with `module purge`, which would strip
# the dfnworks module loaded above and force pydfnworks to fall back to
# ~/.dfnworksrc -> unbuilt DFNTrans -> gcc-15 compile failure. DFN generation
# needs no PFLOTRAN paths.
# prepare_dfn.py builds OUTPUT_ROOT from cwd at IMPORT time, so cwd must be the
# deposit root before extend_matrix.py imports it.
BLOCK=${BLOCK:?set BLOCK=A|B|D}

echo "Block $BLOCK task $SLURM_ARRAY_TASK_ID on $(hostname) | $(date)"
# NOTE: dfnWorks spawns its own mpirun -- do NOT wrap in srun
python "$REVDIR/src/extend_matrix.py" --block "$BLOCK" --index "$SLURM_ARRAY_TASK_ID"
echo "Finished $(date)"
