#!/bin/bash
# hpc01/hpc02 version of run_dirs.sh. Runs pre-staged run directories: the duration series (F), the
# shut-in ensemble (E) and the coupled porosity-permeability runs. Generalized
# from run_pulse.sh. Like it, this does not restage from pflotran_results/ or
# reapply apply_corrections.py: these decks already carry every correction and
# their own injection schedule, and restaging would discard the schedule or
# rescale the volume fractions a second time.
#
# Each line of LIST names one directory: a run ID under runs/ (as in
# pulse_runs.txt), or a path relative to revision/ (runs_aperture/A_p32_100_s1181).
#
#   ls -d runs/E_* | sed 's#^runs/##' > lists/E.txt
#   N=$(wc -l < lists/E.txt)
#   sbatch --array=0-$((N-1))%4 --ntasks=16 --mem=48G --export=ALL,LIST=lists/E.txt slurm/run_dirs_hpc01.sh
#SBATCH --job-name=rtD
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --ntasks=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=rtD-%A_%a.out

set -uo pipefail
# Project root, found the way run_rt.sh finds it, so this runs on hpc01 and hpc02:
# DEPOSIT_ROOT if given, else the submit directory or its parent holding prepare_dfn.py.
SUB="${SLURM_SUBMIT_DIR:-$PWD}"
if   [[ -n "${DEPOSIT_ROOT:-}" ]]; then :
elif [[ -f "$SUB/prepare_dfn.py" ]];    then DEPOSIT_ROOT="$SUB"
elif [[ -f "$SUB/../prepare_dfn.py" ]]; then DEPOSIT_ROOT="$(cd "$SUB/.." && pwd)"
else echo "FATAL: prepare_dfn.py not found from $SUB; set DEPOSIT_ROOT"; exit 1; fi
export DEPOSIT_ROOT
REVDIR="$DEPOSIT_ROOT/revision"
source "$REVDIR/slurm/env.sh"

LIST=${LIST:?set LIST=<file listing one run directory per line>}
[[ "$LIST" = /* ]] || LIST="$REVDIR/$LIST"
[[ -f "$LIST" ]] || { echo "FATAL: $LIST not found"; exit 1; }
RID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
[[ -n "$RID" ]] || { echo "FATAL: no entry at index $SLURM_ARRAY_TASK_ID in $LIST"; exit 1; }
if [[ "$RID" == */* ]]; then D="$REVDIR/$RID"; else D="$REVDIR/runs/$RID"; fi
cd "$D" || { echo "FATAL: $D missing"; exit 1; }

echo "=================================================="
echo " run     : $RID"
echo " list    : $LIST (index $SLURM_ARRAY_TASK_ID)"
echo " started : $(date)"
echo "=================================================="
grep -A9 "RATE LIST" pflotran_co2.in
echo " cells   : $(awk 'NR==1{print $2}' full_mesh.uge)"
grep -q -i "UPDATE_POROSITY" pflotran_co2.in && echo " porosity: coupled (UPDATE_POROSITY)" || echo " porosity: fixed"

# Gravity off (24 Sep 2026): in Richards mode gravity drives no flow, and with it
# on the uniform-pressure outflow face drives circulation. Every run, every block.
# The same step makes every database the deck names reachable from this directory.
"${PY:-python3}" "$REVDIR/src/fix_gravity_off.py" pflotran_co2.in || { echo "FATAL: gravity correction failed"; exit 1; }

# A continuing run (RESTART block, set by make_restart.py) keeps its files:
# the restart file is among pflotran_co2-*.h5, and the earlier part is saved as *_partN.
if grep -qi "^[[:space:]]*RESTART[[:space:]]*$" pflotran_co2.in; then
    echo " restart : continuing from $(awk '/^[[:space:]]*RESTART[[:space:]]*$/{r=1} r&&/FILENAME/{print $2; exit}' pflotran_co2.in)"
else
    rm -f pflotran_co2.h5 pflotran_co2-*.h5 pflotran_co2-mas.dat pflotran.log
fi

SECONDS=0
# Slurm records cpus-per-task twice; with --hint=nomultithread on 2-thread cores
# the two disagree (2 != 1) and srun refuses to start. Clear both, as on Setonix.
unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK
srun --mpi=pmix -n "${SLURM_NTASKS:-4}" pflotran -input_prefix pflotran_co2 \
     > pflotran.log 2>&1
RC=$?
echo " pflotran exit=$RC  wall=${SECONDS}s"
tail -6 pflotran.log
exit $RC
