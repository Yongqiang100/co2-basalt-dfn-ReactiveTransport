#!/bin/bash
# =========================================================================
# Runs the shortened-injection decks built by build_pulse_runs.py.
#
# These are launched outside run_rt.sh on purpose. That script restages the
# deck from pflotran_results/ and reapplies apply_corrections.py on every
# run, either of which would discard the appended stopping time or rescale
# the rate. The decks here already carry the baseline rate and the stop, so
# they are run as they stand.
#
#   python3 src/build_pulse_runs.py          # writes pulse_runs.txt
#   sbatch --array=0-$(($(wc -l < pulse_runs.txt)-1)) slurm/run_pulse.sh
# =========================================================================
#SBATCH --job-name=rtP
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --ntasks=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=rtP-%A_%a.out

set -uo pipefail

export DEPOSIT_ROOT="$HOME/co2-basalt"
REVDIR="$DEPOSIT_ROOT/revision"
source "$REVDIR/slurm/env.sh"

cd "$REVDIR" || exit 1
LIST="$REVDIR/pulse_runs.txt"
[[ -f "$LIST" ]] || { echo "FATAL: $LIST not found"; exit 1; }

RID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
[[ -n "$RID" ]] || { echo "FATAL: no entry at index $SLURM_ARRAY_TASK_ID"; exit 1; }

D="$REVDIR/runs/$RID"
cd "$D" || { echo "FATAL: $D missing"; exit 1; }

echo "=================================================="
echo " run     : $RID"
echo " dir     : $D"
echo " started : $(date)"
echo "=================================================="
grep -A9 "RATE LIST" pflotran_co2.in
NCELLS=$(awk 'NR==1{print $2}' full_mesh.uge)
echo " cells   : $NCELLS"
echo

rm -f pflotran_co2.h5 pflotran_co2-*.h5 pflotran_co2-mas.dat pflotran.log

SECONDS=0
srun --mpi=pmix -n "${SLURM_NTASKS:-4}" pflotran -input_prefix pflotran_co2 \
     > pflotran.log 2>&1
RC=$?
echo " pflotran exit=$RC  wall=${SECONDS}s"
tail -6 pflotran.log

echo
echo "--- validation ---"
"$PY" "$REVDIR/src/validate_run.py" pflotran_co2.h5
echo " -> validate exit=$?"
echo " finished: $(date)"
