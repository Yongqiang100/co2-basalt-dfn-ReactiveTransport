#!/bin/bash
# =============================================================================
# STEP ZERO — throughput calibration on hpc01
#
# Reproduces ONE already-published realisation and measures wall time, peak
# memory and parallel scaling. Every scoping decision for the revision run set
# depends on the number this produces. Do not queue the run set before it runs.
#
#   sbatch --export=ALL,CASE=p32_100_s383,RANKS=4  00_benchmark.sh
#   sbatch --export=ALL,CASE=p32_100_s383,RANKS=8  00_benchmark.sh
#   sbatch --export=ALL,CASE=p32_100_s383,RANKS=16 00_benchmark.sh
#   sbatch --export=ALL,CASE=p32_100_s383,RANKS=24 00_benchmark.sh
#
# Run all four. The scaling curve decides how to pack the node: PFLOTRAN GIRT
# on ~80k cells will stop scaling well below 24 ranks, and once it does you get
# more throughput from N concurrent small jobs than one wide job.
# =============================================================================
#SBATCH --job-name=bench
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=12G
# Explicit --mem is REQUIRED for concurrency: the default is the whole node,
# which would serialise the run set. PFLOTRAN memory scales with CELLS, not
# ranks, so 12G covers up to ~253k cells. Block D (~650k) needs ~40G.
# Override:  sbatch --mem=40G ...
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#
# --hint=nomultithread on the ALLOCATION is what makes one task == one physical
# core. Without it a "CPU" is a thread (ThreadsPerCore=2) and srun's own
# --hint=nomultithread asks for twice what was allocated:
#     srun: error: Unable to create step: More processors requested than permitted
#
# --ntasks above is only a default. OVERRIDE IT ON THE COMMAND LINE and keep it
# equal to RANKS:
#     sbatch --ntasks=24 --export=ALL,CASE=p32_100_s383,RANKS=24 slurm/benchmark.sh
#SBATCH --output=bench-%j.out

set -uo pipefail
source /etc/profile.d/z99-local-modules.sh 2>/dev/null || true
# modules, PFLOTRAN paths and python are all set here -- see slurm/env.sh for
# why the PFLOTRAN_* exports are mandatory on this machine.
_ENV_PENDING=1

CASE=${CASE:-p32_100_s383}
# Trust the ALLOCATION over the RANKS variable; they diverge the moment someone
# forgets --ntasks, and a silent mismatch is how you get a meaningless timing.
RANKS=${SLURM_NTASKS:-${RANKS:-8}}
if [[ -n "${RANKS_REQ:-}" && "$RANKS_REQ" != "$RANKS" ]]; then
    echo "WARNING: RANKS=$RANKS_REQ requested but allocation has $RANKS task(s)."
    echo "         Using $RANKS. Pass --ntasks=$RANKS_REQ to sbatch to fix."
fi

# ---- stage via the published code path -------------------------------------
# ---- resolve paths, whether sbatch ran from revision/ or the deposit root ----
if   [[ -n "${DEPOSIT_ROOT:-}" ]]; then :
elif [[ -f "$SLURM_SUBMIT_DIR/prepare_dfn.py" ]];    then DEPOSIT_ROOT="$SLURM_SUBMIT_DIR"
elif [[ -f "$SLURM_SUBMIT_DIR/../prepare_dfn.py" ]]; then DEPOSIT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
else echo "FATAL: prepare_dfn.py not found from $SLURM_SUBMIT_DIR; set DEPOSIT_ROOT"; exit 1; fi
export DEPOSIT_ROOT
REVDIR="$DEPOSIT_ROOT/revision"
cd "$DEPOSIT_ROOT"
echo "DEPOSIT_ROOT=$DEPOSIT_ROOT"
source "$REVDIR/slurm/env.sh"
[[ -d "dfn_library/$CASE" ]] || { echo "FATAL: dfn_library/$CASE not found"; exit 1; }

# --write_only writes the deck and copies mesh, databases and *.ex, then returns
"$PY" run_pflotran.py --dfn "$CASE" --write_only
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "FATAL: run_pflotran.py --write_only exited $rc"
    echo "       (imports? missing database? check the traceback above)"
    exit 1
fi

# Each rank count gets its OWN directory. Four concurrent sweeps sharing
# pflotran_results/$CASE would overwrite each other's HDF5, .out and timing
# files and every measurement would be worthless.
STAGED="$DEPOSIT_ROOT/pflotran_results/$CASE"
RUNDIR="$DEPOSIT_ROOT/bench/${CASE}_n${RANKS}"
mkdir -p "$RUNDIR"
for f in "$STAGED"/*.in "$STAGED"/*.uge "$STAGED"/*.inp "$STAGED"/*.ex "$STAGED"/*.dat; do
    [[ -e "$f" ]] && { [[ -e "$RUNDIR/$(basename "$f")" ]] || cp "$f" "$RUNDIR"/; }
done
cd "$RUNDIR" || { echo "FATAL: $RUNDIR missing"; exit 1; }
INPUT=$(ls *.in | head -1)
[[ -n "$INPUT" ]] || { echo "FATAL: no deck emitted"; exit 1; }
NCELLS=$(awk 'NR==1{print $2}' full_mesh.uge)   # line 1 is "CELLS <n>"
[[ -f boundary_right_e.ex ]] && echo "outflow boundary: present" \
                            || echo "outflow boundary: ABSENT -> CLOSED SYSTEM case"

echo "=================================================="
echo " case       : $CASE"
echo " cells      : $NCELLS"
echo " ranks      : $RANKS  (SLURM_NTASKS=${SLURM_NTASKS:-unset}, \
CPUs on node=${SLURM_CPUS_ON_NODE:-?})"
echo " deck       : $INPUT"
echo " started    : $(date)"
echo "=================================================="

# ---- run under /usr/bin/time to capture peak RSS ----------------------------
# NOTE on memory: `/usr/bin/time -v` reports the RSS of the srun launcher, not
# the peak across MPI ranks, so it under-reports badly. sacct MaxRSS is the
# number to plan concurrency from.
SECONDS=0
# no --hint here: it is inherited from the allocation (see header)
srun --mpi=pmix -n "$RANKS" \
    pflotran -input_prefix "${INPUT%.in}" > pflotran.log 2>&1 || true
WALL_S=$SECONDS
tail -20 pflotran.log

# ---- extract the numbers that matter ---------------------------------------
sleep 5   # let slurmd flush accounting
MAXRSS=$( { sacct -j "$SLURM_JOB_ID" --noheader --format=MaxRSS 2>/dev/null \
           | tr -d ' ' | grep -v '^$' | sort -h | tail -1; } || true )
if [[ -z "$MAXRSS" ]]; then
    # No slurmdbd -> sacct is empty. Read the cgroup peak instead (v2 then v1).
    for _c in /sys/fs/cgroup/memory.peak \
              /sys/fs/cgroup/memory/memory.max_usage_in_bytes; do
        if [[ -r "$_c" ]]; then
            MAXRSS="$(awk '{printf "%.1fG(cgroup)", $1/1073741824}' "$_c")"
            break
        fi
    done
fi
[[ -n "$MAXRSS" ]] || MAXRSS="n/a(no sacct, no cgroup)"
LASTYR=$( { "$PY" "$REVDIR/src/validate_run.py" "${INPUT%.in}.h5" 2>/dev/null \
            | awk '/time groups/{print $NF" y"}'; } || true )
[[ -n "$LASTYR" ]] || LASTYR=$( { grep -oiE "time *=? *[0-9.eEdD+-]+ *(y|yr|d)" pflotran.log \
                                  2>/dev/null | tail -1; } || true )

echo
echo "=== BENCHMARK RESULT ==============================" | tee -a "$DEPOSIT_ROOT/benchmark_results.txt"
printf "%-18s cells=%-8s ranks=%-3s wall=%-8s maxrss=%-10s last=%-14s %s\n" \
    "$CASE" "$NCELLS" "$RANKS" "${WALL_S}s" "${MAXRSS:-n/a}" "${LASTYR:-UNKNOWN}" "$(date +%F)" \
    | tee -a "$DEPOSIT_ROOT/benchmark_results.txt"
echo "  cells/rank = $((NCELLS/RANKS))"
echo
echo "--- validation -------------------------------------------------------"
"$PY" "$REVDIR/src/validate_run.py" "${INPUT%.in}.h5" --case "$CASE" || \
    echo "  VALIDATOR REPORTED A PROBLEM -- see above"
echo
echo "  VALIDATE BEFORE TRUSTING THIS NUMBER:"
echo "   1. last reported time above must be ~50 y. A deck that stalled at"
echo "      t=3 yr can still exit 0 and looks fast."
echo "   2. carbonate must match the published value for this case"
echo "      (p32_100_s383 -> carb_per_cell_final = 2.202e-5)."
