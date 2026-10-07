#!/bin/bash
# run_local.sh -- Run PFLOTRAN jobs locally on WSL2 (GNU/Linux)
#
# Runs DFNs either in parallel (background) or sequentially.
# Uses mpirun instead of srun.
#
# Usage:
#   1. Generate all input files:
#        python run_pflotran.py --dfn all --write_only
#
#   2. Run all DFNs sequentially (safe, low memory):
#        bash run_local.sh
#
#   3. Run all DFNs in parallel (fast, needs RAM):
#        bash run_local.sh --parallel
#
#   4. Run a single DFN:
#        bash run_local.sh --dfn p32_075_s259
#
#   Options:
#        bash run_local.sh --nprocs 4        # MPI ranks (default: auto-detect)
#        bash run_local.sh --max-parallel 3   # Max concurrent jobs in parallel mode
#        bash run_local.sh --dry-run          # Print commands without running
#        bash run_local.sh --post-only        # Only run postprocessing

set -euo pipefail

# =============================================================
# DEFAULTS -- edit to match your WSL2 setup
# =============================================================
TOTAL_CORES=$(nproc 2>/dev/null || echo 4)
NPROCS=$((TOTAL_CORES > 1 ? TOTAL_CORES - 1 : 1))

PFLOTRAN_EXE="${PFLOTRAN_EXE:-$(which pflotran 2>/dev/null || echo "pflotran")}"
MPIRUN="${MPIRUN:-$(which mpirun 2>/dev/null || echo "mpirun")}"

RESULTS_ROOT="$(pwd)/pflotran_results"
PARALLEL=false
MAX_PARALLEL=2
DRY_RUN=false
POST_ONLY=false
SINGLE_DFN=""

# =============================================================
# Parse arguments
# =============================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)        DRY_RUN=true; shift ;;
        --parallel)       PARALLEL=true; shift ;;
        --post-only)      POST_ONLY=true; shift ;;
        --nprocs)         NPROCS="$2"; shift 2 ;;
        --max-parallel)   MAX_PARALLEL="$2"; shift 2 ;;
        --dfn)            SINGLE_DFN="$2"; shift 2 ;;
        --results)        RESULTS_ROOT="$2"; shift 2 ;;
        --exe)            PFLOTRAN_EXE="$2"; shift 2 ;;
        --mpirun)         MPIRUN="$2"; shift 2 ;;
        -h|--help)
            cat <<EOF
Usage: bash run_local.sh [OPTIONS]

Options:
  --nprocs N         MPI ranks per job (default: $NPROCS, auto-detected)
  --parallel         Run DFNs concurrently in background
  --max-parallel N   Max concurrent jobs in parallel mode (default: $MAX_PARALLEL)
  --dfn NAME         Run only this DFN
  --post-only        Skip PFLOTRAN, only run Python postprocessing
  --dry-run          Print commands without executing
  --results PATH     Results directory (default: ./pflotran_results)
  --exe PATH         PFLOTRAN executable (default: auto-detect)
  --mpirun PATH      mpirun executable (default: auto-detect)
  -h, --help         Show this help

Examples:
  bash run_local.sh                                     # Sequential, all cores
  bash run_local.sh --nprocs 16                         # Sequential, 16 ranks
  bash run_local.sh --parallel --nprocs 16 --max-parallel 4  # 4 DFNs x 16 ranks
  bash run_local.sh --dfn p32_075_s259 --nprocs 8       # Single DFN
  bash run_local.sh --post-only                         # Postprocess only
EOF
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# =============================================================
# Validate
# =============================================================
if [[ ! -d "$RESULTS_ROOT" ]]; then
    echo "ERROR: Results directory not found: $RESULTS_ROOT"
    echo "Run 'python run_pflotran.py --dfn all --write_only' first."
    exit 1
fi

if [[ "$POST_ONLY" == false ]]; then
    if ! command -v "$MPIRUN" &>/dev/null; then
        echo "ERROR: mpirun not found."
        echo "Install OpenMPI:  sudo apt install openmpi-bin libopenmpi-dev"
        echo "  or MPICH:       sudo apt install mpich"
        exit 1
    fi

    if ! command -v "$PFLOTRAN_EXE" &>/dev/null && [[ ! -x "$PFLOTRAN_EXE" ]]; then
        echo "ERROR: PFLOTRAN executable not found: $PFLOTRAN_EXE"
        echo "Set with: --exe /path/to/pflotran"
        echo "  or:     export PFLOTRAN_EXE=/path/to/pflotran"
        exit 1
    fi
fi

# =============================================================
# Build DFN list
# =============================================================
DFN_DIRS=()
if [[ -n "$SINGLE_DFN" ]]; then
    DFN_DIRS=("$RESULTS_ROOT/$SINGLE_DFN")
    if [[ ! -d "${DFN_DIRS[0]}" ]]; then
        echo "ERROR: DFN directory not found: ${DFN_DIRS[0]}"
        exit 1
    fi
else
    for d in "$RESULTS_ROOT"/*/; do
        [[ -d "$d" ]] && DFN_DIRS+=("$d")
    done
fi

# =============================================================
# Header
# =============================================================
echo "=============================================="
echo "PFLOTRAN Local Run -- WSL2"
echo "=============================================="
echo "PFLOTRAN:      $PFLOTRAN_EXE"
echo "mpirun:        $MPIRUN"
echo "MPI ranks:     $NPROCS"
echo "CPU cores:     $TOTAL_CORES"
echo "Mode:          $(if $PARALLEL; then echo "parallel (max $MAX_PARALLEL)"; else echo "sequential"; fi)"
echo "Post-only:     $POST_ONLY"
echo "Dry run:       $DRY_RUN"
echo "DFNs:          ${#DFN_DIRS[@]}"
echo "Results:       $RESULTS_ROOT"
echo "=============================================="
echo ""

# =============================================================
# Run a single DFN
# =============================================================
run_dfn() {
    local DFN_DIR="$1"
    local DFN_NAME
    DFN_NAME=$(basename "$DFN_DIR")
    local INPUT_FILE="$DFN_DIR/pflotran_co2.in"
    local OUT_FILE="$DFN_DIR/pflotran_co2.out"
    local LOG_FILE="$DFN_DIR/run.log"

    # Check input exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "[$DFN_NAME] SKIP: no pflotran_co2.in"
        return 0
    fi

    # Skip if already completed (check both possible completion markers)
    if [[ -f "$OUT_FILE" ]]; then
        if grep -q "Wall Clock Time" "$OUT_FILE" 2>/dev/null && \
           ! grep -q "Simulation failed" "$OUT_FILE" 2>/dev/null; then
            echo "[$DFN_NAME] SKIP: already completed"
            return 0
        fi
    fi

    if [[ "$POST_ONLY" == true ]]; then
        echo "[$DFN_NAME] Postprocessing..."
        if [[ "$DRY_RUN" == false ]]; then
            python run_pflotran.py --dfn "$DFN_NAME" --post_only 2>&1 | tee "$DFN_DIR/postprocess.log"
        fi
        return 0
    fi

    # Clean stale output
    rm -f "$OUT_FILE"
    rm -f "$DFN_DIR"/pflotran_co2*.h5

    echo "[$DFN_NAME] Starting ($NPROCS ranks)..."
    local START_TIME
    START_TIME=$(date +%s)

    local CMD="cd $DFN_DIR && $MPIRUN -n $NPROCS $PFLOTRAN_EXE -pflotranin pflotran_co2.in"

    if [[ "$DRY_RUN" == true ]]; then
        echo "[$DFN_NAME] WOULD RUN: $CMD"
        return 0
    fi

    # Run PFLOTRAN, capture output
    (
        cd "$DFN_DIR"
        $MPIRUN -n $NPROCS $PFLOTRAN_EXE -pflotranin pflotran_co2.in \
            > "$LOG_FILE" 2>&1
        local EXIT_CODE=$?

        local END_TIME
        END_TIME=$(date +%s)
        local ELAPSED=$(( END_TIME - START_TIME ))
        local HOURS=$(( ELAPSED / 3600 ))
        local MINS=$(( (ELAPSED % 3600) / 60 ))
        local SECS=$(( ELAPSED % 60 ))

        if [[ $EXIT_CODE -eq 0 ]]; then
            echo "[$DFN_NAME] DONE (${HOURS}h ${MINS}m ${SECS}s)"
            # Run postprocessing
            cd - > /dev/null
            python run_pflotran.py --dfn "$DFN_NAME" --post_only 2>&1 \
                | tee "$DFN_DIR/postprocess.log"
        else
            echo "[$DFN_NAME] FAILED exit=$EXIT_CODE (${HOURS}h ${MINS}m ${SECS}s)"
            echo "[$DFN_NAME] Last 5 lines of log:"
            tail -5 "$LOG_FILE" | sed "s/^/  /"
        fi
    ) &

    # In sequential mode, wait for each job
    if [[ "$PARALLEL" == false ]]; then
        wait
    fi
}

# =============================================================
# Main loop
# =============================================================
RUNNING=0

for DFN_DIR in "${DFN_DIRS[@]}"; do
    run_dfn "$DFN_DIR"

    if [[ "$PARALLEL" == true && "$DRY_RUN" == false ]]; then
        RUNNING=$((RUNNING + 1))

        # Throttle: wait if we hit the max concurrent limit
        if [[ $RUNNING -ge $MAX_PARALLEL ]]; then
            wait -n 2>/dev/null || wait
            RUNNING=$((RUNNING - 1))
        fi
    fi
done

# Wait for all remaining background jobs
if [[ "$PARALLEL" == true && "$DRY_RUN" == false ]]; then
    echo ""
    echo "Waiting for remaining jobs to finish..."
    wait
fi

echo ""
echo "=============================================="
echo "All done. Summary:"
echo "=============================================="

COMPLETED=0
FAILED=0
PENDING=0

for DFN_DIR in "${DFN_DIRS[@]}"; do
    DFN_NAME=$(basename "$DFN_DIR")
    OUT_FILE="$DFN_DIR/pflotran_co2.out"

    if [[ -f "$OUT_FILE" ]]; then
        if grep -q "Simulation failed" "$OUT_FILE" 2>/dev/null; then
            STEP=$(grep "^ Step" "$OUT_FILE" | tail -1 | awk '{print $2, $3, $4}')
            echo "  FAILED: $DFN_NAME (last: $STEP)"
            FAILED=$((FAILED + 1))
        elif grep -q "Wall Clock Time" "$OUT_FILE" 2>/dev/null; then
            STEP=$(grep "^ Step" "$OUT_FILE" | tail -1 | awk '{print $2, $3, $4}')
            WTIME=$(grep "Wall Clock Time" "$OUT_FILE" | awk '{print $7}')
            echo "  DONE: $DFN_NAME (${WTIME} min, last: $STEP)"
            COMPLETED=$((COMPLETED + 1))
        else
            STEP=$(grep "^ Step" "$OUT_FILE" | tail -1 | awk '{print $2, $3, $4}')
            echo "  RUNNING: $DFN_NAME (last: $STEP)"
            PENDING=$((PENDING + 1))
        fi
    else
        PENDING=$((PENDING + 1))
    fi
done

echo ""
echo "  Completed: $COMPLETED"
echo "  Failed:    $FAILED"
echo "  Pending:   $PENDING"
echo "=============================================="