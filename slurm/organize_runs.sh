#!/bin/bash
# Put this machine's runs into one consistent layout:
#
#   runs/<RUN_ID>                                   every NEW run (gravity off)
#   runs_dirichlet_20260924/runs/<RUN_ID>           superseded runs that lived in runs/
#   runs_dirichlet_20260924/<root>/<RUN_ID>         superseded runs from runs_C050, runs_aperture
#   runs_dirichlet_20260924/tests/<root>            the option-A hydrostatic test runs
#
# A run is NEW when its deck has GRAVITY 0.d0 0.d0 0.d0 (the hook sets it before
# PFLOTRAN starts); anything else in runs/ is superseded. Directories touched in
# the last 10 minutes are left alone, so a job being staged is never moved.
# Nothing is deleted or overwritten. Dry run unless --apply.
#
#   bash slurm/organize_runs.sh            # show the plan
#   bash slurm/organize_runs.sh --apply    # carry it out
set -uo pipefail
REV="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REV" || exit 1
APPLY=0; [[ "${1:-}" == "--apply" ]] && APPLY=1
ARCH=runs_dirichlet_20260924
ROOTS=(runs_C050 runs_aperture)                      # superseded roots of this revision
TESTS=(runs_hydrotest runs_hydrotest_aperture)       # option-A tests
RUNID='^[A-Z][A-Za-z0-9]*_.*p32_[0-9]+_s[0-9]+$'
SRC=(); DST=()          # initialised: an empty array is "unbound" under set -u in older bash
n_new=0; n_recent=0

plan () {   # $1 source, $2 destination
    if [[ -e "$2" ]]; then echo "  SKIP: $1 not moved, $2 already exists; compare them by hand"; return; fi
    SRC+=("$1"); DST+=("$2")
}
gravity_off () {
    grep -i "^[[:space:]]*GRAVITY" "$1/pflotran_co2.in" 2>/dev/null | grep -q "0.d0 0.d0 0.d0"
}

echo "revision: $REV"
# 1. superseded runs still in runs/
if [[ -d runs ]]; then
    for d in runs/*/; do
        d=${d%/}; n=$(basename "$d")
        [[ "$n" =~ $RUNID ]] || continue
        if gravity_off "$d"; then ((n_new++)); continue; fi
        if [[ -n "$(find "$d" -maxdepth 1 -mmin -10 -print -quit 2>/dev/null)" ]]; then ((n_recent++)); continue; fi
        plan "$d" "$ARCH/runs/$n"
    done
fi
# 2. runs moved earlier straight into the archive, without the runs/ level
if [[ -d "$ARCH" ]]; then
    for d in "$ARCH"/*/; do
        d=${d%/}; n=$(basename "$d")
        [[ "$n" =~ $RUNID ]] && plan "$d" "$ARCH/runs/$n"
    done
fi
# 3. superseded roots, and 4. the option-A tests
for r in "${ROOTS[@]}"; do [[ -d "$r" ]] && plan "$r" "$ARCH/$r"; done
for r in "${TESTS[@]}"; do [[ -d "$r" ]] && plan "$r" "$ARCH/tests/$r"; done

echo "  new runs left in runs/ (gravity off): $n_new"
(( n_recent )) && echo "  left alone, touched in the last 10 minutes: $n_recent"
if (( ${#SRC[@]} == 0 )); then echo "nothing to move"; exit 0; fi
declare -A count
for i in "${!SRC[@]}"; do k=$(dirname "${DST[$i]}"); count[$k]=$(( ${count[$k]:-0} + 1 )); done
for k in "${!count[@]}"; do echo "  ${count[$k]} -> $k/"; done
if (( ! APPLY )); then
    echo "dry run: nothing moved. Re-run with --apply to carry out the plan above."; exit 0
fi
mkdir -p "$ARCH/runs" "$ARCH/tests"
for i in "${!SRC[@]}"; do
    mv -n "${SRC[$i]}" "${DST[$i]}" || { echo "FATAL: could not move ${SRC[$i]}"; exit 1; }
done
echo "moved ${#SRC[@]}. runs/ now holds only new runs; superseded runs are under $ARCH/."
