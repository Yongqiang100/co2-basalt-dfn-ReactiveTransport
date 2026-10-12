#!/bin/bash
# Build pre-staged run directories for run_dirs.sh, then list them.
#
#   --variant feedback (default)  coupled (evolving aperture):  runs/<BLOCK>_feedback__<network>
#   --variant fixed               fixed porosity, no coupling:  runs/<BLOCK>_<network>
#                                 (needs --networks: the same name run_rt.sh uses, so
#                                 networks are named explicitly and never enumerated)
#   --name NAME                   directory prefix instead of BLOCK (networks still come
#                                 from BLOCK):  runs/<NAME>_<network> or <NAME>__<network>
#                                 when NAME has a tag (F_10d), and <NAME>_feedback__<network>
#   --shutin YEARS                stop the injection at YEARS (add_shutin.py)
#   --offset REL                  stop offset (default 1e-6, the E decks; the old F decks use 1e-4)
#
# Each network is staged exactly as run_rt.sh stages a plain run (template deck,
# mesh, region and property files from pflotran_results/ and dfn_library/),
# then:  apply_corrections.py (RESCALE_VF=1, so rescaled once)
#        fix_gravity_off.py   (gravity off, databases)
#        add_shutin.py        (with --shutin: the injection stops at the given time)
#        add_coupling.py      (feedback only: porosity-permeability coupling, 5 lines)
# and checked: VF 0.5000, gravity off, coupling present (feedback) or absent
# (fixed), SCALED_MASS_RATE, real mesh.
# No simulation is run. Existing directories are never touched.
#
#   bash slurm/build_dirs.sh --skip 43 --extra p32_200_s117 --apply                  # coupled Block A
#   bash slurm/build_dirs.sh --variant fixed --networks p32_200_s117 --apply         # fixed s117 only
#   bash slurm/build_dirs.sh --name E --shutin 10 --skip 43 --extra p32_200_s117 --apply                  # E, coupled
#   bash slurm/build_dirs.sh --name E --shutin 10 --variant fixed --skip 43 --extra p32_200_s117 --apply  # E, fixed
# Without --apply, only the plan is shown.
set -uo pipefail
REV="$(cd "$(dirname "$0")/.." && pwd)"; ROOT="$(dirname "$REV")"
PY="${PY:-$(command -v python3)}"
BLOCK=A; SKIP=""; EXTRA=""; NETWORKS=""; VARIANT=feedback; NAME=""; SHUTIN=""; OFFSET=1e-6; APPLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --block) BLOCK="$2"; shift 2 ;;
        --skip)  SKIP="$2";  shift 2 ;;
        --extra) EXTRA="$2"; shift 2 ;;
        --networks) NETWORKS="$2"; shift 2 ;;
        --variant)  VARIANT="$2";  shift 2 ;;
        --name)     NAME="$2";     shift 2 ;;
        --shutin)   SHUTIN="$2";   shift 2 ;;
        --offset)   OFFSET="$2";   shift 2 ;;
        --apply) APPLY=1;    shift ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done
cd "$REV" || exit 1
NAME="${NAME:-$BLOCK}"
case "$VARIANT" in
    feedback) SUFFIX="_feedback__" ;;
    fixed)    SUFFIX="_"; [[ "$NAME" == *_* ]] && SUFFIX="__"
              if [[ "$NAME" == "$BLOCK" && -z "$NETWORKS" ]]; then
                  echo "FATAL: --variant fixed with the block's own name needs --networks (run_rt.sh uses the same names)"; exit 1; fi ;;
    *) echo "FATAL: unknown --variant $VARIANT (feedback or fixed)"; exit 1 ;;
esac
nets=()
if [[ -n "$NETWORKS" ]]; then
    for e in ${NETWORKS//,/ }; do nets+=("$e"); done
else
    N=$("$PY" src/blocks.py --block "$BLOCK" --count) || { echo "FATAL: blocks.py --count failed"; exit 1; }
    for ((i = 0; i < N; i++)); do
        [[ ",$SKIP," == *",$i,"* ]] && continue
        nets+=("$("$PY" src/blocks.py --block "$BLOCK" --index "$i" | cut -f2)")
    done
    for e in ${EXTRA//,/ }; do nets+=("$e"); done
fi

if [[ "$VARIANT" == fixed && "$NAME" == "$BLOCK" ]]; then   # a network run_rt.sh runs for this block is never built here
    N=$("$PY" src/blocks.py --block "$BLOCK" --count) || { echo "FATAL: blocks.py --count failed"; exit 1; }
    own=" "; for ((i = 0; i < N; i++)); do own+="$("$PY" src/blocks.py --block "$BLOCK" --index "$i" | cut -f2) "; done
    for n in "${nets[@]}"; do
        [[ "$own" == *" $n "* ]] && { echo "FATAL: $n is in block $BLOCK of blocks.py; run_rt.sh runs it as ${BLOCK}_$n"; exit 1; }
    done
fi

LIST="lists/${NAME}_${VARIANT}.txt"
vfsum () { awk '$1 ~ /^(Anorthite|Albite|Diopside|Forsterite|Fayalite|Enstatite)$/ &&
                $2 ~ /^[0-9]*\.?[0-9]+([dDeE][+-]?[0-9]+)?$/ {gsub(/[dD]/,"e",$2); t+=$2} END {printf "%.4f", t}' "$1"; }

check () {   # $1 = run directory; prints what is wrong, nothing if it is ready
    local rd="$1" in d why=""
    in=$(cd "$rd" 2>/dev/null && ls *.in 2>/dev/null | head -1)
    [[ -n "$in" ]] || { echo " no-deck"; return; }
    d="$rd/$in"
    [[ "$(vfsum "$d")" == "0.5000" ]] || why+=" VF=$(vfsum "$d")"
    grep -qi "^[[:space:]]*GRAVITY[[:space:]]\+0.d0 0.d0 0.d0" "$d" || why+=" gravity-on"
    if [[ "$VARIANT" == feedback ]]; then
        grep -qi "^[[:space:]]*UPDATE_POROSITY" "$d"  || why+=" not-coupled"
        grep -qi "^[[:space:]]*MINIMUM_POROSITY" "$d" || why+=" no-porosity-floor"
    else
        grep -qi "^[[:space:]]*UPDATE_POROSITY" "$d" && why+=" coupled-but-should-be-fixed"
    fi
    grep -q  "SCALED_MASS_RATE" "$d"                              || why+=" rate-uncorrected"
    if [[ -n "$SHUTIN" ]]; then
        "$PY" "$REV/src/add_shutin.py" --check --offset "$OFFSET" "$d" "$SHUTIN" || why+=" no-stop-at-${SHUTIN}y"
    fi
    [[ -f "$rd/full_mesh.uge" && ! -L "$rd/full_mesh.uge" ]]      || why+=" no-mesh"
    echo "$why"
}

ok=(); todo=0; skipped=0; missing=0; failed=0; invalid=0
for n in "${nets[@]}"; do
    id="${NAME}${SUFFIX}$n"; rd="runs/$id"
    st="$ROOT/pflotran_results/$n"; lib="$ROOT/dfn_library/$n"
    if [[ ! -f "$st/pflotran_co2.in" || ! -f "$lib/full_mesh.uge" ]]; then
        echo "  MISSING inputs: $n"; ((missing++)); continue; fi
    if [[ -e "$rd" ]]; then
        why=$(check "$rd")
        if [[ -z "$why" ]]; then echo "  already built and valid: $rd"; ok+=("$id"); ((skipped++))
        else echo "  EXISTS BUT NOT VALID, left alone and not listed: $rd:$why"; ((invalid++)); fi
        continue
    fi
    ((todo++))
    (( APPLY )) || continue
    mkdir -p "$rd"
    for f in "$st"/*.in "$st"/*.uge "$st"/*.inp "$st"/*.ex "$st"/*.dat \
             "$lib"/*.uge "$lib"/*.inp "$lib"/*.ex "$lib"/dfn_properties.h5; do
        [[ -e "$f" && ! -e "$rd/$(basename "$f")" ]] && cp -L "$f" "$rd"/
    done
    in=$(cd "$rd" && ls *.in 2>/dev/null | head -1)
    log="$rd/build.log"
    if ! ( cd "$rd" && RESCALE_VF=1 "$PY" "$REV/src/apply_corrections.py" "$in" \
                    && "$PY" "$REV/src/fix_gravity_off.py" "$in" \
                    && { [[ -z "$SHUTIN" ]] || "$PY" "$REV/src/add_shutin.py" --offset "$OFFSET" "$in" "$SHUTIN"; } \
                    && { [[ "$VARIANT" != feedback ]] || "$PY" "$REV/src/add_coupling.py" "$in"; } ) > "$log" 2>&1; then
        echo "  BUILD FAILED: $n (see below)"; tail -3 "$log" | sed 's/^/      /'
        rm -rf "$rd"; ((failed++)); continue
    fi
    why=$(check "$rd")
    if [[ -n "$why" ]]; then
        echo "  CHECK FAILED: $n:$why"; rm -rf "$rd"; ((failed++)); continue; fi
    ok+=("$id")
done

echo "$NAME from block $BLOCK ($VARIANT${SHUTIN:+, stop at ${SHUTIN} y}): ${#nets[@]} networks | to build $todo | already built $skipped | existing but invalid $invalid | missing inputs $missing | failed $failed"
if (( ! APPLY )); then echo "dry run: nothing built. Re-run with --apply."; exit 0; fi
mkdir -p lists
if (( ${#ok[@]} )); then printf "%s\n" "${ok[@]}" > "$LIST"; else : > "$LIST"; fi
echo "list: $LIST ($(wc -l < "$LIST") directories, each checked: VF 0.5000, gravity off, $VARIANT)"
(( failed == 0 && invalid == 0 ))
