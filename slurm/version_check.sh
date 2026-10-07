#!/bin/bash
# =============================================================================
# dfnWorks version audit + seed-reproducibility test
#
# Answers three questions, in order:
#   1. What version is ACTUALLY installed on hpc01?
#   2. What version does the archive say produced the published 25?
#   3. Does the installed version reproduce an archived realisation from its seed?
#
# Question 3 is the one that matters. Run from ~/co2-basalt.
#
#   sbatch 03_version_check.sh          # or run interactively, it's ~2 min
# =============================================================================
#SBATCH --job-name=vercheck
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --output=vercheck-%j.out

set -uo pipefail
CASE=${CASE:-p32_100_s383}
SRC=/opt/sw/dfnworks/dfnWorks

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


source /etc/profile.d/z99-local-modules.sh 2>/dev/null || true
module purge
module load dfnworks/2.7
source /opt/miniforge3/etc/profile.d/conda.sh
conda activate /opt/sw/conda/envs/dfnworks
export PMIX_MCA_psec=native

echo "############ 1. WHAT IS INSTALLED ############"
# The setup doc clones with no tag checkout, so the module label "2.7" is an
# admin label, not an upstream version. The commit hash is the real identity.
if [[ -d "$SRC/.git" ]]; then
    echo "  commit    : $(git -C "$SRC" rev-parse HEAD 2>/dev/null)"
    echo "  date      : $(git -C "$SRC" log -1 --format=%ad --date=iso 2>/dev/null)"
    echo "  describe  : $(git -C "$SRC" describe --tags --always --dirty 2>/dev/null)"
    echo "  branch    : $(git -C "$SRC" rev-parse --abbrev-ref HEAD 2>/dev/null)"
    echo "  --- does it contain the 2026-08-01 seeding change? ---"
    for c in 210a910 fb5f0c2 fce27ef; do
        if git -C "$SRC" merge-base --is-ancestor $c HEAD 2>/dev/null; then
            echo "    $c  PRESENT   <-- changes seed -> network mapping"
        else
            echo "    $c  absent"
        fi
    done
else
    echo "  WARNING: $SRC is not a git checkout; version not recoverable from source"
fi
python -c "import pydfnworks,os;print('  pydfnworks:',getattr(pydfnworks,'__version__','no __version__ attr'),'at',os.path.dirname(pydfnworks.__file__))" 2>/dev/null

echo
echo "############ 2. WHAT THE ARCHIVE RECORDS ############"
python3 - << 'PY'
import json, glob, os
f = sorted(glob.glob("dfn_library/*/dfn_summary.json"))
if not f:
    print("  no dfn_summary.json found — is dfn_library/ unpacked?")
else:
    print(f"  {len(f)} summaries present")
    d = json.load(open(f[0]))
    print(f"  keys in {os.path.basename(os.path.dirname(f[0]))}: {sorted(d)}")
    hits = {k: v for k, v in d.items()
            if any(t in k.lower() for t in ("version", "dfnworks", "commit", "git", "seed"))}
    print("  version/seed fields:", hits if hits else "NONE RECORDED")
PY

echo
echo "############ 3. SEED-REPRODUCIBILITY TEST ############"
REF="dfn_library/$CASE"
[[ -d "$REF" ]] || { echo "FATAL: $REF missing"; exit 1; }

SEED=$(echo "$CASE" | sed 's/.*_s//')
MULT=$(echo "$CASE" | sed 's/p32_0*//; s/_s.*//' | awk '{printf "%.2f", $1/100}')
echo "  case=$CASE  seed=$SEED  P32 multiplier=$MULT"

# Regenerate into a scratch library so the archive is never touched
mkdir -p vercheck && cd vercheck
cp "$DEPOSIT_ROOT/prepare_dfn.py" . 2>/dev/null
python3 - "$MULT" "$SEED" << 'PY'
import sys, os
sys.path.insert(0, ".")
import prepare_dfn as P
P.OUTPUT_ROOT = os.path.join(os.getcwd(), "dfn_library")
os.makedirs(P.OUTPUT_ROOT, exist_ok=True)
mult, seed = float(sys.argv[1]), int(sys.argv[2])
P.generate_single_dfn(P.make_config(mult), "REGEN", seed=seed)
PY
cd ..

NEW="vercheck/dfn_library/REGEN"
[[ -d "$NEW" ]] || { echo "  regeneration FAILED — see log above"; exit 1; }

echo
printf "  %-26s %14s %14s  %s\n" "quantity" "archived" "regenerated" "match"
cmp_num () {  # label  refval  newval
    local m="DIFFER"; [[ "$2" == "$3" ]] && m="ok"
    printf "  %-26s %14s %14s  %s\n" "$1" "$2" "$3" "$m"
}
# .uge line 1 is "CELLS <n>" -- field 2, not 1. Comparing field 1 compared
# the literal string CELLS against itself and could never fail.
cmp_num "mesh cells (CELLS field)" \
    "$(awk 'NR==1{print $2}' "$REF/full_mesh.uge" 2>/dev/null)" \
    "$(awk 'NR==1{print $2}' "$NEW/full_mesh.uge" 2>/dev/null)"
cmp_num "uge total lines" \
    "$(wc -l < "$REF/full_mesh.uge" 2>/dev/null)" \
    "$(wc -l < "$NEW/full_mesh.uge" 2>/dev/null)"
for ex in "$REF"/*.ex; do
    b=$(basename "$ex")
    cmp_num "  ${b%.ex} (.ex lines)" \
        "$(wc -l < "$ex" 2>/dev/null)" \
        "$(wc -l < "$NEW/$b" 2>/dev/null || echo MISSING)"
done
# Raw dfnGen geometry, if the archive kept it — this is the sharpest test
for g in radii_Final.dat translations.dat normal_vectors.dat; do
    if [[ -f "$REF/$g" && -f "$NEW/$g" ]]; then
        cmp_num "$g (md5)" \
            "$(md5sum < "$REF/$g" | cut -c1-12)" \
            "$(md5sum < "$NEW/$g" | cut -c1-12)"
    fi
done

echo
echo "############ VERDICT ############"
REF_N=$(awk 'NR==1{print $2}' "$REF/full_mesh.uge")
NEW_N=$(awk 'NR==1{print $2}' "$NEW/full_mesh.uge")
REF_L=$(wc -l < "$REF/full_mesh.uge"); NEW_L=$(wc -l < "$NEW/full_mesh.uge")

# Geometry is authoritative, not cell count. Upstream commit fce27ef records
# that Poisson-disc MESHING drew from the unseeded global NumPy generator, so a
# build lacking set_seed() can yield a different cell count from an IDENTICAL
# fracture network. Prefer radii/translations/normals where archived; the .uge
# line count is the next-best structural fingerprint.
GEO_FILES=0; GEO_MATCH=0
for g in radii_Final.dat translations.dat normal_vectors.dat; do
    if [[ -f "$REF/$g" && -f "$NEW/$g" ]]; then
        GEO_FILES=$((GEO_FILES+1))
        [[ "$(md5sum < "$REF/$g")" == "$(md5sum < "$NEW/$g")" ]] && GEO_MATCH=$((GEO_MATCH+1))
    fi
done

# Refuse to compare empty values. An unreadable file makes both sides empty,
# and "" == "" would pass vacuously -- the same failure as comparing the literal
# word CELLS against itself.
for v in REF_N NEW_N REF_L NEW_L; do
    [[ -n "${!v}" ]] || { echo "  FATAL: $v is empty -- could not read full_mesh.uge"; exit 1; }
done
echo "  CELLS:     archived=$REF_N  regenerated=$NEW_N"
echo "  uge lines: archived=$REF_L  regenerated=$NEW_L"
if [[ "$GEO_FILES" -gt 0 ]]; then
    echo "  dfnGen geometry files matching: $GEO_MATCH / $GEO_FILES"
    if [[ "$GEO_MATCH" -eq "$GEO_FILES" ]]; then
        echo "  -> FRACTURE NETWORK IS SEED-REPRODUCIBLE (geometry identical)."
        [[ "$REF_N" != "$NEW_N" ]] && echo "     Cell counts differ only because meshing is unseeded in this build."
    else
        echo "  -> GEOMETRY DIFFERS. Do NOT pool new realisations with the published 25"
        echo "     until distributional equivalence is tested (README section 6)."
    fi
elif [[ "$REF_N" == "$NEW_N" && "$REF_L" == "$NEW_L" ]]; then
    echo "  No dfnGen geometry archived, but CELLS and .uge line count both match"
    echo "  exactly -> the mesh is structurally identical. Strong evidence of"
    echo "  seed-reproducibility."
else
    echo "  No dfnGen geometry archived and the mesh differs."
    echo "  -> INDETERMINATE. If this build lacks set_seed() (see section 1),"
    echo "     unseeded meshing alone can explain a cell-count difference."
    echo "     Compare fracture counts and radii distributions before concluding."
fi
echo
echo "  Record in the methods: pydfnworks __version__ from section 1"
echo "  (the manuscript's \"v2.10\" corresponds to pydfnworks 2.10.x)."
