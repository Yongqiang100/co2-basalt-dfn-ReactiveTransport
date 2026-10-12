# =============================================================================
# Shared environment setup. Source this from every Slurm script.
#
#   source "$REVDIR/slurm/env.sh"          # after DEPOSIT_ROOT/REVDIR are set
#
# WHY THIS EXISTS
# ---------------
# run_pflotran.py:37-57 decides it is on Setonix if ANY of these hold:
#     /software/projects/pawsey1284 exists
#     PAWSEY_PROJECT in environ
#     SLURM_JOB_ID in environ        <-- true inside ANY Slurm job
#     HOSTNAME starts with setonix / nid
#
# So every sbatch job on hpc01 is misdetected as Setonix and resolves
# /software/projects/pawsey1284/ychen6/pflotran/... which does not exist here.
# It works interactively and fails under Slurm for that reason alone.
#
# The 'local' branch would not help either: it searches ~/pflotran, /usr/local
# and /opt/pflotran, never /opt/sw/pflotran.
#
# Both branches wrap the paths in os.environ.get(), so exporting the three
# variables below fixes it WITHOUT modifying the deposit. The misdetection also
# leaves launcher="srun", which is what we want under Slurm anyway.
# =============================================================================

# ---- modules ---------------------------------------------------------------
module purge
module load pflotran/6.0 || { echo "FATAL: cannot load pflotran/6.0"; exit 1; }

# ---- PFLOTRAN executable ---------------------------------------------------
if [[ -z "${PFLOTRAN_EXE:-}" || ! -f "${PFLOTRAN_EXE:-}" ]]; then
    PFLOTRAN_EXE="$(command -v pflotran || true)"
fi
[[ -f "$PFLOTRAN_EXE" ]] || {
    echo "FATAL: pflotran executable not found."
    echo "       Set PFLOTRAN_EXE explicitly and resubmit."
    exit 1; }
export PFLOTRAN_EXE

# ---- thermodynamic databases ----------------------------------------------
# The source tree is <root>/src/pflotran/pflotran with <root>/database alongside,
# so derive <root> from the executable path.
_pf_root="${PFLOTRAN_EXE%/src/pflotran/pflotran}"
# hpc01: /opt/sw/pflotran/src/pflotran/pflotran, with database/ alongside.
# hpc02: flat /opt/sw/pflotran/pflotran and NO database directory, so the
# md5-identical copy shipped in $DEPOSIT_ROOT/database is used instead.
_pf_dir="$(dirname "$PFLOTRAN_EXE")"
_db_dirs=("$_pf_root/database" "$_pf_root/../database" "$_pf_dir/database"
          "$_pf_dir/../database" "/opt/sw/pflotran/database"
          "$DEPOSIT_ROOT/database")

_find_db () {   # $1 = filename, $2 = existing value
    local want="$1" have="$2" d
    if [[ -n "$have" && -f "$have" ]]; then echo "$have"; return 0; fi
    for d in "${_db_dirs[@]}"; do
        [[ -f "$d/$want" ]] && { echo "$d/$want"; return 0; }
    done
    # last resort: bounded search under the pflotran tree
    d="$(find "$_pf_root" -maxdepth 4 -name "$want" -type f 2>/dev/null | head -1)"
    [[ -n "$d" ]] && { echo "$d"; return 0; }
    return 1
}

PFLOTRAN_DB="$(_find_db hanford.dat "${PFLOTRAN_DB:-}")" || {
    echo "FATAL: hanford.dat not found under $_pf_root"
    echo "       Locate it and export PFLOTRAN_DB, e.g.:"
    echo "         find /opt/sw -name hanford.dat 2>/dev/null"
    exit 1; }
CO2_DB="$(_find_db co2_sw.dat "${CO2_DB:-}")" || {
    echo "FATAL: co2_sw.dat not found under $_pf_root"
    echo "         find /opt/sw -name co2_sw.dat 2>/dev/null"
    exit 1; }
export PFLOTRAN_DB CO2_DB

# ---- python ----------------------------------------------------------------
# Pin the interpreter rather than inheriting whichever conda env happened to be
# active at submission time; --export=ALL means that varies between jobs.
if [[ -n "${REVISION_PYTHON:-}" ]]; then
    PY="$REVISION_PYTHON"
elif [[ -x /opt/sw/conda/envs/dfnworks/bin/python3 ]]; then
    PY=/opt/sw/conda/envs/dfnworks/bin/python3
else
    PY="$(command -v python3)"
fi
export PY

echo "--- environment ------------------------------------------------------"
echo "  PFLOTRAN_EXE : $PFLOTRAN_EXE"
echo "  PFLOTRAN_DB  : $PFLOTRAN_DB"
echo "  CO2_DB       : $CO2_DB"
echo "  python       : $PY ($("$PY" -c 'import sys;print(sys.version.split()[0])' 2>/dev/null || echo unusable))"
echo "  modules      : ${LOADEDMODULES:-none}"
echo "----------------------------------------------------------------------"
