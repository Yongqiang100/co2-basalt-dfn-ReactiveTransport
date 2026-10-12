#!/bin/bash
# Regenerate all manuscript and letter figures from the corrected coupled runs (revision R2).
# Each figure comes from its original script, so layouts are unchanged.
# Usage:  bash src/make_figures.sh              (from the revision directory)
#         bash src/make_figures.sh --skip-3d    # skip slow 3D figures (1, 6)
set -e

RUNS=runs_gravityoff                    # corrected runs: gravity omitted, uniform outflow pressure
LIBRARY=../dfn_library
OUTDIR=figures_R2                       # review here before copying into figures/
# Co-location pair: evolving-aperture Block A at P32 x 1.00, highest and lowest carbonate
# per cell at 50 years (net of the seed), the measure _pick_pair uses.
read HIGH LOW < <(python - "$RUNS" << 'PYEOF'
import glob, os, re, sys
import numpy as np, h5py
runs = sys.argv[1]; res = []
for d in sorted(glob.glob(os.path.join(runs, "A_feedback__p32_100_s*"))):
    h5 = os.path.join(d, "pflotran_co2.h5")
    if not os.path.isfile(h5):
        continue
    with h5py.File(h5, "r") as f:
        tg = sorted((float(k.split("Time")[1].split()[0]), k) for k in f if "Time" in k)
        if not tg or tg[-1][0] < 49.99:
            continue
        g = f[tg[-1][1]]; n = g["pH"].shape[0]
        c = sum(np.clip(np.asarray(g[k][:], float).ravel() - 1e-6, 0, None).sum()
                for m in ("Calcite", "Magnesite", "Siderite", "Dawsonite") for k in g if k.startswith(m + " VF"))
    res.append((c / n, os.path.basename(d)))
res.sort()
for v, nm in res:
    print(f"   {nm}: {v:.3e} per cell", file=sys.stderr)
print(res[-1][1], res[0][1])
PYEOF
)
echo "  co-location pair (evolving aperture): high $HIGH, low $LOW"

echo "============================================================"
echo "  Generating all figures -> ${OUTDIR}/"
echo "============================================================"

# --- Manuscript figures (generate_figures.py): 1-9, and 16 (dissolution, carbonate and pH over time) ---
echo ""
echo "--- Manuscript figures (generate_figures.py) ---"
python src/generate_figures.py \
    --results-dir "$RUNS" \
    --dfn-dir "$LIBRARY" \
    --output-dir "$OUTDIR" \
    --only 1 2 3 4 5 6 7 8 9 16 \
    "$@"

# --- Study design (fig_study_design.py) ---
echo ""
echo "--- Study design ---"
python src/fig_study_design.py --out "$OUTDIR"

# --- Co-location pair (fig_cations.py): two realisations, 2 x 6 (dissolution in panel c) ---
# --no-age drops the Age x Flux panel, as in the original figure.
echo ""
echo "--- Co-location: pair ---"
python src/fig_cations.py \
    --runs "$RUNS" \
    --high "$HIGH" \
    --low "$LOW" \
    --library "$LIBRARY" \
    --no-age --time 5 --time-end 7 \
    --out "${OUTDIR}/fig_cations_pair.pdf"

# --- Ensemble statistics figures (ensemble_statistics_figures.py): variability and sensitivity ---
echo ""
echo "--- Ensemble statistics figures ---"
python src/prepare_figure_inputs.py
python src/ensemble_statistics_figures.py \
    --runs "$RUNS" \
    --out "$OUTDIR" \
    --blocka results/figure_input_blockA.csv \
    --sens results/figure_input_sensitivity.csv

echo ""
echo "============================================================"
echo "  Done. All figures in ${OUTDIR}/"
echo "============================================================"
