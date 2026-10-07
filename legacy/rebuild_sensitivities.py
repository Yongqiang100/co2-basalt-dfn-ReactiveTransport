#!/usr/bin/env python3
"""
Rebuild sensitivities_corrected.csv from the rescaled Block C runs.

The existing file was written from decks whose mineral volume fractions summed
to 0.85. runs_C050/ holds the same variants at 1 - phi = 0.50, so the ratios
change: anor_as30 at p32_150_s117 was 2.56 and is now about 1.95.

The schema is preserved so src/letter_figures.py needs no change.

    python3 src/rebuild_sensitivities.py
    python3 src/rebuild_sensitivities.py --runs runs_aperture --out sensitivities_aperture.csv
"""
from __future__ import annotations
import argparse
import csv
import os
import sys

import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6


def read(run_dir):
    """Total carbonate per cell, per precipitating cell, and dominant phase."""
    path = os.path.join(run_dir, "pflotran_co2.h5")
    with h5py.File(path, "r") as f:
        g = f[sorted(f.keys())[-1]]
        per_phase, total, ncell = {}, None, None
        for p in CARB:
            k = f"{p} VF [m^3 mnrl_m^3 bulk]"
            if k in g:
                d = np.maximum(g[k][:].flatten() - SEED, 0.0)
                ncell = len(d)
                per_phase[p] = d.sum()
                total = d if total is None else total + d
    if total is None or ncell is None:
        return None
    nprec = int((total > 0).sum())
    s = sum(per_phase.values())
    dominant = max(per_phase, key=per_phase.get) if s > 0 else None
    short = {"Calcite": "cal", "Magnesite": "mag",
             "Siderite": "sid", "Dawsonite": "daw"}
    return {
        "per_cell": total.sum() / ncell,
        "intensity": total.sum() / nprec if nprec else 0.0,
        "dominant": short.get(dominant, "none"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_C050")
    ap.add_argument("--out", default="sensitivities_corrected.csv")
    a = ap.parse_args()

    cases = sorted(d for d in os.listdir(a.runs) if d.startswith("C_"))
    rows, skipped = [], []

    for case in cases:
        if case.startswith("C_baseline"):
            continue
        variant, _, tail = case.partition("__")
        variant = variant[2:]                     # drop the C_ prefix
        base = f"C_baseline__{tail}"
        if not os.path.isfile(f"{a.runs}/{base}/pflotran_co2.h5"):
            skipped.append((case, "no same-seed baseline"))
            continue
        v = read(f"{a.runs}/{case}")
        b = read(f"{a.runs}/{base}")
        if v is None or b is None:
            skipped.append((case, "unreadable"))
            continue
        if b["per_cell"] <= 0:
            skipped.append((case, "baseline has no carbonate"))
            continue
        rows.append({
            "variant": variant,
            "case": tail,
            "base_per_cell": b["per_cell"],
            "var_per_cell": v["per_cell"],
            "ratio": v["per_cell"] / b["per_cell"],
            "base_intensity": b["intensity"],
            "var_intensity": v["intensity"],
            "intensity_ratio": (v["intensity"] / b["intensity"]
                                if b["intensity"] > 0 else float("nan")),
            "base_dominant": b["dominant"],
            "var_dominant": v["dominant"],
        })

    fields = ["variant", "case", "base_per_cell", "var_per_cell", "ratio",
              "base_intensity", "var_intensity", "intensity_ratio",
              "base_dominant", "var_dominant"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {a.out}: {len(rows)} variant comparison(s) from {a.runs}")
    print(f"variant types: {len({r['variant'] for r in rows})}")
    if skipped:
        print(f"\nskipped {len(skipped)}:")
        for c, why in skipped:
            print(f"  {c}: {why}")

    print("\nmedian ratio by variant:")
    by = {}
    for r in rows:
        by.setdefault(r["variant"], []).append(r["ratio"])
    for v in sorted(by):
        x = by[v]
        print(f"  {v:20s} n={len(x)}  median {np.median(x):.2f}  "
              f"range {min(x):.2f} to {max(x):.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
