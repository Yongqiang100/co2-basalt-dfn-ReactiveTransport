#!/usr/bin/env python3
"""
Where is the carbonate lost under aperture evolution?

Ca++, Mg++ and porosity were not written to the output files, but the primary
mineral volume fractions were, so the cation release per cell is recoverable.
pH was also written. Between them these separate two possibilities:

  supply fell      -- less primary mineral dissolved, so fewer cations were
                      released into the cells that precipitated
  retention fell   -- the same cations were released but less was captured,
                      which points to the chemical conditions rather than the
                      supply

Cation release is measured as the drop in the volume fraction of the
Ca-bearing and Mg-bearing primary minerals from their initial values.

Usage
-----
    python3 src/supply_vs_retention.py
    python3 src/supply_vs_retention.py --block E_p32_
"""
from __future__ import annotations
import argparse, os
import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6

# initial volume fractions after rescaling to sum to 1 - phi = 0.50, and the
# cations each primary mineral releases
PRIMARY = {
    "Anorthite":  (0.1765, "Ca"),
    "Diopside":   (0.1471, "CaMg"),
    "Forsterite": (0.0176, "Mg"),
    "Enstatite":  (0.0294, "Mg"),
}


def read(run_dir):
    """Final-time carbonate, pH, and the released fraction of each primary."""
    out = {}
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[sorted(f.keys())[-1]]

        carb = None
        for phase in CARB:
            key = f"{phase} VF [m^3 mnrl_m^3 bulk]"
            if key in g:
                d = np.maximum(g[key][:].flatten() - SEED, 0.0)
                carb = d if carb is None else carb + d
        out["carb"] = carb

        out["pH"] = g["pH"][:].flatten() if "pH" in g else None

        released = None
        for mineral, (vf0, _) in PRIMARY.items():
            key = f"{mineral} VF [m^3 mnrl_m^3 bulk]"
            if key in g:
                lost = np.maximum(vf0 - g[key][:].flatten(), 0.0)
                released = lost if released is None else released + lost
        out["released"] = released
    return out


def med(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", default="A_p32_")
    ap.add_argument("--fixed", default="runs_gravityoff")
    ap.add_argument("--aperture", default="runs_aperture")
    ap.add_argument("--min-cells", type=int, default=5)
    a = ap.parse_args()

    vol_r, sup_r, eff_r, ph_f, ph_a, dph = [], [], [], [], [], []

    print("%-24s %7s %7s %7s %7s %7s"
          % ("network", "vol_r", "supp_r", "capt_r", "pH_fix", "pH_ap"))
    print("-" * 66)

    for case in sorted(d for d in os.listdir(a.aperture)
                       if d.startswith(a.block)):
        fx = os.path.join(a.fixed, case)
        apr = os.path.join(a.aperture, case)
        if not (os.path.isfile(os.path.join(fx, "pflotran_co2.h5"))
                and os.path.isfile(os.path.join(apr, "pflotran_co2.h5"))):
            continue
        try:
            F, A = read(fx), read(apr)
        except (OSError, KeyError):
            continue
        if F["carb"] is None or A["carb"] is None:
            continue
        if F["released"] is None or A["released"] is None:
            continue

        mask = F["carb"] > 0          # cells that precipitated with a fixed aperture
        if mask.sum() < a.min_cells:
            continue

        # supply into those cells, and the fraction of it captured as carbonate
        sf, sa = F["released"][mask].sum(), A["released"][mask].sum()
        cf, ca = F["carb"][mask].sum(), A["carb"][mask].sum()
        s_r = sa / sf if sf > 0 else np.nan
        capt_f = cf / sf if sf > 0 else np.nan
        capt_a = ca / sa if sa > 0 else np.nan
        e_r = capt_a / capt_f if capt_f and np.isfinite(capt_f) else np.nan

        pf = med(F["pH"][mask]) if F["pH"] is not None else np.nan
        pa = med(A["pH"][mask]) if A["pH"] is not None else np.nan

        v_r = A["carb"].sum() / F["carb"].sum()

        print("%-24s %7.2f %7.2f %7.2f %7.2f %7.2f"
              % (case, v_r, s_r, e_r, pf, pa))

        vol_r.append(v_r); sup_r.append(s_r); eff_r.append(e_r)
        ph_f.append(pf); ph_a.append(pa); dph.append(pa - pf)

    print("-" * 66)
    print(f"networks: {len(vol_r)}\n")
    print(f"carbonate volume ratio      : {med(vol_r):.3f}")
    print(f"cation supply ratio         : {med(sup_r):.3f}   "
          f"(release into the same cells)")
    print(f"capture efficiency ratio    : {med(eff_r):.3f}   "
          f"(carbonate formed per unit released)")
    print(f"pH, fixed aperture          : {med(ph_f):.2f}")
    print(f"pH, evolving aperture       : {med(ph_a):.2f}")
    print(f"pH change                   : {med(dph):+.2f}")
    print()
    print("A supply ratio near 1 with a capture ratio below 1 means the")
    print("cations arrived but were not captured. A supply ratio below 1")
    print("means less primary mineral dissolved in those cells.")


if __name__ == "__main__":
    main()
