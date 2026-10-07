#!/usr/bin/env python3
"""
Did the flow field actually reorganise?

The cell-count comparison implies that coupling relocates precipitation, but
that is an inference. This measures the permeability field directly: how far it
departs from the fixed-aperture case, and whether the cells that gained or lost
permeability are the cells that gained or lost carbonate.

Usage
-----
    python3 src/flowfield_change.py
    python3 src/flowfield_change.py --block E_p32_
"""
from __future__ import annotations
import argparse, os
import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6
PERM_KEY = "Permeability [m^2]"


def fields(run_dir):
    """Final-time permeability and net carbonate per cell."""
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[sorted(f.keys())[-1]]
        perm = g[PERM_KEY][:].flatten() if PERM_KEY in g else None
        carb = None
        for phase in CARB:
            key = f"{phase} VF [m^3 mnrl_m^3 bulk]"
            if key in g:
                d = np.maximum(g[key][:].flatten() - SEED, 0.0)
                carb = d if carb is None else carb + d
    return perm, carb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", default="A_p32_")
    ap.add_argument("--fixed", default="runs_gravityoff")
    ap.add_argument("--aperture", default="runs_aperture")
    ap.add_argument("--min-cells", type=int, default=5)
    a = ap.parse_args()

    print("%-24s %8s %9s %9s %8s %8s"
          % ("network", "k_ratio", "k_spread", "frac_chg", "vol_r", "rho"))
    print("-" * 74)

    kr, spread, changed, rho_all = [], [], [], []

    for case in sorted(d for d in os.listdir(a.aperture)
                       if d.startswith(a.block)):
        fx = os.path.join(a.fixed, case)
        apr = os.path.join(a.aperture, case)
        if not (os.path.isfile(os.path.join(fx, "pflotran_co2.h5"))
                and os.path.isfile(os.path.join(apr, "pflotran_co2.h5"))):
            continue
        try:
            pf, cf = fields(fx)
            pa, ca = fields(apr)
        except (OSError, KeyError):
            continue
        if pf is None or pa is None or cf is None or ca is None:
            continue
        if len(pf) != len(pa) or int((cf > 0).sum()) < a.min_cells:
            continue

        # how far the permeability field moved
        ratio = pa / np.maximum(pf, 1e-30)
        med = float(np.median(ratio))
        # spread of the per-cell permeability change, in orders of magnitude
        lg = np.log10(np.maximum(ratio, 1e-30))
        iqr = float(np.percentile(lg, 75) - np.percentile(lg, 25))
        # fraction of cells whose permeability moved by more than 1 %
        frac = float((np.abs(ratio - 1.0) > 0.01).mean())

        vol_r = ca.sum() / cf.sum() if cf.sum() > 0 else np.nan

        # do the cells that gained permeability gain carbonate?
        dk = lg
        dc = np.log10(np.maximum(ca, 1e-30)) - np.log10(np.maximum(cf, 1e-30))
        ok = (cf > 0) | (ca > 0)
        if ok.sum() > 10:
            x, y = dk[ok], dc[ok]
            if x.std() > 0 and y.std() > 0:
                rho = float(np.corrcoef(x, y)[0, 1])
            else:
                rho = np.nan
        else:
            rho = np.nan

        print("%-24s %8.3f %9.2f %9.3f %8.2f %8.2f"
              % (case, med, iqr, frac, vol_r, rho))
        kr.append(med); spread.append(iqr); changed.append(frac)
        if not np.isnan(rho):
            rho_all.append(rho)

    print("-" * 74)
    if not kr:
        print("no comparable pairs")
        return
    print(f"networks: {len(kr)}")
    print(f"median permeability ratio:        {np.median(kr):.3f}")
    print(f"spread of per-cell change (IQR of log10): "
          f"{np.median(spread):.2f} orders")
    print(f"cells with permeability change >1%: "
          f"{np.median(changed) * 100:.1f}%")
    if rho_all:
        print(f"correlation, per-cell log permeability change against "
              f"log carbonate change: median rho = {np.median(rho_all):+.2f} "
              f"(n={len(rho_all)})")
    print()
    print("Reading: k_ratio near 1 with a wide spread means the field was")
    print("reorganised without a net opening or closing. A positive rho means")
    print("cells that gained permeability gained carbonate.")


if __name__ == "__main__":
    main()
