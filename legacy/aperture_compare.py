#!/usr/bin/env python3
"""
Compare the fixed-aperture and evolving-aperture ensembles.

All networks are included. Ratios are unstable where the fixed-aperture run
precipitated in very few cells, so the summary reports the median, which is
insensitive to those values, alongside the count of networks in each direction.
The per-network table carries the cell counts so the unstable cases are
identifiable.

Usage
-----
    python3 src/aperture_compare.py
    python3 src/aperture_compare.py --block E_p32_
"""
from __future__ import annotations
import argparse, os
import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6


def carbonate_field(run_dir):
    """Net carbonate volume fraction per cell, seed value removed."""
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[sorted(f.keys())[-1]]
        total = None
        for phase in CARB:
            key = f"{phase} VF [m^3 mnrl_m^3 bulk]"
            if key in g:
                d = np.maximum(g[key][:].flatten() - SEED, 0.0)
                total = d if total is None else total + d
    return total


def summarise(x, label):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        print(f"  {label}: no cases")
        return
    print(f"  {label}: n={len(x)}  median={np.median(x):.2f}  "
          f"IQR={np.percentile(x, 25):.2f}-{np.percentile(x, 75):.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", default="A_p32_")
    ap.add_argument("--fixed", default="runs_gravityoff")
    ap.add_argument("--aperture", default="runs_aperture")
    a = ap.parse_args()

    print("%-24s %10s %7s %10s %7s %8s %7s"
          % ("network", "vol_fix", "n_fix", "vol_ap", "n_ap", "vol_r", "n_r"))
    print("-" * 80)

    vol_r, cell_r, zero_fix = [], [], []

    for case in sorted(d for d in os.listdir(a.aperture)
                       if d.startswith(a.block)):
        fx = os.path.join(a.fixed, case)
        apr = os.path.join(a.aperture, case)
        if not (os.path.isfile(os.path.join(fx, "pflotran_co2.h5"))
                and os.path.isfile(os.path.join(apr, "pflotran_co2.h5"))):
            continue
        try:
            f_fld = carbonate_field(fx)
            a_fld = carbonate_field(apr)
        except (OSError, KeyError):
            continue
        if f_fld is None or a_fld is None:
            continue

        vf, va = f_fld.sum(), a_fld.sum()
        nf, na = int((f_fld > 0).sum()), int((a_fld > 0).sum())

        if vf <= 0:
            # no carbonate under a fixed aperture: the ratio is undefined
            zero_fix.append((case, va, na))
            print("%-24s %10.2e %7d %10.2e %7d %8s %7s"
                  % (case, vf, nf, va, na, "--", "--"))
            continue

        v_r = va / vf
        c_r = na / nf if nf > 0 else np.nan
        vol_r.append(v_r)
        cell_r.append(c_r)
        print("%-24s %10.2e %7d %10.2e %7d %8.2f %7.2f"
              % (case, vf, nf, va, na, v_r, c_r))

    vr = np.asarray(vol_r)
    cr = np.asarray(cell_r)

    print("-" * 80)
    print(f"networks compared: {len(vr)}")
    if zero_fix:
        print(f"networks with no carbonate under a fixed aperture: "
              f"{len(zero_fix)} (ratio undefined, listed above)")
    print()
    summarise(vr, "carbonate volume ratio")
    summarise(cr, "precipitating cell ratio")
    print()

    down, up = vr <= 1.0, vr > 1.0
    print(f"  volume falls in {down.sum()} of {len(vr)} "
          f"({down.sum() / len(vr) * 100:.0f}%)")
    summarise(vr[down], "    volume ratio")
    summarise(cr[down], "    cell ratio")
    print(f"  volume rises in {up.sum()} of {len(vr)} "
          f"({up.sum() / len(vr) * 100:.0f}%)")
    summarise(vr[up], "    volume ratio")
    summarise(cr[up], "    cell ratio")
    print()

    ok = np.isfinite(cr) & (cr > 0)
    inten = np.full_like(vr, np.nan)
    inten[ok] = vr[ok] / cr[ok]
    print("  carbonate per precipitating cell:")
    summarise(inten[down & ok], "    where the volume falls")
    summarise(inten[up & ok], "    where the volume rises")
    print()
    print("  sensitivity of the headline figures to a minimum-cell filter:")
    print("  %-14s %6s %8s %9s" % ("threshold", "n", "median", "falls"))
    for thr in (0, 3, 5, 10, 20):
        keep = []
        for case in sorted(d for d in os.listdir(a.aperture)
                           if d.startswith(a.block)):
            fx = os.path.join(a.fixed, case)
            apr = os.path.join(a.aperture, case)
            if not (os.path.isfile(os.path.join(fx, "pflotran_co2.h5"))
                    and os.path.isfile(os.path.join(apr, "pflotran_co2.h5"))):
                continue
            try:
                f_fld, a_fld = carbonate_field(fx), carbonate_field(apr)
            except (OSError, KeyError):
                continue
            if f_fld is None or a_fld is None or f_fld.sum() <= 0:
                continue
            if int((f_fld > 0).sum()) < thr:
                continue
            keep.append(a_fld.sum() / f_fld.sum())
        k = np.asarray(keep)
        print("  %-14s %6d %8.2f %8.0f%%"
              % (f">= {thr} cells", len(k), np.median(k),
                 (k <= 1).sum() / len(k) * 100))


if __name__ == "__main__":
    main()
