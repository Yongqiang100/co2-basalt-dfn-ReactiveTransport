#!/usr/bin/env python3
"""
Why does coupling lower carbonate?

Injection is a prescribed mass rate, so the volumetric throughput is set by the
boundary condition and does not respond to permeability. Dissolution raises the
porosity, so the same cation release occupies a larger water volume. This tests
the two consequences that follow:

  outflow ratio near 1      -- throughput is fixed, as the boundary condition
                               implies, so residence time rises with porosity
  concentration ratio < 1   -- aqueous species are diluted, lowering the
                               saturation index and with it the precipitate

Usage
-----
    python3 src/dilution_test.py
    python3 src/dilution_test.py --block E_p32_
"""
from __future__ import annotations
import argparse, os
import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6
# aqueous species that carry the cations and the carbonate ion
AQ = ["Ca++ [M]", "Mg++ [M]", "CO3-- [M]", "HCO3- [M]"]


def mas_column(run_dir, want):
    """Final value of a named column in the mass-balance file."""
    path = os.path.join(run_dir, "pflotran_co2-mas.dat")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        header = f.readline()
        last = None
        for line in f:
            if line.strip():
                last = line.split()
    if last is None:
        return None
    cols = [c.strip().strip('"') for c in header.split(",")]
    for i, c in enumerate(cols):
        if want.lower() in c.lower():
            try:
                return abs(float(last[i]))
            except (IndexError, ValueError):
                return None
    return None


def h5_fields(run_dir):
    """Final-time carbonate, porosity if present, and aqueous means."""
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
        for k in g.keys():
            if "orosity" in k:
                out["poros"] = g[k][:].flatten()
                break
        for sp in AQ:
            if sp in g:
                out[sp] = g[sp][:].flatten()
    return out


def med(x):
    return float(np.median(x)) if len(x) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", default="A_p32_")
    ap.add_argument("--fixed", default="runs_gravityoff")
    ap.add_argument("--aperture", default="runs_aperture")
    ap.add_argument("--min-cells", type=int, default=5)
    a = ap.parse_args()

    outq, porq, volq = [], [], []
    conc = {sp: [] for sp in AQ}

    print("%-24s %8s %8s %8s %9s"
          % ("network", "outflow", "poros", "vol_r", "conc_Ca"))
    print("-" * 64)

    for case in sorted(d for d in os.listdir(a.aperture)
                       if d.startswith(a.block)):
        fx = os.path.join(a.fixed, case)
        apr = os.path.join(a.aperture, case)
        if not (os.path.isfile(os.path.join(fx, "pflotran_co2.h5"))
                and os.path.isfile(os.path.join(apr, "pflotran_co2.h5"))):
            continue
        try:
            F = h5_fields(fx)
            A = h5_fields(apr)
        except (OSError, KeyError):
            continue
        if F.get("carb") is None or A.get("carb") is None:
            continue
        if int((F["carb"] > 0).sum()) < a.min_cells:
            continue

        of = mas_column(fx, "outflow Water Mass [kg/y]")
        oa = mas_column(apr, "outflow Water Mass [kg/y]")
        o_r = oa / of if (of and oa and of > 0) else np.nan

        p_r = np.nan
        if "poros" in F and "poros" in A:
            p_r = float(np.median(A["poros"]) / np.median(F["poros"]))
        elif "poros" in A:
            p_r = float(np.median(A["poros"]) / 0.50)

        v_r = A["carb"].sum() / F["carb"].sum()

        # concentrations where the fixed run precipitated: the sites that matter
        mask = F["carb"] > 0
        c_r = {}
        for sp in AQ:
            if sp in F and sp in A and mask.sum():
                fv = float(np.median(F[sp][mask]))
                av = float(np.median(A[sp][mask]))
                c_r[sp] = av / fv if fv > 0 else np.nan

        print("%-24s %8.2f %8.2f %8.2f %9.2f"
              % (case, o_r, p_r, v_r, c_r.get("Ca++ [M]", np.nan)))

        if not np.isnan(o_r):
            outq.append(o_r)
        if not np.isnan(p_r):
            porq.append(p_r)
        volq.append(v_r)
        for sp, v in c_r.items():
            if not np.isnan(v):
                conc[sp].append(v)

    print("-" * 64)
    print(f"networks: {len(volq)}\n")
    print(f"outflow ratio           : {med(outq):.3f}   "
          f"(near 1 confirms the throughput is set by the injection rate)")
    print(f"porosity ratio          : {med(porq):.3f}   "
          f"(above 1 means dissolution opened the pore space)")
    print(f"carbonate volume ratio  : {med(volq):.3f}")
    print()
    print("aqueous concentration at the cells that precipitated with a")
    print("fixed aperture, coupled value over fixed value:")
    for sp in AQ:
        if conc[sp]:
            print(f"  {sp:<12s} {med(conc[sp]):.3f}  (n={len(conc[sp])})")
    print()
    print("If the porosity ratio and the inverse concentration ratio agree,")
    print("the loss is dilution: the same cation release spread through more")
    print("water, lowering the saturation index.")


if __name__ == "__main__":
    main()
