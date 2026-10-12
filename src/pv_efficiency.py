#!/usr/bin/env python3
"""
Trapping efficiency measured against one pore volume of injectate.

Efficiency reported as carbon fixed over carbon injected is dominated by how
much was injected rather than by how much the rock can hold. These simulations
flush about four thousand pore volumes through the domain in fifty years, so
the ratio is small by construction: even a domain that precipitated everywhere
it thermodynamically could would return a figure near zero, because the
denominator is four thousand times the fluid the domain contains.

This computes the alternative denominator. One pore volume of injectate carries

    C_pv = (pore water mass) x (injectate concentration)

moles of carbon, and the ratio of the carbonate formed to that quantity says
what fraction of a single charge of injectate ends up as mineral. It does not
change with the duration of injection, so it describes the rock and the network
rather than the operating schedule.

The carbonate is taken from the volume fractions in the HDF5 output, converted
through the molar volumes, with the seeded volume fraction removed. That route
agrees with the efficiency column in finalised.csv to one part in a thousand.
The mineral totals in the mass-balance file give a lower figure, by a factor
that differs per phase, and are not used.

    python3 src/pv_efficiency.py
    python3 src/pv_efficiency.py --block A --top 10
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py

# Molar volumes, cubic centimetres per mole, read from the second field of
# each mineral record in hanford.dat. Calcite 36.9340 and dawsonite 58.5366
# are confirmed against that file; an earlier version of this script used
# 59.30 for dawsonite, which overstated its contribution by 1.3 per cent and
# the total by about 0.1 per cent.
MOLAR_VOL_CM3 = {
    "Calcite": 36.9340,
    "Magnesite": 28.018,
    "Siderite": 29.378,
    "Dawsonite": 58.5366,
}
CARB = list(MOLAR_VOL_CM3)
SEED = 1.0e-6


def carbonate_moles(h5_path, vol):
    with h5py.File(h5_path, "r") as f:
        timed = [(float(m.group(1)), k) for k in f.keys()
                 if (m := re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k))]
        if not timed:
            return None, None
        t, grp = max(timed)
        g = f[grp]
        tot = 0.0
        ok = False
        for p in CARB:
            kk = next((x for x in g if x.startswith(p + " VF")), None)
            if kk is None:
                continue
            ok = True
            v = np.maximum(np.asarray(g[kk][:], float).flatten() - SEED, 0.0)
            tot += float((v * vol).sum()) / (MOLAR_VOL_CM3[p] * 1e-6)
    return (t, tot) if ok else (None, None)


def read_volumes(uge):
    with open(uge) as f:
        n = int(f.readline().split()[1])
        v = np.empty(n)
        for i in range(n):
            v[i] = float(f.readline().split()[4])
    return v


def pore_water_kg(mas):
    """Water mass in the domain, from the final mass-balance row."""
    with open(mas) as f:
        header = [h.strip().strip('"') for h in f.readline().split(",")]
        last = None
        for line in f:
            if line.strip():
                last = line
    vals = [float(x) for x in last.split()]
    for n, v in zip(header, vals):
        if "global water mass" in n.lower():
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--block", default="A")
    ap.add_argument("--conc", type=float, default=0.82,
                    help="injectate concentration, mol per kg of water")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--csv")
    a = ap.parse_args()

    dirs = sorted(glob.glob(os.path.join(a.runs, f"{a.block}_p32_*")))
    if not dirs:
        print(f"no {a.block}_p32_* directories under {a.runs}")
        return 1

    rows = []
    for d in dirs:
        h5 = os.path.join(d, "pflotran_co2.h5")
        uge = os.path.join(d, "full_mesh.uge")
        mas = os.path.join(d, "pflotran_co2-mas.dat")
        if not all(os.path.isfile(x) for x in (h5, uge, mas)):
            continue
        try:
            vol = read_volumes(uge)
            t, carb = carbonate_moles(h5, vol)
            if carb is None:
                continue
            w = pore_water_kg(mas)
            if not w:
                continue
            c_pv = w * a.conc
            rows.append((os.path.basename(d), carb, c_pv, w, len(vol),
                         100.0 * carb / c_pv))
        except (OSError, KeyError, ValueError):
            continue

    if not rows:
        print("nothing computed")
        return 1

    rows.sort(key=lambda r: -r[5])
    print(f"carbonate formed against the carbon in one pore volume of "
          f"injectate, at {a.conc} mol/kg")
    print(f"{len(rows)} realizations in block {a.block}")
    print()
    print("%-22s %12s %12s %10s" %
          ("network", "carbonate", "C in 1 PV", "per cent"))
    print("%-22s %12s %12s %10s" % ("", "mol", "mol", ""))
    print("-" * 60)
    for nm, carb, cpv, w, n, pct in rows[:a.top]:
        print("%-22s %12.2f %12.1f %10.3f" % (nm, carb, cpv, pct))
    if len(rows) > a.top:
        print(f"... {len(rows)-a.top} more")

    p = np.array([r[5] for r in rows])
    print()
    print(f"highest {p.max():.3f} %   in {rows[0][0]}")
    print(f"median  {np.median(p):.3f} %")
    print(f"lowest  {p.min():.3f} %")
    print()
    print("for comparison, the same carbonate against all the carbon injected:")
    print("  that denominator is about four thousand pore volumes, so the "
          "figure is")
    print(f"  smaller by that factor")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w_ = _csv.writer(fh)
            w_.writerow(["run_id", "carbonate_mol", "carbon_in_one_pv_mol",
                         "pore_water_kg", "n_cells", "percent_of_one_pv"])
            w_.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
