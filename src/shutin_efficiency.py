#!/usr/bin/env python3
"""
Compares continuous injection against shut-in, on carbonate amount and on
trapping efficiency separately.

The letter reports that a shut-in ensemble produces no systematic increase in
carbonate, with a median ratio of 1.1 across the matched pairs. That comparison
is of the carbonate amount at fifty years, and it is correct.

Efficiency is a different comparison, because the two cases inject different
amounts of carbon. The continuous case injects for fifty years; the shut-in
case injects for ten and then monitors for forty. The shut-in denominator is
therefore about a fifth of the continuous one, so two cases that form the same
carbonate do not have the same efficiency.

This matters because the efficiency of the continuous case peaks at three to
five years and declines thereafter: the carbonate continues to accumulate while
the injected carbon grows faster. Stopping injection at ten years freezes the
denominator while the numerator keeps rising, so the shut-in efficiency should
continue to improve after the continuous efficiency has begun to fall.

The script reports both ratios for every matched pair, so the amount and the
efficiency can be stated separately rather than one being taken for the other.

    python3 src/shutin_efficiency.py
    python3 src/shutin_efficiency.py --csv shutin.csv
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py

MOLAR_VOL_CM3 = {
    "Calcite": 36.9340,
    "Magnesite": 28.0180,
    "Siderite": 29.3780,
    "Dawsonite": 58.5366,
}
SEED = 1.0e-6
YEAR = 3.156e7
RATE_KG_S = 1.06725e-2
CONC = 0.82


def read_volumes(uge):
    with open(uge) as f:
        n = int(f.readline().split()[1])
        v = np.empty(n)
        for i in range(n):
            v[i] = float(f.readline().split()[4])
    return v


def carbonate_at_end(h5, vol):
    """Carbonate in moles at the last snapshot, and that time."""
    with h5py.File(h5, "r") as f:
        tg = sorted((float(m.group(1)), k) for k in f.keys()
                    if (m := re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)))
        if not tg:
            return None, None
        t, k = tg[-1]
        g = f[k]
        tot = 0.0
        seen = False
        for p, mv in MOLAR_VOL_CM3.items():
            kk = next((x for x in g if x.startswith(p + " VF")), None)
            if kk is None:
                continue
            seen = True
            v = np.maximum(np.asarray(g[kk][:], float).flatten() - SEED, 0.0)
            tot += float((v * vol).sum()) / (mv * 1e-6)
    return (t, tot) if seen else (None, None)


def injected_kg(path):
    """Water injected, from the final row of the mass-balance file."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        header = [h.strip().strip('"') for h in f.readline().split(",")]
        last = None
        for line in f:
            if line.strip():
                last = line
    if last is None:
        return None
    vals = [float(x) for x in last.split()]
    for n, v in zip(header, vals):
        if "injector water [kg]" in n.lower():
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--csv")
    a = ap.parse_args()

    rows = []
    for e in sorted(glob.glob(os.path.join(a.runs, "E_p32_*"))):
        name = os.path.basename(e)
        c = os.path.join(a.runs, "A" + name[1:])
        if not os.path.isdir(c):
            continue
        try:
            vol = read_volumes(os.path.join(c, "full_mesh.uge"))
            t_c, carb_c = carbonate_at_end(
                os.path.join(c, "pflotran_co2.h5"), vol)
            t_e, carb_e = carbonate_at_end(
                os.path.join(e, "pflotran_co2.h5"), vol)
        except (OSError, KeyError, ValueError):
            continue
        if carb_c is None or carb_e is None:
            continue

        # injected carbon: from the mass balance where available, otherwise
        # from the rate and the injection period
        kg_c = injected_kg(os.path.join(c, "pflotran_co2-mas.dat"))
        kg_e = injected_kg(os.path.join(e, "pflotran_co2-mas.dat"))
        if kg_c is None:
            kg_c = RATE_KG_S * t_c * YEAR
        if kg_e is None:
            kg_e = RATE_KG_S * 10.0 * YEAR
        C_c, C_e = kg_c * CONC, kg_e * CONC
        eff_c = 100.0 * carb_c / C_c if C_c > 0 else np.nan
        eff_e = 100.0 * carb_e / C_e if C_e > 0 else np.nan
        rows.append((name[2:], carb_c, carb_e, C_c, C_e, eff_c, eff_e))

    if not rows:
        print(f"no matched pairs found under {a.runs}")
        return 1

    print(f"{len(rows)} matched pairs")
    print()
    print("carbon injected, median over the pairs")
    print(f"  continuous {np.median([r[3] for r in rows]):.4e} mol")
    print(f"  shut-in    {np.median([r[4] for r in rows]):.4e} mol")
    print(f"  ratio      {np.median([r[4]/r[3] for r in rows if r[3]>0]):.3f}")
    print()

    ok = [r for r in rows if r[1] > 0 and r[2] > 0]
    if not ok:
        print("no pair has carbonate in both cases")
        return 1
    amt = np.array([r[2] / r[1] for r in ok])
    eff = np.array([r[6] / r[5] for r in ok if np.isfinite(r[5]) and r[5] > 0])

    print(f"{len(ok)} pairs with carbonate in both cases")
    print()
    print("shut-in over continuous")
    print(f"  carbonate amount : median {np.median(amt):.3f}, "
          f"quartiles {np.percentile(amt,25):.3f} to "
          f"{np.percentile(amt,75):.3f}")
    print(f"  efficiency       : median {np.median(eff):.3f}, "
          f"quartiles {np.percentile(eff,25):.3f} to "
          f"{np.percentile(eff,75):.3f}")
    print()
    print(f"  amount higher under shut-in in {int((amt>1).sum())} of "
          f"{len(amt)} pairs")
    print(f"  efficiency higher under shut-in in {int((eff>1).sum())} of "
          f"{len(eff)} pairs")
    print()
    print("the ten highest by continuous efficiency")
    print("%-18s %10s %10s %11s %11s" %
          ("network", "carb cont", "carb shut", "eff cont %", "eff shut %"))
    for nm, cc, ce, Cc, Ce, ec, ee in sorted(ok, key=lambda r: -r[5])[:10]:
        print("%-18s %10.2f %10.2f %11.6f %11.6f" % (nm, cc, ce, ec, ee))

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["network", "carbonate_continuous_mol",
                        "carbonate_shutin_mol", "carbon_continuous_mol",
                        "carbon_shutin_mol", "efficiency_continuous_pct",
                        "efficiency_shutin_pct"])
            w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
