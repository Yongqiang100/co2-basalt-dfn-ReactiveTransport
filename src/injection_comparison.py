#!/usr/bin/env python3
"""
How much carbon is mineralized under each injection schedule.

Four quantities are needed to compare schedules, and reporting any one alone
misleads:

  carbonate formed        moles of carbon fixed as calcite, magnesite,
                          siderite and dawsonite
  carbon injected         the injected water mass from the PFLOTRAN mass
                          balance, times the injectate concentration
  efficiency              the first divided by the second
  carbon left over        injected minus fixed, which is the carbon that
                          passed through the domain without reacting

A short pulse traps less carbonate in absolute terms and yet reaches a much
higher efficiency, because it injects far less carbon. The fill series on one
network shows why: raising the injected mass fifty-fold leaves the carbonate
unchanged at about 2.26 moles, so the carbonate is set by something other than
the carbon supply and every additional mole injected only lowers the
efficiency.

Baselines are the matching C_baseline__<case> runs, which inject continuously
for fifty years. Molar volumes are those in hanford.dat.

    python3 src/injection_comparison.py
    python3 src/injection_comparison.py --csv injection_comparison.csv
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
CONC = 0.82


def carbonate_mol(run_dir):
    uge = os.path.join(run_dir, "full_mesh.uge")
    h5 = os.path.join(run_dir, "pflotran_co2.h5")
    if not (os.path.isfile(uge) and os.path.isfile(h5)):
        return None, None
    with open(uge) as f:
        n = int(f.readline().split()[1])
        vol = np.empty(n)
        for i in range(n):
            vol[i] = float(f.readline().split()[4])
    with h5py.File(h5, "r") as f:
        timed = [(float(m.group(1)), k) for k in f.keys()
                 if (m := re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k))]
        if not timed:
            return None, None
        t, grp = max(timed)
        g = f[grp]
        tot = 0.0
        seen = False
        for p, mv in MOLAR_VOL_CM3.items():
            kk = next((x for x in g if x.startswith(p + " VF")), None)
            if kk is None:
                continue
            seen = True
            v = np.maximum(np.asarray(g[kk][:], float).flatten() - SEED, 0.0)
            tot += float((v * vol).sum()) / (mv * 1e-6)
    return (tot, t) if seen else (None, None)


def injected_kg(run_dir):
    mas = os.path.join(run_dir, "pflotran_co2-mas.dat")
    if not os.path.isfile(mas):
        return None
    with open(mas) as f:
        header = [h.strip().strip('"') for h in f.readline().split(",")]
        last = None
        for line in f:
            if line.strip():
                last = line
    if last is None:
        return None
    vals = [float(x) for x in last.split()]
    for n_, v_ in zip(header, vals):
        if "injector water [kg]" in n_.lower():
            return v_
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--csv")
    a = ap.parse_args()

    # every variant run whose name carries a schedule, plus the baselines
    want = []
    for d in sorted(glob.glob(os.path.join(a.runs, "C_*"))):
        nm = os.path.basename(d)
        if "__" not in nm:
            continue
        variant, case = nm.split("__", 1)
        variant = variant[2:]
        if variant.startswith(("shutin", "fieldscale", "baseline")):
            want.append((nm, variant, case, d))
    if not want:
        print(f"no schedule runs under {a.runs}")
        return 1

    rows = []
    for nm, variant, case, d in want:
        carb, t = carbonate_mol(d)
        kg = injected_kg(d)
        if carb is None or kg is None or kg <= 0:
            rows.append((nm, variant, case, np.nan, np.nan, np.nan, np.nan))
            continue
        C_in = kg * CONC
        rows.append((nm, variant, case, carb, C_in,
                     100.0 * carb / C_in, C_in - carb))

    base = {r[2]: r for r in rows if r[1] == "baseline" and np.isfinite(r[3])}

    print("%-36s %10s %12s %12s %9s" %
          ("run", "carb mol", "C in mol", "efficiency %", "vs base"))
    print("-" * 84)
    for nm, variant, case, carb, C_in, eff, left in rows:
        if variant == "baseline":
            continue
        if not np.isfinite(carb):
            print("%-36s %10s" % (nm, "no output"))
            continue
        b = base.get(case)
        rel = f"x{carb/b[3]:.3f}" if b and b[3] > 0 else "-"
        print("%-36s %10.3f %12.4e %12.6f %9s"
              % (nm, carb, C_in, eff, rel))
    print()
    print("the continuous baselines, for reference")
    for case, r in sorted(base.items()):
        print("%-36s %10.3f %12.4e %12.6f"
              % (r[0], r[3], r[4], r[5]))
    print()

    # the fill series, where only the injected mass changes
    fs = [r for r in rows if r[1].startswith("fieldscale")
          and np.isfinite(r[3])]
    if len(fs) > 1:
        fs.sort(key=lambda r: r[4])
        print("the fill series, one network, injected mass varied")
        print("%-36s %10s %12s %12s" %
              ("run", "carb mol", "C in mol", "efficiency %"))
        for nm, v, c, carb, C_in, eff, left in fs:
            print("%-36s %10.3f %12.4e %12.6f" % (nm, carb, C_in, eff))
        span_c = fs[-1][4] / fs[0][4]
        span_m = fs[-1][3] / fs[0][3]
        print()
        print(f"  injected carbon spans {span_c:.0f} times")
        print(f"  carbonate spans {span_m:.3f} times")
        if span_c > 5 and abs(span_m - 1) < 0.2:
            print("  The carbonate does not follow the injected carbon, so the")
            print("  carbon is not the limiting reactant and the additional")
            print("  injection only lowers the efficiency.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["run", "variant", "case", "carbonate_mol",
                        "carbon_injected_mol", "efficiency_pct",
                        "carbon_unreacted_mol"])
            w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
