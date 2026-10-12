#!/usr/bin/env python3
"""
Reconciles the trapping-efficiency figures computed three different ways.

Three values for the same quantity appeared during the revision and they do not
agree. The manuscript reports 0.0002 per cent. The efficiency_pct column that
finalise.py writes into finalised.csv gives 0.000947 per cent for the
highest-carbonate network. A calculation from the PFLOTRAN mass-balance file
gave 0.00058 per cent for that same network. A factor of 1.6 separates the last
two and they should be identical.

This script computes the efficiency by every available route on one network and
prints the intermediate quantities, so the step where they diverge is visible
rather than inferred.

  route A   the mass-balance file, pflotran_co2-mas.dat, using the mineral
            totals in moles that PFLOTRAN itself reports, and the injected
            water mass times the injectate concentration

  route B   the HDF5 output, converting each carbonate volume fraction to moles
            through its molar volume, and the injected water from the same
            mass-balance file

  route C   the columns finalise.py already wrote into finalised.csv

Where A and B differ, the molar volumes or the seeded volume fraction are the
cause. Where they differ from C, the definition in finalise.py differs.

    python3 check_efficiency.py --run runs/A_p32_100_s1181
    python3 check_efficiency.py --run runs/A_p32_100_s1181 --csv finalised.csv
"""
from __future__ import annotations
import argparse
import csv as _csv
import os
import re
import sys

import numpy as np
import h5py

# molar volumes, cubic centimetres per mole, from the PFLOTRAN database
MOLAR_VOL_CM3 = {
    "Calcite": 36.934,
    "Magnesite": 28.018,
    "Siderite": 29.378,
    "Dawsonite": 59.30,
}
CARB = list(MOLAR_VOL_CM3)
SEED = 1.0e-6


def mas_final_row(path):
    """Column names and the final row of the mass-balance file."""
    with open(path) as f:
        header = f.readline()
        last = None
        for line in f:
            if line.strip():
                last = line
    names = [h.strip().strip('"') for h in header.split(",")]
    vals = [float(x) for x in last.split()]
    return names, vals


def pick(names, vals, needle):
    for n, v in zip(names, vals):
        if needle.lower() in n.lower():
            return n, v
    return None, None


def read_h5(path):
    with h5py.File(path, "r") as f:
        timed = [(float(m.group(1)), k) for k in f.keys()
                 if (m := re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k))]
        t, grp = max(timed)
        g = f[grp]
        vf = {}
        for p in CARB:
            kk = next((x for x in g if x.startswith(p + " VF")), None)
            if kk:
                vf[p] = np.asarray(g[kk][:], float).flatten()
    return t, vf


def read_uge_volumes(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        vol = np.empty(n)
        for i in range(n):
            vol[i] = float(f.readline().split()[4])
    return vol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--csv", default="finalised.csv")
    ap.add_argument("--conc", type=float, default=0.82,
                    help="injectate concentration, mol per kg of water")
    a = ap.parse_args()

    d = a.run
    name = os.path.basename(d.rstrip("/"))
    mas = os.path.join(d, "pflotran_co2-mas.dat")
    h5 = os.path.join(d, "pflotran_co2.h5")
    uge = os.path.join(d, "full_mesh.uge")

    print(f"network: {name}")
    print()

    # ---------------------------------------------------------------- route A
    print("route A, from the mass-balance file")
    names, vals = mas_final_row(mas)
    tA = vals[0]
    carbA = {}
    for p in CARB:
        n_, v_ = pick(names, vals, f"{p} Total Mass")
        if n_:
            carbA[p] = v_
    _, w_inj = pick(names, vals, "injector Water [kg]")
    _, w_dom = pick(names, vals, "Global Water Mass")
    totA = sum(carbA.values())
    C_in = w_inj * a.conc
    print(f"  final time              {tA:g} y")
    for p, v in carbA.items():
        print(f"  {p:10s}              {v:.4f} mol")
    print(f"  carbonate total         {totA:.3f} mol")
    print(f"  water injected          {w_inj:.5e} kg")
    print(f"  carbon injected         {C_in:.5e} mol   "
          f"(water x {a.conc} mol/kg)")
    effA = 100.0 * totA / C_in
    print(f"  efficiency              {effA:.6f} %")
    print()

    # ---------------------------------------------------------------- route B
    print("route B, from the HDF5 volume fractions and the molar volumes")
    tB, vf = read_h5(h5)
    vol = read_uge_volumes(uge)
    ncell = len(vol)
    print(f"  snapshot                {tB:g} y")
    print(f"  cells                   {ncell}")
    totB = 0.0
    totB_seeded = 0.0
    for p, f_ in vf.items():
        mv = MOLAR_VOL_CM3[p] * 1e-6           # m3 per mol
        raw = float((f_ * vol).sum()) / mv
        net = float((np.maximum(f_ - SEED, 0.0) * vol).sum()) / mv
        totB += net
        totB_seeded += raw
        print(f"  {p:10s} net {net:10.3f} mol   with the seed included "
              f"{raw:10.3f} mol")
    print(f"  carbonate total, net    {totB:.3f} mol")
    print(f"  carbonate total, raw    {totB_seeded:.3f} mol")
    effB = 100.0 * totB / C_in
    print(f"  efficiency, net         {effB:.6f} %")
    print()

    # ---------------------------------------------------------------- route C
    print("route C, the columns finalise.py wrote")
    row = None
    if os.path.isfile(a.csv):
        for r in _csv.DictReader(open(a.csv)):
            if r["run_id"] == name:
                row = r
                break
    if row is None:
        print(f"  {name} not found in {a.csv}")
    else:
        for k in ("carb_CO2_mol", "injected_CO2_mol", "injected_kg",
                  "efficiency_pct", "carb_per_cell", "pore_volumes_flushed"):
            if k in row:
                print(f"  {k:22s} {row[k]}")
        effC = float(row["efficiency_pct"])
        print()

        # ------------------------------------------------------- reconciling
        print("=== where the routes diverge ===")
        cc = float(row.get("carb_CO2_mol") or "nan")
        ic = float(row.get("injected_CO2_mol") or "nan")
        print(f"  carbonate, A {totA:10.3f}   B {totB:10.3f}   "
              f"C {cc:10.3f} mol")
        print(f"  injected,   A {C_in:10.4e}                  "
              f"C {ic:10.4e} mol")
        print()
        if np.isfinite(cc) and totA > 0:
            print(f"  carbonate C over A      {cc/totA:.4f}")
        if np.isfinite(cc) and totB > 0:
            print(f"  carbonate C over B      {cc/totB:.4f}")
        if np.isfinite(ic) and C_in > 0:
            print(f"  injected  C over A      {ic/C_in:.4f}")
        print()
        print(f"  efficiency  A {effA:.6f} %   B {effB:.6f} %   "
              f"C {effC:.6f} %")
        print()
        print("  A carbonate close to B means the molar volumes are right.")
        print("  C injected differing from A points at the injectate "
              "concentration")
        print("  or at whether finalise.py counts water or carbon.")
        print("  C carbonate differing from B points at the seed subtraction.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
