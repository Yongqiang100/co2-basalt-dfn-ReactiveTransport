#!/usr/bin/env python3
"""
Dissolved carbon at shut-in, measured in the verification runs, against the injectate-equivalent
estimate used for the 108 % ratio.

    python3 src/verify_shutin_dic.py

For each run in runs_verify_shutin/:
  DIC_actual      sum over cells of Total CO2(aq) [M] x 1000 L/m3 x cell volume x porosity at 10 y
                  (porosity = 1 - sum of mineral volume fractions, liquid saturation 1)
  inj_equivalent  pore-water mass at 10 y (mass balance) x injectate concentration of the deck
  carbonate_50    carbonate at 50 y of the original shut-in run (results/shutin_carbon_balance.csv)
Reports DIC_actual / inj_equivalent and carbonate_50 / DIC_actual. Writes
results/verify_shutin_dic.csv.
"""
import csv, glob, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import shutin_carbon_balance as sb


def dic(d, t=10.0):
    snaps = cb.snapshots(d)
    if not snaps:
        raise ValueError("no snapshot output yet: run PFLOTRAN on this deck first")
    tt, p, k = min(snaps, key=lambda s: abs(s[0] - t))
    if abs(tt - t) > 0.01:
        raise ValueError(f"no snapshot at {t} y (nearest {tt})")
    vol = cb.cell_volumes(d)
    with h5py.File(p, "r") as f:
        g = f[k]
        key = next((n for n in g if re.match(r"Total[ _]CO2\(aq\)", n)), None)
        if key is None:
            raise ValueError("Total CO2(aq) not in the output: " + ", ".join(sorted(n for n in g if "CO2" in n)))
        c = np.asarray(g[key][:], float).ravel()
        phi = 1.0 - sum(np.asarray(g[n][:], float).ravel() for n in g if " VF" in n)
    return float((c * 1000.0 * vol * phi).sum()), key


def main():
    ref = {}
    p = os.path.join(REV, "results", "shutin_carbon_balance.csv")
    if os.path.isfile(p):
        ref = {r["network"]: r for r in csv.DictReader(open(p))}
    out = []
    for d in sorted(glob.glob(os.path.join(REV, "runs_verify_shutin", "E_feedback__p32_*"))):
        nm = os.path.basename(d)
        try:
            actual, key = dic(d)
        except Exception as e:
            print(f"  {nm}: {e}"); continue
        cols, rows = sb.series(d); w10 = sb.at(rows, sb.col(cols, r"Global Water Mass"), 10.0)
        ieq = w10 * sb.conc(d)
        c50 = float(ref[nm]["carbonate_50"]) if nm in ref else float("nan")
        r = dict(network=nm, dic_actual=actual, inj_equivalent=ieq, dic_over_injeq=actual / ieq,
                 carbonate_50=c50, carb50_over_dic=c50 / actual, carb50_over_injeq=c50 / ieq, dataset=key)
        out.append(r)
        print(f"  {nm}: DIC at shut-in {actual:.5g} mol, injectate-equivalent {ieq:.5g} mol, ratio {actual / ieq:.3f}; "
              f"carbonate at 50 y / DIC {c50 / actual:.3f} (against the estimate {c50 / ieq:.3f})")
    if not out:
        sys.exit("no verification runs with output at 10 y")
    q = os.path.join(REV, "results", "verify_shutin_dic.csv")
    with open(q, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print("written: results/verify_shutin_dic.csv")


if __name__ == "__main__":
    main()
