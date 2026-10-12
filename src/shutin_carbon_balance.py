#!/usr/bin/env python3
"""
Carbon balance of the shut-in runs (evolving aperture, injection shut in at 10 years),
from the PFLOTRAN mass-balance files (all parts of continued runs).

    python3 src/shutin_carbon_balance.py --root runs_gravityoff --prefix E_feedback__ [--t-shut 10]

Per network, in mol of carbon:
  injected          injected water up to shut-in x injectate CO2(aq) concentration of the deck
  carbonate_shut    carbonate formed by shut-in (calcite + magnesite + siderite + dawsonite,
                    net of the value at t = 0)
  carbonate_50      carbonate formed by 50 years
  after_shut        carbonate formed after shut-in (carbonate_50 - carbonate_shut)
  inj_equivalent    pore-water mass at shut-in x injectate concentration
  dissolved_shut    dissolved carbon in the domain at shut-in, if a total is written
  exported          injected - carbonate_shut - dissolved_shut (carbon that left before shut-in)
Ratios: carbonate_50 / inj_equivalent (the 108 % of the letter), after_shut / inj_equivalent,
after_shut / dissolved_shut, carbonate_50 / injected (efficiency).
Prints the carbon-related mass-balance columns found. Writes results/shutin_carbon_balance.csv.
"""
import argparse, csv, glob, os, re, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")


def series(d):
    """(columns, rows) over all mass-balance parts, rows sorted by time, later parts win."""
    cols, byt = None, {}
    for p in cb._parts(d, "pflotran_co2-mas", ".dat"):
        c, rows = cb._mas_rows(p)
        cols = cols or c
        for r in rows:
            vals = [float(x) for x in r]
            byt[round(vals[0], 8)] = dict(zip(c, vals))
    return cols, [byt[t] for t in sorted(byt)]


def col(cols, *pats):
    for p in pats:
        k = next((c for c in cols if re.search(p, c, re.I)), None)
        if k:
            return k
    return None


def at(rows, key, t):
    r = min(rows, key=lambda r: abs(list(r.values())[0] - t))
    return r.get(key, float("nan")) if key else float("nan")


def conc(d):
    m = re.search(r"CONSTRAINT\s+co2_rich_water.*?CO2\(aq\)\s+([0-9.eEdD+-]+)", open(os.path.join(d, "pflotran_co2.in")).read(), re.S | re.I)
    return float(m.group(1).lower().replace("d", "e"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff"); ap.add_argument("--prefix", default="E_feedback__")
    ap.add_argument("--t-shut", type=float, default=10.0)
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    out, shown = [], False
    for d in sorted(glob.glob(os.path.join(root, a.prefix + "p32_*"))):
        try:
            cols, rows = series(d)
        except Exception as e:
            print(f"  skipped {os.path.basename(d)}: {e}"); continue
        if not rows:
            continue
        if not shown:
            print("carbon-related mass-balance columns:")
            for c in cols:
                if re.search(r"CO2|HCO3|CO3|carbon|Calcite|Magnesite|Siderite|Dawsonite|injector|Water Mass", c, re.I):
                    print("   " + c)
            shown = True
        k_inj = col(cols, r"injector.*Water"); k_w = col(cols, r"Global Water Mass")
        k_diss = col(cols, r"Global.*Total CO2\(aq\)", r"Global.*CO2\(aq\).*\[mol\]", r"Global.*HCO3-.*\[mol\]")
        k_carb = [col(cols, rf"Region all {m} Total Mass") for m in CARB]
        c0 = conc(d); t0 = list(rows[0].values())[0]
        carb = lambda t: sum(at(rows, k, t) - at(rows, k, t0) for k in k_carb if k)
        inj = at(rows, k_inj, a.t_shut) * c0
        cs, c50 = carb(a.t_shut), carb(50.0)
        ieq = at(rows, k_w, a.t_shut) * c0
        diss = at(rows, k_diss, a.t_shut)
        row = dict(network=os.path.basename(d), injected=inj, carbonate_shut=cs, carbonate_50=c50, after_shut=c50 - cs,
                   inj_equivalent=ieq, dissolved_shut=diss, exported=inj - cs - diss if np.isfinite(diss) else float("nan"),
                   carb50_over_injeq=c50 / ieq, after_over_injeq=(c50 - cs) / ieq,
                   after_over_dissolved=(c50 - cs) / diss if np.isfinite(diss) and diss > 0 else float("nan"),
                   efficiency_pct=100 * c50 / inj)
        out.append(row)
    if not out:
        sys.exit("no runs")
    med = lambda k: np.nanmedian([r[k] for r in out]); rng = lambda k: (np.nanmin([r[k] for r in out]), np.nanmax([r[k] for r in out]))
    print(f"\n{len(out)} networks (medians, ranges)")
    for k, lab in (("injected", "carbon injected to shut-in (mol)"), ("carbonate_shut", "carbonate at shut-in (mol)"), ("carbonate_50", "carbonate at 50 years (mol)"),
                   ("after_shut", "carbonate formed after shut-in (mol)"), ("inj_equivalent", "injectate-equivalent carbon of the pore water (mol)"),
                   ("dissolved_shut", "dissolved carbon at shut-in (mol)"), ("exported", "carbon exported before shut-in (mol)"),
                   ("carb50_over_injeq", "carbonate at 50 y / injectate-equivalent"), ("after_over_injeq", "carbonate after shut-in / injectate-equivalent"),
                   ("after_over_dissolved", "carbonate after shut-in / dissolved carbon at shut-in"), ("efficiency_pct", "efficiency, carbonate at 50 y / injected (%)")):
        lo, hi = rng(k)
        print(f"  {lab:54s} {med(k):12.5g}   ({lo:.4g} to {hi:.4g})")
    p = os.path.join(REV, "results", "shutin_carbon_balance.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print("written: results/shutin_carbon_balance.csv")


if __name__ == "__main__":
    main()
