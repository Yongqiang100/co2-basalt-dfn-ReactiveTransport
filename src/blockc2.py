#!/usr/bin/env python3
"""Block C sensitivity report against the CORRECTED baselines.

Supersedes blockc.py, which compared corrected variants against the published
values in betweenness_results.csv. Those were computed on uncorrected geometry
(cell volumes were fracture areas, ~1000x too large) with MASS_RATE applied per
cell (236x to 782,061x the stated rate), so that comparison is invalid.

Every variant runs on an archived case for which C_baseline__<case> now exists
with identical mesh and identical corrections, differing only in the variant op.
That is the only valid reference.

Reports per-cell average (as the manuscript does) and volume-weighted intensity,
plus the dominant carbonate phase in each, since the corrected baselines show
phase dominance varying case to case rather than with fracture intensity.
"""
from __future__ import annotations
import sys, os, glob, csv, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED, PHI = 1e-6, 0.50

COMMENTS = {
    "anor_as30":      "R3-6, R1-2c  anorthite A_s 10 -> 30",
    "anor_as50":      "R3-6, R1-2c  anorthite A_s 10 -> 50",
    "anor_as100":     "R3-6, R1-13  anorthite A_s 10 -> 100",
    "global_as_x0.1": "R1-2c        all A_s x0.1",
    "global_as_x10":  "R1-2c, R1-13 all A_s x10",
    "vf_consistent":  "R2-1         primary VF rescaled to 1-phi",
    "no_dawsonite":   "R1-2a, R1-8  dawsonite removed",
    "seed_vf_lo":     "R3-5, R1-2b  secondary seed VF 1e-6 -> 1e-8",
    "seed_vf_hi":     "R3-5, R1-2b  secondary seed VF 1e-6 -> 1e-4",
    "sec_as_lo":      "R3-5, R1-2b  secondary A_s 1 -> 0.1",
    "sec_as_hi":      "R3-5, R1-2b  secondary A_s 1 -> 10",
    "feedback":       "R2-4, R3-4   porosity-permeability feedback",
    "feedback_vfc":   "R2-4, R3-4, R2-1  feedback + corrected VF",
}

def read_run(d, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(d, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(d, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "no output"
    with open(uge) as f:
        n = int(f.readline().split()[1])
        vol = np.array([float(f.readline().split()[4]) for _ in range(n)])
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, f"incomplete ({tg[-1][0]:g} yr)"
        g = f[tg[-1][1]]
        per_vf, per_vol = {}, {}
        for m in CARB:
            k = next((x for x in g if x.startswith(f"{m} VF")), None)
            if k is None:
                per_vf[m] = per_vol[m] = 0.0; continue
            a = np.clip(g[k][:] - SEED, 0.0, None)
            if a.size != vol.size:
                return None, f"length mismatch {a.size} vs {vol.size}"
            per_vf[m] = float(a.sum()); per_vol[m] = float((a * vol).sum())
    pore = float((vol * PHI).sum())
    return dict(n_cells=n, per_cell=sum(per_vf.values()) / n,
                intensity=sum(per_vol.values()) / pore,
                per_vol=per_vol, total_vol=sum(per_vol.values())), None

def dominant(per_vol):
    t = sum(per_vol.values())
    if t <= 0: return "-", 0.0
    m = max(per_vol, key=per_vol.get)
    return m[:3].lower(), 100 * per_vol[m] / t

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--csv")
    a = ap.parse_args()

    base = {}
    for d in sorted(glob.glob(os.path.join(a.runs, "C_baseline__*"))):
        case = os.path.basename(d).split("__")[1]
        r, why = read_run(d)
        if r: base[case] = r
        else: print(f"  baseline {case}: {why}", file=sys.stderr)
    if not base:
        sys.exit("no corrected baselines found -- run C_baseline first")
    print(f"reference: {len(base)} corrected baselines\n")

    variants = {}
    for d in sorted(glob.glob(os.path.join(a.runs, "C_*__*"))):
        v, case = os.path.basename(d)[2:].split("__", 1)
        if v == "baseline" or v.startswith(("corrected", "shutin", "fieldscale")):
            continue
        variants.setdefault(v, []).append((case, d))

    rows, order = [], list(COMMENTS)
    for v in sorted(variants, key=lambda k: order.index(k) if k in order else 99):
        print(f"{v}\n  {COMMENTS.get(v, '(unmapped)')}")
        print(f"    {'case':<16}{'base':>11}{'variant':>11}{'ratio':>8}"
              f"{'base_int':>11}{'var_int':>11}{'i_ratio':>8}  dominant")
        ratios = []
        for case, d in sorted(variants[v]):
            r, why = read_run(d)
            b = base.get(case)
            if r is None:
                print(f"    {case:<16}{'':>11}{'FAILED':>11}  {why}"); continue
            if b is None:
                print(f"    {case:<16}  no corrected baseline"); continue
            pr = r["per_cell"]/b["per_cell"] if b["per_cell"] > 0 else float("nan")
            ir = r["intensity"]/b["intensity"] if b["intensity"] > 0 else float("nan")
            if np.isfinite(pr): ratios.append(pr)
            dm, dp = dominant(r["per_vol"]); db, dbp = dominant(b["per_vol"])
            print(f"    {case:<16}{b['per_cell']:>11.3e}{r['per_cell']:>11.3e}"
                  f"{(f'x{pr:.2f}' if np.isfinite(pr) else 'base=0'):>8}"
                  f"{b['intensity']:>11.3e}{r['intensity']:>11.3e}"
                  f"{(f'x{ir:.2f}' if np.isfinite(ir) else '-'):>8}"
                  f"  {db} {dbp:.0f}% -> {dm} {dp:.0f}%")
            rows.append(dict(variant=v, case=case,
                             base_per_cell=b["per_cell"], var_per_cell=r["per_cell"],
                             ratio=pr, base_intensity=b["intensity"],
                             var_intensity=r["intensity"], intensity_ratio=ir,
                             base_dominant=db, var_dominant=dm))
        if ratios:
            lo, hi = min(ratios), max(ratios)
            flag = "within a factor of 2" if (hi <= 2.0 and lo >= 0.5) else "EXCEEDS a factor of 2"
            print(f"    -> effect: x{lo:.2f} to x{hi:.2f}   ({flag})")
        print()

    missing = [v for v in COMMENTS if v not in variants]
    if missing:
        print("=" * 78)
        print("NO USABLE OUTPUT -- these comments remain unanswered:")
        for v in missing: print(f"  {v:<16} {COMMENTS[v]}")

    if a.csv and rows:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {a.csv} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
