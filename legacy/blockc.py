#!/usr/bin/env python3
"""Block C sensitivity report: variant vs baseline, mapped to reviewer comments.

Each variant ran on an ARCHIVED case, so the baseline is that case's published
value and the mesh is identical between the two -- the per-cell normalisation
cancels and both metrics are directly comparable.
"""
from __future__ import annotations
import sys, csv, os
import numpy as np

COMMENTS = {
    "anor_as30":      "R3-6, R1-2c  anorthite reactive surface area 10 -> 30",
    "anor_as50":      "R3-6, R1-2c  anorthite reactive surface area 10 -> 50",
    "anor_as100":     "R3-6, R1-13  anorthite 10 -> 100 (convergence limit test)",
    "global_as_x0.1": "R1-2c        ALL surface areas x0.1 (Damkohler down)",
    "global_as_x10":  "R1-2c, R1-13 ALL surface areas x10 (Damkohler up)",
    "vf_consistent":  "R2-1         primary volume fractions rescaled to sum 1-phi",
    "no_dawsonite":   "R1-2a, R1-8  dawsonite removed from the secondary set",
    "seed_vf_lo":     "R3-5, R1-2b  secondary seed VF 1e-6 -> 1e-8",
    "seed_vf_hi":     "R3-5, R1-2b  secondary seed VF 1e-6 -> 1e-4",
    "sec_as_lo":      "R3-5, R1-2b  secondary surface area 1 -> 0.1",
    "sec_as_hi":      "R3-5, R1-2b  secondary surface area 1 -> 10",
    "feedback":       "R2-4, R3-4   porosity-permeability feedback enabled",
}

def load(path, key, val):
    out = {}
    if not os.path.exists(path): return out
    for r in csv.DictReader(open(path)):
        try: out[r[key]] = float(r[val])
        except (KeyError, ValueError, TypeError): pass
    return out

def main():
    vw = sys.argv[1] if len(sys.argv) > 1 else "volweighted.csv"
    bl = sys.argv[2] if len(sys.argv) > 2 else "../betweenness_results.csv"
    base_pc = load(bl, "name", "carb_per_cell_final")
    rows = list(csv.DictReader(open(vw)))
    vwi = {r["run_id"]: float(r["carb_intensity"]) for r in rows}
    vwp = {r["run_id"]: float(r["carb_per_cell"]) for r in rows}

    def base_intensity(case, run_id):
        if run_id in vwp and vwp[run_id] > 0 and case in base_pc:
            return base_pc[case] * (vwi[run_id] / vwp[run_id])
        return None

    variants = {}
    for r in rows:
        rid = r["run_id"]
        if not rid.startswith("C_"): continue
        v, case = rid[2:].split("__", 1)
        variants.setdefault(v, []).append((case, rid))

    print("=" * 96)
    print("BLOCK C SENSITIVITIES -- variant vs published baseline, same mesh")
    print("=" * 96)
    for v in sorted(variants, key=lambda k: list(COMMENTS).index(k) if k in COMMENTS else 99):
        print(f"\n{v}\n  {COMMENTS.get(v,'unmapped')}")
        print(f"    {'case':<16}{'baseline':>12}{'variant':>12}{'ratio':>9}"
              f"{'baseline_int':>14}{'variant_int':>13}")
        ratios = []
        for case, rid in sorted(variants[v]):
            b, x = base_pc.get(case), vwp.get(rid)
            bi, xi = base_intensity(case, rid), vwi.get(rid)
            if b is None or x is None:
                print(f"    {case:<16}{'n/a':>12}{'n/a':>12}"); continue
            if b > 0: ratios.append(x/b); rs = f"x{x/b:.2f}"
            else: rs = "base=0"
            print(f"    {case:<16}{b:>12.3e}{x:>12.3e}{rs:>9}"
                  f"{(f'{bi:.3e}' if bi else 'n/a'):>14}"
                  f"{(f'{xi:.3e}' if xi is not None else 'n/a'):>13}")
        if ratios:
            lo, hi = min(ratios), max(ratios)
            print(f"    -> effect on trapping: x{lo:.2f} to x{hi:.2f}"
                  f"   ({'within a factor of 2' if hi < 2 and lo > 0.5 else 'EXCEEDS a factor of 2'})")
    missing = [v for v in COMMENTS if v not in variants]
    if missing:
        print(f"\n{'='*96}\nNO USABLE OUTPUT -- these comments remain unanswered:")
        for v in missing: print(f"  {v:<16} {COMMENTS[v]}")

if __name__ == "__main__":
    main()
