#!/usr/bin/env python3
"""Audit every realisation for precipitation exceeding the available pore space.

The fixed-porosity model never updates porosity, so nothing prevents a cell's
secondary volume fraction from exceeding unity. Where precipitation is spatially
confined this happens, and the reported trapping in those cells is not physical.

  seal_share  fraction of carbonate VOLUME in cells where closure exceeds
              opening -- cells a coupled calculation would seal
  over_share  fraction in cells whose secondary VF exceeds 1.0 -- cells holding
              more than their own volume

seal_share bounds how much a coupled model would change. over_share is the
harder limit: that part of the result is already invalid.
"""
from __future__ import annotations
import sys, os, glob, csv, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PRIM = ("Anorthite","Albite","Diopside","Forsterite","Fayalite","Enstatite")
CARB = ("Calcite","Magnesite","Siderite","Dawsonite")
SEC  = CARB + ("Kaolinite","Chalcedony")

def audit(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir,"*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge): return None, "no output"
    with open(uge) as f:
        n = int(f.readline().split()[1])
        vol = np.array([float(f.readline().split()[4]) for _ in range(n)])
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg: return None, "no time groups"
        if tg[-1][0] < final_year*0.999: return None, f"incomplete ({tg[-1][0]:g} yr)"
        g0, g1 = f[tg[0][1]], f[tg[-1][1]]
        fld = lambda g,p: next((np.asarray(g[k][:],float) for k in g
                                if k.startswith(p)), None)
        opened = np.zeros(n)
        for m in PRIM:
            a0, a1 = fld(g0,f"{m} VF"), fld(g1,f"{m} VF")
            if a0 is not None and a1 is not None: opened += a0 - a1
        carb = np.zeros(n); sec = np.zeros(n)
        for m in SEC:
            a = fld(g1, f"{m} VF")
            if a is None: continue
            net = np.clip(a-1e-6, 0.0, None); sec += net
            if m in CARB: carb += net
    tot = float((carb*vol).sum())
    if tot <= 0:
        return dict(n_cells=n, n_carb=0, carb_volume=0.0, seal_share=None,
                    over_share=None, n_seal=0, n_over=0,
                    max_sec_vf=float(sec.max()), max_carb_vf=float(carb.max())), None
    seal = (carb - opened) > 0
    over = sec > 1.0
    return dict(n_cells=n, n_carb=int((carb>0).sum()), carb_volume=tot,
                n_seal=int(seal.sum()), n_over=int(over.sum()),
                seal_share=100*float((carb[seal]*vol[seal]).sum())/tot,
                over_share=100*float((carb[over]*vol[over]).sum())/tot,
                max_sec_vf=float(sec.max()), max_carb_vf=float(carb.max())), None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_domainscaled")
    ap.add_argument("--prefix", default="C_baseline__")
    ap.add_argument("--csv")
    a = ap.parse_args()
    dirs = sorted(d for d in glob.glob(os.path.join(a.runs, a.prefix+"*"))
                  if os.path.isdir(d))
    rows, skipped = [], []
    for d in dirs:
        r, why = audit(d)
        if r is None: skipped.append((os.path.basename(d), why)); continue
        r["run_id"] = os.path.basename(d); rows.append(r)
    if skipped:
        print(f"skipped {len(skipped)}:")
        for nm, why in skipped[:6]: print(f"  {nm}: {why}")
    if not rows: sys.exit("nothing audited")
    print(f"\n{'case':<26}{'cells':>9}{'carb':>6}{'seal':>6}{'>1.0':>6}"
          f"{'seal %':>9}{'over %':>9}{'max sec VF':>12}")
    for r in sorted(rows, key=lambda x: -(x["over_share"] or 0)):
        ss = f"{r['seal_share']:.1f}" if r["seal_share"] is not None else "-"
        os_ = f"{r['over_share']:.1f}" if r["over_share"] is not None else "-"
        print(f"{r['run_id'].replace(a.prefix,''):<26}{r['n_cells']:>9,}"
              f"{r['n_carb']:>6}{r['n_seal']:>6}{r['n_over']:>6}"
              f"{ss:>9}{os_:>9}{r['max_sec_vf']:>12.3f}")
    v = [r for r in rows if r["over_share"] is not None]
    if v:
        ov = np.array([r["over_share"] for r in v])
        sl = np.array([r["seal_share"] for r in v])
        print(f"\n{'='*74}")
        print(f"  realisations with carbonate         {len(v)} of {len(rows)}")
        print(f"  with any cell above VF 1.0          {sum(1 for r in v if r['n_over']>0)}")
        print(f"  invalid share  median {np.median(ov):.1f}%   range {ov.min():.1f}-{ov.max():.1f}%")
        print(f"  affected share median {np.median(sl):.1f}%   range {sl.min():.1f}-{sl.max():.1f}%")
        bad = [r["run_id"].replace(a.prefix,"") for r in v if r["over_share"] > 50]
        print(f"\n  over half the carbonate in cells holding more than their")
        print(f"  own volume: {len(bad)} realisation(s)")
        for b in bad: print(f"    {b}")
        print(f"\n  median over_share < 10%  -> report as a caveat")
        print(f"  median over_share > 50%  -> absolute amounts unsupported")
    if a.csv and rows:
        keys = sorted({k for r in rows for k in r})
        with open(a.csv,"w",newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
        print(f"\n  wrote {a.csv}")

if __name__ == "__main__":
    main()
