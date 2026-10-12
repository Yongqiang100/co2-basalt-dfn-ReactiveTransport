#!/usr/bin/env python3
"""Two checks on the porosity-permeability feedback question (R2-4, R3-4).

CHECK A -- is the converged feedback run physically valid?
  feedback_vfc__p32_150_s42 completed and reports 2143x the carbonate of its
  feedback-off twin. Before citing that, its final volume fraction field must
  be inspected: the fixed-porosity runs reach carbonate VF of 2.75 in some
  cells, which exceeds unity, so a large ratio may reflect a different
  unphysical state rather than a physical effect.

CHECK B -- how much of the reported trapping is in cells that would seal?
  Closure exceeds opening in the highest-carbonate cells, so a coupled model
  would seal them. If those cells carry a small share of the domain total the
  upper-bound caveat is minor; if they carry most of it, the reported trapping
  in the high-trapping realisations is substantially overstated.
"""
from __future__ import annotations
import sys, os, glob, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PRIM = ("Anorthite","Albite","Diopside","Forsterite","Fayalite","Enstatite")
CARB = ("Calcite","Magnesite","Siderite","Dawsonite")
SEC  = CARB + ("Kaolinite","Chalcedony")

def load(run_dir):
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
        t_end = tg[-1][0]; g0, g1 = f[tg[0][1]], f[tg[-1][1]]
        fld = lambda g,p: next((np.asarray(g[k][:],float) for k in g
                                if k.startswith(p)), None)
        return dict(n=n, vol=vol, t_end=t_end,
                    prim0={m: fld(g0,f"{m} VF") for m in PRIM},
                    prim1={m: fld(g1,f"{m} VF") for m in PRIM},
                    sec1={m: fld(g1,f"{m} VF") for m in SEC}), None

def check_a(run_dir, label):
    d, why = load(run_dir)
    if d is None: print(f"  {label}: {why}"); return
    print(f"\n  {label}   (t_end = {d['t_end']:g} yr, {d['n']:,} cells)")
    tot_sec = sum(np.clip(v-1e-6,0,None) for v in d["sec1"].values() if v is not None)
    tot_prim = sum(v for v in d["prim1"].values() if v is not None)
    solid = tot_prim + tot_sec
    print(f"    secondary VF   max {tot_sec.max():.4f}   "
          f"cells > 0.5: {int((tot_sec>0.5).sum())}   > 1.0: {int((tot_sec>1.0).sum())}")
    print(f"    total solid VF max {solid.max():.4f}   cells > 1.0: {int((solid>1.0).sum())}")
    if (tot_sec > 1.0).any():
        print("    *** UNPHYSICAL: secondary volume fraction exceeds unity")
        print("        the result cannot be cited as a physical effect")
    elif (solid > 1.0).any():
        print("    *** UNPHYSICAL: solid volume fraction exceeds unity")
    else:
        print("    volume fractions remain physical")
    for m in CARB:
        v = d["sec1"].get(m)
        if v is not None and v.max() > 1e-3:
            print(f"      {m:<10} max VF {v.max():.4f}")

def check_b(run_dir, label):
    d, why = load(run_dir)
    if d is None: print(f"  {label}: {why}"); return None
    vol = d["vol"]
    opened = sum((d["prim0"][m]-d["prim1"][m]) for m in PRIM
                 if d["prim0"][m] is not None and d["prim1"][m] is not None)
    carb = sum(np.clip(d["sec1"][m]-1e-6,0,None) for m in CARB
               if d["sec1"][m] is not None)
    net = opened - carb
    seal = net < 0
    tot = float((carb*vol).sum())
    sv = float((carb[seal]*vol[seal]).sum()) if seal.any() else 0.0
    frac = 100*sv/tot if tot > 0 else 0.0
    print(f"\n  {label}")
    print(f"    cells                              {d['n']:>12,}")
    print(f"    cells with any carbonate           {int((carb>0).sum()):>12,}")
    print(f"    cells where closure > opening      {int(seal.sum()):>12,}"
          f"   ({100*seal.mean():.4f}% of domain)")
    print(f"    carbonate volume, total            {tot:>12.4e} m3")
    print(f"    carbonate volume in sealing cells  {sv:>12.4e} m3")
    print(f"    -> share of trapping that would seal  {frac:>9.1f}%")
    if carb.max() > 1.0:
        n_over = int((carb>1.0).sum())
        v_over = float((carb[carb>1.0]*vol[carb>1.0]).sum())
        print(f"    cells with carbonate VF > 1        {n_over:>12,}"
              f"   carrying {100*v_over/tot:.1f}% of the total")
    return frac

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--alt", default="runs_domainscaled")
    a = ap.parse_args()
    def find(name):
        for root in (a.runs, a.alt, "runs_uncorrected"):
            p = os.path.join(root, name)
            if os.path.isdir(p): return p
        return None
    print("="*74); print("CHECK A  is the converged feedback run physically valid?"); print("="*74)
    for nm in ("C_feedback_vfc__p32_150_s42","C_vf_consistent__p32_150_s42",
               "C_baseline__p32_150_s42"):
        p = find(nm)
        check_a(p, nm) if p else print(f"  {nm}: not found")
    print("\n"+"="*74); print("CHECK B  what share of trapping is in cells that would seal?"); print("="*74)
    fracs = []
    for nm in ("C_baseline__p32_150_s117","C_baseline__p32_150_s42",
               "C_baseline__p32_200_s117","C_baseline__p32_100_s383"):
        p = find(nm)
        if p:
            f = check_b(p, nm)
            if f is not None: fracs.append(f)
        else: print(f"  {nm}: not found")
    if fracs:
        print(f"\n  share of trapping in sealing cells: {min(fracs):.1f}% to "
              f"{max(fracs):.1f}% across {len(fracs)} cases")

if __name__ == "__main__":
    main()
