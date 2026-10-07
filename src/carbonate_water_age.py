#!/usr/bin/env python3
"""
Cell-level check of the mechanism behind the network index: does the carbonate
sit in the oldest water of each network?

    python3 src/carbonate_water_age.py --root runs_gravityoff            # Blocks B and A

For each coupled run (B_feedback__*, A_feedback__*), the steady-state mean
groundwater age of every cell comes from flowfield.solve() on the run's own mesh
(uniform permeability, as for the index), and the carbonate from the 50 y output.
Per network:
  age_pct_carb     carbonate-weighted mean of each cell's age percentile, where the
                   percentile is volume-weighted within the network: 50 = no preference,
                   above 50 = carbonate sits in older water
  share_oldest10   share of the carbonate in the oldest 10% of the pore volume (10 = no preference)
  share_stagnant   share of the carbonate in stagnant cells (no outflow). solve() sets
                   their age to zero; here they are treated as the oldest water.
Writes results/carbonate_water_age.csv.
"""
import argparse, csv, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import flowfield


def carbonate_mol(d, vol):
    t, p, k = cb.snapshots(d)[-1]
    mv = cb.molar_volumes(d)
    with h5py.File(p, "r") as f:
        g = f[k]; mol = np.zeros(len(vol))
        for m in cb.CARB:
            n = next((n for n in g if " VF" in n and n.split(" VF")[0] == m), None)
            if n is not None:
                mol += np.clip(np.asarray(g[n][:], float).ravel() - cb.SEED, 0, None) * vol / mv[m]
    return mol, t


def measures(age, stag, vol, mol):
    old = np.where(stag | ~np.isfinite(age), np.inf, age)       # stagnant or undefined: oldest
    o = np.argsort(old, kind="stable")
    pv = vol / vol.sum()
    pct = np.empty(len(vol)); pct[o] = 100 * (np.cumsum(pv[o]) - pv[o] / 2)
    w = mol / mol.sum()
    oldest = np.zeros(len(vol), bool); oldest[o[np.searchsorted(np.cumsum(pv[o]), 0.9):]] = True
    return float((w * pct).sum()), float(100 * w[oldest].sum()), float(100 * w[stag].sum()), float(100 * pv[stag].sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--blocks", nargs="+", default=["B", "A"])
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    out = []
    for blk in a.blocks:
        rows = []
        for name in sorted(x for x in os.listdir(root) if x.startswith(f"{blk}_feedback__p32_")):
            d = os.path.join(root, name)
            try:
                sol = flowfield.solve(d, verbose=False)
                vol = np.asarray(sol["vol"], float)
                mol, t = carbonate_mol(d, vol)
                if t < 49.99 or mol.sum() <= 0:
                    print(f"  skipped {name}: output at {t:g} y, carbonate {mol.sum():.3g} mol"); continue
                ap_, s10, sst, vst = measures(np.asarray(sol["age"], float), np.asarray(sol["stagnant"], bool), vol, mol)
            except Exception as e:
                print(f"  skipped {name}: {e}"); continue
            m = re.search(r"p32_(\d{3})_s(\d+)", name)
            rows.append(dict(block=blk, network=name, p32_factor=int(m.group(1)) / 100, age_pct_carb=ap_,
                             share_oldest10=s10, share_stagnant=sst, stagnant_vol_pct=vst))
            print(f"  {name:30s} carbonate at age percentile {ap_:5.1f}   in oldest 10%: {s10:5.1f}%   "
                  f"in stagnant cells: {sst:5.1f}% (of {vst:4.1f}% of the volume)")
        if rows:
            med = lambda k: np.median([r[k] for r in rows])
            print(f"\nBlock {blk}: {len(rows)} networks, medians: carbonate-weighted age percentile {med('age_pct_carb'):.1f} "
                  f"(50 = none), share in the oldest 10% of the volume {med('share_oldest10'):.1f}% (10 = none), "
                  f"in stagnant cells {med('share_stagnant'):.1f}% of the carbonate in {med('stagnant_vol_pct'):.1f}% of the volume\n")
        out += rows
    if out:
        p = os.path.join(REV, "results", "carbonate_water_age.csv")
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
        print(f"written: results/{os.path.basename(p)}")


if __name__ == "__main__":
    main()
