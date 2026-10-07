#!/usr/bin/env python3
"""
Flow-index test specified in advance, on the corrected coupled runs.

    python3 src/flow_index_test.py --root runs_gravityoff          # Block B (specified in advance) and Block A

For each coupled run directory (B_feedback__*, A_feedback__*):
  shape99        the a priori index: p99 / median of the mean groundwater age field,
                 from flowfield.solve() on that run's own mesh (full_mesh.uge and
                 boundary_right_e.ex), using shape99() from inlet_sensitivity.py
  carb_per_cell  the response specified in advance: mean over cells of the net carbonate
                 volume fraction at 50 y (seed subtracted, clipped at zero)
  efficiency     mineralized / injected CO2, % (from results/<Block>_coupled_networks.csv),
                 which does not depend on network volume
Spearman correlations: shape99 with each response, and shape99 with P32 (confounding).
Block B: one-sided, H1 rho > 0 (the direction specified in advance). Block A: two-sided, supporting.
Writes results/flow_index_test.csv.
"""
import argparse, csv, inspect, os, re, sys
import numpy as np, h5py
from scipy.stats import spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import flowfield

try:
    from inlet_sensitivity import shape99
    SOURCE = "inlet_sensitivity.shape99"
except Exception as e:                      # only if that module cannot be imported
    def shape99(age):
        a = np.asarray(age, float); a = a[np.isfinite(a) & (a > 0)]
        return float(np.percentile(a, 99) / np.median(a))
    SOURCE = f"local p99/median of finite ages (inlet_sensitivity not importable: {e})"


def carb_per_cell(d):
    t, p, k = cb.snapshots(d)[-1]
    with h5py.File(p, "r") as f:
        g = f[k]
        tot = 0.0
        for m in cb.CARB:
            n = next((n for n in g if " VF" in n and n.split(" VF")[0] == m), None)
            if n is not None:
                tot = tot + np.clip(np.asarray(g[n][:], float).ravel() - cb.SEED, 0, None)
    return float(np.mean(tot)), t


def one_sided(rho, p):
    return p / 2 if rho > 0 else 1 - p / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--blocks", nargs="+", default=["B", "A"])
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    print(f"index function: {SOURCE}")
    if SOURCE.startswith("inlet"):
        print("  " + inspect.getsource(shape99).strip().replace("\n", "\n  "))
    out = []
    for blk in a.blocks:
        effs = {}
        pcsv = os.path.join(REV, "results", f"{blk}_coupled_networks.csv")
        if os.path.isfile(pcsv):
            effs = {r["network"]: float(r["efficiency_pct"]) for r in csv.DictReader(open(pcsv))}
        else:
            print(f"  note: {pcsv} missing; run summarize_block.py for {blk} first (efficiency left out)")
        dirs = sorted(x for x in os.listdir(root) if x.startswith(f"{blk}_feedback__p32_"))
        rows = []
        for name in dirs:
            d = os.path.join(root, name)
            if not os.path.isfile(os.path.join(d, "boundary_right_e.ex")):
                print(f"  skipped {name}: no boundary_right_e.ex"); continue
            try:
                sol = flowfield.solve(d, verbose=False)
                s = float(shape99(sol["age"]))
                cpc, t = carb_per_cell(d)
            except Exception as e:
                print(f"  skipped {name}: {e}"); continue
            if t < 49.99:
                print(f"  skipped {name}: output reaches only {t:g} y"); continue
            m = re.search(r"p32_(\d{3})_s(\d+)", name)
            rows.append(dict(block=blk, network=name, p32_factor=int(m.group(1)) / 100, seed=int(m.group(2)),
                             shape99=s, carb_per_cell=cpc, efficiency_pct=effs.get(name, np.nan)))
            print(f"  {name:30s} shape99 {s:8.3f}   carbonate per cell {cpc:.3e}   efficiency {effs.get(name, np.nan):.4g}%")
        if len(rows) < 4:
            print(f"{blk}: too few networks ({len(rows)})"); continue
        s = [r["shape99"] for r in rows]
        print(f"\nBlock {blk}: {len(rows)} networks  ({'test specified in advance, one-sided, H1 rho > 0' if blk == 'B' else 'supporting, two-sided'})")
        for lab, key in (("carbonate per cell (specified in advance)", "carb_per_cell"), ("efficiency", "efficiency_pct"),
                         ("P32 (confounding check)", "p32_factor")):
            v = np.array([r[key] for r in rows], float); ok = np.isfinite(v)
            if ok.sum() < 4:
                continue
            rho = spearmanr(np.array(s)[ok], v[ok])
            pv = one_sided(rho.correlation, rho.pvalue) if blk == "B" and key != "p32_factor" else rho.pvalue
            print(f"  shape99 vs {lab:38s} rho = {rho.correlation:+.3f}   p = {pv:.3g}{' (one-sided)' if blk == 'B' and key != 'p32_factor' else ''}   n = {ok.sum()}")
        print()
        out += rows
    if out:
        p = os.path.join(REV, "results", "flow_index_test.csv")
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
        print(f"written: results/{os.path.basename(p)}")


if __name__ == "__main__":
    main()
