#!/usr/bin/env python3
"""Pore-volume-weighted carbonate metrics, replacing the per-cell average.

compute_betweenness.py:169,172 reports sum(VF)/n_cells -- an UNWEIGHTED cell
average. DFN meshes grade cell size by distance from fracture intersections, so
cells span orders of magnitude in volume and that average weights a sliver
equally with a large interior cell. It moves when the mesh moves.

Evidence: Block B's x0.90 runs (~71k cells) report HIGHER carbonate per cell than
x1.75 (~206k cells) -- the reverse of the manuscript's P32 trend. That is the
divisor, not chemistry.

Computes:
  carb_per_cell   sum(VF)/N                      as published, for comparison
  carb_volume_m3  sum(VF_i * V_i)                absolute
  carb_intensity  sum(VF_i*V_i)/sum(phi*V_i)     mesh-independent
"""
from __future__ import annotations
import os, sys, csv, glob, argparse
import numpy as np

CARBONATES = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED_VF = 1e-6
POROSITY = 0.50

def uge_volumes(path):
    with open(path) as f:
        hdr = f.readline().split()
        if hdr[0].upper() != "CELLS":
            raise ValueError(f"{path}: expected 'CELLS n'")
        n = int(hdr[1]); vol = np.empty(n)
        for i in range(n):
            vol[i] = float(f.readline().split()[4])
    return vol

def time_groups(f):
    out = []
    for k in f.keys():
        if "Time" not in k or "failure" in k.lower() or "cut" in k.lower():
            continue
        parts = k.strip().split()
        try:
            ti = parts.index("Time"); t = float(parts[ti + 1])
            t_yr = t if (len(parts) > ti + 2 and parts[ti + 2] == "y") else t / 3.156e7
            out.append((t_yr, k))
        except Exception:
            continue
    out.sort(); return out

def analyse(run_dir, h5py):
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5"))) if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "missing .h5 or full_mesh.uge"
    vol = uge_volumes(uge)
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, "no parseable Time groups"
        t_end, key = tg[-1]; grp = f[key]
        per_min, total_vf, total_vol, ncells = {}, 0.0, 0.0, None
        for m in CARBONATES:
            ds = next((k for k in grp.keys() if k.startswith(f"{m} VF")), None)
            if ds is None:
                per_min[m] = 0.0; continue
            arr = np.asarray(grp[ds][:], dtype=float)
            if ncells is None:
                ncells = arr.size
                if arr.size != vol.size:
                    return None, (f"length mismatch: h5 {arr.size} cells vs "
                                  f".uge {vol.size} -- cannot pair volumes")
            net = np.clip(arr - SEED_VF, 0.0, None)
            per_min[m] = float((net * vol).sum())
            total_vf += float(net.sum()); total_vol += per_min[m]
    pore = float((vol * POROSITY).sum())
    return dict(run_id=os.path.basename(run_dir.rstrip("/")),
                t_end_yr=round(t_end, 3), n_cells=ncells,
                mesh_volume_m3=float(vol.sum()), pore_volume_m3=pore,
                carb_per_cell=total_vf / ncells,
                carb_volume_m3=total_vol,
                carb_intensity=total_vol / pore,
                **{f"vol_{m.lower()}_m3": per_min[m] for m in CARBONATES}), None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--out", default="volweighted.csv")
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    try:
        import h5py
    except ImportError:
        sys.exit("FATAL: h5py not available")
    dirs = sorted(d for d in glob.glob(os.path.join(a.runs, "*")) if os.path.isdir(d))
    rows, skipped = [], []
    for d in dirs:
        r, why = analyse(d, h5py)
        (rows.append(r) if r else skipped.append((os.path.basename(d), why)))
        if r:
            print(f"  {r['run_id']:<34} cells={r['n_cells']:>7,}  "
                  f"per_cell={r['carb_per_cell']:.3e}  "
                  f"intensity={r['carb_intensity']:.3e}", flush=True)
    if skipped:
        print(f"\n  skipped {len(skipped)}:")
        for n, w in skipped[:10]: print(f"    {n}: {w}")
    if not rows: sys.exit("no runs analysed")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {a.out} ({len(rows)} runs)")
    if a.compare:
        from scipy import stats
        pc = np.array([r["carb_per_cell"] for r in rows])
        ci = np.array([r["carb_intensity"] for r in rows])
        nc = np.array([r["n_cells"] for r in rows], dtype=float)
        print("\n=== do the two metrics rank runs the same way? ===")
        print(f"  per_cell vs intensity   rho={stats.spearmanr(pc,ci).statistic:+.3f}")
        print(f"  per_cell vs n_cells     rho={stats.spearmanr(pc,nc).statistic:+.3f}")
        print(f"  intensity vs n_cells    rho={stats.spearmanr(ci,nc).statistic:+.3f}")

if __name__ == "__main__":
    main()
