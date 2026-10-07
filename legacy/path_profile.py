#!/usr/bin/env python3
"""
Does the dissolved cation content rise as water travels along a flow path?

Comparing cells at different positions at one instant does not answer this: a
cell far from the inlet may be fed by a wide fracture that released little
material, and precipitation removes cations wherever it occurs. The comparison
across cells gave no relationship, median Spearman +0.011 over ten
realisations.

This follows individual paths instead. Starting from each inlet cell, the walk
steps to whichever neighbour receives the largest outgoing flow, and continues
until it reaches an outlet, revisits a cell, or runs out of onward flow. Along
each path it records the distance from the start, the dissolved Ca and Mg
content, the pH, and whether the cell holds carbonate.

The question is whether the cation content rises with distance along a single
path. It reports the Spearman coefficient for each path, the fraction of paths
with a rising trend, and the position along the path at which carbonate first
appears.

    python3 src/path_profile.py --dfn A_p32_100_s1181
    python3 src/path_profile.py --dfn A_p32_100_s1181 --paths 40
    python3 src/path_profile.py --csv paths.csv
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py
from scipy.stats import spearmanr

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
CATION = ["CaCO3(aq) [M]", "CaHCO3+ [M]", "MgCO3(aq) [M]", "MgHCO3+ [M]"]
SEED = 1.0e-6


def final_group(f):
    timed = []
    for k in f.keys():
        m = re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)
        if m:
            timed.append((float(m.group(1)), k))
    if not timed:
        raise KeyError("no time groups")
    return max(timed)[1]


def read_fields(run_dir):
    """Cation content, pH and the carbonate mask, one value per cell."""
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[final_group(f)]
        cat = None
        for s in CATION:
            k = next((x for x in g if x == s or x.startswith(s.split(" [")[0])), None)
            if k:
                v = np.asarray(g[k][:], float).flatten()
                cat = v if cat is None else cat + v
        ph = None
        k = next((x for x in g if x == "pH" or x.startswith("pH")), None)
        if k:
            ph = np.asarray(g[k][:], float).flatten()
        carb = None
        for p in CARB:
            k = next((x for x in g if x.startswith(p + " VF")), None)
            if k:
                v = np.maximum(np.asarray(g[k][:], float).flatten() - SEED, 0.0)
                carb = v if carb is None else carb + v
    return cat, ph, (None if carb is None else carb > 0)


def build_outflow(up, dn, q, ncell):
    """For each cell, the neighbours it sends water to and how much."""
    out = [[] for _ in range(ncell)]
    for e in range(len(q)):
        i, j, f = int(up[e]), int(dn[e]), float(q[e])
        if f >= 0:
            out[i].append((j, f))
        else:
            out[j].append((i, -f))
    return out


def walk(start, out, centroids, max_steps=4000):
    """Step downstream, always taking the largest outgoing flow."""
    path = [start]
    seen = {start}
    dist = [0.0]
    total = 0.0
    cur = start
    for _ in range(max_steps):
        nxt = [(j, f) for j, f in out[cur] if j not in seen]
        if not nxt:
            break
        j, _ = max(nxt, key=lambda x: x[1])
        if centroids is not None:
            total += float(np.linalg.norm(centroids[j] - centroids[cur]))
        else:
            total += 1.0
        cur = j
        seen.add(cur)
        path.append(cur)
        dist.append(total)
    return np.array(path), np.array(dist)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn")
    ap.add_argument("--paths", type=int, default=25, help="how many paths to walk")
    ap.add_argument("--min-len", type=int, default=20, help="shortest path to keep")
    ap.add_argument("--csv")
    a = ap.parse_args()

    sys.path.insert(0, "src")
    import flowfield

    cases = ([a.dfn] if a.dfn else
             sorted(os.path.basename(p) for p in glob.glob("runs/A_p32_*")))

    summary = []
    for name in cases:
        d = os.path.join("runs", name)
        if not os.path.isfile(os.path.join(d, "pflotran_co2.h5")):
            continue
        try:
            cat, ph, has_carb = read_fields(d)
            if cat is None or has_carb is None:
                print(f"  {name:24s} fields missing")
                continue
            sol = flowfield.solve(d, verbose=False)
            up = np.asarray(sol["up"], int)
            dn = np.asarray(sol["dn"], int)
            q = np.asarray(sol["flux"], float)
            ncell = len(cat)
            cen = None
            for key in ("centroids", "xyz", "cc"):
                if key in sol:
                    c = np.asarray(sol[key], float)
                    if c.ndim == 2 and len(c) == ncell:
                        cen = c
                    break
            out = build_outflow(up, dn, q, ncell)

            # start from the cells with the largest outgoing flow
            send = np.array([sum(f for _, f in out[i]) for i in range(ncell)])
            starts = np.argsort(send)[::-1][:a.paths * 4]

            rhos, onsets, kept = [], [], 0
            for s in starts:
                if kept >= a.paths:
                    break
                path, dist = walk(int(s), out, cen)
                if len(path) < a.min_len:
                    continue
                kept += 1
                c = cat[path]
                ok = np.isfinite(c) & np.isfinite(dist)
                if ok.sum() > 5:
                    r, _ = spearmanr(dist[ok], c[ok])
                    if np.isfinite(r):
                        rhos.append(r)
                cp = has_carb[path]
                if cp.any():
                    onsets.append(float(np.argmax(cp)) / len(path))

            if not rhos:
                print(f"  {name:24s} no usable paths")
                continue
            r = np.array(rhos)
            frac_up = float((r > 0).mean())
            on = np.array(onsets) if onsets else np.array([np.nan])
            print(f"  {name:24s} paths={kept:3d}  median rho={np.median(r):+.3f}"
                  f"  rising in {frac_up*100:3.0f}%"
                  f"  carbonate onset at {np.nanmedian(on):.2f} of path length")
            summary.append((name, kept, float(np.median(r)), frac_up,
                            float(np.nanmedian(on))))
        except Exception as e:
            print(f"  {name:24s} failed: {str(e)[:50]}")

    if not summary:
        return 1
    r = np.array([s[2] for s in summary])
    f = np.array([s[3] for s in summary])
    print()
    print(f"realisations: {len(summary)}")
    print(f"median of the per-realisation median rho: {np.median(r):+.3f}")
    print(f"paths with a rising cation trend: {np.median(f)*100:.0f}% (median "
          f"across realisations)")
    print()
    if np.median(r) > 0.3:
        print("Cation content rises along the flow path, which supports the")
        print("proposed mechanism.")
    elif np.median(r) < -0.3:
        print("Cation content falls along the flow path, contradicting the")
        print("proposed mechanism.")
    else:
        print("Cation content shows no consistent trend along the flow path, so")
        print("the proposed mechanism is not supported by this test.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["name", "paths", "median_rho", "fraction_rising",
                        "median_onset_fraction"])
            w.writerows(summary)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
