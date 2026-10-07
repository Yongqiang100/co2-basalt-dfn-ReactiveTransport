#!/usr/bin/env python3
"""
Traces the route water took to reach each carbonate cell, and asks whether the
dissolved cation content rises along that route.

An earlier version walked downstream from the inlet, taking the largest
outgoing flow at every junction. That traces the trunk of the flow network and
never enters the branches, so it passed through no carbonate cells at all: the
cells that precipitate carry little flow. The result described the walk, not
the rock.

This walks the other way. It starts at a cell holding carbonate and steps
upstream, each time to whichever neighbour supplies the most inflow, until it
reaches an inlet cell or has nowhere further to go. The path is then read
forwards, from the inlet end to the carbonate cell, and the cation content is
plotted against position along it.

If the cation content climbs steadily towards the carbonate cell, the route is
loading as the water travels and the high concentration is the end of a
gradient. If the content is flat until close to the cell, the concentration is
produced locally rather than along the path. If the content is high throughout,
the water arrived already loaded.

    python3 src/upstream_profile.py --dfn A_p32_100_s1181
    python3 src/upstream_profile.py --dfn A_p32_100_s1181 --targets 40 --show 3
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
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[final_group(f)]
        cat = None
        for s in CATION:
            k = next((x for x in g if x == s or x.startswith(s.split(" [")[0])), None)
            if k:
                v = np.asarray(g[k][:], float).flatten()
                cat = v if cat is None else cat + v
        k = next((x for x in g if x == "pH" or x.startswith("pH")), None)
        ph = np.asarray(g[k][:], float).flatten() if k else None
        carb = None
        for p in CARB:
            k = next((x for x in g if x.startswith(p + " VF")), None)
            if k:
                v = np.maximum(np.asarray(g[k][:], float).flatten() - SEED, 0.0)
                carb = v if carb is None else carb + v
    return cat, ph, carb


def build_inflow(up, dn, q, ncell):
    """For each cell, the neighbours that send water into it, and how much."""
    inflow = [[] for _ in range(ncell)]
    for e in range(len(q)):
        i, j, f = int(up[e]), int(dn[e]), float(q[e])
        if f >= 0:
            inflow[j].append((i, f))
        else:
            inflow[i].append((j, -f))
    return inflow


def walk_upstream(target, inflow, max_steps=6000):
    """Step to the largest supplier until the water has no further source."""
    path = [target]
    seen = {target}
    cur = target
    for _ in range(max_steps):
        src = [(i, f) for i, f in inflow[cur] if i not in seen]
        if not src:
            break
        i, _ = max(src, key=lambda x: x[1])
        cur = i
        seen.add(cur)
        path.append(cur)
    return np.array(path[::-1])          # inlet end first


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn")
    ap.add_argument("--targets", type=int, default=25,
                    help="how many carbonate cells to trace back from")
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--show", type=int, default=0,
                    help="print the profile of this many individual paths")
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
            cat, ph, carb = read_fields(d)
            if cat is None or carb is None:
                print(f"  {name:24s} fields missing")
                continue
            prec = np.where(carb > 0)[0]
            if len(prec) == 0:
                print(f"  {name:24s} no carbonate")
                continue
            sol = flowfield.solve(d, verbose=False)
            up = np.asarray(sol["up"], int)
            dn = np.asarray(sol["dn"], int)
            q = np.asarray(sol["flux"], float)
            inflow = build_inflow(up, dn, q, len(cat))

            # trace back from the cells holding the most carbonate
            order = prec[np.argsort(carb[prec])[::-1]]
            targets = order[:a.targets]

            rhos, lens, ratios, shown = [], [], [], 0
            for tcell in targets:
                path = walk_upstream(int(tcell), inflow)
                if len(path) < a.min_len:
                    continue
                c = cat[path]
                pos = np.arange(len(path)) / (len(path) - 1)
                ok = np.isfinite(c)
                if ok.sum() > 5:
                    r, _ = spearmanr(pos[ok], c[ok])
                    if np.isfinite(r):
                        rhos.append(r)
                lens.append(len(path))
                head = np.nanmedian(c[:max(1, len(c) // 4)])
                if head > 0:
                    ratios.append(float(c[-1] / head))
                if shown < a.show:
                    shown += 1
                    idx = np.linspace(0, len(path) - 1, 8).astype(int)
                    print(f"    path to cell {tcell}, {len(path)} cells")
                    print("      position " +
                          " ".join(f"{pos[i]:7.2f}" for i in idx))
                    print("      cation   " +
                          " ".join(f"{c[i]:7.1e}" for i in idx))
                    if ph is not None:
                        print("      pH       " +
                              " ".join(f"{ph[path[i]]:7.2f}" for i in idx))

            if not rhos:
                print(f"  {name:24s} no usable paths")
                continue
            r = np.array(rhos)
            print(f"  {name:24s} paths={len(rhos):3d}  "
                  f"median length={int(np.median(lens)):4d} cells  "
                  f"median rho={np.median(r):+.3f}  "
                  f"rising in {(r > 0).mean()*100:3.0f}%  "
                  f"end/start ratio={np.median(ratios) if ratios else float('nan'):.1f}")
            summary.append((name, len(rhos), float(np.median(lens)),
                            float(np.median(r)), float((r > 0).mean()),
                            float(np.median(ratios)) if ratios else float("nan")))
        except Exception as e:
            print(f"  {name:24s} failed: {str(e)[:60]}")

    if not summary:
        return 1
    r = np.array([s[3] for s in summary])
    fr = np.array([s[4] for s in summary])
    ra = np.array([s[5] for s in summary])
    print()
    print(f"realisations: {len(summary)}")
    print(f"median rho along the route: {np.median(r):+.3f}")
    print(f"routes with a rising trend: {np.median(fr)*100:.0f}%")
    print(f"cation content at the carbonate cell over the inlet quarter: "
          f"{np.nanmedian(ra):.1f} times")
    print()
    if np.median(r) > 0.3:
        print("The cation content rises along the route to the carbonate cell,")
        print("so the high concentration is the end of a gradient rather than")
        print("an isolated spike.")
    elif np.median(r) < -0.3:
        print("The cation content falls along the route, so the water arrives")
        print("already loaded and loses cations on the way.")
    else:
        print("The cation content shows no consistent trend along the route.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["name", "paths", "median_length", "median_rho",
                        "fraction_rising", "end_over_start"])
            w.writerows(summary)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
