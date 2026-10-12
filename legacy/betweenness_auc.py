#!/usr/bin/env python3
"""
Within-realisation separation of precipitating cells in node betweenness.

The reduction reported in Table 8, the mean node betweenness over precipitating
cells, is a between-network average correlated against a between-network total.
Reviewer 3 objected that this is a property measured on the dependent
variable's own support, and recomputation on the corrected runs bears that out:
the correlation falls from 0.865 to 0.357 across the 46 precipitating
realisations, and to 0.279 (p = 0.067) once the two trace cases are also
dropped. The association with the precipitating-cell count is 0.313 (p = 0.027),
which is the mechanical effect the reviewer predicted.

This computes a different statistic that does not have that defect. Within each
realisation, it asks how completely the betweenness of precipitating cells
separates from the betweenness of the remainder,

    AUC = P(betweenness_precipitating > betweenness_other),

the same quantity reported for the chemical fields. The selection is still on
the outcome, as it is for pH or cation concentration, but nothing is averaged
and nothing is correlated against a network total, so cluster size does not
enter. A value near 0.5 means precipitation is topologically indifferent, in
which case the connectivity statement should be withdrawn rather than restated.

Requires the cell graph, so it repeats the mesh read and the betweenness solve
from compute_betweenness.py. Roughly 20 s per realisation at k = 50.

    python3 src/betweenness_auc.py                 # all Block A
    python3 src/betweenness_auc.py --dfn A_p32_100_s1181
    python3 src/betweenness_auc.py --k 100         # more source samples
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py
import networkx as nx
from scipy.stats import mannwhitneyu

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
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


def carbonate_mask(run_dir):
    """Boolean mask of cells holding carbonate above the seed."""
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[final_group(f)]
        tot = None
        for p in CARB:
            k = next((x for x in g if x.startswith(p + " VF")), None)
            if k:
                v = np.maximum(np.asarray(g[k][:], float).flatten() - SEED, 0.0)
                tot = v if tot is None else tot + v
    if tot is None:
        return None
    return tot > 0


def read_uge_graph(path):
    """Cell connectivity graph from a PFLOTRAN .uge file."""
    with open(path) as fh:
        first = fh.readline().split()
        ncell = int(first[1])
        for _ in range(ncell):
            fh.readline()
        conn_line = fh.readline().split()
        nconn = int(conn_line[1])
        g = nx.Graph()
        g.add_nodes_from(range(ncell))
        for _ in range(nconn):
            parts = fh.readline().split()
            if len(parts) < 2:
                continue
            i, j = int(parts[0]) - 1, int(parts[1]) - 1
            g.add_edge(i, j)
    return g, ncell


def auc(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    try:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return np.nan, np.nan
    return u / (len(a) * len(b)), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn", help="single realisation")
    ap.add_argument("--k", type=int, default=50, help="source samples")
    ap.add_argument("--csv", help="write per-realisation values here")
    a = ap.parse_args()

    cases = ([a.dfn] if a.dfn else
             sorted(os.path.basename(p) for p in glob.glob("runs/A_p32_*")))

    rows = []
    for name in cases:
        d = os.path.join("runs", name)
        uge = os.path.join(d, "full_mesh.uge")
        if not os.path.isfile(uge) or not os.path.isfile(
                os.path.join(d, "pflotran_co2.h5")):
            continue
        try:
            prec = carbonate_mask(d)
            if prec is None or prec.sum() == 0 or prec.sum() == len(prec):
                print(f"  {name:24s} no usable precipitation mask")
                continue
            g, ncell = read_uge_graph(uge)
            if ncell != len(prec):
                print(f"  {name:24s} mesh {ncell} != field {len(prec)}")
                continue
            k = min(a.k, ncell)
            bet = nx.betweenness_centrality(g, k=k, seed=0, normalized=True)
            b = np.array([bet.get(i, 0.0) for i in range(ncell)])
            val, p = auc(b[prec], b[~prec])
            rows.append((name, val, p, int(prec.sum()), ncell))
            print(f"  {name:24s} AUC = {val:.3f}  p = {p:.4f}  "
                  f"{prec.sum()} of {ncell} cells")
        except Exception as e:
            print(f"  {name:24s} failed: {str(e)[:50]}")

    if not rows:
        print("no realisations processed")
        return 1

    v = np.array([r[1] for r in rows])
    s = sum(1 for r in rows if r[2] < 0.05)
    print()
    print(f"realisations: {len(rows)}")
    print(f"median AUC:   {np.median(v):.3f}")
    print(f"range:        {v.min():.3f} to {v.max():.3f}")
    print(f"significant:  {s} of {len(rows)} at p < 0.05")
    print()
    if abs(np.median(v) - 0.5) < 0.1:
        print("The median is close to 0.5, so precipitation is topologically")
        print("indifferent and the connectivity statement should be withdrawn.")
    elif np.median(v) > 0.5:
        print("Precipitating cells sit at higher betweenness than the remainder,")
        print("which supports the statement in Section 4.1 that precipitation")
        print("localizes at well-connected cells.")
    else:
        print("Precipitating cells sit at LOWER betweenness than the remainder.")
        print("This contradicts the statement in Section 4.1 that precipitation")
        print("localizes at well-connected cells, and agrees with the specific")
        print("discharge result: precipitation occurs off the main flow corridors.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["name", "auc", "p", "n_precip", "n_cells"])
            w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
