#!/usr/bin/env python3
"""
Are the cells that precipitate carbonate more or less connected than the rest?

For each realisation this compares two groups of cells: those holding carbonate
at 50 years, and all the others. Three properties are compared.

  neighbour count      how many cells share a face with this cell. A cell at
                       the end of a branch has one neighbour; a cell in the
                       middle of a fracture has two; a cell at a fracture
                       intersection has more. This is counted exactly from the
                       mesh, with no approximation.

  share of total flow  the water passing through the cell each second, divided
                       by the total for the whole domain. Taken from the steady
                       Darcy solution.

  local residence time the time the cell takes to replace its own water, equal
                       to its pore volume divided by its throughput.

Each comparison is reported as the probability that a cell holding carbonate
has a higher value than a cell that does not. A probability of 0.5 means the
two groups are indistinguishable. Values above 0.5 mean the carbonate cells
have higher values, and below 0.5 that they have lower values.

The comparison is made within each realisation and then summarised across
realisations, so a realisation with many carbonate cells does not count for
more than one with few. Realisations with no carbonate at all are skipped,
because there is nothing to compare.

    python3 src/connectivity_auc.py
    python3 src/connectivity_auc.py --dfn A_p32_100_s1181
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py
from scipy.stats import mannwhitneyu

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6


def final_group(f):
    """The group for the last simulated time, found by reading the time out of
    each group name rather than sorting the names as text."""
    timed = []
    for k in f.keys():
        m = re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)
        if m:
            timed.append((float(m.group(1)), k))
    if not timed:
        raise KeyError("no time groups")
    return max(timed)[1]


def carbonate_cells(run_dir):
    """True for each cell holding carbonate above the seeded amount."""
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[final_group(f)]
        total = None
        for phase in CARB:
            key = next((x for x in g if x.startswith(phase + " VF")), None)
            if key:
                v = np.maximum(np.asarray(g[key][:], float).flatten() - SEED, 0.0)
                total = v if total is None else total + v
    return None if total is None else total > 0


def neighbour_counts(uge_path, n_expected):
    """How many cells share a face with each cell, read from the mesh file."""
    with open(uge_path) as fh:
        header = fh.readline().split()
        ncell = int(header[1])
        for _ in range(ncell):
            fh.readline()
        conn_header = fh.readline().split()
        nconn = int(conn_header[1])
        count = np.zeros(ncell, dtype=int)
        for _ in range(nconn):
            parts = fh.readline().split()
            if len(parts) < 2:
                continue
            i, j = int(parts[0]) - 1, int(parts[1]) - 1
            count[i] += 1
            count[j] += 1
    if ncell != n_expected:
        raise ValueError(f"mesh has {ncell} cells, field has {n_expected}")
    return count


def probability_higher(a, b):
    """Probability that a value drawn from a exceeds one drawn from b."""
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    try:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return np.nan, np.nan
    return u / (len(a) * len(b)), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn", help="one realisation instead of all of them")
    ap.add_argument("--csv", help="write the per-realisation values here")
    a = ap.parse_args()

    sys.path.insert(0, "src")
    import flowfield

    labels = ["neighbour count", "share of total flow", "local residence time"]
    values = {k: [] for k in labels}
    n_signif = {k: 0 for k in labels}
    rows = []

    cases = ([a.dfn] if a.dfn else
             sorted(os.path.basename(p) for p in glob.glob("runs/A_p32_*")))

    for name in cases:
        d = os.path.join("runs", name)
        if not os.path.isfile(os.path.join(d, "pflotran_co2.h5")):
            continue
        try:
            has_carb = carbonate_cells(d)
            if has_carb is None or has_carb.sum() == 0 \
                    or has_carb.sum() == len(has_carb):
                continue

            nbr = neighbour_counts(os.path.join(d, "full_mesh.uge"),
                                   len(has_carb))
            sol = flowfield.solve(d, verbose=False)
            q = np.abs(np.asarray(sol["flux"], float))
            up = np.asarray(sol["up"], int)
            dn = np.asarray(sol["dn"], int)
            vol = np.asarray(sol["vol"], float)

            throughput = np.zeros(len(has_carb))
            np.add.at(throughput, up, q)
            np.add.at(throughput, dn, q)
            throughput *= 0.5

            flow_share = throughput / max(throughput.sum(), 1e-30)
            residence = vol / np.maximum(throughput, 1e-30)

            row = [name]
            for label, field in zip(labels, [nbr, flow_share, residence]):
                ok = np.isfinite(field)
                pr, p = probability_higher(field[has_carb & ok],
                                           field[(~has_carb) & ok])
                if np.isfinite(pr):
                    values[label].append(pr)
                    if p < 0.05:
                        n_signif[label] += 1
                row.append(pr)
            rows.append(row)
            print(f"  {name:24s} " + "  ".join(
                f"{label.split()[0]}={v:.3f}" for label, v in zip(labels, row[1:])))
        except Exception as e:
            print(f"  {name:24s} skipped: {str(e)[:50]}")

    if not rows:
        print("nothing processed")
        return 1

    print()
    print(f"realisations compared: {len(rows)}")
    print()
    print("%-24s %26s %s" % ("property of the cell",
                             "P(carbonate cell higher)", "significant in"))
    print("-" * 74)
    for label in labels:
        v = np.array(values[label])
        if len(v) == 0:
            continue
        print("%-24s %26.3f %9d of %-3d"
              % (label, np.median(v), n_signif[label], len(v)))

    nbr_med = np.median(np.array(values["neighbour count"]))
    print()
    if nbr_med > 0.55:
        print("Cells holding carbonate have more neighbours than cells without,")
        print("so carbonate forms at the better connected parts of the network.")
    elif nbr_med < 0.45:
        print("Cells holding carbonate have fewer neighbours than cells without,")
        print("so carbonate forms at the less connected parts of the network.")
    else:
        print("The two groups have similar neighbour counts, so how well a cell")
        print("is connected does not by itself say whether carbonate forms there.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["name"] + labels)
            w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
