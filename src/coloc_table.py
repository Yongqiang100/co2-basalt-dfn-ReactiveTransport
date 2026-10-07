#!/usr/bin/env python3
"""
Effect size of each field at precipitating cells, on one common set.

Replaces the median-percentile-rank summary used earlier. Median rank
saturates: once several fields sit at 99.9 per cent it cannot separate a field
whose precipitating cells are always at the very top from one whose cells are
merely usually high, and it discards magnitude entirely.

The primary statistic here is the rank-biserial correlation, the effect size
that belongs with the Mann-Whitney test,

    r_rb = 2 * P(x_prec > x_other) - 1,

which runs from -1 to +1, is dimensionless, comparable across fields, and does
not saturate. A value near +1 means a precipitating cell almost always exceeds
a non-precipitating one. The median ratio is reported alongside for the
concentration fields, where magnitude is the point; it is omitted for pH, which
is already logarithmic.

The denominator is the 46 production realisations with at least one
precipitating cell. Four of the 50 produce none, so the comparison is undefined
for those. The carbonate mineral rate is deliberately absent: it is the
precipitation itself, so a statistic computed on precipitating cells is
definitional rather than evidential. Forsterite is reported but flagged,
because it dissolves to 97-100 per cent everywhere and so little remains
anywhere that its rank is near-random by construction.

    python3 src/coloc_table.py
    python3 src/coloc_table.py --tex
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

# label, fields summed, report a median ratio, footnote marker
FIELDS = [
    ("dissolved Ca and Mg",
     ["CaCO3(aq) [M]", "CaHCO3+ [M]", "MgCO3(aq) [M]", "MgHCO3+ [M]"], True, ""),
    ("dissolved carbonate ion",
     ["CO3-- [M]"], True, ""),
    ("dissolved inorganic carbon",
     ["HCO3- [M]", "CO3-- [M]", "CaCO3(aq) [M]", "MgCO3(aq) [M]"], True, ""),
    ("pH", ["pH"], False, ""),
    ("primary volume remaining, total",
     ["Anorthite VF", "Albite VF", "Diopside VF", "Forsterite VF",
      "Fayalite VF", "Enstatite VF"], True, ""),
    ("  anorthite", ["Anorthite VF"], True, ""),
    ("  albite", ["Albite VF"], True, ""),
    ("  diopside", ["Diopside VF"], True, ""),
    ("  enstatite", ["Enstatite VF"], True, ""),
    ("  fayalite", ["Fayalite VF"], True, ""),
    ("  forsterite", ["Forsterite VF"], True, "*"),
]


def final_group(f):
    """Last timestep by time, not by the alphabetical order of the names."""
    timed = []
    for k in f.keys():
        m = re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)
        if m:
            timed.append((float(m.group(1)), k))
    if not timed:
        raise KeyError("no time groups")
    return max(timed)[1]


def field(g, names):
    out = None
    for n in names:
        key = next((x for x in g if x == n or x.startswith(n.split(" [")[0])), None)
        if key is None:
            continue
        v = np.asarray(g[key][:], float).flatten()
        out = v if out is None else out + v
    return out


def rank_biserial(a, b):
    """2 * P(a > b) - 1, from the Mann-Whitney U statistic."""
    try:
        u, _ = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return np.nan, np.nan
    n1, n2 = len(a), len(b)
    auc = u / (n1 * n2)
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return 2.0 * auc - 1.0, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", action="store_true")
    a = ap.parse_args()

    rb = {lab: [] for lab, _, _, _ in FIELDS}
    ratio = {lab: [] for lab, _, _, _ in FIELDS}
    sig = {lab: 0 for lab, _, _, _ in FIELDS}
    n_read = 0

    for d in sorted(glob.glob("runs/A_p32_*")):
        path = os.path.join(d, "pflotran_co2.h5")
        if not os.path.isfile(path):
            continue
        try:
            with h5py.File(path, "r") as f:
                g = f[final_group(f)]
                carb = None
                for p in CARB:
                    k = next((x for x in g if x.startswith(p + " VF")), None)
                    if k:
                        v = np.maximum(np.asarray(g[k][:], float).flatten() - SEED, 0.0)
                        carb = v if carb is None else carb + v
                if carb is None or carb.sum() == 0:
                    continue
                prec = carb > 0
                if prec.sum() == 0 or prec.sum() == len(carb):
                    continue
                n_read += 1
                for lab, names, want_ratio, _ in FIELDS:
                    v = field(g, names)
                    if v is None or len(v) != len(carb):
                        continue
                    r, p = rank_biserial(v[prec], v[~prec])
                    if np.isfinite(r):
                        rb[lab].append(r)
                        if p < 0.05:
                            sig[lab] += 1
                    if want_ratio:
                        m2 = float(np.median(v[~prec]))
                        if m2 > 0:
                            ratio[lab].append(float(np.median(v[prec])) / m2)
        except (OSError, KeyError):
            continue

    print(f"precipitating realisations: {n_read} of 50 simulated\n")
    print("%-32s %9s %14s %14s" % ("field", "r_rb", "median ratio", "significant in"))
    print("-" * 74)
    rows = []
    for lab, _, want_ratio, mark in FIELDS:
        v = np.array(rb[lab])
        if len(v) == 0:
            print("%-32s %9s" % (lab, "no data"))
            continue
        med = float(np.median(v))
        rr = np.array(ratio[lab])
        rtxt = f"{np.median(rr):.1f}x" if want_ratio and len(rr) else "--"
        print("%-32s %+9.3f %14s %9d of %-3d" % (lab + mark, med, rtxt, sig[lab], len(v)))
        rows.append((lab + mark, med, rtxt, sig[lab], len(v)))

    if a.tex:
        print("\n% --- SI table body ---")
        for lab, med, rtxt, s, n in rows:
            print("%s & %+.3f & %s & %d of %d \\\\" % (lab.strip(), med, rtxt, s, n))

    print("\n* forsterite dissolves to 97--100 per cent throughout the domain, so"
          "\n  little remains anywhere and the comparison is near-random.")
    print("\nThe carbonate mineral rate is excluded: it is the precipitation itself.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
