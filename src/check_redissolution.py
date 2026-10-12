#!/usr/bin/env python3
"""
Did carbonate form and then re-dissolve in the low-yield realizations?

The submitted manuscript claimed three cases reached a peak in the first decade
and declined to zero by 50 years as acid delivery depressed the pH. That claim
was withdrawn, but the replacement text asserts no realization re-dissolves,
which needs measuring rather than asserting.

Reads every saved time group, sums the carbonate above the 1e-6 seed, and
reports the trajectory: peak value, time of peak, final value, and whether the
final value falls materially below the peak.

    python3 src/check_redissolution.py
    python3 src/check_redissolution.py --all      # every Block A realization
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6
LOW_YIELD = ["A_p32_075_s1063", "A_p32_075_s727", "A_p32_125_s727",
             "A_p32_150_s1289", "A_p32_100_s839", "A_p32_125_s1397"]


def trajectory(run_dir):
    """Total carbonate above seed at each saved time, as (times, totals)."""
    path = os.path.join(run_dir, "pflotran_co2.h5")
    times, totals = [], []
    with h5py.File(path, "r") as f:
        groups = []
        for k in f.keys():
            m = re.search(r"([-\d.eE+]+)\s*y", k)
            if m:
                try:
                    groups.append((float(m.group(1)), k))
                except ValueError:
                    continue
        for t, k in sorted(groups):
            g = f[k]
            tot = 0.0
            for p in CARB:
                kk = next((x for x in g if x.startswith(p + " VF")), None)
                if kk:
                    tot += np.maximum(np.asarray(g[kk][:], float).flatten()
                                      - SEED, 0.0).sum()
            times.append(t)
            totals.append(tot)
    return np.array(times), np.array(totals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="check every Block A realization, not just the six")
    a = ap.parse_args()

    cases = (sorted(os.path.basename(p) for p in glob.glob("runs/A_p32_*"))
             if a.all else LOW_YIELD)

    print("%-22s %4s %11s %7s %11s %s"
          % ("case", "nt", "peak", "at yr", "final", "verdict"))
    print("-" * 78)
    redissolved = []
    for c in cases:
        d = os.path.join("runs", c)
        if not os.path.isfile(os.path.join(d, "pflotran_co2.h5")):
            print("%-22s  no output" % c)
            continue
        try:
            t, v = trajectory(d)
        except (OSError, KeyError) as e:
            print("%-22s  unreadable: %s" % (c, e))
            continue
        if len(v) == 0:
            print("%-22s  no time groups" % c)
            continue
        ipk = int(np.argmax(v))
        peak, final = v[ipk], v[-1]
        if peak <= 0:
            verdict = "never precipitated"
        elif final < 0.5 * peak:
            verdict = "PEAK THEN DECLINE"
            redissolved.append(c)
        elif final < 0.95 * peak:
            verdict = "slight decline"
            redissolved.append(c)
        else:
            verdict = "monotonic"
        print("%-22s %4d %11.3e %7.1f %11.3e %s"
              % (c, len(v), peak, t[ipk], final, verdict))

    print()
    if redissolved:
        print(f"{len(redissolved)} case(s) lose carbonate after a peak:")
        for c in redissolved:
            print(f"  {c}")
        print("\nThe response must report this rather than saying no case "
              "re-dissolves.")
    else:
        print("No case loses carbonate after a peak, so the statement that no "
              "realization\nre-dissolves is supported.")
    print("\nNote: only the saved snapshots are visible. A peak between two "
          "snapshots\nwould not appear here.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
