#!/usr/bin/env python3
"""
Carbonate assemblage by fracture intensity, reported as medians across runs.

Volume weighting answers "of all the carbonate formed, what fraction is
calcite", which one high-carbonate realisation can dominate: at P32 x0.75 the
run A_p32_075_s1481 holds more carbonate than the next two combined and is
99.85% dawsonite, so the volume-weighted calcite share for that level reads 37%
while eight of the ten runs are calcite-dominated.

The median across runs describes the typical realisation instead. Medians of
separate percentages do not sum to 100, and the caption must say so.

Values come from finalised.csv, which reads the fixed-aperture ensemble in
runs/ with mineral volume fractions summing to 1 - phi = 0.50.

    python3 src/assemblage_median.py
"""
import csv
import sys

import numpy as np

MINERALS = ["calcite", "magnesite", "siderite", "dawsonite"]
LEVELS = ["075", "100", "125", "150", "200"]


def main(path="finalised.csv"):
    rows = [r for r in csv.DictReader(open(path))
            if r["run_id"].startswith("A_p32_")]
    if not rows:
        print("no Block A rows found")
        return 1

    print(f"Block A, fixed aperture, {len(rows)} runs\n")
    print("%-8s %4s %4s %9s %9s %9s %9s"
          % ("P32", "n", "zero", "calcite", "magnesite", "siderite", "dawsonite"))
    print("-" * 62)

    allvals = {m: [] for m in MINERALS}
    for lev in LEVELS:
        sub = [r for r in rows if f"_p32_{lev}_" in r["run_id"]]
        # a run with no carbonate has no assemblage to report
        prec = [r for r in sub if float(r["carb_per_cell"]) > 0]
        zero = len(sub) - len(prec)
        med = {}
        for m in MINERALS:
            v = np.array([float(r["pct_" + m]) for r in prec])
            med[m] = np.median(v) if len(v) else np.nan
            allvals[m].extend(v.tolist())
        print("x%-7s %4d %4d %8.0f%% %8.0f%% %8.0f%% %8.0f%%"
              % (lev, len(prec), zero, med["calcite"], med["magnesite"],
                 med["siderite"], med["dawsonite"]))

    print("-" * 62)
    print("%-8s %4d %4s %8.0f%% %8.0f%% %8.0f%% %8.0f%%"
          % ("all", len(allvals["calcite"]), "",
             *[np.median(allvals[m]) for m in MINERALS]))

    print("\ncalcite-dominated runs per level, for the caption:")
    for lev in LEVELS:
        sub = [r for r in rows if f"_p32_{lev}_" in r["run_id"]]
        prec = [r for r in sub if float(r["carb_per_cell"]) > 0]
        n = sum(1 for r in prec if float(r["pct_calcite"]) > 50)
        print(f"  x{lev}: {n} of {len(prec)}")

    print("\nMedians of separate percentages do not sum to 100.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "finalised.csv"))
