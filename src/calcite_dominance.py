#!/usr/bin/env python3
"""
How often is calcite the largest carbonate phase?

The medians in Table 6 are 51 to 85 % calcite, so the typical realization is
calcite-dominated. That is not the same as calcite always exceeding 50 %, and
the distinction matters because the production runs use a dawsonite rate
constant of 1e-7, which exceeds the measured value by more than three orders of
magnitude.

    python3 src/calcite_dominance.py
"""
import csv
import sys

LEVELS = ["075", "100", "125", "150", "200"]


def main(path="finalised.csv"):
    rows = [r for r in csv.DictReader(open(path))
            if r["run_id"].startswith("A_p32_")
            and float(r["carb_per_cell"]) > 0]

    print("Block A, realizations with carbonate\n")
    print("%-8s %4s %9s %9s %9s" % ("P32", "n", ">50% cal", "daw > cal", "cal < 10%"))
    print("-" * 46)
    tot = {"n": 0, "cal": 0, "daw": 0, "low": 0}
    for lev in LEVELS:
        sub = [r for r in rows if f"_p32_{lev}_" in r["run_id"]]
        cal = sum(1 for r in sub if float(r["pct_calcite"]) > 50)
        daw = sum(1 for r in sub
                  if float(r["pct_dawsonite"]) > float(r["pct_calcite"]))
        low = sum(1 for r in sub if float(r["pct_calcite"]) < 10)
        print("x%-7s %4d %9d %9d %9d" % (lev, len(sub), cal, daw, low))
        tot["n"] += len(sub); tot["cal"] += cal
        tot["daw"] += daw; tot["low"] += low
    print("-" * 46)
    print("%-8s %4d %9d %9d %9d"
          % ("all", tot["n"], tot["cal"], tot["daw"], tot["low"]))

    print("\nthe dawsonite-dominated realizations:")
    for r in sorted(rows, key=lambda x: -float(x["pct_dawsonite"])):
        if float(r["pct_dawsonite"]) > float(r["pct_calcite"]):
            print("  %-24s calcite %5.1f%%  dawsonite %5.1f%%  carb/cell %.2e"
                  % (r["run_id"], float(r["pct_calcite"]),
                     float(r["pct_dawsonite"]), float(r["carb_per_cell"])))

    print(f"\ncalcite exceeds 50 % in {tot['cal']} of {tot['n']} realizations; "
          f"dawsonite is the larger phase in {tot['daw']}.")
    print("These runs use k_dawsonite = 1e-7. The sensitivity variants give 5 to "
          "28 % dawsonite\nat 1e-9 and 1 % or less at 1e-11, so the "
          "dawsonite-dominated cases above reflect\nthe rate constant rather "
          "than the network.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "finalised.csv"))
