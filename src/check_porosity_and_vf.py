#!/usr/bin/env python3
"""
Check the porosity setting and the mineral volume fraction sum in the
PFLOTRAN decks behind the reported numbers.

  porosity treatment      fixed at 0.50, or recomputed from the mineral
                          inventory at each step (UPDATE_POROSITY)
  volume fraction sum     the primary minerals must sum to 1 - phi = 0.50,
                          not the 0.85 used before the rescaling

runs/ holds the fixed-porosity cases and runs_aperture/ the coupled ones, on
the same networks. Every number in the results sections except the aperture
comparison must come from runs/. Neither setting is recorded in the output
files, so a case built in the wrong directory would not show up downstream.

Four checks:

  1. finalise.py reads runs/ and not runs_aperture/
  2. finalised.csv contains no case that exists only in runs_aperture/
  3. every deck in runs/ sums to 0.50
  4. no deck in runs/ sets UPDATE_POROSITY, and every deck in runs_aperture/
     does

The aperture comparison reads both directories, which is intended.

    python3 src/check_porosity_and_vf.py
"""
import csv
import glob
import os
import re
import sys

TOL = 0.005
TARGET_VF = 0.50

# The first attempts at coupling sit in runs/ rather than runs_aperture/. They
# are incomplete and finalise.py skips them, so they carry UPDATE_POROSITY
# without affecting any reported number. Check 2 confirms they are absent from
# finalised.csv.
EARLY_COUPLED_TESTS = {
    "C_feedback__p32_150_s117", "C_feedback__p32_150_s42",
    "C_feedback__p32_200_s117", "C_feedback__p32_200_s383",
    "test_feedback",
}


def vf_sum(deck):
    """Sum the primary mineral volume fractions.

    In the CONSTRAINT block each mineral is one line,
        <name>  <volume_fraction>  <specific_area>  cm^2/cm^3
    so the second field is the volume fraction. Secondary phases are seeded at
    1e-6 and are excluded, since the sum is over the primary assemblage.
    """
    total = 0.0
    for line in open(deck):
        if "cm^2/cm^3" not in line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            v = float(parts[1].lower().replace("d", "e"))
        except ValueError:
            continue
        if v > 1e-3:
            total += v
    return total


def main():
    fails = []

    print("1. directory read by finalise.py")
    src = open("src/finalise.py").read()
    print(f"   mentions runs_aperture: {'runs_aperture' in src}")
    if "runs_aperture" in src:
        fails.append("finalise.py references runs_aperture")

    print("\n2. cases listed in finalised.csv")
    ids = [r["run_id"] for r in csv.DictReader(open("finalised.csv"))]
    in_runs = {os.path.basename(p) for p in glob.glob("runs/*") if os.path.isdir(p)}
    in_ap = {os.path.basename(p) for p in glob.glob("runs_aperture/*")
             if os.path.isdir(p)}
    strays = [i for i in ids if i not in in_runs and i in in_ap]
    print(f"   {len(ids)} rows, {len(strays)} present only in runs_aperture")
    if strays:
        fails.append(f"aperture-only cases in finalised.csv: {strays[:5]}")

    print("\n3. mineral volume fractions in runs/")
    bad, checked = [], 0
    for d in sorted(in_runs):
        deck = f"runs/{d}/pflotran_co2.in"
        if not os.path.isfile(deck):
            continue
        checked += 1
        s = vf_sum(deck)
        if abs(s - TARGET_VF) > TOL:
            bad.append((d, s))
    print(f"   {checked} decks, {len(bad)} away from {TARGET_VF}")
    for d, s in bad[:6]:
        print(f"     {d}: sum = {s:.4f}")
    if bad:
        fails.append(f"{len(bad)} deck(s) in runs/ not at {TARGET_VF}")

    print("\n4. porosity coupling flag")
    listed = [i for i in ids if i in EARLY_COUPLED_TESTS]
    print(f"   {len(EARLY_COUPLED_TESTS)} early coupled test(s) in runs/ exempted; "
          f"{len(listed)} of them appear in finalised.csv")
    if listed:
        fails.append(f"early coupled test(s) in finalised.csv: {listed}")
    for label, dirs in (("runs", in_runs), ("runs_aperture", in_ap)):
        want = label == "runs_aperture"
        wrong = [d for d in sorted(dirs)
                 if d not in EARLY_COUPLED_TESTS
                 and os.path.isfile(f"{label}/{d}/pflotran_co2.in")
                 and ("UPDATE_POROSITY" in open(f"{label}/{d}/pflotran_co2.in").read())
                 != want]
        state = "with" if want else "without"
        print(f"   {label}: {len(dirs)} cases, expected {state}, {len(wrong)} wrong")
        for d in wrong[:6]:
            print(f"     {d}")
        if wrong:
            fails.append(f"{len(wrong)} case(s) in {label}/ set the wrong "
                         f"porosity option")

    print()
    if fails:
        print(f"{len(fails)} problem(s):")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("reported results come from the fixed-porosity, corrected-VF runs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
