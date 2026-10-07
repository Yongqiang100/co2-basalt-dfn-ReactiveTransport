#!/usr/bin/env python3
"""
Validate a completed PFLOTRAN run from its HDF5 output.

Answers the two questions that matter for every run in the run set:
  1. did it actually reach the final simulation time?
  2. does the chemistry reproduce the published value (for archived cases)?

A run that stalls at t = 3 yr exits 0 and looks fast. A stalled run entering the
ensemble would corrupt the variance statistics that are the paper's principal
result, so this reads the output rather than trusting the exit code.

Usage
-----
    python3 src/validate_run.py bench/p32_100_s383_n24/pflotran_co2.h5
    python3 src/validate_run.py <h5> --expect 2.202e-5
    python3 src/validate_run.py <h5> --case p32_100_s383   # look up the reference
"""
from __future__ import annotations
import sys, os, re, csv, argparse

CARBONATES = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED_VF = 1e-6          # secondary phases are initialised at this value


def time_groups(f):
    """Mirror compute_betweenness.py:105-120. Group names carry LEADING
    whitespace ('  Time  5.00000E+01 y'), so startswith('Time') finds nothing."""
    out = []
    for k in f.keys():
        if "Time" not in k or "failure" in k.lower() or "cut" in k.lower():
            continue
        parts = k.strip().split()
        try:
            ti = parts.index("Time")
            t = float(parts[ti + 1])
            t_yr = t if (len(parts) > ti + 2 and parts[ti + 2] == "y") else t / 3.156e7
            out.append((t_yr, k))
        except Exception:
            continue
    out.sort()
    return out


def find_var(grp, prefix):
    for k in grp.keys():
        if k.startswith(prefix):
            return k
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--expect", type=float, help="reference carbonate per cell")
    ap.add_argument("--case", help="look the reference up in betweenness_results.csv")
    ap.add_argument("--final-year", type=float, default=50.0)
    ap.add_argument("--tol", type=float, default=0.05, help="relative tolerance")
    a = ap.parse_args()

    try:
        import h5py
    except ImportError:
        sys.exit("FATAL: h5py not available in this interpreter")
    if not os.path.isfile(a.h5):
        sys.exit(f"FATAL: {a.h5} not found")

    expect = a.expect
    if expect is None and a.case:
        here = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(here, "..", "..", "betweenness_results.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                if r.get("name") == a.case:
                    try:
                        expect = float(r["carb_per_cell_final"])
                    except (KeyError, ValueError):
                        pass
                    break

    with h5py.File(a.h5, "r") as f:
        tg = time_groups(f)
        if not tg:
            print("keys:", list(f.keys())[:20]); sys.exit("naming mismatch")
        last_yr, last_key = tg[-1]
        print(f"  file         : {a.h5}")
        print(f"  time groups  : {len(tg)}  (first {tg[0][0]:g} y, last {last_yr:g} y)")

        reached = last_yr >= a.final_year * 0.999
        print(f"  reached {a.final_year:g} y : {'YES' if reached else 'NO -- STALLED'}")

        grp = f[last_key]
        total, per_min, ncells = 0.0, {}, None
        for m in CARBONATES:
            k = find_var(grp, f"{m} VF")
            if k is None:
                per_min[m] = None
                continue
            arr = grp[k][:]
            ncells = len(arr)
            net = (arr - SEED_VF).clip(min=0.0)   # net of the nucleation seed
            per_min[m] = float(net.sum())
            total += per_min[m]

        print(f"  cells        : {ncells:,}" if ncells else "  cells        : ?")
        for m, v in per_min.items():
            print(f"    {m:<10}   {'absent' if v is None else f'{v:.6e}'}")
        if ncells:
            pc = total / ncells
            print(f"  TOTAL carbonate VF per cell : {pc:.4e}")
            if expect is not None:
                rel = abs(pc - expect) / expect if expect else float("inf")
                ok = rel <= a.tol
                print(f"  published reference         : {expect:.4e}")
                print(f"  relative difference         : {rel*100:.2f}%  "
                      f"{'MATCH' if ok else 'MISMATCH'}")
                if not ok:
                    print("\n  A mismatch here means the toolchain does not reproduce the")
                    print("  archive. Do not start the run set until it is understood --")
                    print("  candidates: database version, solver tolerances, mesh, or a")
                    print("  deck difference.")
                    sys.exit(2)

    if not reached:
        sys.exit(3)
    print("\n  VALID")


if __name__ == "__main__":
    main()
