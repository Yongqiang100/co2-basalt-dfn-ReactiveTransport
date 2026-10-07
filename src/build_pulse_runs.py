#!/usr/bin/env python3
"""
Builds shortened-injection runs from the corrected baseline decks.

The variant path in run_rt.sh applies the deck edit before the deck
corrections, so set_shutin reads the injection rate from a deck that has not
yet been scaled and copies the template value of 0.01 kg/s. The correction
then sees a shut-in schedule and skips the scaling by design. Every run built
that way injects at 0.01 kg/s while its baseline injects at 7.35e-03 to
1.95e-02, so the comparison mixes the schedule with a 2.65-fold spread in dose
that correlates with network size.

The E-block runs avoided this because they were built the other way round:
their decks are the already-corrected baseline decks with two rows appended,
the network's own rate at the stopping time and zero just after. This script
does the same, so the flux matches the baseline and only the duration differs.
The delivered volume follows from the duration, 4000 x T / 50 pore volumes.

Times below the deck's 0.01-year ramp need the ramp shortened, or the
tabulated times would not increase monotonically and PFLOTRAN would reject the
schedule. For those the plateau is placed at the stopping time and the
intermediate row dropped, so the rate ramps over the whole pulse.

Mesh and boundary files are symlinked to the baseline directory, as the
E-block runs do, so the pore volume resolves to the same mesh the baseline
used. The thermodynamic databases are copied. No output is copied.

    python3 src/build_pulse_runs.py --dry-run
    python3 src/build_pulse_runs.py
    python3 src/build_pulse_runs.py --durations 2.0:2yr 0.123:45d
"""
from __future__ import annotations
import argparse
import os
import re
import shutil
import sys

CASES = ["p32_075_s117", "p32_100_s383", "p32_125_s383",
         "p32_150_s117", "p32_150_s42", "p32_200_s117"]

DURATIONS = [("2.0", "2yr"), ("0.123", "45d"), ("0.0822", "30d"),
             ("0.0274", "10d")]

LINK = ("full_mesh.uge", "full_mesh.inp", "dfn_properties.h5",
        "boundary_back_s.ex", "boundary_bottom.ex", "boundary_front_n.ex",
        "boundary_left_w.ex", "boundary_right_e.ex", "boundary_top.ex")
COPY = ("hanford.dat", "co2_sw.dat")
RAMP_END = 0.01          # the deck's ramp reaches the plateau here


def plateau_rate(deck_path):
    """The largest positive rate in the RATE LIST block."""
    txt = open(deck_path).read()
    m = re.search(r"RATE LIST(.*?)\n\s*/", txt, re.S)
    if not m:
        return None
    best = 0.0
    for line in m.group(1).splitlines():
        p = line.split()
        if len(p) == 2:
            try:
                v = float(p[1].replace("d", "e"))
            except ValueError:
                continue
            best = max(best, v)
    return best or None


def build_deck(src, dst, years, rate):
    """The baseline deck with the injection stopped at `years`."""
    txt = open(src).read()
    m = re.search(r"(RATE LIST\s*\n\s*TIME_UNITS[^\n]*\n\s*DATA_UNITS[^\n]*\n)"
                  r"(.*?)(\n\s*/)", txt, re.S)
    if not m:
        return False, "RATE LIST block not matched"

    if years > RAMP_END:
        rows = [(0.0, 0.0), (1e-4, rate / 10.0), (RAMP_END, rate),
                (years, rate), (years * 1.0001, 0.0)]
    else:
        # the pulse ends inside the deck's ramp, so ramp over the pulse itself
        rows = [(0.0, 0.0), (1e-4, rate / 10.0),
                (years, rate), (years * 1.0001, 0.0)]

    for a, b in zip(rows, rows[1:]):
        if b[0] <= a[0]:
            return False, f"times not increasing: {a[0]:g} then {b[0]:g}"

    # The opening row of a baseline deck reads 0.000000e+00 and the closing
    # zero reads 0.d0. Both are valid Fortran, but keeping each in its
    # customary place means a diff against the baseline shows only the
    # appended rows.
    body = "\n".join(
        f"    {t:.6e}   {r:.6e}" if r
        else f"    {t:.6e}   {'0.000000e+00' if i == 0 else '0.d0'}"
        for i, (t, r) in enumerate(rows))
    open(dst, "w").write(txt[:m.start(2)] + body + txt[m.end(2):])
    return True, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--prefix", default="F")
    ap.add_argument("--cases", nargs="*", default=CASES)
    ap.add_argument("--durations", nargs="*",
                    default=[f"{y}:{t}" for y, t in DURATIONS],
                    help="years:tag pairs, e.g. 0.0274:10d")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    durations = []
    for d in a.durations:
        y, _, t = d.partition(":")
        durations.append((float(y), t or y))

    made, failed = [], []
    print("%-34s %12s %10s %9s" % ("run", "rate kg/s", "stop y", "PV inj"))
    print("-" * 70)

    for case in a.cases:
        base = os.path.join(a.runs, f"C_baseline__{case}")
        src = os.path.join(base, "pflotran_co2.in")
        if not os.path.isfile(src):
            failed.append((case, "no baseline deck"))
            continue
        rate = plateau_rate(src)
        if rate is None:
            failed.append((case, "no rate in the baseline deck"))
            continue

        for years, tag in durations:
            rid = f"{a.prefix}_{tag}__{case}"
            wd = os.path.join(a.runs, rid)
            pv_inj = 4000.0 * years / 50.0
            print("%-34s %12.6e %10.5f %9.3f" % (rid, rate, years, pv_inj))
            if a.dry_run:
                continue

            os.makedirs(wd, exist_ok=True)
            for f in LINK:
                s = os.path.abspath(os.path.join(base, f))
                d = os.path.join(wd, f)
                if os.path.islink(d) or os.path.exists(d):
                    os.remove(d)
                if os.path.exists(s):
                    os.symlink(s, d)
            for f in COPY:
                s = os.path.join(base, f)
                if os.path.isfile(s):
                    shutil.copy2(s, wd)
            # never carry output across
            for f in ("pflotran_co2.h5", "pflotran_co2-mas.dat",
                      "pflotran_co2.out", "pflotran.log"):
                p = os.path.join(wd, f)
                if os.path.exists(p):
                    os.remove(p)

            ok, info = build_deck(src, os.path.join(wd, "pflotran_co2.in"),
                                  years, rate)
            if ok:
                made.append(rid)
            else:
                failed.append((rid, info))

    print()
    if a.dry_run:
        print(f"would build {len(a.cases) * len(durations)} runs")
        return 0

    print(f"built {len(made)} run(s)")
    for rid, why in failed:
        print(f"  FAILED {rid}: {why}")

    if made:
        print()
        print("Check one deck before submitting:")
        print(f"  grep -A9 'RATE LIST' {a.runs}/{made[0]}/pflotran_co2.in")
        print()
        print("The rate must equal the baseline's and the times must increase.")
        print()
        print(f"Runs are named {a.prefix}_<tag>__<case>. Submit them with an")
        print("array over the list written to pulse_runs.txt.")
        # The list sits beside the runs directory, because run_pulse.sh
        # reads it from the revision root rather than from the working
        # directory. Existing entries are kept, so building in more than one
        # pass does not discard the earlier runs, and an array sized from the
        # file covers everything built.
        listing = os.path.join(os.path.dirname(os.path.abspath(a.runs)),
                               "pulse_runs.txt")
        have = []
        if os.path.isfile(listing):
            have = [l.strip() for l in open(listing) if l.strip()]
        merged = have + [r for r in made if r not in have]
        with open(listing, "w") as f:
            f.write("\n".join(merged) + "\n")
        print(f"  {listing}: {len(merged)} run(s), "
              f"{len(merged) - len(have)} added")
        print(f"  submit: sbatch --array=0-{len(merged)-1} "
              f"slurm/run_pulse.sh")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
