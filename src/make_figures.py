#!/usr/bin/env python3
"""
Generate every figure for the manuscript and the response letter in one run.

Three scripts currently produce figures, each with its own arguments, its own
input paths and its own idea of where output belongs. This calls all three and
writes everything into one directory, which is what both LaTeX documents expect
under `figures/`.

    manuscript      1  fig_dfn_geometry            generate_figures.py
                    2  fig_timeseries
                    3  fig_carbonate_budget
                    4  fig_connectivity_trapping
                    5  fig_dissolution
                    6  fig_p32_spatial_comparison
                    8  fig_trapping_efficiency
                   10  fig_study_design            fig_study_design.py
                       fig_cations                 fig_cations.py

    letter             letter_index                letter_figures.py
                       letter_variability
                       letter_sensitivity

Figures 7 and 9 are omitted deliberately. fig_topology_trapping repeats
fig_connectivity_trapping now that neither shows a trend, and
fig_stagnation_zones averages a quantity that varies from 0 to 11 times
enrichment across realisations, so a single mean describes no case. Pass
--include-dropped to produce them anyway.

The co-location schematic is not here. It was drawn by hand and has no script.

Usage
-----
    python3 src/make_figures.py                     # everything, into figures/
    python3 src/make_figures.py --out figures_new   # somewhere else
    python3 src/make_figures.py --only manuscript   # or letter, or design
    python3 src/make_figures.py --list              # what would run
"""
from __future__ import annotations
import argparse, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # the revision directory

# name, group, argv relative to ROOT, expected output files
JOBS = [
    ("generate_figures (manuscript 1-6, 8)", "manuscript",
     ["python3", "src/generate_figures.py",
      "--results-dir", "runs", "--dfn-dir", "../dfn_library",
      "--only", "1", "2", "3", "4", "5", "6", "8"],
     ["fig_dfn_geometry.pdf", "fig_timeseries.pdf",
      "fig_carbonate_budget.pdf", "fig_connectivity_trapping.pdf",
      "fig_dissolution.pdf", "fig_p32_spatial_comparison.pdf",
      "fig_trapping_efficiency.pdf"]),

    ("generate_figures (dropped 7, 9)", "dropped",
     ["python3", "src/generate_figures.py",
      "--results-dir", "runs", "--dfn-dir", "../dfn_library",
      "--only", "7", "9"],
     ["fig_topology_trapping.pdf", "fig_stagnation_zones.pdf"]),

    ("fig_study_design", "design",
     ["python3", "src/fig_study_design.py"],
     ["fig_study_design.pdf"]),

    # fig_cations solves a steady flow field per realisation for the age panel,
    # so it is the slowest of these by some margin.
    ("fig_cations", "manuscript",
     ["python3", "src/fig_cations.py", "--runs", "runs",
      "--library", "../dfn_library", "--high", "C_baseline__p32_100_s383"],
     ["fig_cations.pdf"]),

    ("letter_figures", "letter",
     ["python3", "src/letter_figures.py", "--runs", "runs"],
     ["letter_index.pdf", "letter_variability.pdf",
      "letter_sensitivity.pdf"]),
]

# how each script is told where to write
OUTFLAG = {
    "src/generate_figures.py": "--output-dir",
    "src/fig_study_design.py": "--out",
    "src/letter_figures.py":   "--out",
    "src/fig_cations.py":      None,      # takes a full path via --out
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures",
                    help="directory for every figure (default: figures)")
    ap.add_argument("--only", choices=["manuscript", "letter", "design"],
                    action="append", default=None,
                    help="restrict to one group; repeatable")
    ap.add_argument("--include-dropped", action="store_true",
                    help="also produce figures 7 and 9")
    ap.add_argument("--list", action="store_true",
                    help="print what would run and stop")
    a = ap.parse_args()

    # Every one of these scripts needs numpy, matplotlib and h5py. Running
    # under the wrong conda environment fails each of them in turn with the
    # same import error, which is easier to read once here than four times.
    if not a.list:
        missing_mods = []
        for m in ("numpy", "matplotlib", "h5py", "scipy"):
            try:
                __import__(m)
            except ImportError:
                missing_mods.append(m)
        if missing_mods:
            env = os.environ.get("CONDA_DEFAULT_ENV", "?")
            sys.exit(f"missing modules in the active environment ({env}): "
                     f"{', '.join(missing_mods)}\n"
                     f"  try:  conda activate geochem")

    groups = set(a.only) if a.only else {"manuscript", "letter", "design"}
    if a.include_dropped:
        groups.add("dropped")
    todo = [j for j in JOBS if j[1] in groups]
    if not todo:
        sys.exit("nothing selected")
    absent = sorted({j[2][1] for j in todo
                     if not os.path.exists(os.path.join(ROOT, j[2][1]))})
    if absent and not a.list:
        sys.exit("script(s) not found: " + ", ".join(absent))

    out_abs = os.path.abspath(os.path.join(ROOT, a.out))
    if a.list:
        print(f"would write to {out_abs}\n")
        for name, grp, argv, outs in todo:
            print(f"  [{grp}] {name}")
            for o in outs:
                print(f"      {o}")
        return 0

    os.makedirs(out_abs, exist_ok=True)
    print("=" * 66)
    print(f"  writing to {out_abs}")
    print("=" * 66)

    # A figure left over from an earlier run would otherwise be reported as a
    # success when the script that should have replaced it has failed. Every
    # output must therefore be newer than the moment this run started. The
    # margin allows for a filesystem whose timestamps lag slightly.
    t_start = time.time() - 2.0

    def is_fresh(path):
        try:
            return os.path.getmtime(path) >= t_start
        except OSError:
            return False

    failed, produced, stale = [], [], []
    for name, grp, argv, outs in todo:
        script = argv[1]
        cmd = list(argv)
        flag = OUTFLAG.get(script, "--out")
        if script == "src/fig_cations.py":
            cmd += ["--out", os.path.join(out_abs, "fig_cations.pdf")]
        elif flag:
            cmd += [flag, out_abs]
        print(f"\n--- {name}")
        t0 = time.time()
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        dt = time.time() - t0
        # the child scripts already report per-figure progress; surface only
        # the lines that matter here
        for ln in (r.stdout or "").splitlines():
            if any(k in ln for k in ("Saved", "wrote", "skipped", "WARNING",
                                     "Loaded", "ERROR", "not found")):
                print("   ", ln.strip())
        if r.returncode != 0:
            failed.append(name)
            tail = (r.stderr or "").strip().splitlines()[-4:]
            for ln in tail:
                print("    !", ln)
        for o in outs:
            f = os.path.join(out_abs, o)
            if is_fresh(f):
                produced.append(o)
            elif os.path.exists(f):
                stale.append(o)
        print(f"    ({dt:.0f} s)")

    print("\n" + "=" * 66)
    have = sorted(set(produced))
    old_files = sorted(set(stale) - set(have))
    want = sorted({o for _, _, _, outs in todo for o in outs})
    missing = [o for o in want if o not in have and o not in old_files]
    print(f"  {len(have)} of {len(want)} figures written by this run")
    for o in have:
        sz = os.path.getsize(os.path.join(out_abs, o)) / 1024
        print(f"    ok       {o:<34} {sz:>7.0f} kB")
    for o in old_files:
        age = (time.time() - os.path.getmtime(os.path.join(out_abs, o))) / 3600
        print(f"    STALE    {o:<34} {age:>6.1f} h old, not replaced")
    for o in missing:
        print(f"    MISSING  {o}")
    if failed:
        print(f"\n  scripts that exited non-zero: {', '.join(failed)}")
    print("=" * 66)

    # The schematic has no generating script, so it cannot appear here. Say so
    # rather than let its absence look like a failure.
    if not os.path.exists(os.path.join(out_abs, "fig_colocation_schematic.pdf")):
        print("\n  fig_colocation_schematic.pdf is not generated by any script")
        print("  and must be copied in. It also needs redrawing: it shows the")
        print("  flow path intersecting the dissolution zones, whereas the")
        print("  precipitating cells lie at percentile 1.4 for local")
        print("  dissolution in all 19 realisations.")

    return 1 if (missing or old_files or failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
