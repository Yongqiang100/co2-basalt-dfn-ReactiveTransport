#!/usr/bin/env python3
"""
Does the a priori index depend on the assumed inlet width? (AE-1)

THE ISSUE
---------
flowfield.py imposes the age boundary condition on the left 20% of the domain
(INJ_FRACTION = 0.20, 13,471 cells in the reference realisation). The reactive
decks inject into REGION injection_zone, which is x <= -9 and covers 2,722
cells, or 3.2%. The two were consistent when the region was first written and
diverged when the deck's region was narrowed.

The index shape99 = p99/median of the mean groundwater age field was computed
with the 20% inlet, and the out-of-sample predictions were frozen and hashed on
that basis before the test simulations ran. The pre-registration therefore
stands on the 20% definition, and recomputing it now would remove the property
that makes the test worth reporting.

The defensible position is that the inlet width is part of the index
definition rather than a mistake: the index has to be computable before the
chemistry and applied identically to every realisation, and both hold. It does
not have to share the reactive deck's boundary condition to predict the
outcome. What the revision must not do is leave a reader to assume the two
match.

This script supplies the accompanying sensitivity. It recomputes the index with
the inlet matched to the deck and re-tests the correlation, so the response can
state whether the result depends on the choice.

    registered   inlet = left 20% of the domain     (as frozen)
    matched      inlet = x <= -9, the deck's region (0.05 of a 20 m domain)

INTERPRETATION
--------------
  both correlations similar   the choice does not matter; report the registered
                              value and cite this as a robustness check.
  matched much weaker         the registered result depends on a wider inlet
                              than the simulations use, and the index should be
                              presented with that caveat.
  matched much stronger       worth reporting, but the registered value remains
                              the one that was predicted in advance.

Usage
-----
    python3 src/inlet_sensitivity.py --runs runs --library ../dfn_library
    python3 src/inlet_sensitivity.py --runs runs --library ../dfn_library \\
        --prefix B_ --csv inlet_sensitivity.csv
"""
from __future__ import annotations
import argparse, csv, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEC = CARB + ("Kaolinite", "Chalcedony")
SEED_VF = 1e-6

REGISTERED_FRACTION = 0.20      # as frozen, flowfield.INJ_FRACTION
DECK_X_MAX = -9.0               # REGION injection_zone upper x bound


def shape99(age):
    """The index: ratio of the 99th percentile of the age field to its median."""
    a = age[np.isfinite(age) & (age > 0)]
    if a.size < 100:
        return None
    med = float(np.median(a))
    return float(np.percentile(a, 99) / med) if med > 0 else None


def carbonate_total(run_dir, final_year=50.0):
    """Volume-weighted carbonate, excluding cells that hold more than their own."""
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "no output"
    with open(uge) as f:
        n = int(f.readline().split()[1])
        vol = np.array([float(f.readline().split()[4]) for _ in range(n)])
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, f"incomplete ({tg[-1][0]:g} yr)"
        g = f[tg[-1][1]]
        fld = lambda pre: next((np.asarray(g[k][:], float) for k in g
                                if k.startswith(pre)), None)
        carb = sec = None
        for m in SEC:
            a = fld(f"{m} VF")
            if a is None:
                continue
            net = np.clip(a - SEED_VF, 0.0, None)
            sec = net if sec is None else sec + net
            if m in CARB:
                carb = net if carb is None else carb + net
    if carb is None:
        return None, "no carbonate fields"
    keep = sec <= 1.0
    return float((carb[keep] * vol[keep]).sum() / vol.sum()), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--library", default="../dfn_library")
    ap.add_argument("--prefix", default="B_")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    import flowfield
    from scipy import stats

    rows = []
    print(f"{'case':<20}{'carbonate':>12}{'reg. index':>12}{'matched':>10}"
          f"{'reg. inlet':>12}{'matched inlet':>15}")
    for d in sorted(glob.glob(os.path.join(a.runs, a.prefix + "*"))):
        if not os.path.isdir(d):
            continue
        case = os.path.basename(d)[len(a.prefix):]
        y, why = carbonate_total(d)
        if y is None:
            print(f"{case:<20}  {why}")
            continue
        mesh = os.path.join(a.library, case)
        if not os.path.isdir(mesh):
            print(f"{case:<20}  mesh not found")
            continue

        # the deck's region expressed as a fraction of this mesh's extent
        with open(os.path.join(d, "full_mesh.uge")) as f:
            n = int(f.readline().split()[1])
            xs = np.array([float(f.readline().split()[1]) for _ in range(n)])
        xmin, xmax = float(xs.min()), float(xs.max())
        matched_frac = max((DECK_X_MAX - xmin) / (xmax - xmin), 1e-6)
        n_reg = int((xs <= xmin + REGISTERED_FRACTION * (xmax - xmin)).sum())
        n_mat = int((xs <= DECK_X_MAX).sum())

        try:
            s_reg = shape99(np.asarray(
                flowfield.solve(mesh, verbose=False,
                                inj_fraction=REGISTERED_FRACTION)["age"], float))
            s_mat = shape99(np.asarray(
                flowfield.solve(mesh, verbose=False,
                                inj_fraction=matched_frac)["age"], float))
        except Exception as e:
            print(f"{case:<20}  flow solve failed: {str(e)[:40]}")
            continue
        if s_reg is None or s_mat is None:
            print(f"{case:<20}  index undefined")
            continue

        print(f"{case:<20}{y:>12.4e}{s_reg:>12.2f}{s_mat:>10.2f}"
              f"{n_reg:>12,}{n_mat:>15,}")
        rows.append(dict(case=case, carbonate=y, index_registered=s_reg,
                         index_matched=s_mat, cells_registered=n_reg,
                         cells_matched=n_mat, matched_fraction=matched_frac))

    if len(rows) < 5:
        sys.exit(f"only {len(rows)} realisations; too few to correlate")

    y = np.array([r["carbonate"] for r in rows])
    print("\n" + "=" * 74)
    print(f"  Spearman correlation with carbonate, n = {len(rows)}")
    print("=" * 74)
    out = {}
    for key, lab in (("index_registered", "registered inlet (left 20%)"),
                     ("index_matched", "matched inlet (deck region)")):
        x = np.array([r[key] for r in rows])
        rho, p2 = stats.spearmanr(x, y)
        one = p2 / 2 if rho > 0 else 1 - p2 / 2
        out[key] = (rho, one)
        print(f"  {lab:<32} rho = {rho:+.3f}   one-sided p = {one:.4f}")

    rr, _ = out["index_registered"]
    rm, _ = out["index_matched"]
    ix = np.array([r["index_registered"] for r in rows])
    im = np.array([r["index_matched"] for r in rows])
    agree, _ = stats.spearmanr(ix, im)
    print(f"\n  the two index values rank the realisations at rho = {agree:+.3f}")
    print()
    if abs(rr - rm) < 0.10:
        print("  The choice of inlet width does not materially affect the result.")
        print("  Report the registered value, which was predicted in advance, and")
        print("  cite this comparison as a robustness check.")
    elif rm < rr - 0.10:
        print("  The matched inlet gives a weaker correlation. The registered")
        print("  result depends in part on an inlet wider than the simulations")
        print("  use, and the index should be presented with that stated.")
    else:
        print("  The matched inlet gives a stronger correlation. The registered")
        print("  value remains the one predicted in advance and should be the")
        print("  headline; this comparison belongs alongside it.")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\n  wrote {a.csv}")


if __name__ == "__main__":
    main()
