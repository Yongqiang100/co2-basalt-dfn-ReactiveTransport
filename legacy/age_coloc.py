#!/usr/bin/env python3
"""
Do the precipitating cells sit on old water? (AE-1, AE-2)

WHY THIS TEST
-------------
The a priori index shape99 is the ratio of the 99th percentile to the median of
the mean groundwater age field. It correlates with total carbonate at
rho = +0.398 across realisations, which is a BETWEEN-network result.

The co-location analysis is a WITHIN-network result: precipitating cells lie at
a median percentile of 99.4 in dissolved Ca and Mg, 91.6 in pH, and 1.4 in
local dissolution.

The response letter currently bridges the two by asserting that slow-moving
water has longer to neutralise, and that precipitation therefore favours slow
paths. That bridge is untested, and one measured result points against it:
precipitating cells have ELEVATED flux, median percentile 62.0, significant in
11 of 19 realisations. High flux implies short residence, which is the opposite
of the asserted mechanism at the cell level.

The two can coexist -- a network may carry a heavy slow tail while precipitation
within it occurs on faster-than-median cells -- but the letter asserts the cell
level version. This script tests it directly: it computes the percentile rank of
the precipitating cells within their own network's age field, exactly as was
done for flux, pH and the cations.

INTERPRETATION
--------------
  age percentile high (>70)  the asserted mechanism holds; precipitation
                             favours old water and the paragraph stands.
  age percentile near 50     age does not discriminate; shape99 works for some
                             other reason and the paragraph must be rewritten.
  age percentile low  (<30)  precipitation favours YOUNG water; the paragraph
                             is wrong in sign and must be replaced.

Usage
-----
    python3 src/age_coloc.py --runs runs --library ../dfn_library
    python3 src/age_coloc.py --runs runs --library ../dfn_library --csv age_coloc.csv
"""
from __future__ import annotations
import argparse, csv, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEC = CARB + ("Kaolinite", "Chalcedony")
SEED_VF = 1e-6


def pct_rank(x, sel):
    """Mean percentile rank of the selected cells within the full field."""
    ok = np.isfinite(x)
    r = np.full(x.size, np.nan)
    r[ok] = np.argsort(np.argsort(x[ok])) / max(ok.sum() - 1, 1) * 100.0
    v = r[sel & ok]
    return float(np.mean(v)) if v.size else float("nan")


def rank_test(x, sel):
    """Mann-Whitney U, one-sided (selected greater), with rank-biserial r."""
    from scipy import stats
    ok = np.isfinite(x)
    a, b = x[sel & ok], x[~sel & ok]
    if a.size < 3 or b.size < 3:
        return None
    u, p = stats.mannwhitneyu(a, b, alternative="greater")
    auc = u / (a.size * b.size)
    return dict(p=float(p), r=float(2 * auc - 1))


def carbonate(run_dir, final_year=50.0):
    """Net carbonate per cell at the final time, and the unphysical mask."""
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    if not h5:
        return None, None, "no output"
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, None, f"incomplete ({tg[-1][0]:g} yr)"
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
        return None, None, "no carbonate fields"
    return carb, (sec > 1.0), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--library", default="../dfn_library",
                    help="directory holding the DFN meshes")
    ap.add_argument("--prefix", default="C_baseline__")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--min-precip", type=int, default=5)
    a = ap.parse_args()

    import flowfield

    rows = []
    print(f"{'case':<20}{'precip':>7}{'age pct':>9}{'r_age':>8}"
          f"{'p_age':>10}{'age med (yr)':>14}")
    for d in sorted(glob.glob(os.path.join(a.runs, a.prefix + "*"))):
        if not os.path.isdir(d):
            continue
        case = os.path.basename(d).replace(a.prefix, "")
        carb, unphys, why = carbonate(d)
        if carb is None:
            print(f"{case:<20}  {why}")
            continue
        sel = (carb > 0) & ~unphys
        if sel.sum() < a.min_precip:
            print(f"{case:<20}{int(sel.sum()):>7}   fewer than "
                  f"{a.min_precip} precipitating cells")
            continue

        mesh = os.path.join(a.library, case)
        if not os.path.isdir(mesh):
            print(f"{case:<20}  mesh not found at {mesh}")
            continue
        try:
            sol = flowfield.solve(mesh, verbose=False)
        except Exception as e:
            print(f"{case:<20}  flow solve failed: {str(e)[:44]}")
            continue
        age = np.asarray(sol["age"], float)
        if age.size != carb.size:
            print(f"{case:<20}  age field {age.size} != carbonate "
                  f"{carb.size}; cell ordering differs")
            continue

        keep = ~unphys
        ap_ = pct_rank(age[keep], sel[keep])
        t = rank_test(age[keep], sel[keep])
        med = float(np.nanmedian(age)) / 3.156e7
        print(f"{case:<20}{int(sel.sum()):>7}{ap_:>9.1f}"
              f"{(t['r'] if t else float('nan')):>8.3f}"
              f"{(t['p'] if t else float('nan')):>10.2e}{med:>14.3f}")
        rows.append(dict(case=case, n_precip=int(sel.sum()), age_pct=ap_,
                         r_age=(t["r"] if t else None),
                         p_age=(t["p"] if t else None),
                         age_median_yr=med))

    if not rows:
        sys.exit("no realisation could be tested")

    v = np.array([r["age_pct"] for r in rows])
    alpha = 0.05 / len(rows)
    hi = sum(1 for r in rows if r["p_age"] is not None and r["p_age"] < alpha)
    print("\n" + "=" * 68)
    print(f"  realisations tested            {len(rows)}")
    print(f"  age percentile of precipitating cells: median {np.median(v):.1f}, "
          f"range {v.min():.1f}-{v.max():.1f}")
    print(f"  significantly elevated (Bonferroni, alpha = {alpha:.2e}): "
          f"{hi} of {len(rows)}")
    print()
    print("  For comparison, the same statistic on the other fields:")
    print("    dissolved Ca and Mg   99.4   (17 of 17)")
    print("    pH                    91.6   (15 of 19)")
    print("    cell flux             62.0   (11 of 19)")
    print("    local dissolution      1.4   ( 0 of 19)")
    print()
    m = np.median(v)
    if m > 70:
        print("  Age is ELEVATED in the precipitating cells. The mechanism as")
        print("  stated in the response letter holds, and shape99 measures the")
        print("  abundance of the slow paths on which precipitation occurs.")
    elif m < 30:
        print("  Age is DEPRESSED in the precipitating cells. Precipitation")
        print("  favours young water, and the mechanism as stated in the")
        print("  response letter is wrong in sign. shape99 predicts trapping")
        print("  for some other reason, which must be identified or the")
        print("  index presented as an empirical relationship.")
    else:
        print("  Age does not discriminate the precipitating cells. shape99")
        print("  works as a between-network descriptor without implying that")
        print("  precipitation occurs on the oldest water, and the response")
        print("  letter should state the relationship empirically rather than")
        print("  asserting the cell-level mechanism.")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\n  wrote {a.csv}")


if __name__ == "__main__":
    main()
