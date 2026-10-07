#!/usr/bin/env python3
"""
Is the domain-size trend in the aperture effect resolvable at n = 8?

The median volume ratio falls from 0.70 at L = 20 m to 0.56 at 30 m and 0.47 at
40 m. With eight realisations per domain size the question is whether that is a
trend or sampling noise, so the difference between each pair is bootstrapped.

    python3 src/domain_aperture_trend.py
"""
import os
import sys

import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
SEED = 1.0e-6
RNG = np.random.default_rng(20260914)


def carbonate(run_dir):
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[sorted(f.keys())[-1]]
        tot = None
        for p in CARB:
            k = f"{p} VF [m^3 mnrl_m^3 bulk]"
            if k in g:
                d = np.maximum(g[k][:].flatten() - SEED, 0.0)
                tot = d if tot is None else tot + d
    return tot


def ratios(prefix):
    out = []
    for case in sorted(os.listdir("runs_aperture")):
        if not case.startswith(prefix):
            continue
        fx, ap = f"runs/{case}", f"runs_aperture/{case}"
        if not (os.path.isfile(fx + "/pflotran_co2.h5")
                and os.path.isfile(ap + "/pflotran_co2.h5")):
            continue
        try:
            a, b = carbonate(fx), carbonate(ap)
        except (OSError, KeyError):
            continue
        if a is None or b is None or a.sum() <= 0:
            continue
        out.append(b.sum() / a.sum())
    return np.array(out)


def boot_median_diff(x, y, n=20000):
    """Bootstrap interval on median(x) - median(y)."""
    d = [np.median(RNG.choice(x, len(x), replace=True))
         - np.median(RNG.choice(y, len(y), replace=True)) for _ in range(n)]
    return np.percentile(d, [2.5, 97.5])


def main():
    # A is the 20 m continuous ensemble; D30 and D40 are the same schedule at
    # larger domains, so the three are comparable.
    groups = [("20 m (block A)", "A_p32_"), ("30 m (D30)", "D30_"),
              ("40 m (D40)", "D40_")]
    data = {}
    print("%-16s %4s %8s %14s %8s" % ("domain", "n", "median", "IQR", "falls"))
    print("-" * 58)
    for label, prefix in groups:
        r = ratios(prefix)
        data[label] = r
        if len(r) == 0:
            print(f"  {label}: no pairs")
            continue
        print("%-16s %4d %8.2f  %5.2f--%-6.2f %6.0f%%"
              % (label, len(r), np.median(r),
                 np.percentile(r, 25), np.percentile(r, 75),
                 (r <= 1).sum() / len(r) * 100))

    print("\nbootstrap of the difference in median ratio:")
    labels = [l for l, _ in groups if len(data.get(l, [])) > 1]
    for i in range(len(labels) - 1):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            lo, hi = boot_median_diff(data[a], data[b])
            sign = "resolvable" if lo * hi > 0 else "not resolvable"
            print(f"  {a} minus {b}: "
                  f"{np.median(data[a]) - np.median(data[b]):+.2f} "
                  f"95% CI [{lo:+.2f}, {hi:+.2f}]  {sign}")

    print("\nNote: D30 and D40 have eight realisations each, so the intervals "
          "are wide.\nA monotonic trend across three sizes is weaker evidence "
          "than a resolvable\npairwise difference, and the manuscript should "
          "say which it has.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
