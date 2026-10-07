#!/usr/bin/env python3
"""Bootstrap confidence intervals on coefficients of variation.

R3-3: CVs from small samples "are themselves highly sensitive to a single
outlier". Two CVs cannot be compared on point estimates alone. Answers:
  1. Is there a resolvable CV trend across P32 levels?
  2. Does the variance reduction with domain size plateau from L=30 to 40 m?
Also reports log10 IQR, a more stable spread measure for a right-skewed
variable bounded below by zero.
"""
from __future__ import annotations
import sys, csv, itertools
import numpy as np

NBOOT = 20000
RNG = np.random.default_rng(20260818)

def cv(x):
    m = np.mean(x)
    return 100.0 * np.std(x, ddof=1) / m if m > 0 else np.nan

def boot_ci(x, stat, nboot=NBOOT, lo=2.5, hi=97.5):
    x = np.asarray(x, float); n = len(x); vals = np.empty(nboot)
    for i in range(nboot):
        vals[i] = stat(x[RNG.integers(0, n, n)])
    vals = vals[np.isfinite(vals)]
    return stat(x), np.percentile(vals, lo), np.percentile(vals, hi), vals

def log_iqr(x):
    x = np.asarray(x, float); x = x[x > 0]
    if len(x) < 3: return np.nan
    l = np.log10(x)
    return float(np.percentile(l, 75) - np.percentile(l, 25))

def overlap(a, b):
    return not (a[1] < b[0] or b[1] < a[0])

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "volweighted.csv"
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r["carb_intensity"] = float(r["carb_intensity"])
    def grp(pref):
        return [r["carb_intensity"] for r in rows if r["run_id"].startswith(pref)]

    print("="*74); print("Q1  CV of trapping intensity across P32 levels"); print("="*74)
    levels = [("0.75","A_p32_075"),("1.00","A_p32_100"),("1.25","A_p32_125"),
              ("1.50","A_p32_150"),("2.00","A_p32_200")]
    res = {}
    print(f"  {'P32':>6}{'n':>4}{'CV %':>8}{'95% CI':>22}{'log10 IQR':>12}")
    for lab, pref in levels:
        x = grp(pref)
        if len(x) < 3: continue
        pt, lo, hi, _ = boot_ci(x, cv); res[lab] = (lo, hi, pt, len(x))
        print(f"  {lab:>6}{len(x):>4}{pt:>8.0f}   [{lo:>6.0f}, {hi:>6.0f}]{log_iqr(x):>12.2f}")
    print("\n  pairwise: do the intervals overlap?")
    n_over = n_pair = 0
    for (a,_),(b,_) in itertools.combinations(levels, 2):
        if a not in res or b not in res: continue
        n_pair += 1; ov = overlap(res[a][:2], res[b][:2]); n_over += ov
        if not ov: print(f"    x{a} vs x{b}: SEPARATED")
    print(f"    {n_over} of {n_pair} pairs overlap")
    if n_over == n_pair:
        print("\n  -> No CV difference between ANY pair is resolvable. Supported claim:")
        print("     'variability is large at every fracture intensity tested, with no")
        print("     resolvable trend'. The published ranking is not supported.")

    print(); print("="*74)
    print("Q2  Does variance fall with domain size?  (R3-7, sub-REV)"); print("="*74)
    dom = [("20 m","A_p32_100"),("30 m","D30"),("40 m","D40")]
    dres = {}
    print(f"  {'L':>6}{'n':>4}{'CV %':>8}{'95% CI':>22}{'median intensity':>20}")
    for lab, pref in dom:
        x = grp(pref)
        if len(x) < 3: continue
        pt, lo, hi, _ = boot_ci(x, cv); dres[lab] = (lo, hi, pt, len(x))
        print(f"  {lab:>6}{len(x):>4}{pt:>8.0f}   [{lo:>6.0f}, {hi:>6.0f}]{np.median(x):>20.3e}")
    if "20 m" in dres and "30 m" in dres:
        print(f"\n  20 m vs 30 m: {'overlap' if overlap(dres['20 m'][:2], dres['30 m'][:2]) else 'SEPARATED'}")
    if "30 m" in dres and "40 m" in dres:
        sep = not overlap(dres["30 m"][:2], dres["40 m"][:2])
        print(f"  30 m vs 40 m: {'SEPARATED' if sep else 'overlap'}")
        if not sep:
            print("\n  -> The 30->40 m plateau is NOT resolvable at n=8. Supported:")
            print("     'variance is lower at L>=30 m than at 20 m; no further")
            print("     change is resolvable'.")
    print("\n  bootstrap of the CV difference (30 m - 40 m), the direct test:")
    a, b = np.asarray(grp("D30"), float), np.asarray(grp("D40"), float)
    if len(a) >= 3 and len(b) >= 3:
        d = np.empty(NBOOT)
        for i in range(NBOOT):
            d[i] = cv(a[RNG.integers(0,len(a),len(a))]) - cv(b[RNG.integers(0,len(b),len(b))])
        d = d[np.isfinite(d)]; lo, hi = np.percentile(d, [2.5, 97.5])
        print(f"    difference = {cv(a)-cv(b):+.0f} pp,  95% CI [{lo:+.0f}, {hi:+.0f}] pp")
        print(f"    interval {'excludes' if (lo>0 or hi<0) else 'INCLUDES'} zero"
              f" -> {'resolvable' if (lo>0 or hi<0) else 'not resolvable'}")

    print(); print("="*74)
    print("How much can a CV be trusted at this sample size?"); print("="*74)
    print("  Sampling from a lognormal with true CV = 150%:")
    for n in (5, 8, 10, 15, 25):
        est = [cv(RNG.lognormal(0, np.sqrt(np.log(1+1.5**2)), n)) for _ in range(4000)]
        lo, hi = np.percentile(est, [2.5, 97.5])
        print(f"    n={n:>3}: 95% of estimates fall in [{lo:>5.0f}, {hi:>5.0f}] %")

if __name__ == "__main__":
    main()
