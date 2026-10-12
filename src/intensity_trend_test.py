#!/usr/bin/env python3
"""
Tests of an efficiency trend with fracture intensity that respect the grouping of
networks by P32 level (Block A: five levels, ten networks each; Block B: two levels).

    python3 src/intensity_trend_test.py            # reads results/A_coupled_networks.csv, results/B_coupled_networks.csv

  Kruskal-Wallis     any difference between levels
  Jonckheere-Terpstra  monotonic trend across ordered levels; p from 20000 permutations
                     of the level labels over the networks
  level medians      Spearman rho of the five level medians against P32; exact p over
                     all 120 orderings of the medians
  Block B            Mann-Whitney U between its two levels
Also the carbonate (co2_kg), for comparison. Writes results/intensity_trend_tests.csv.
"""
import csv, itertools, os, sys
import numpy as np
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE


def load(name):
    p = os.path.join(REV, "results", f"{name}_networks.csv")
    if not os.path.isfile(p):
        sys.exit(f"missing {p}")
    return list(csv.DictReader(open(p)))


def jt_stat(levels, values):
    """Jonckheere-Terpstra statistic: pairs (i<j levels) with x_i < x_j, ties counted 1/2."""
    u = sorted(set(levels)); g = [values[levels == L] for L in u]; s = 0.0
    for a in range(len(g)):
        for b in range(a + 1, len(g)):
            d = g[b][None, :] - g[a][:, None]
            s += (d > 0).sum() + 0.5 * (d == 0).sum()
    return s


def tests(rows, key, rng, nperm=20000):
    lev = np.array([float(r["p32_factor"]) for r in rows]); val = np.array([float(r[key]) for r in rows])
    u = sorted(set(lev)); groups = [val[lev == L] for L in u]
    kw = stats.kruskal(*groups)
    jt = jt_stat(lev, val); perm = np.array([jt_stat(rng.permutation(lev), val) for _ in range(nperm)])
    mean_jt = perm.mean(); p_jt = (np.abs(perm - mean_jt) >= abs(jt - mean_jt)).mean()
    med = np.array([np.median(g) for g in groups]); rho_med = stats.spearmanr(u, med).correlation
    rhos = [stats.spearmanr(u, perm_m).correlation for perm_m in itertools.permutations(med)]
    p_med = np.mean([abs(r) >= abs(rho_med) - 1e-12 for r in rhos])
    return dict(quantity=key, n=len(val), levels=len(u), kruskal_H=kw.statistic, kruskal_p=kw.pvalue,
                jt=jt, jt_p_two_sided=p_jt, median_rho=rho_med, median_p_exact=p_med,
                level_medians=" ".join(f"{m:.4g}" for m in med))


def main():
    rng = np.random.default_rng(1); out = []
    A = load("A_coupled")
    for key in ("efficiency_pct", "co2_kg"):
        r = tests(A, key, rng); out.append(dict(block="A", **r))
        print(f"Block A, {key}: level medians {r['level_medians']}")
        print(f"   Kruskal-Wallis H = {r['kruskal_H']:.2f}, p = {r['kruskal_p']:.3g}")
        pj = r['jt_p_two_sided']; pj_txt = f"< {1 / 20000:.0e}" if pj == 0 else f"= {pj:.3g}"
        print(f"   Jonckheere-Terpstra (monotonic trend), permutation p {pj_txt}")
        print(f"   level medians vs P32: rho = {r['median_rho']:+.2f}, exact p = {r['median_p_exact']:.3g} (120 orderings)")
    try:
        B = load("B_coupled")
        lev = np.array([float(x["p32_factor"]) for x in B]); u = sorted(set(lev))
        for key in ("efficiency_pct", "co2_kg"):
            v = np.array([float(x[key]) for x in B]); mw = stats.mannwhitneyu(v[lev == u[0]], v[lev == u[1]], alternative="two-sided")
            print(f"Block B, {key}: Mann-Whitney U = {mw.statistic:.1f}, p = {mw.pvalue:.3g} (levels {[float(x) for x in u]})")
            out.append(dict(block="B", quantity=key, n=len(v), levels=2, mannwhitney_U=mw.statistic, mannwhitney_p=mw.pvalue))
    except SystemExit as e:
        print(f"Block B skipped: {e}")
    p = os.path.join(REV, "results", "intensity_trend_tests.csv"); keys = list(dict.fromkeys(k for r in out for k in r))
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(out)
    print(f"written: results/intensity_trend_tests.csv")


if __name__ == "__main__":
    main()
