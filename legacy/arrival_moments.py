#!/usr/bin/env python3
"""
Does the spread of arrival times predict where carbonate forms?

WHY
---
The mean groundwater age is elevated in the precipitating cells (median
percentile 89.3, in 17 of 19 realisations) but old water occurs without
precipitation, so age alone does not locate the precipitation zone.

Age cannot in fact do better than it does, for a reason worth stating. The
natural proxy for cation supply is the reactive surface the water has passed
over, and that field obeys the same advection equation as age with the reactive
surface area per volume as its source instead of the porosity:

    div(q a) = phi          div(q s) = sigma

Both sources are uniform in this model, phi = 0.50 and sigma = 100 cm2/cm3
everywhere, so s = (sigma/phi) a and the two fields are identical up to a
constant. Within a homogeneous mineralogy, time in the domain and rock
contacted are the same quantity.

What the mean age omits is DILUTION. Two cells can share a mean age and differ
entirely in concentration: one fed by a single old path carries concentrated
cations, another fed by old and young paths converging is diluted at the
junction. The mean cannot distinguish them; the spread can.

WHAT THIS COMPUTES
------------------
Goode's formulation gives moments of the arrival-time distribution, and the
second moment satisfies the same operator with the mean age as its source:

    div(q a)  = phi                     (first moment,  already in flowfield)
    div(q m2) = 2 phi a                 (second moment, added here)
    var = m2 - a^2

The candidate discriminator is

    C = a / sqrt(var)

which is large only where old water arrives without having been mixed with
young water. It is dimensionless, and it is the reciprocal of the coefficient
of variation of the arrival time.

flowfield.py is NOT modified. It produced the predictions that were frozen and
hashed before the test simulations were run, and changing it would invalidate
that record. The matrix is instead rebuilt here from the fields solve()
returns.

CAUTION -- the numerical floor
------------------------------
A first-order upwind scheme makes every cell a mixing tank, so a single path
with no physical mixing already carries the variance of n tanks in series:

    a = n tau,      var = n tau^2,      C = a/sqrt(var) = sqrt(n)

Verified on a 1-D chain of 40 cells at tau = 0.25 s, where the scheme returns
var = 2.500 against a true value of zero. The discriminator therefore measures
the number of cells along the path before it measures anything physical.

Physical mixing adds a term on top of that floor. For a cell fed by two paths
of mean age a1 and a2 in equal proportion,

    var  ~  n tau^2  +  ((a1 - a2)/2)^2

so the mixing signal is detectable only where the paths differ in age by more
than about sqrt(n) tau. Whether that holds in these networks is the question
this script answers, and the answer may be no.

The script therefore reports the rank correlation between C and the mean age
across cells. If it is near unity, C is a monotone function of age, carries no
information the mean age does not, and the test has failed for numerical rather
than physical reasons.

The scheme can also return m2 < a^2 on a strongly bimodal distribution, giving
a negative variance. That fraction is reported rather than clipped.

This quantity was NOT pre-registered. Any result from it is exploratory and
must be reported as such, alongside the frozen result rather than in place of
it.

Usage
-----
    python3 src/arrival_moments.py --runs runs --library ../dfn_library
    python3 src/arrival_moments.py --runs runs --library ../dfn_library \\
        --prefix A_ --csv arrival_moments.csv
"""
from __future__ import annotations
import argparse, csv, glob, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6
POROSITY = 0.50          # matches flowfield.POROSITY


def second_moment(sol):
    """Second moment of the arrival-time distribution, and its variance.

    Rebuilds the operator flowfield.solve() used for the first moment. For a
    non-stagnant cell k that operator is

        (total outflow)_k a_k  -  sum over upstream u of  f_u a_u  =  phi V_k

    and the second moment obeys the same left-hand side with 2 phi V a on the
    right. Stagnant cells were pinned to zero in the first solve and are pinned
    again here.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu

    up, dn = np.asarray(sol["up"]), np.asarray(sol["dn"])
    f = np.asarray(sol["flux"], float)
    vol = np.asarray(sol["vol"], float)
    stag = np.asarray(sol["stagnant"], bool)
    bnd = np.asarray(sol["bnd_out"], float)
    age = np.asarray(sol["age"], float)
    n = vol.size

    i_out = np.zeros(n)
    np.add.at(i_out, up, f)
    denom = i_out + bnd

    r, c, v = [], [], []
    for k in np.where(~stag)[0]:
        r.append(k); c.append(k); v.append(denom[k])
    for u, d_, ff in zip(up, dn, f):
        if not stag[d_]:
            r.append(d_); c.append(u); v.append(-ff)
    for k in np.where(stag)[0]:
        r.append(k); c.append(k); v.append(1.0)

    A = coo_matrix((v, (r, c)), shape=(n, n)).tocsc()
    rhs = 2.0 * vol * POROSITY * np.nan_to_num(age)
    rhs[stag] = 0.0
    try:
        m2 = splu(A).solve(rhs)
    except Exception as e:
        return None, None, f"second-moment solve failed: {str(e)[:50]}"

    var = m2 - age ** 2
    neg = float(np.mean(var[np.isfinite(var)] < 0)) if np.isfinite(var).any() else 1.0
    return m2, var, None if neg < 0.30 else f"variance negative in {100*neg:.0f}% of cells"


def pct_rank(x, sel):
    ok = np.isfinite(x)
    r = np.full(x.size, np.nan)
    r[ok] = np.argsort(np.argsort(x[ok])) / max(ok.sum() - 1, 1) * 100.0
    v = r[sel & ok]
    return float(np.mean(v)) if v.size else float("nan")


def rank_test(x, sel):
    from scipy import stats
    ok = np.isfinite(x)
    a, b = x[sel & ok], x[~sel & ok]
    if a.size < 3 or b.size < 3:
        return None
    u, p = stats.mannwhitneyu(a, b, alternative="greater")
    return dict(p=float(p), r=float(2 * u / (a.size * b.size) - 1))


def carbonate(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    if not h5:
        return None
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg or tg[-1][0] < final_year * 0.999:
            return None
        g = f[tg[-1][1]]
        tot = None
        for m in CARB:
            k = next((x for x in g if x.startswith(f"{m} VF")), None)
            if k:
                a = np.clip(np.asarray(g[k][:], float) - SEED_VF, 0.0, None)
                tot = a if tot is None else tot + a
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--library", default="../dfn_library")
    ap.add_argument("--prefix", default="A_")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--min-precip", type=int, default=5)
    a = ap.parse_args()

    import flowfield

    rows = []
    print(f"{'case':<20}{'precip':>7}{'age pct':>9}{'C pct':>8}"
          f"{'r_C':>8}{'rho(C,age)':>12}{'neg var':>9}")
    for d in sorted(glob.glob(os.path.join(a.runs, a.prefix + "*"))):
        if not os.path.isdir(d):
            continue
        nm = os.path.basename(d)
        mm = re.search(r"(L\d+_)?p32_\d+_s\d+", nm)
        case = mm.group(0) if mm else nm
        carb = carbonate(d)
        if carb is None:
            continue
        sel = carb > 0
        if sel.sum() < a.min_precip:
            continue
        mesh = os.path.join(a.library, case)
        if not os.path.isdir(mesh):
            print(f"{case:<20}  no mesh at {mesh}")
            continue
        try:
            sol = flowfield.solve(mesh, verbose=False)
        except Exception as e:
            print(f"{case:<20}  flow solve failed: {str(e)[:38]}")
            continue
        age = np.asarray(sol["age"], float)
        if age.size != carb.size:
            print(f"{case:<20}  field {age.size} != carbonate {carb.size}")
            continue
        m2, var, why = second_moment(sol)
        if m2 is None:
            print(f"{case:<20}  {why}")
            continue
        neg = float(np.mean(var[np.isfinite(var)] < 0))
        with np.errstate(invalid="ignore", divide="ignore"):
            C = np.where(var > 0, age / np.sqrt(var), np.nan)

        # If C is a monotone function of age it adds nothing, whatever its
        # percentile. The rank correlation across cells settles that.
        from scipy import stats as _st
        ok = np.isfinite(C) & np.isfinite(age)
        rho_ca = (float(_st.spearmanr(C[ok], age[ok])[0])
                  if ok.sum() > 100 else float("nan"))

        ap_ = pct_rank(age, sel)
        cp_ = pct_rank(C, sel)
        t = rank_test(C, sel)
        print(f"{case:<20}{int(sel.sum()):>7}{ap_:>9.1f}{cp_:>8.1f}"
              f"{(t['r'] if t else float('nan')):>8.3f}{rho_ca:>12.3f}"
              f"{100*neg:>8.0f}%")
        rows.append(dict(case=case, n_precip=int(sel.sum()), age_pct=ap_,
                         C_pct=cp_, r_C=(t["r"] if t else None),
                         p_C=(t["p"] if t else None), rho_C_age=rho_ca,
                         neg_var_frac=neg))

    if not rows:
        sys.exit("no realisation could be tested")

    ageq = np.array([r["age_pct"] for r in rows])
    cq = np.array([r["C_pct"] for r in rows])
    negs = np.array([r["neg_var_frac"] for r in rows])
    alpha = 0.05 / len(rows)
    hi = sum(1 for r in rows if r["p_C"] is not None and r["p_C"] < alpha)
    print("\n" + "=" * 70)
    print(f"  realisations tested        {len(rows)}")
    print(f"  negative variance          median {100*np.median(negs):.0f}% of cells")
    print(f"  mean age percentile        median {np.median(ageq):.1f}")
    print(f"  discriminator percentile   median {np.median(cq):.1f}")
    print(f"  significantly elevated     {hi} of {len(rows)} "
          f"(Bonferroni, alpha = {alpha:.1e})")
    rca = np.array([r["rho_C_age"] for r in rows if r["rho_C_age"] is not None])
    if rca.size:
        print(f"  rank corr C against age    median {np.median(rca):+.3f}")
    print()
    if rca.size and abs(np.median(rca)) > 0.90:
        print("  C is very nearly a monotone function of the mean age, which is")
        print("  what the tanks-in-series floor predicts: on this discretisation")
        print("  C ~ sqrt(number of cells along the path). It carries no")
        print("  information beyond the mean age, and the numerical dispersion of")
        print("  the first-order scheme is the reason. A higher-order transport")
        print("  scheme, or particle tracking, would be needed to separate")
        print("  physical mixing from the discretisation.")
        print()
    if np.median(negs) > 0.10:
        print("  The variance is negative in more than a tenth of cells, which is")
        print("  numerical dispersion in the first-order upwind scheme rather than")
        print("  a property of the flow. The discriminator is not usable on this")
        print("  discretisation and the comparison below should be disregarded.")
    elif np.median(cq) > np.median(ageq) + 3:
        print("  The discriminator locates the precipitating cells better than the")
        print("  mean age does. Worth reporting as an exploratory result alongside")
        print("  the frozen index, not in place of it.")
    elif np.median(cq) < np.median(ageq) - 3:
        print("  The discriminator locates them less well than the mean age. The")
        print("  dilution correction does not help, and the mean age remains the")
        print("  better of the two flow-only quantities.")
    else:
        print("  The discriminator and the mean age locate the precipitating cells")
        print("  equally well, so the arrival-time spread adds nothing beyond the")
        print("  mean. Reported as a negative result.")
    print()
    print("  This quantity was not pre-registered. Whatever it shows is")
    print("  exploratory and belongs alongside the frozen result.")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\n  wrote {a.csv}")


if __name__ == "__main__":
    main()
