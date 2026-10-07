#!/usr/bin/env python3
"""
Do the age and the flux together locate the precipitation zone?

WHY
---
Neither flow-only field is sufficient on its own. Across the 38 realisations of
the primary ensemble that produced measurable carbonate, the precipitating
cells lie at a median percentile of

    mean groundwater age    82.6    significant in 26 of 38
    cell throughflow        61.9    significant in 18 of 38

and old water occurs where no carbonate forms.

The two measure opposing things. A cell with high age and low flux is a
stagnant pocket: the water has had time to neutralise but arrives too slowly to
deliver cations. A cell with high flux and low age sits near the inlet, where
the acid has not yet been spent. Precipitation should require both conditions
at once, and neither field alone expresses that.

Three combinations are therefore tested, all computed from percentile ranks
within each network so that no scaling is assumed:

    product     r_age * r_flux           rewards both being high
    minimum     min(r_age, r_flux)       demands both, ignores by how much
    geometric   sqrt(r_age * r_flux)     as the product, on the original scale

The minimum is the strictest reading of "both conditions at once". The product
allows one to compensate for the other.

WHAT WOULD COUNT
----------------
The comparison is against the better of the two single fields, which is the age
at 82.6. A combination worth reporting must place the precipitating cells
materially higher than that, and must do so in more realisations than the 26 of
38 the age achieves.

If none does, then within a homogeneous mineralogy the flow field alone cannot
locate the precipitation zone. That is a result in itself: it points at the
spatially uniform mineralogy as the limiting assumption, which is the substance
of Reviewer 1's first comment.

None of these was pre-registered. Any positive result is exploratory and
belongs alongside the frozen index rather than in place of it.

Usage
-----
    python3 src/joint_predictor.py --runs runs --library ../dfn_library
    python3 src/joint_predictor.py --runs runs --library ../dfn_library \\
        --prefix A_ --csv joint_predictor.csv
"""
from __future__ import annotations
import argparse, csv, glob, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6

# the age percentile achieved by the age field alone, for comparison
AGE_ALONE_PCT, AGE_ALONE_SIG, AGE_ALONE_N = 82.6, 26, 38



def run_flux(run_dir):
    """Cell throughflow from the simulated pressure field.

    flowfield.solve() imposes a uniform permeability of 8.33e-8, whereas the
    meshes carry aperture-dependent permeability and the aperture variation is
    the dominant control on flow in a fracture network. Its flux field is
    therefore not the one the reaction saw: measured against it the
    precipitating cells appear at percentile 7 rather than the 62 that
    coloc_test.py reports from the run itself.

    This reads the pressure the simulation solved and forms the same quantity
    coloc_test.py forms: half the sum of the connection fluxes at each cell,
    with the connection areas taken from the .uge.
    """
    import h5py
    from validate_run import time_groups
    uge = os.path.join(run_dir, "full_mesh.uge")
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    if not (os.path.exists(uge) and h5):
        return None
    xyz, i_idx, j_idx, area = [], [], [], []
    with open(uge) as f:
        n = int(f.readline().split()[1])
        for _ in range(n):
            t = f.readline().split()
            xyz.append((float(t[1]), float(t[2]), float(t[3])))
        m = int(f.readline().split()[1])
        for _ in range(m):
            t = f.readline().split()
            i_idx.append(int(t[0]) - 1); j_idx.append(int(t[1]) - 1)
            area.append(float(t[5]))
    xyz = np.asarray(xyz); i_idx = np.asarray(i_idx); j_idx = np.asarray(j_idx)
    area = np.asarray(area)
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None
        g = f[tg[-1][1]]
        k = next((x for x in g if x.startswith("Liquid Pressure")), None)
        if k is None:
            return None
        press = np.asarray(g[k][:], float).flatten()[:n]
    d = np.linalg.norm(xyz[i_idx] - xyz[j_idx], axis=1)
    d[d <= 0] = np.finfo(float).eps
    # k/mu is absorbed into the scaling; only ratios are used
    q = area / d * (press[i_idx] - press[j_idx])
    flux = np.zeros(n)
    np.add.at(flux, i_idx, np.abs(q))
    np.add.at(flux, j_idx, np.abs(q))
    return flux * 0.5


def ranks(x):
    """Percentile rank of each finite element, 0 to 100."""
    r = np.full(x.size, np.nan)
    ok = np.isfinite(x)
    if ok.sum() > 1:
        r[ok] = np.argsort(np.argsort(x[ok])) / (ok.sum() - 1) * 100.0
    return r


def pct_rank(x, sel):
    r = ranks(x)
    v = r[sel & np.isfinite(r)]
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

    COMBOS = ("age", "flux", "product", "minimum", "geometric")
    rows = []
    print(f"{'case':<20}{'precip':>7}" + "".join(f"{c:>10}" for c in COMBOS))
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
            print(f"{case:<20}  flow solve failed: {str(e)[:36]}")
            continue
        age = np.asarray(sol["age"], float)
        if age.size != carb.size:
            print(f"{case:<20}  field {age.size} != carbonate {carb.size}")
            continue

        # The flux must come from the run's own pressure field, not from the
        # uniform-permeability solve that produced the age.
        flux = run_flux(d)
        if flux is None or flux.size != age.size:
            print(f"{case:<20}  no usable pressure field")
            continue

        r_age, r_flux = ranks(age), ranks(flux)
        fields = {
            "age": age,
            "flux": flux,
            "product": r_age * r_flux,
            "minimum": np.minimum(r_age, r_flux),
            "geometric": np.sqrt(np.clip(r_age * r_flux, 0, None)),
        }
        rec = dict(case=case, n_precip=int(sel.sum()))
        line = f"{case:<20}{int(sel.sum()):>7}"
        for c in COMBOS:
            pc = pct_rank(fields[c], sel)
            t = rank_test(fields[c], sel)
            rec[f"{c}_pct"] = pc
            rec[f"{c}_p"] = t["p"] if t else None
            line += f"{pc:>10.1f}"
        print(line)
        rows.append(rec)

    if not rows:
        sys.exit("no realisation could be tested")

    alpha = 0.05 / (len(rows) * len(COMBOS))
    print("\n" + "=" * 72)
    print(f"  {len(rows)} realisations, Bonferroni alpha = {alpha:.2e}")
    print(f"\n  {'quantity':<12}{'median pct':>12}{'range':>16}{'significant':>14}")
    summary = {}
    for c in COMBOS:
        v = np.array([r[f"{c}_pct"] for r in rows])
        v = v[np.isfinite(v)]
        sig = sum(1 for r in rows
                  if r[f"{c}_p"] is not None and r[f"{c}_p"] < alpha)
        summary[c] = (float(np.median(v)), sig)
        print(f"  {c:<12}{np.median(v):>12.1f}"
              f"{f'{v.min():.1f}-{v.max():.1f}':>16}{f'{sig} of {len(rows)}':>14}")

    print()
    best = max((c for c in COMBOS if c not in ("age", "flux")),
               key=lambda c: summary[c][0])
    b_pct, b_sig = summary[best]
    a_pct, a_sig = summary["age"]
    print(f"  best combination: {best}, at {b_pct:.1f} in {b_sig} of {len(rows)}")
    print(f"  age alone:                 {a_pct:.1f} in {a_sig} of {len(rows)}")
    print()
    if b_pct > a_pct + 4 and b_sig >= a_sig:
        print("  The combination locates the precipitating cells better than the")
        print("  age alone, on both the percentile and the count. Worth reporting")
        print("  as an exploratory predictor alongside the frozen index.")
    elif b_pct > a_pct + 4:
        print("  The combination raises the percentile but not the number of")
        print("  realisations in which the elevation is significant. That is a")
        print("  weaker result than it appears and should be reported with both")
        print("  figures.")
    else:
        print("  No combination of the age and the flux improves on the age")
        print("  alone. Within a spatially uniform mineralogy the flow field")
        print("  does not locate the precipitation zone: the missing quantity is")
        print("  where the reactive minerals are, not how the water moves. That")
        print("  is the substance of Reviewer 1's first comment, and it is the")
        print("  honest answer to give the Associate Editor.")
    print()
    print("  Not pre-registered. Any positive result is exploratory.")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\n  wrote {a.csv}")


if __name__ == "__main__":
    main()
