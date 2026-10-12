#!/usr/bin/env python3
"""
Is precipitation co-located with fast flow AND active dissolution? (AE-2, R2-6)

The associate editor asked whether the flow-reactive co-location index is more
than a description of where mineralisation happened. That is a statistical
question about the precipitating cells relative to the population, and it can be
answered directly rather than through a map.

WHY NOT A MAP OR A PROFILE
--------------------------
Two presentations were tried first and both fail on these data, for reasons that
are themselves the paper's finding:

  * A projection of the network collapses the fracture planes into overlapping
    streaks. Precipitation occupies 0.1% of cells and the dissolution field
    spans only 0.13 to 0.18 volume fraction, so the backdrop carries no
    information.
  * Binning by distance from the injection face destroys the signal entirely:
    the two realisations are indistinguishable in dissolution rate, cation
    complexes and pH, and a median over ~2000 cells per bin renders 87
    precipitating cells as exactly zero. Precipitation is not organised by
    distance. It is organised by which cells the flow connects to dissolution,
    which is precisely the claim under test.

WHAT THIS SCRIPT DOES
---------------------
    --mode rank     For every realisation, report where the precipitating cells
                    sit in the flux and dissolution distributions of their own
                    domain, as percentiles. If co-location holds they should be
                    high in both. A Mann-Whitney test against the non-
                    precipitating population gives a p-value, and the rank-
                    biserial correlation gives an effect size that does not
                    depend on the mesh size.

    --mode joint    Two-panel figure for a chosen realisation: flux against
                    dissolution for all cells as a density, with the
                    precipitating cells overlaid; and the marginal
                    distributions. If precipitation requires both conditions
                    the red points occupy the upper-right corner.

    --mode pick     Rank all realisations by how many cells precipitate, so the
                    clearest case can be chosen for the figure rather than
                    guessed at.

Cells whose secondary volume fraction exceeds unity are excluded throughout:
the fixed-porosity model does not prevent precipitation beyond the available
pore space, and those cells are not physical.

Usage
-----
    python3 src/coloc_test.py --runs runs --mode pick
    python3 src/coloc_test.py --runs runs --mode rank --csv coloc_rank.csv
    python3 src/coloc_test.py --runs runs --mode joint \\
        --case C_baseline__p32_200_s117 --out figures/fig_coloc_joint.pdf
"""
from __future__ import annotations
import argparse, csv, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEC = CARB + ("Kaolinite", "Chalcedony")
SEED_VF = 1e-6


def read_uge(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        xyz = np.empty((n, 3)); vol = np.empty(n)
        for i in range(n):
            t = f.readline().split()
            xyz[i] = (float(t[1]), float(t[2]), float(t[3])); vol[i] = float(t[4])
        m = int(f.readline().split()[1])
        ids = np.empty((m, 2), dtype=np.int64); area = np.empty(m)
        for j in range(m):
            t = f.readline().split()
            ids[j] = (int(t[0]) - 1, int(t[1]) - 1); area[j] = float(t[5])
    return xyz, vol, ids, area


def gather(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "no output"
    xyz, vol, ids, area = read_uge(uge)

    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, f"incomplete ({tg[-1][0]:g} yr)"
        g0, g1 = f[tg[0][1]], f[tg[-1][1]]

        def fld(grp, *pre):
            for p in pre:
                k = next((x for x in grp.keys() if x.startswith(p)), None)
                if k is not None:
                    return np.asarray(grp[k][:], dtype=float)
            return None

        press = fld(g1, "Liquid Pressure", "Pressure")
        anor0, anor1 = fld(g0, "Anorthite VF"), fld(g1, "Anorthite VF")
        diss_rate = fld(g1, "Anorthite Rate")
        ph = fld(g1, "pH")

        # Saturation state. PFLOTRAN does not write a saturation index, but the
        # transition-state rate carries the factor (1 - Q/K), so the SIGN of the
        # carbonate reaction rate is the sign of the saturation index: positive
        # where the fluid is supersaturated and the mineral is growing, negative
        # where it is undersaturated and dissolving. Summing over the four
        # carbonates gives the net carbonate saturation state of each cell,
        # which is the quantity that actually determines whether precipitation
        # can occur -- pH and cation supply are only proxies for it.
        sat_rate = None
        for m in CARB:
            v = fld(g1, f"{m} Rate")
            if v is not None:
                sat_rate = v if sat_rate is None else sat_rate + v
        # CO3(2-) activity, the other half of the ion activity product
        co3 = fld(g1, "CO3--")
        carb = np.zeros(len(xyz)); sec = np.zeros(len(xyz))
        for m in SEC:
            a = fld(g1, f"{m} VF")
            if a is None:
                continue
            net = np.clip(a - SEED_VF, 0.0, None)
            sec += net
            if m in CARB:
                carb += net

    if press is None:
        return None, "no pressure field in the output"
    if anor0 is None or anor1 is None:
        return None, "no anorthite field"

    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1); d[d <= 0] = np.nan
    q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                      posinf=0.0, neginf=0.0)
    flux = np.zeros(len(xyz))
    np.add.at(flux, i, np.abs(q)); np.add.at(flux, j, np.abs(q)); flux *= 0.5

    unphys = sec > 1.0
    return dict(name=os.path.basename(run_dir),
                short=os.path.basename(run_dir).replace("C_baseline__", ""),
                n=len(xyz), vol=vol, xyz=xyz, flux=flux, ph=ph,
                diss=anor0 - anor1,
                diss_rate=np.abs(diss_rate) if diss_rate is not None else None,
                sat_rate=sat_rate, co3=co3,
                carb=carb, unphys=unphys,
                n_carb=int(((carb > 0) & ~unphys).sum()),
                n_unphys=int(unphys.sum()),
                carb_per_cell=float(carb[~unphys].sum() / len(xyz))), None


def pct_of(x, sel):
    """Mean percentile rank of the selected cells within the full distribution."""
    r = np.argsort(np.argsort(x)) / max(len(x) - 1, 1) * 100.0
    return float(np.mean(r[sel]))


def rank_test(x, sel):
    """Mann-Whitney U with the rank-biserial correlation as effect size.

    Rank-biserial r = 2*AUC - 1, so r = 0 is no separation and r = 1 means
    every precipitating cell exceeds every other cell. Unlike a difference in
    means it does not depend on the units or on the mesh size.
    """
    from scipy import stats
    a, b = x[sel], x[~sel]
    if len(a) < 3 or len(b) < 3:
        return None
    u, p = stats.mannwhitneyu(a, b, alternative="greater")
    auc = u / (len(a) * len(b))
    return dict(p=float(p), r=float(2 * auc - 1), auc=float(auc))


def mode_pick(runs, prefix):
    rows = []
    for d in sorted(glob.glob(os.path.join(runs, prefix + "*"))):
        if not os.path.isdir(d):
            continue
        c, why = gather(d)
        if c is None:
            print(f"  {os.path.basename(d):<34} {why}")
            continue
        rows.append(c)
    # rank by the FRACTION of cells precipitating, not the raw count:
    # a 250k-cell mesh will out-count an 18k-cell mesh regardless of behaviour
    rows.sort(key=lambda c: -(c["n_carb"] / max(c["n"], 1)))
    print(f"\n{'case':<20}{'cells':>10}{'precip':>8}{'frac %':>9}"
          f"{'carb/cell':>12}{'VF>1':>6}")
    for c in rows:
        print(f"{c['short']:<20}{c['n']:>10,}{c['n_carb']:>8}"
              f"{100*c['n_carb']/c['n']:>9.3f}{c['carb_per_cell']:>12.3e}"
              f"{c['n_unphys']:>6}")
    if rows:
        best = rows[0]
        print(f"\n  highest precipitating fraction: {best['short']} "
              f"({best['n_carb']} of {best['n']:,} = "
              f"{100*best['n_carb']/best['n']:.3f}%)")
        by_count = max(rows, key=lambda c: c["n_carb"])
        if by_count is not best:
            print(f"  most precipitating cells:       {by_count['short']} "
                  f"({by_count['n_carb']})")
        print(f"  use:  --mode joint --case {best['name']}")
    return rows


def mode_rank(runs, prefix, csv_path=None):
    out = []
    print(f"{'case':<20}{'precip':>7}{'flux':>7}{'diss':>7}{'pH':>7}"
          f"{'satrate':>9}{'CO3':>7}   {'% of cells SI>0':>15}")
    for d in sorted(glob.glob(os.path.join(runs, prefix + "*"))):
        if not os.path.isdir(d):
            continue
        c, why = gather(d)
        if c is None:
            continue
        sel = (c["carb"] > 0) & ~c["unphys"]
        if sel.sum() < 5:
            print(f"{c['short']:<20}{int(sel.sum()):>7}   too few precipitating cells")
            continue
        keep = ~c["unphys"]
        s = sel[keep]
        rec = dict(case=c["short"], n_cells=c["n"], n_precip=int(sel.sum()),
                   carb_per_cell=c["carb_per_cell"])
        cells = []
        for key in ("flux", "diss", "ph", "sat_rate", "co3"):
            x = c[key]
            if x is None:
                rec[f"{key}_pct"] = rec[f"r_{key}"] = rec[f"p_{key}"] = None
                cells.append("     -")
                continue
            x = x[keep]
            t = rank_test(x, s)
            rec[f"{key}_pct"] = pct_of(x, s)
            rec[f"r_{key}"] = t["r"] if t else None
            rec[f"p_{key}"] = t["p"] if t else None
            cells.append(f"{rec[f'{key}_pct']:>6.1f}"
                         if key != "sat_rate" else
                         f"{rec[f'{key}_pct']:>8.1f}")
        # fraction of the domain that is supersaturated at 50 yr, and the
        # fraction of THOSE cells that actually precipitated
        if c["sat_rate"] is not None:
            sr = c["sat_rate"][keep]
            frac_sup = 100.0 * float((sr > 0).mean())
            rec["pct_supersat"] = frac_sup
            rec["precip_of_supersat"] = (100.0 * float((s & (sr > 0)).sum()
                                                       / max((sr > 0).sum(), 1)))
            tail = f"{frac_sup:>10.3f}"
        else:
            rec["pct_supersat"] = rec["precip_of_supersat"] = None
            tail = f"{'-':>10}"
        print(f"{c['short']:<20}{int(sel.sum()):>7}" + "".join(cells) + tail)
        out.append(rec)
    if out:
        print("\n  Percentile 50 means the precipitating cells are typical of "
              "the domain;")
        print("  above 50 means they favour higher values. r is the rank-"
              "biserial correlation:")
        print("  0 = no separation, 1 = every precipitating cell exceeds every "
              "other cell.")
        print(f"\n  across {len(out)} realisations, percentile rank of the "
              f"precipitating cells:")
        for key, lab in (("flux", "cell flux"), ("diss", "local dissolution"),
                         ("ph", "pH"), ("sat_rate", "carbonate rate (SI sign)"),
                         ("co3", "CO3(2-) concentration")):
            v = np.array([r[f"{key}_pct"] for r in out
                          if r.get(f"{key}_pct") is not None])
            if v.size:
                print(f"    {lab:<26} median {np.median(v):>5.1f}  "
                      f"range {v.min():>5.1f}-{v.max():>5.1f}")
        sup = np.array([r["pct_supersat"] for r in out
                        if r.get("pct_supersat") is not None])
        if sup.size:
            print(f"\n    supersaturated cells, % of domain: "
                  f"median {np.median(sup):.3f}  "
                  f"range {sup.min():.3f}-{sup.max():.3f}")
        # A percentile marginally above 50 is not evidence: use the test.
        # Bonferroni across the realisations, two tests each.
        alpha = 0.05 / (2 * len(out))
        def sig(r, k):
            p_ = r.get(f"p_{k}")
            return p_ is not None and p_ < alpha
        both = sum(1 for r in out if sig(r, "flux") and sig(r, "diss"))
        fonly = sum(1 for r in out if sig(r, "flux") and not sig(r, "diss"))
        donly = sum(1 for r in out if sig(r, "diss") and not sig(r, "flux"))
        neither = len(out) - both - fonly - donly
        sat = sum(1 for r in out if sig(r, "sat_rate"))
        phn = sum(1 for r in out if sig(r, "ph"))
        print(f"\n  at alpha = 0.05 Bonferroni-corrected for "
              f"{2*len(out)} tests ({alpha:.2e}):")
        print(f"    elevated in BOTH flux and dissolution : {both} of {len(out)}")
        print(f"    elevated in flux only                 : {fonly}")
        print(f"    elevated in dissolution only          : {donly}")
        print(f"    neither                               : {neither}")
        print(f"    elevated in pH                        : {phn} of {len(out)}")
        print(f"    elevated in carbonate rate (SI>0)     : {sat} of {len(out)}")
        print(f"\n  'both' supports co-location as stated. 'flux only' would")
        print(f"  mean precipitation follows the flow field irrespective of")
        print(f"  where dissolution occurs, which is a weaker claim.")
    if csv_path and out:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader()
            w.writerows(out)
        print(f"\n  wrote {csv_path}")
    return out


def mode_joint(runs, case, out_path):
    c, why = gather(os.path.join(runs, case))
    if c is None:
        sys.exit(f"{case}: {why}")
    keep = ~c["unphys"]
    sel = ((c["carb"] > 0) & keep)[keep]
    fl, ds = c["flux"][keep], c["diss"][keep]
    if sel.sum() < 5:
        sys.exit(f"{case}: only {int(sel.sum())} precipitating cells")

    pf, pd_ = pct_of(fl, sel), pct_of(ds, sel)
    tf, td = rank_test(fl, sel), rank_test(ds, sel)
    print(f"  {c['short']}: {int(sel.sum())} of {keep.sum():,} cells precipitate")
    print(f"    flux        percentile {pf:.1f}, r = {tf['r']:+.3f}, "
          f"p = {tf['p']:.2e}")
    print(f"    dissolution percentile {pd_:.1f}, r = {td['r']:+.3f}, "
          f"p = {td['p']:.2e}")

    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                                "xtick.labelsize": 7, "ytick.labelsize": 7})
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.2, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.30)

    # (a) joint distribution
    ax = fig.add_subplot(gs[0, 0])
    pos = fl > 0
    ax.hexbin(np.log10(fl[pos]), ds[pos], gridsize=55, cmap="Greys",
              bins="log", mincnt=1, linewidths=0)
    ax.plot(np.log10(fl[sel & pos]), ds[sel & pos], "o", ms=4.0,
            mfc="none", mec="#c1272d", mew=0.9,
            label=f"precipitating (n={int(sel.sum())})")
    ax.axvline(np.log10(np.percentile(fl[pos], 90)), color="#12406e",
               lw=0.7, ls="--")
    ax.axhline(np.percentile(ds, 90), color="#c98a3a", lw=0.7, ls="--")
    ax.text(0.985, 0.03, "dashed: 90th percentile", transform=ax.transAxes,
            fontsize=6, ha="right", color="0.4")
    ax.set_xlabel("log$_{10}$ cell flux magnitude (arb. units)", fontsize=8)
    ax.set_ylabel("anorthite volume fraction consumed", fontsize=8)
    ax.set_title(f"(a) {c['short']}: where precipitation sits", fontsize=8)
    ax.legend(fontsize=6.5, loc="upper left", framealpha=0.9)
    ax.grid(alpha=0.2, lw=0.4)

    # (b) marginals as empirical CDFs, which show the shift without binning
    ax2 = fig.add_subplot(gs[0, 1])
    for arr, lab, col, ls in ((fl, "flux", "#12406e", "-"),
                              (ds, "dissolution", "#c98a3a", "-")):
        r = np.argsort(np.argsort(arr)) / max(len(arr) - 1, 1)
        xs = np.sort(r); ax2.plot(xs, np.linspace(0, 1, len(xs)),
                                  color="0.7", lw=0.8, ls=ls)
        rs = np.sort(r[sel])
        ax2.plot(rs, np.linspace(0, 1, len(rs)), color=col, lw=1.4, ls=ls,
                 label=f"{lab}: precipitating cells")
    ax2.plot([0, 1], [0, 1], color="0.4", lw=0.6, ls=":",
             label="all cells (uniform by construction)")
    ax2.set_xlabel("percentile rank within the domain", fontsize=8)
    ax2.set_ylabel("cumulative fraction", fontsize=8)
    ax2.set_title("(b) are precipitating cells atypical?", fontsize=8)
    ax2.legend(fontsize=6.2, loc="upper left", framealpha=0.9)
    ax2.grid(alpha=0.2, lw=0.4)
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
    ax2.text(0.98, 0.06,
             f"flux: pct {pf:.0f}, r={tf['r']:+.2f}\n"
             f"diss: pct {pd_:.0f}, r={td['r']:+.2f}",
             transform=ax2.transAxes, fontsize=6.5, ha="right",
             bbox=dict(fc="white", ec="0.8", lw=0.4, pad=2.5))

    fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.145)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=300)
    fig.savefig(out_path.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\n  wrote {out_path}")
    print(f"        {out_path.rsplit('.',1)[0]}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--prefix", default="C_baseline__")
    ap.add_argument("--mode", choices=("pick", "rank", "joint"), default="pick")
    ap.add_argument("--case", default=None)
    ap.add_argument("--out", default="figures/fig_coloc_joint.pdf")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    if a.mode == "pick":
        mode_pick(a.runs, a.prefix)
    elif a.mode == "rank":
        mode_rank(a.runs, a.prefix, a.csv)
    else:
        if not a.case:
            sys.exit("--mode joint needs --case; run --mode pick first")
        mode_joint(a.runs, a.case, a.out)


if __name__ == "__main__":
    main()
