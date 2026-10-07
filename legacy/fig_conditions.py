#!/usr/bin/env python3
"""
The conditions under which carbonate precipitates (AE-2, R2-6).

THE COMMENT
-----------
The associate editor asked that the flow paths, dissolution zones and
precipitation patterns be shown together, that the co-location index be either
replaced by something computable before the reactive simulation or reframed as
diagnostic, and that the handling of negligible-carbonate cases be explained.

THE ANSWER THIS FIGURE GIVES
----------------------------
Rather than a map, the figure reports where the precipitating cells sit within
the distributions of their own domain, for the three quantities that determine
the saturation index. Across the realisations that precipitate at all:

    cell flux             median percentile  62   elevated
    local anorthite loss  median percentile   2   strongly depressed
    pH                    median percentile  90   strongly elevated

The second is the informative one. Precipitation does not occur where anorthite
is dissolving; it occurs downstream of there. Anorthite dissolution consumes
protons, so a cell with high local dissolution is one where acid is still
arriving and being neutralised, and the fluid is still too acidic for CO3(2-)
to be present at the concentration calcite requires. A cell with low local
dissolution but high pH is one where the proton load was spent upstream, and
cations advected from those upstream surfaces can precipitate.

That is a co-location statement, but of flow with the *downstream* side of the
dissolution front rather than with the front itself.

PANELS
------
  (a) three empirical CDFs of percentile rank, precipitating cells against the
      uniform expectation for a random subset. Displacement above the diagonal
      means depressed, below means elevated. One panel, three curves, no
      binning and no spatial projection.
  (b) per-realisation percentiles for the three quantities, so the reader sees
      the result is not carried by one case.
  (c) the joint distribution of pH against flux for one realisation, with the
      precipitating cells overlaid, showing that they occupy the corner where
      both conditions are met.

Realisations with no measurable carbonate are reported explicitly in panel (b)
rather than omitted: they are the cases in which the conditions were never met,
which is the point.

Usage
-----
    python3 src/fig_conditions.py --runs runs --out figures/fig_conditions.pdf
    python3 src/fig_conditions.py --runs runs --case C_baseline__p32_150_s501 \\
        --out figures/fig_conditions.pdf
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# The four quantities that bear on whether carbonate can form. sat_rate is the
# net carbonate reaction rate: PFLOTRAN writes no saturation index, but the
# transition-state rate carries the factor (1 - Q/K), so its sign is the sign of
# the saturation index. It is the quantity that actually decides precipitation;
# the other three are the conditions that produce it.
QTY = (("sat_rate", "net carbonate rate (SI sign)", "#8c3b8c"),
       ("flux", "cell flux magnitude", "#12406e"),
       ("ph", "pH", "#4a8c5c"),
       ("diss", "anorthite consumed locally", "#c98a3a"))


def collect(runs, prefix, min_precip=5):
    from coloc_test import gather, pct_of, rank_test
    rows, empty, skipped = [], [], []
    for d in sorted(glob.glob(os.path.join(runs, prefix + "*"))):
        if not os.path.isdir(d):
            continue
        c, why = gather(d)
        if c is None:
            skipped.append((os.path.basename(d), why)); continue
        if c["ph"] is None:
            skipped.append((c["short"], "no pH field")); continue
        keep = ~c["unphys"]
        sel = ((c["carb"] > 0) & keep)[keep]
        if sel.sum() < min_precip:
            empty.append((c["short"], int(sel.sum()), c["n"])); continue
        rec = dict(short=c["short"], n=c["n"], n_precip=int(sel.sum()),
                   ranks={}, pct={}, test={})
        for key, _, _ in QTY:
            if c.get(key) is None:
                continue
            x = c[key][keep]
            r = np.argsort(np.argsort(x)) / max(len(x) - 1, 1)
            rec["ranks"][key] = r[sel]
            rec["pct"][key] = pct_of(x, sel)
            rec["test"][key] = rank_test(x, sel)
        rows.append(rec)
    return rows, empty, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--prefix", default="C_baseline__")
    ap.add_argument("--case", default=None,
                    help="realisation for panel (c); default is the one with "
                         "the most precipitating cells")
    ap.add_argument("--out", default="figures/fig_conditions.pdf")
    a = ap.parse_args()

    rows, empty, skipped = collect(a.runs, a.prefix)
    if not rows:
        sys.exit("no realisation has enough precipitating cells")
    for nm, why in skipped:
        print(f"  skipped {nm}: {why}")
    print(f"\n  {len(rows)} realisations with >=5 precipitating cells, "
          f"{len(empty)} with fewer")
    keys = [k for k, _, _ in QTY if all(k in r["pct"] for r in rows)]
    print(f"\n  {'case':<20}{'precip':>7}" +
          "".join(f"{lab[:13]:>15}" for k, lab, _ in QTY if k in keys))
    for r in rows:
        print(f"  {r['short']:<20}{r['n_precip']:>7}" +
              "".join(f"{r['pct'][k]:>15.1f}" for k in keys))
    print(f"  {'-'*20}{'-'*7}" + "-" * (15 * len(keys)))
    print(f"  {'median':<20}{'':>7}" +
          "".join(f"{np.median([r['pct'][k] for r in rows]):>15.1f}"
                  for k in keys))

    case = a.case or "C_baseline__" + max(rows, key=lambda r: r["n_precip"])["short"]

    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                                "xtick.labelsize": 7, "ytick.labelsize": 7})
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0],
                          width_ratios=[1.0, 1.0], hspace=0.34, wspace=0.28)

    # ---- (a) pooled CDFs of percentile rank -----------------------------
    ax = fig.add_subplot(gs[0, 0])
    for key, lab, col in QTY:
        if not all(key in r["ranks"] for r in rows):
            continue
        pooled = np.concatenate([r["ranks"][key] for r in rows])
        xs = np.sort(pooled)
        ax.plot(xs, np.linspace(0, 1, len(xs)), color=col, lw=1.5, label=lab)
    ax.plot([0, 1], [0, 1], color="0.45", lw=0.8, ls=":",
            label="a random subset of cells")
    ax.set_xlabel("percentile rank within the realisation", fontsize=8)
    ax.set_ylabel("cumulative fraction of precipitating cells", fontsize=8)
    ax.set_title("(a) precipitating cells are atypical", fontsize=8.5)
    ax.legend(fontsize=6.4, loc="upper left", framealpha=0.92)
    ax.grid(alpha=0.22, lw=0.4); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.97, 0.30, "above the diagonal:\ndepressed",
            transform=ax.transAxes, fontsize=6, ha="right", color="0.4")
    ax.text(0.97, 0.10, "below: elevated", transform=ax.transAxes,
            fontsize=6, ha="right", color="0.4")

    # ---- (b) per-realisation percentiles --------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    skey = "ph" if all("ph" in r["pct"] for r in rows) else "flux"
    order = sorted(rows, key=lambda r: r["pct"][skey])
    y = np.arange(len(order))
    for key, lab, col in QTY:
        if not all(key in r["pct"] for r in order):
            continue
        ax2.plot([r["pct"][key] for r in order], y, "o", ms=4.2, color=col,
                 mec="none", label=lab)
    ax2.axvline(50, color="0.45", lw=0.8, ls=":")
    ax2.set_yticks(y)
    ax2.set_yticklabels([r["short"].replace("p32_", "") for r in order],
                        fontsize=5.6)
    ax2.set_xlabel("percentile rank of the precipitating cells", fontsize=8)
    ax2.set_title(f"(b) each realisation separately (n={len(order)})",
                  fontsize=8.5)
    ax2.legend(fontsize=6.2, loc="upper left", ncol=1, framealpha=0.92,
               borderpad=0.35, handletextpad=0.4)
    ax2.grid(alpha=0.22, lw=0.4, axis="x"); ax2.set_xlim(-3, 103)
    if empty:
        # inside the axes, lower right, so it does not collide with the label
        ax2.text(0.985, 0.02,
                 f"{len(empty)} further realisations produced\n"
                 f"no measurable carbonate and so\nappear in neither panel",
                 transform=ax2.transAxes, fontsize=5.8, ha="right", va="bottom",
                 color="0.35", linespacing=1.35,
                 bbox=dict(fc="white", ec="0.85", lw=0.4, pad=2.0))

    # ---- (c) joint pH-flux for one realisation --------------------------
    from coloc_test import gather, pct_of, rank_test
    c, why = gather(os.path.join(a.runs, case))
    if c is None:
        print(f"  panel (c): {case}: {why}")
    else:
        keep = ~c["unphys"]
        sel = ((c["carb"] > 0) & keep)[keep]
        fl, ph, ds = c["flux"][keep], c["ph"][keep], c["diss"][keep]
        pos = fl > 0
        ax3 = fig.add_subplot(gs[1, :])
        if c.get("sat_rate") is not None:
            sr = c["sat_rate"][keep]
            # fraction of cells in each bin that are supersaturated: this shows
            # WHERE in the flux-pH plane carbonate can grow at all
            hb = ax3.hexbin(np.log10(fl[pos]), ph[pos], C=(sr[pos] > 0),
                            reduce_C_function=np.mean, gridsize=(70, 34),
                            cmap="BuPu", mincnt=1, linewidths=0, vmin=0, vmax=1)
            cb_lab = "fraction of cells supersaturated"
        else:
            hb = ax3.hexbin(np.log10(fl[pos]), ph[pos], gridsize=(70, 34),
                            cmap="Greys", bins="log", mincnt=1, linewidths=0)
            cb_lab = "cells per bin"
        ax3.plot(np.log10(fl[sel & pos]), ph[sel & pos], "o", ms=4.4,
                 mfc="none", mec="#c1272d", mew=0.9,
                 label=f"precipitating cells (n={int(sel.sum())})")
        ax3.axvline(np.log10(np.percentile(fl[pos], 50)), color="#12406e",
                    lw=0.7, ls="--")
        ax3.axhline(np.percentile(ph, 90), color="#4a8c5c", lw=0.7, ls="--")
        ax3.set_xlabel("log$_{10}$ cell flux magnitude (arbitrary units)",
                       fontsize=8)
        ax3.set_ylabel("pH at 50 yr", fontsize=8)
        ax3.set_title(f"(c) {c['short']}: precipitation requires flow "
                      f"(dashed blue: median) and a recovered pH "
                      f"(dashed green: 90th percentile)", fontsize=8.5)
        ax3.legend(fontsize=6.6, loc="lower left", framealpha=0.92)
        ax3.grid(alpha=0.2, lw=0.4)
        cb = fig.colorbar(hb, ax=ax3, pad=0.012, fraction=0.035)
        cb.set_label(cb_lab, fontsize=7); cb.ax.tick_params(labelsize=6)
        msg = (f"\n  panel (c) {c['short']}: flux pct {pct_of(fl, sel):.1f}, "
               f"diss pct {pct_of(ds, sel):.1f}, pH pct {pct_of(ph, sel):.1f}")
        if c.get("sat_rate") is not None:
            sr = c["sat_rate"][keep]
            msg += (f", SI>0 in {100*float((sr>0).mean()):.3f}% of cells; "
                    f"{100*float((sel & (sr>0)).sum()/max((sr>0).sum(),1)):.1f}% "
                    f"of supersaturated cells precipitated")
        print(msg)

    fig.subplots_adjust(left=0.085, right=0.985, top=0.945, bottom=0.075)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300)
    fig.savefig(a.out.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\n  wrote {a.out}")
    print(f"        {a.out.rsplit('.',1)[0]}.png")


if __name__ == "__main__":
    main()
