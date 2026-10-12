#!/usr/bin/env python3
"""
Ensemble statistics figures: flow-index test, carbonate and dissolution variability, and parameter sensitivity.

Three claims in the letter are currently carried by numbers in prose, and each
is better shown than stated:

  fig_flow_index_test.pdf    predictions recorded in advance against the carbonate
                      those simulations produced. n = 20, rho = +0.398.
                      Answers AE-1, R2-5 and R3-2.

  fig_carbonate_dissolution_variability.pdf  the central result. Carbonate per cell for every
                      realisation grouped by fracture intensity, beside the
                      primary dissolution on the same axis. The contrast is the
                      finding: dissolution varies by under four percentage
                      points while carbonate spans orders of magnitude with no
                      ordering. Answers AE-3, R2-3, R3-3 and R1-12.

  fig_parameter_sensitivity.pdf  the single-parameter variants as ratios to their
                      own baselines, with the factor-of-two band marked.
                      Answers AE-5, and supports R1-2c, R3-6 and R2-1.

The scripts reads what already exists: prereg_blockB.json and the Block B runs
for the first, final_blockA.csv for the second, sens_final.csv for the third.
Any input that is missing is reported and that figure skipped, so a partial run
still produces what it can.

Usage
-----
    python3 src/ensemble_statistics_figures.py --runs runs --out figures
"""
from __future__ import annotations
import argparse, csv, glob, json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6
RED, BLUE, GREY = "#b2182b", "#2166ac", "0.55"

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42,
})


def per_cell(run_dir, final_year=50.0):
    """Net carbonate volume fraction per cell, the metric named in the frozen file."""
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
        tot, n = 0.0, None
        for m in CARB:
            k = next((x for x in g if x.startswith(f"{m} VF")), None)
            if k:
                a = np.clip(np.asarray(g[k][:], float) - SEED_VF, 0.0, None)
                n = a.size
                tot += float(a.sum())
    return tot / n if n else None


# ---------------------------------------------------------------- figure 1
def fig_index(runs, out):
    """Predictions recorded in advance, against what the simulations produced."""
    from scipy import stats
    pre_path = "prereg_blockB.json"
    if not os.path.isfile(pre_path):
        return f"{pre_path} not found"
    pre = json.load(open(pre_path))
    rows = []
    for r in pre["predictions"]:
        v = per_cell(os.path.join(runs, "B_" + r["name"]))
        if v is None:
            continue
        rows.append((r["name"], float(r["shape99"]), v, float(r["p32_mult"])))
    if len(rows) < 5:
        return f"only {len(rows)} complete Block B runs"

    x = np.array([r[1] for r in rows])
    y = np.array([r[2] for r in rows])
    p32 = np.array([r[3] for r in rows])
    rho, p2 = stats.spearmanr(x, y)
    one = p2 / 2 if rho > 0 else 1 - p2 / 2

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))

    # (a) the values as measured
    ax = axes[0]
    zero = y <= 0
    floor = y[~zero].min() / 3 if (~zero).any() else 1e-8
    for lv, mk, c in ((0.90, "o", BLUE), (1.75, "s", RED)):
        m = (p32 == lv) & ~zero
        ax.scatter(x[m], y[m], s=26, marker=mk, facecolor=c, edgecolor="w",
                   linewidth=0.5, label=f"P$_{{32}}\\times${lv:g}", zorder=3)
    if zero.any():
        ax.scatter(x[zero], np.full(zero.sum(), floor), s=26, marker="v",
                   facecolor="none", edgecolor=GREY, linewidth=0.8,
                   label="no carbonate", zorder=3)
        ax.axhline(floor, lw=0.5, ls=":", color=GREY)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("index value, recorded before the simulations")
    ax.set_ylabel("carbonate per cell at 50 yr")
    ax.set_title("(a) prediction against outcome", fontsize=8.5, pad=4)
    ax.legend(fontsize=6.4, frameon=False, loc="lower right")

    # (b) ranks, since the test is a rank correlation
    ax = axes[1]
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    for lv, mk, c in ((0.90, "o", BLUE), (1.75, "s", RED)):
        m = p32 == lv
        ax.scatter(rx[m], ry[m], s=26, marker=mk, facecolor=c, edgecolor="w",
                   linewidth=0.5, zorder=3)
    lim = [0, len(rows) + 1]
    ax.plot(lim, lim, lw=0.6, ls="--", color=GREY, zorder=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("rank of the index"); ax.set_ylabel("rank of carbonate")
    ax.set_title("(b) the same data as ranks", fontsize=8.5, pad=4)
    ax.text(0.04, 0.96, f"$\\rho$ = {rho:+.3f}\none-sided $p$ = {one:.3f}\n"
            f"$n$ = {len(rows)}", transform=ax.transAxes, fontsize=7.2,
            va="top", ha="left")

    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_flow_index_test.pdf"))
    plt.close(fig)
    return None


# ---------------------------------------------------------------- figure 2
def fig_variability(out, csv_path="final_blockA.csv"):
    """The contrast between a tight dissolution response and a wide one."""
    if not os.path.isfile(csv_path):
        return f"{csv_path} not found"
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return f"{csv_path} is empty"
    cols = rows[0].keys()

    def pick(*cands):
        for c in cands:
            if c in cols:
                return c
        return None

    c_carb = pick("carb_per_cell", "carbonate_per_cell", "carb", "per_cell")
    c_diss = pick("diss_anorthite_pct", "anorthite_pct", "anorthite",
                  "diss_anorthite", "anorthite_dissolved_pct")
    # finalise.py records no intensity column, so it is read from the run
    # identifier as level_of() does there: p32_100_s613 -> 1.00.
    c_p32 = pick("p32_mult", "p32", "level", "mult")
    c_id = pick("run_id", "case", "name")
    if not c_carb or not (c_p32 or c_id):
        return ("could not identify the columns; found: "
                + ", ".join(sorted(cols)))

    def num(r, c):
        try:
            return float(r[c])
        except (TypeError, ValueError):
            return np.nan

    def p32_of(r):
        if c_p32:
            return num(r, c_p32)
        m = re.search(r"p32_(\d+)", r.get(c_id, "") or "")
        return int(m.group(1)) / 100.0 if m else np.nan

    lv = sorted({p32_of(r) for r in rows if np.isfinite(p32_of(r))})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))

    # (a) carbonate: wide, unordered
    ax = axes[0]
    allv = np.array([num(r, c_carb) for r in rows])
    pos = allv[np.isfinite(allv) & (allv > 0)]
    floor = pos.min() / 3 if pos.size else 1e-8
    for i, L in enumerate(lv):
        v = np.array([num(r, c_carb) for r in rows
                      if np.isclose(p32_of(r), L)])
        v = v[np.isfinite(v)]
        jit = (np.random.default_rng(i).random(v.size) - 0.5) * 0.30
        z = v <= 0
        ax.scatter(i + jit[~z], v[~z], s=20, facecolor=BLUE, edgecolor="w",
                   linewidth=0.4, zorder=3)
        if z.any():
            ax.scatter(i + jit[z], np.full(z.sum(), floor), s=20, marker="v",
                       facecolor="none", edgecolor=GREY, linewidth=0.7, zorder=3)
        if (~z).any():
            ax.plot([i - 0.28, i + 0.28], [np.median(v[~z])] * 2, lw=1.6,
                    color=RED, zorder=4)
    if pos.size:
        ax.axhline(np.median(pos), lw=0.6, ls="--", color=GREY, zorder=1)
    ax.set_yscale("log")
    ax.set_xticks(range(len(lv)))
    ax.set_xticklabels([f"$\\times${L:g}" for L in lv])
    ax.set_xlabel("fracture intensity P$_{32}$")
    ax.set_ylabel("carbonate per cell at 50 yr")
    ax.set_title("(a) carbonate: no ordering", fontsize=8.5, pad=4)

    # (b) dissolution of the three cation-supplying minerals: narrow, on a linear axis
    ax = axes[1]
    MIN = (("forsterite", "#1a9850"), ("diopside", "#8073ac"), ("anorthite", BLUE))
    c_min = {m: pick(f"{m}_dissolved_pct", f"diss_{m}_pct", f"{m}_pct") for m, _ in MIN}
    if not any(c_min.values()) and c_diss:
        c_min = {"anorthite": c_diss}
    shown = [(m, c) for m, c in MIN if c_min.get(m)]
    if shown:
        off = np.linspace(-0.24, 0.24, len(shown)) if len(shown) > 1 else [0.0]
        for i, L in enumerate(lv):
            for (m, col), dx in zip(shown, off):
                v = np.array([num(r, c_min[m]) for r in rows if np.isclose(p32_of(r), L)])
                v = v[np.isfinite(v)]
                jit = (np.random.default_rng(100 + i).random(v.size) - 0.5) * 0.10
                ax.scatter(i + dx + jit, v, s=12, facecolor=col, edgecolor="w", linewidth=0.3, zorder=3,
                           label=m if i == 0 else None)
                if v.size:
                    ax.plot([i + dx - 0.09, i + dx + 0.09], [np.median(v)] * 2, lw=1.4, color=RED, zorder=4)
        lo = min(np.nanmin([num(r, c_min[m]) for r in rows]) for m, _ in shown)
        ax.set_ylim(np.floor(lo) - 1.0, 101.4)
        ax.set_yticks(np.arange(np.floor(lo) - 1.0, 100.5, 1.0))     # no ticks above 100 %; the band above holds the legend
        ax.set_xticks(range(len(lv)))
        ax.set_xticklabels([f"$\\times${L:g}" for L in lv])
        ax.set_xlabel("fracture intensity P$_{32}$")
        ax.set_ylabel("primary mineral dissolved (\\%)")
        ax.set_title("(b) dissolution: nearly uniform", fontsize=8.5, pad=4)
        ax.legend(fontsize=6.2, frameon=False, loc="upper center", ncol=len(shown), handletextpad=0.2, columnspacing=1.0)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no dissolution column in\n" + csv_path,
                ha="center", va="center", fontsize=7, color=GREY)

    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_carbonate_dissolution_variability.pdf"))
    plt.close(fig)
    return None


# ---------------------------------------------------------------- figure 3
def fig_sensitivity(out, csv_path="sens_final.csv"):
    """Every variant as a ratio to its own baseline, on one axis."""
    if not os.path.isfile(csv_path):
        return f"{csv_path} not found"
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return f"{csv_path} is empty"
    cols = rows[0].keys()

    def pick(*cands):
        for c in cands:
            if c in cols:
                return c
        return None

    c_var = pick("variant", "name", "tag")
    c_rat = pick("ratio", "r", "factor")
    if not (c_var and c_rat):
        return "could not identify the columns; found: " + ", ".join(sorted(cols))

    agg = {}
    for r in rows:
        try:
            v = float(r[c_rat])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(v) or v <= 0:
            continue
        agg.setdefault(r[c_var], []).append(v)
    if not agg:
        return "no finite positive ratios"

    names = sorted(agg, key=lambda k: np.median(agg[k]))
    fig, ax = plt.subplots(figsize=(6.4, 0.30 * len(names) + 1.3))
    ax.axvspan(0.5, 2.0, color=BLUE, alpha=0.07, lw=0, zorder=0)
    ax.axvline(1.0, lw=0.8, color="0.35", zorder=1)
    for i, nm in enumerate(names):
        v = np.array(agg[nm])
        lo, hi, md = v.min(), v.max(), float(np.median(v))
        big = (md < 0.5) or (md > 2.0)
        c = RED if big else BLUE
        ax.plot([lo, hi], [i, i], lw=1.1, color=c, zorder=2,
                solid_capstyle="round")
        ax.scatter([md], [i], s=24, facecolor=c, edgecolor="w",
                   linewidth=0.5, zorder=3)
        ax.text(hi * 1.18, i, f"n={v.size}", fontsize=5.8, va="center",
                color="0.45")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=6.6)
    ax.set_xscale("log")
    # Keep every dot and its n= label inside the axes: margin below the lowest
    # value and room above the highest for the label written at 1.18 x.
    lo_all = min(min(agg[k]) for k in names); hi_all = max(max(agg[k]) for k in names)
    ax.set_xlim(min(lo_all, 0.5) / 1.8, max(hi_all, 2.0) * 2.4)
    ax.set_xlabel("carbonate relative to the same realisation's baseline")
    ax.set_ylim(-0.8, len(names) - 0.2)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
    ax.set_title("Single-parameter sensitivity; shaded band is a factor of two",
                 fontsize=8.5, pad=5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_parameter_sensitivity.pdf"))
    plt.close(fig)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--blocka", default="final_blockA.csv")
    ap.add_argument("--sens", default="sens_final.csv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    for lab, fn in (("fig_flow_index_test.pdf", lambda: fig_index(a.runs, a.out)),
                    ("fig_carbonate_dissolution_variability.pdf",
                     lambda: fig_variability(a.out, a.blocka)),
                    ("fig_parameter_sensitivity.pdf",
                     lambda: fig_sensitivity(a.out, a.sens))):
        try:
            why = fn()
        except Exception as e:
            why = f"{type(e).__name__}: {e}"
        if why:
            print(f"  skipped {lab}: {why}")
        else:
            print(f"  wrote   {os.path.join(a.out, lab)}")


if __name__ == "__main__":
    main()
