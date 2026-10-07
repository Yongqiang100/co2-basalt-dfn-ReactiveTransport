#!/usr/bin/env python3
"""
Quantitative results cited in the revision, computed from the simulation outputs.

    python3 src/revision_statistics.py                                   # all sections
    python3 src/revision_statistics.py --only injection_duration shutin_porewater_carbon

Sections, each writing results/<section>.csv:
  volume_effect              carbonate, efficiency, CO2 per pore volume and pore volume against P32
  domain_size_variability    CV at 20, 30 and 40 m, independent bootstrap and permutation test
  sensitivity_ratios         each coupled Block C variant against the baseline, by network: carbonate
                             per cell and intensity as in blockc2.py, and total CO2
  injection_duration         network s1181: water, carbon, carbonate and efficiency from the mass
                             balance, dissolution, calcite share, calcium share of released Ca + Mg
  shutin_porewater_carbon    shut-in (E) against continuous injection (A), and carbonate against
                             the dissolved carbon in the pore water at shut-in
  porosity_limit_pairs       runs without the porosity limit against their reruns with it
  platform_comparison        identical runs on Setonix, hpc02 and hpc01
  flow_index_surface_area    flow index against efficiency with constant and shrinking surface areas
  clogging_cells             cells at porosity 0.02 or less at 50 y, and their new mineral volume
  shrinking_surface_ratios   carbonate with shrinking surface areas against constant, by network
  anorthite_surface_area     anorthite constraint line in each Block C deck
  silicate_regrowth          diopside and albite growth after their minimum (shut-in and duration runs)
Uses results/*_networks.csv written by summarize_block.py, and the run folders on hpc02.
"""
import argparse, csv, glob, os, re, statistics as st, sys
import numpy as np, h5py
from scipy.stats import spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import summarize_block as sb

A = argparse.Namespace()
R = lambda p: p if os.path.isabs(p) else os.path.join(REV, p)
NET = re.compile(r"p32_(\d{3})_s(\d+)")


def csv_rows(name, folder="results"):
    p = R(os.path.join(folder, f"{name}_networks.csv"))
    return list(csv.DictReader(open(p))) if os.path.isfile(p) else []


def write(key, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    p = R(f"results/{key}.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"  written: results/{os.path.basename(p)}")


def rng_text(v, f="{:.3g}"):
    return f"median {f.format(st.median(v))}, range {f.format(min(v))} to {f.format(max(v))}" if v else "no values"


def mas_value_at(d, t_want, col):
    best = None
    for p in cb._parts(d, "pflotran_co2-mas", ".dat"):
        cols, rows = cb._mas_rows(p)
        i = next((i for i, c in enumerate(cols) if col in c), None)
        if i is None:
            continue
        for r in rows:
            if best is None or abs(float(r[0]) - t_want) < abs(best[0] - t_want):
                best = (float(r[0]), float(r[i]))
    return best[1] if best else None


def co2_conc(d):
    m = re.search(r"CONSTRAINT\s+co2_rich_water.*?CO2\(aq\)\s+([0-9.eEdD+-]+)", open(os.path.join(d, "pflotran_co2.in")).read(), re.S | re.I)
    return float(m.group(1).lower().replace("d", "e"))


# ---------------------------------------------------------------- a
def volume_effect():
    print("\nVolume effect: Spearman correlation with P32")
    out = []
    for name in ("A_coupled", "B_coupled"):
        rows = csv_rows(name)
        if not rows:
            print(f"  {name}: results file missing"); continue
        p32, co2, eff, pv = [], [], [], []
        for r in rows:
            w0 = mas_value_at(os.path.join(R(A.root), r["network"]), 0.0, "Global Water Mass")
            if w0 is None:
                continue
            p32.append(float(r["p32_factor"])); co2.append(float(r["co2_kg"])); eff.append(float(r["efficiency_pct"])); pv.append(w0 / 1000.0)
        for lab, v in (("carbonate (kg CO2)", co2), ("efficiency", eff), ("CO2 per m3 of pore space", [c / x for c, x in zip(co2, pv)]), ("pore volume", pv)):
            rho = spearmanr(p32, v)
            print(f"  {name}: {lab:26s} rho = {rho.correlation:+.3f}  p = {rho.pvalue:.2g}  n = {len(v)}")
            out.append(dict(set=name, quantity=lab, rho=rho.correlation, p=rho.pvalue, n=len(v)))
    write("volume_effect", out)


# ---------------------------------------------------------------- b
def domain_size_variability():
    print("\nDomain size, coupled, P32 x 1.00 (independent samples)")
    rng = np.random.default_rng(1)
    def eff(name):
        return np.array([float(r["efficiency_pct"]) for r in csv_rows(name)
                         if abs(float(r["p32_factor"]) - 1.0) < 1e-6 and float(r["co2_kg"]) > 0])
    X = {L: eff(n) for L, n in (("20", "A_coupled"), ("30", "D30_coupled"), ("40", "D40_coupled"))}
    short = [f"{L} m: {len(x)}" for L, x in X.items() if len(x) < 3]
    if short:
        print("  fewer than 3 networks at P32 x 1.00 (" + ", ".join(short) + "): check the results files"); return
    cv = lambda x: x.std(ddof=1) / x.mean(); out = []
    for L, x in X.items():
        b = [cv(rng.choice(x, len(x))) for _ in range(10000)]
        lo, hi = np.percentile(b, [2.5, 97.5])
        print(f"  {L} m: n = {len(x)}, median {np.median(x):.3g}%, CV {100 * cv(x):.0f}% (95% interval {100 * lo:.0f} to {100 * hi:.0f}%)")
        out.append(dict(item=f"{L} m", n=len(x), median_eff_pct=np.median(x), cv_pct=100 * cv(x), cv_lo=100 * lo, cv_hi=100 * hi))
    for a_, b_ in (("20", "40"), ("30", "40"), ("20", "30")):
        r = np.array([cv(rng.choice(X[b_], len(X[b_]))) / cv(rng.choice(X[a_], len(X[a_]))) for _ in range(10000)])
        pooled = np.concatenate([X[a_], X[b_]]); obs = cv(X[b_]) / cv(X[a_]); n_a = len(X[a_])
        perm = np.array([cv(q[n_a:]) / cv(q[:n_a]) for q in (rng.permutation(pooled) for _ in range(10000))])
        lo, hi = np.percentile(r, [2.5, 97.5]); pp = (perm <= obs).mean()
        print(f"  CV {b_} m / CV {a_} m: {obs:.2f} (bootstrap 95% {lo:.2f} to {hi:.2f}), permutation p (lower at {b_} m) = {pp:.3f}")
        out.append(dict(item=f"CV {b_} m / CV {a_} m", ratio=obs, ratio_lo=lo, ratio_hi=hi, permutation_p=pp))
    write("domain_size_variability", out)


# ---------------------------------------------------------------- c
PHI = 0.50


def per_cell_intensity(d):
    """blockc2.py definitions at 50 y: carbonate per cell (sum over carbonates and cells of the
    volume fraction above the seed, divided by the cell count) and intensity (volume-weighted
    carbonate divided by the pore volume at porosity 0.50), plus the dominant carbonate."""
    t, p, k = cb.snapshots(d)[-1]; vol = cb.cell_volumes(d); per_vf, per_vol = {}, {}
    with h5py.File(p, "r") as f:
        g = f[k]
        for m in cb.CARB:
            n = next((x for x in g if x.startswith(f"{m} VF")), None)
            a = np.clip(np.asarray(g[n][:], float).ravel() - cb.SEED, 0, None) if n else np.zeros(len(vol))
            per_vf[m] = float(a.sum()); per_vol[m] = float((a * vol).sum())
    tot = sum(per_vol.values()); dom = max(per_vol, key=per_vol.get) if tot > 0 else "-"
    return dict(years=t, per_cell=sum(per_vf.values()) / len(vol), intensity=tot / float((vol * PHI).sum()),
                dominant=dom, dominant_pct=100 * per_vol.get(dom, 0) / tot if tot > 0 else 0.0)


def sensitivity_ratios():
    print("\nBlock C sensitivity, coupled, variant / baseline by network (definitions of blockc2.py)")
    root = R(A.root); key = lambda d: re.search(r"__(p32_\d{3}_s\d+)$", d).group(1)
    base = {key(d): (per_cell_intensity(d), sb.network(d)) for d in glob.glob(os.path.join(root, "C_baseline_feedback__p32_*"))}
    variants = sorted({re.match(r"C_(.+)_feedback__", os.path.basename(d)).group(1)
                       for d in glob.glob(os.path.join(root, "C_*_feedback__p32_*"))} - {"baseline"})
    out = []
    for v in variants:
        q = []
        for d in sorted(glob.glob(os.path.join(root, f"C_{v}_feedback__p32_*"))):
            n = key(d); b = base.get(n)
            if not b:
                continue
            r = per_cell_intensity(d); rs = sb.network(d); bi, bs = b
            if r["years"] < 49.99 or bi["years"] < 49.99 or bi["per_cell"] <= 0:
                continue
            q.append(r["per_cell"] / bi["per_cell"])
            out.append(dict(variant=v, network=n, ratio=q[-1], intensity_ratio=r["intensity"] / bi["intensity"],
                            co2_ratio=rs["co2_kg"] / bs["co2_kg"] if bs["co2_kg"] > 0 else float("nan"),
                            base_dominant=f"{bi['dominant']} {bi['dominant_pct']:.0f}%", var_dominant=f"{r['dominant']} {r['dominant_pct']:.0f}%",
                            dawsonite_pct=rs["dawsonite_pct"], anorthite_dissolved_pct=rs["anorthite_dissolved_pct"],
                            base_anorthite_pct=bs["anorthite_dissolved_pct"]))
        if q:
            flag = "within a factor of 2" if (max(q) <= 2.0 and min(q) >= 0.5) else "exceeds a factor of 2"
            print(f"  {v:18s} n = {len(q)}  per-cell ratio x{min(q):.2f} to x{max(q):.2f}, median {st.median(q):.2f}  ({flag})")
    write("sensitivity_ratios", out)


# ---------------------------------------------------------------- d
DUR = (("1 day", "F_1d"), ("10 days", "F_10d"), ("30 days", "F_30d"), ("45 days", "F_45d"),
       ("2 years", "F_2yr"), ("5 years", "F_5yr"), ("10 years", "E"), ("50 years", "A"))


def dur_dirs(pre, n="p32_100_s1181"):
    root = R(A.root)
    return (os.path.join(root, f"{pre}{'__' if '_' in pre else '_'}{n}"), os.path.join(root, f"{pre}_feedback__{n}"))


def released(r):
    g = lambda k: r.get(k) if r.get(k) == r.get(k) and r.get(k) is not None else 0.0
    ca = g("anorthite_dissolved_mol") + g("diopside_dissolved_mol")
    mg = 2 * g("forsterite_dissolved_mol") + g("diopside_dissolved_mol") + g("enstatite_dissolved_mol")
    return ca, mg


def injection_duration():
    print("\nInjection duration, network s1181, P32 x 1.00 (mass balance)")
    out = []
    print(f"  {'injection':>9} {'config':>8} {'water kg':>11} {'C in mol':>11} {'carbonate mol':>14} {'eff %':>9} {'fo %':>6} {'di %':>6} {'an %':>6} {'Ca share':>9} {'calcite':>8}")
    for lab, pre in DUR:
        for cfg, d in zip(("fixed", "coupled"), dur_dirs(pre)):
            if not os.path.isdir(d):
                continue
            w = cb.injected_water_mas(d); c = cb.carbonate_mas(d)
            if w is None or c is None:
                continue
            cin = w * co2_conc(d); r = sb.network(d) or {}; ca, mg = released(r)
            row = dict(injection=lab, config=cfg, water_kg=w, carbon_in_mol=cin, carbonate_mol=c, efficiency_pct=100 * c / cin,
                       forsterite_pct=r.get("forsterite_dissolved_pct"), diopside_pct=r.get("diopside_dissolved_pct"),
                       anorthite_pct=r.get("anorthite_dissolved_pct"), ca_share_released_pct=100 * ca / (ca + mg) if ca + mg else float("nan"),
                       calcite_pct=r.get("calcite_pct"))
            out.append(row)
            print(f"  {lab:>9} {cfg:>8} {w:11.5g} {cin:11.5g} {c:14.5g} {row['efficiency_pct']:9.4g} {row['forsterite_pct'] or 0:6.1f} "
                  f"{row['diopside_pct'] or 0:6.1f} {row['anorthite_pct'] or 0:6.1f} {row['ca_share_released_pct']:8.0f}% {row['calcite_pct'] or 0:7.0f}%")
    write("injection_duration", out)


# ---------------------------------------------------------------- e
def shutin_porewater_carbon():
    print("\nShut-in (coupled E) against continuous injection (coupled A), mass balance")
    root = R(A.root); out = []; ratio, rel, grow, p32 = [], [], [], []
    for d in sorted(glob.glob(os.path.join(root, "E_feedback__p32_*"))):
        n = d.split("__")[-1]; a = os.path.join(root, f"A_feedback__{n}")
        ce, ca = cb.carbonate_mas(d), cb.carbonate_mas(a)
        w0, w10 = mas_value_at(d, 0.0, "Global Water Mass"), mas_value_at(d, 10.0, "Global Water Mass")
        if ce is None or w10 is None:
            continue
        q = ce / (w10 * co2_conc(d)); rel.append(q); grow.append(w10 / w0); p32.append(int(NET.search(n).group(1)) / 100)
        rr = ce / ca if ca else float("nan")
        if ca: ratio.append(rr)
        out.append(dict(network=n, carbonate_E_mol=ce, carbonate_A_mol=ca, ratio_E_over_A=rr, carbonate_over_porewater_carbon=q, porewater_10y_over_0y=w10 / w0))
    print(f"  carbonate, shut-in / continuous: {rng_text(ratio, '{:.2f}')}, n = {len(ratio)}")
    print(f"  carbonate / carbon in pore water at shut-in: {rng_text(rel, '{:.3f}')}, n = {len(rel)}")
    if len(rel) > 3:
        rho = spearmanr(p32, rel); print(f"  that ratio vs P32: rho = {rho.correlation:+.2f}, p = {rho.pvalue:.2g}")
    print(f"  pore water at 10 y / at 0 y: {rng_text(grow, '{:.3f}')}")
    write("shutin_porewater_carbon", out)


# ---------------------------------------------------------------- f, g
def pair_change(label, a_root, b_root, names):
    q = []
    for n in names:
        a, b = cb.carbonate_mas(os.path.join(a_root, n)), cb.carbonate_mas(os.path.join(b_root, n))
        if a and b:
            q.append(dict(pair=label, network=n, a_mol=a, b_mol=b, change_pct=100 * (b - a) / a))
    v = [x["change_pct"] for x in q]
    print(f"  {label}: {len(v)} pairs, change {rng_text(v, '{:+.2f}')}%" if v else f"  {label}: no complete pairs")
    return q


def porosity_limit_pairs():
    print("\nPorosity limit: runs without the limit against their reruns with it")
    pre = R(A.prefloor)
    write("porosity_limit_pairs", pair_change("with / without the limit", pre, R(A.root), sorted(os.listdir(pre)) if os.path.isdir(pre) else []))


def platform_comparison():
    print("\nPlatform comparison: identical runs on different systems")
    ov = R(A.overlap); names = sorted(os.listdir(ov)) if os.path.isdir(ov) else []
    out = pair_change("Setonix -> hpc02, coupled A x2.00", ov, R(A.root), [n for n in names if n.startswith("A_feedback__")])
    out += pair_change("Setonix -> hpc02, fixed C", ov, R(A.root), [n for n in names if n.startswith("C_") and "_feedback__" not in n])
    a = {r["network"].split("__")[1]: float(r["co2_kg"]) for r in csv_rows("A_coupled")}
    q = [dict(pair="hpc01 S0 -> Block A", network=r["network"], change_pct=100 * (float(r["co2_kg"]) - a[k]) / a[k])
         for r in csv_rows("S0", A.hpc01) for k in [r["network"].split("__")[1]] if k in a]
    v = [x["change_pct"] for x in q]
    print(f"  hpc01 S0 -> Block A: {len(v)} pairs, change {rng_text(v, '{:+.2f}')}%" if v else "  hpc01 S0 -> Block A: no pairs")
    write("platform_comparison", out + q)


# ---------------------------------------------------------------- h
def flow_index_surface_area():
    print("\nFlow index against efficiency with constant and shrinking surface areas")
    p = R("results/flow_index_test.csv")
    if not os.path.isfile(p):
        print("  results/flow_index_test.csv missing"); return
    idx = {r["network"].split("__")[1]: float(r["shape99"]) for r in csv.DictReader(open(p)) if r["block"] == "A"}
    out = []
    for name in ("S0", "S1"):
        rows = [(idx[k], float(r["efficiency_pct"])) for r in csv_rows(name, A.hpc01) for k in [r["network"].split("__")[1]] if k in idx]
        if len(rows) > 3:
            rho = spearmanr([x for x, _ in rows], [y for _, y in rows])
            print(f"  {name}: rho = {rho.correlation:+.2f}, p = {rho.pvalue:.2g}, n = {len(rows)}")
            out.append(dict(set=name, rho=rho.correlation, p=rho.pvalue, n=len(rows)))
    write("flow_index_surface_area", out)


# ---------------------------------------------------------------- i
def clogging_cells():
    print("\nClogging cells in coupled Block A (porosity at or below 0.02 at 50 y)")
    tot = {"carbonate": 0.0, "Kaolinite": 0.0, "Chalcedony": 0.0}; ncl = nall = 0; out = []
    for d in sorted(glob.glob(os.path.join(R(A.root), "A_feedback__p32_*"))):
        t, p, k = cb.snapshots(d)[-1]; vol = cb.cell_volumes(d)
        with h5py.File(p, "r") as f:
            g = f[k]
            v = lambda m: next((np.asarray(g[n][:], float).ravel() for n in g if n.startswith(m + " VF")), np.zeros(len(vol)))
            phi = 1 - sum(np.asarray(g[n][:], float).ravel() for n in g if " VF" in n)
            cl = phi <= 0.02; ncl += int(cl.sum()); nall += len(phi)
            this = {"carbonate": sum(((v(m) - cb.SEED) * vol)[cl].sum() for m in cb.CARB)}
            for m in ("Kaolinite", "Chalcedony"):
                this[m] = ((v(m) - cb.SEED) * vol)[cl].sum()
        for m in tot:
            tot[m] += this[m]
        out.append(dict(network=os.path.basename(d), clog_cells=int(cl.sum()), **{f"{m}_m3": this[m] for m in this}))
    s = sum(tot.values())
    print(f"  {ncl} of {nall} cells; new mineral volume: " + ", ".join(f"{m} {100 * x / s:.1f}%" for m, x in tot.items()) if s else "  none")
    write("clogging_cells", out)


# ---------------------------------------------------------------- j
def domain_totals(d, minerals):
    vol = cb.cell_volumes(d); res = []
    files = {}
    try:
        for t, p, k in cb.snapshots(d):
            f = files.setdefault(p, h5py.File(p, "r")); g = f[k]
            res.append((t, {m: float(sum((np.asarray(g[n][:], float).ravel() * vol).sum()
                                         for n in g if " VF" in n and n.split(" VF")[0] == m)) for m in minerals}))
    finally:
        for f in files.values():
            f.close()
    return res


def growth(series, m):
    """(final - minimum) / initial, %, and the time of the minimum"""
    v = [x[1][m] for x in series]
    i = int(np.argmin(v))
    return (100 * (v[-1] - v[i]) / v[0] if v[0] else float("nan")), series[i][0]


def shrinking_surface_ratios():
    print("\nShrinking surface areas: carbonate S1 / S0 by network (hpc01)")
    s0 = {r["network"].split("__")[1]: r for r in csv_rows("S0", A.hpc01)}; out = []; q = []
    for r in csv_rows("S1", A.hpc01):
        k = r["network"].split("__")[1]; a = s0.get(k)
        if a and float(a["co2_kg"]) > 0:
            q.append(float(r["co2_kg"]) / float(a["co2_kg"]))
            out.append(dict(network=k, ratio=q[-1], anorthite_S0=float(a["anorthite_dissolved_pct"]), anorthite_S1=float(r["anorthite_dissolved_pct"])))
            print(f"  {k}: {q[-1]:.2f}   anorthite {float(a['anorthite_dissolved_pct']):.1f} -> {float(r['anorthite_dissolved_pct']):.1f}%")
    if q:
        print(f"  ratio {rng_text(q, '{:.2f}')}, n = {len(q)}")
    write("shrinking_surface_ratios", out)


def anorthite_surface_area():
    print("\nAnorthite in the Block C decks (constraint line: volume fraction, surface area, units)")
    seen = {}; out = []
    for d in sorted(glob.glob(os.path.join(R(A.root), "C_*_feedback__p32_*"))):
        v = re.match(r"C_(.+)_feedback__", os.path.basename(d)).group(1)
        if v != "baseline" and not v.startswith(("anor_as", "global_as")):
            continue
        lines = [l.strip() for l in open(os.path.join(d, "pflotran_co2.in")).read().split("\n") if re.match(r"\s*Anorthite\s+[0-9.]", l)]
        seen.setdefault(v, set()).update(lines)
    for v, ls in seen.items():
        print(f"  {v:18s} " + " | ".join(sorted(ls)))
        out.append(dict(variant=v, lines=" | ".join(sorted(ls))))
    write("anorthite_surface_area", out)


def silicate_regrowth():
    print("\nDiopside and albite growth after their minimum (domain totals, % of the initial volume)")
    out = []
    for label, pattern in (("shut-in E", "E_feedback__p32_*"), ("injection duration F", "F_*_feedback__p32_100_s1181")):
        gd, ga = [], []
        for d in sorted(glob.glob(os.path.join(R(A.root), pattern))):
            try:
                ser = domain_totals(d, ("Diopside", "Albite"))
            except Exception as e:
                print(f"  skipped {os.path.basename(d)}: {e}"); continue
            (g1, t1), (g2, t2) = growth(ser, "Diopside"), growth(ser, "Albite")
            gd.append(g1); ga.append(g2)
            out.append(dict(set=label, network=os.path.basename(d), diopside_growth_pct=g1, diopside_min_at_y=t1,
                            albite_growth_pct=g2, albite_min_at_y=t2))
            if label.startswith("injection"):
                print(f"  {os.path.basename(d):32s} diopside +{g1:.2f}% after its minimum at {t1:g} y,  albite +{g2:.2f}% after {t2:g} y")
        if gd:
            print(f"  {label}: diopside growth {rng_text(gd, '{:.2f}')}%, albite growth {rng_text(ga, '{:.2f}')}%, n = {len(gd)}")
    write("silicate_regrowth", out)


SECTIONS = dict(volume_effect=volume_effect, domain_size_variability=domain_size_variability, sensitivity_ratios=sensitivity_ratios,
                injection_duration=injection_duration, shutin_porewater_carbon=shutin_porewater_carbon,
                porosity_limit_pairs=porosity_limit_pairs, platform_comparison=platform_comparison,
                flow_index_surface_area=flow_index_surface_area, clogging_cells=clogging_cells,
                shrinking_surface_ratios=shrinking_surface_ratios, anorthite_surface_area=anorthite_surface_area,
                silicate_regrowth=silicate_regrowth)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--overlap", default="runs_gravityoff_setonix_overlap")
    ap.add_argument("--prefloor", default="runs_gravityoff_prefloor")
    ap.add_argument("--hpc01", default="archive_hpc01/results", help="folder with S0_networks.csv and S1_networks.csv")
    ap.add_argument("--only", nargs="+", choices=list(SECTIONS), default=list(SECTIONS))
    global A
    A = ap.parse_args()
    for k in A.only:
        try:
            SECTIONS[k]()
        except Exception as e:
            print(f"  {k} failed: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
