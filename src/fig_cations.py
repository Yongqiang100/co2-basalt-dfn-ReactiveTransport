#!/usr/bin/env python3
"""
Flow, cation supply, carbon supply and precipitation together (AE-2, R2-6).

THE COMMENT
-----------
The associate editor asked that the flow paths, the dissolution zones and the
precipitation patterns be shown together, and that the co-location argument be
made visible rather than asserted through a single correlation coefficient.

WHAT THIS FIGURE SHOWS
----------------------
One realisation in five columns. A second realisation can be added as a second
row with --low, for instance one that traps nothing.

  (a) Flow        the volumetric throughflow of each cell, obtained from the
                  simulated pressure field and the aperture-corrected
                  connection areas by Darcy's law. Relative units: k/mu is
                  absorbed into the scaling, so only ratios are meaningful.
  (b) Ca/Mg       the sum of the four dissolved Ca- and Mg-bearing carbonate
                  complexes written by PFLOTRAN: CaCO3(aq), CaHCO3+, MgCO3(aq)
                  and MgHCO3+. Free Ca2+ and Mg2+ are not in the output, so
                  these complexes stand for the dissolved cation supply. Log
                  scale: the field spans five orders of magnitude.
  (c) Carbon      HCO3(-) + CO3(2-), the delivered inorganic carbon. CO3(2-)
                  alone was tried and is set almost entirely by pH, so it
                  duplicates panel (d); the sum instead shows delivery.
  (d) pH          Carbonate stability requires the proton load to have been
                  consumed, so this is the third condition alongside cation and
                  carbon supply. Diverging colour scale pivoted at PH_PIVOT.
  (e) Co-location the fastest cells in grey, the cation-rich cells in blue and
                  the precipitating cells ringed in red, so the spatial
                  coincidence is visible in one frame

WHY THESE FIELDS
----------------
Panel (d) shows pH rather than the anorthite consumed: the dissolution extent
is nearly uniform across the domain (0.13 to 0.18 volume fraction) and carries
almost no spatial information, whereas pH varies over four units and is one of
the conditions that must be met.

Testing every field against the precipitating cells across 17 realisations gave:

    Ca complexes        median percentile 99.4   elevated in 17 of 17
    Mg complexes        median percentile 99.8   elevated in 17 of 17
    pH                  median percentile 89.6   elevated in 12 of 16
    cell flux           median percentile 62.5   elevated in 10 of 16
    local dissolution   median percentile  1.9   elevated in  0 of 16

The dissolved cation concentration is the one quantity elevated without
exception, by a median factor of about 200 and up to 2e5. Local dissolution is
anti-correlated, which is the signature of the mechanism rather than a
contradiction of it: dissolved cations arrive by advection from upstream
surfaces, so a precipitating cell has a high cation concentration and a low
local dissolution extent. That is transport-controlled cation supply.

Note the cation fields available are the complexes CaCO3(aq), CaHCO3+,
MgCO3(aq) and MgHCO3+; free Ca2+ and Mg2+ are not written to the output. The
complexes therefore depend on carbonate as well as on the cation, and the
figure labels them accordingly.

Usage
-----
    python3 src/fig_cations.py --runs runs \\
        --high C_baseline__p32_100_s383 --out figures/fig_cations.pdf

    # add --low <case> for a second row, e.g. a realisation that traps nothing
    python3 src/fig_cations.py --runs runs \\
        --high C_baseline__p32_100_s383 --low C_baseline__p32_100_s42 \\
        --out figures/fig_cations_pair.pdf
"""
from __future__ import annotations
import argparse, glob, os, sys
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEC = CARB + ("Kaolinite", "Chalcedony")
SEED_VF = 1e-6

# Colour scheme. Every panel uses the same red-to-blue family with RED at the
# high end, so the key is learned once. For the magnitude fields (flow,
# cations, carbon, age) red marks a large value. For pH red marks the alkaline
# side, on which carbonate is stable. pH uses a diverging map because the
# quantity has a chemically meaningful midpoint, pivoted at 6.3, above which
# CO3(2-) rises steeply enough for calcite to approach saturation at these
# calcium concentrations; a sequential map would imply that low pH is merely
# less of something.
CMAP_SEQ = "RdBu_r"        # magnitude fields
CMAP_PH = "RdBu_r"         # pH, pivoted at PH_PIVOT
PH_PIVOT = 6.3
# Set in main(): pivot of the pH scale ("mean": volume-weighted network mean, the
# definition of a high pH in the co-location mapping; "calcite": PH_PIVOT) and
# the time shown.
PH_MODE = "mean"
PH_LABEL = f"pH  (pivot {PH_PIVOT:g})"
TIME_TEXT = "All fields are at t = 50 yr."
DISS_LABEL = "Ca + Mg released (mol m$^{-3}$)"
OVERLAY_RICH = "#b2182b"   # the dark red end of the ramp, for the overlay

# Overlay thresholds for the co-location panel, named once so that the plotted
# selection and the caption cannot disagree. The age cut is 75 rather than 90
# because the precipitating cells average the 89th percentile of the age field
# and would otherwise fall just outside the highlighted population. The cation
# cut is 99.9 so that the highlighted count is comparable to the number of
# precipitating cells rather than an order of magnitude larger.
OVERLAY_AGE_PCT = 75.0
OVERLAY_CAT_PCT = 99.9

# Cells at or above this pH are drawn larger in panel (d). The value is the
# initial pH of the formation water, so the emphasised population is the part
# of the domain the injected acid has not depressed, which is where carbonate
# remains stable. At the common marker size those cells are lost among the
# acidified majority.
PH_EMPHASIS = 7.5


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


def _snaps(run_dir):
    """(time in years, file, group) of every snapshot over all output files of
    the run, sorted by time; a later file wins at a repeated time (continued runs)."""
    import h5py
    out = {}
    for p in sorted(glob.glob(os.path.join(run_dir, "*.h5"))):
        if "dfn_properties" in p:
            continue
        with h5py.File(p, "r") as f:
            for k in f.keys():
                m = re.search(r"Time\s+([0-9.eE+-]+)\s*y", k)
                if m:
                    out[round(float(m.group(1)), 8)] = (p, k)
    return [(t,) + out[t] for t in sorted(out)]


def _fields(path, key):
    """The fields the figure needs from one snapshot."""
    import h5py
    with h5py.File(path, "r") as f:
        grp = f[key]

        def fld(*pre):
            for p in pre:
                k = next((x for x in grp.keys() if x.startswith(p)), None)
                if k is not None:
                    return np.asarray(grp[k][:], dtype=float).ravel()
            return None

        def summed(*pre):
            t = None
            for p in pre:
                v = fld(p)
                if v is not None:
                    t = v if t is None else t + v
            return t
        vf = {m: fld(f"{m} VF") for m in SEC}
        return dict(press=fld("Liquid Pressure", "Pressure"), ph=fld("pH"),
                    ca=summed("CaCO3(aq)", "CaHCO3+"), mg=summed("MgCO3(aq)", "MgHCO3+"),
                    carbon=summed("HCO3-", "CO3--"), anor=fld("Anorthite VF"),
                    diop=fld("Diopside VF"), fors=fld("Forsterite VF"), vf=vf)


def _net(vf, minerals):
    t = None
    for m in minerals:
        a = vf.get(m)
        if a is None:
            continue
        a = np.clip(a - SEED_VF, 0.0, None)
        t = a if t is None else t + a
    return t


def gather(run_dir, final_year=50.0, want_age=True,
           library="../dfn_library", t_start=None, t_end=None):
    """Fields at t_start (default: the final snapshot) and the carbonate formed
    between t_start and t_end (default: net of the seed at the final snapshot)."""
    snaps = _snaps(run_dir)
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not snaps or not os.path.isfile(uge):
        return None, "no output"
    xyz, vol, ids, area = read_uge(uge)
    if snaps[-1][0] < final_year * 0.999:
        return None, f"incomplete ({snaps[-1][0]:g} yr)"

    def at(t):
        s_ = min(snaps, key=lambda x: abs(x[0] - t))
        if abs(s_[0] - t) > 0.01:
            raise SystemExit(f"{os.path.basename(run_dir)}: no snapshot at {t:g} yr (nearest {s_[0]:g})")
        return _fields(s_[1], s_[2])
    F0 = _fields(snaps[0][1], snaps[0][2])
    if t_start is None:
        F1 = FE = _fields(snaps[-1][1], snaps[-1][2])
        carb = _net(FE["vf"], CARB)
    else:
        F1, FE = at(t_start), at(t_end)
        carb = np.clip(_net(FE["vf"], CARB) - _net(F1["vf"], CARB), 0.0, None)
    carb = np.zeros(len(xyz)) if carb is None else carb
    sec = _net(FE["vf"], SEC)
    sec = np.zeros(len(xyz)) if sec is None else sec
    press, ph, ca, mg, carbon = F1["press"], F1["ph"], F1["ca"], F1["mg"], F1["carbon"]
    anor0, anor1 = F0["anor"], F1["anor"]
    # Ca + Mg released by dissolution over the period in which the ringed carbonate forms
    # (t_start to t_end; whole run without --time), per m3 of bulk rock. Molar volumes and
    # stoichiometry as in local_cation_balance.py: anorthite 1 Ca, diopside 1 Ca + 1 Mg,
    # forsterite 2 Mg. Regrowth counts as negative release and is clipped at zero.
    Fa, Fb = (F0, FE) if t_start is None else (F1, FE)
    release = np.zeros(len(xyz))
    for key, vm, nion in (("anor", 100.79e-6, 1), ("diop", 66.09e-6, 2), ("fors", 43.79e-6, 2)):
        if Fa.get(key) is not None and Fb.get(key) is not None:
            release += nion * (Fa[key] - Fb[key]) / vm
    release = np.clip(release, 0.0, None)

    if press is None:
        return None, "no pressure field in the output"
    for nm, v in (("Ca complexes", ca), ("Mg complexes", mg),
                  ("CO3(2-)", carbon)):
        if v is None:
            return None, f"no {nm} in the output"

    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1); d[d <= 0] = np.nan
    q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                      posinf=0.0, neginf=0.0)
    flux = np.zeros(len(xyz))
    np.add.at(flux, i, np.abs(q)); np.add.at(flux, j, np.abs(q)); flux *= 0.5

    # Mean groundwater age, from a steady Darcy solution on the same mesh. This
    # is the field from which the a priori index shape99 is computed, so
    # showing it connects the index to the precipitation sites. The
    # precipitating cells lie at a median percentile of 89.3 in it, elevated in
    # 17 of 19 realisations, against 62.0 for the local throughflow.
    age = None
    if want_age:
        try:
            import flowfield
            # The mesh directory is named for the realisation alone, while a
            # run directory carries its block prefix: C_baseline__p32_100_s383,
            # A_p32_100_s1181, D30_L30_p32_100_s613. Take the realisation part
            # wherever it sits rather than assuming a particular separator.
            nm = os.path.basename(run_dir.rstrip("/"))
            mm = re.search(r"(L\d+_)?p32_\d+_s\d+", nm)
            case = mm.group(0) if mm else nm.split("__")[-1]
            mesh = os.path.join(library, case)
            if not os.path.isdir(mesh):
                print(f"    no mesh at {mesh}; age panel unavailable")
            else:
                sol = flowfield.solve(mesh, verbose=False)
                a_ = np.asarray(sol["age"], float)
                if a_.size == len(xyz):
                    age = a_ / 3.156e7          # seconds -> years
        except Exception as e:
            print(f"    age field unavailable: {str(e)[:60]}")

    unphys = sec > 1.0
    # Age × flux rank product: the flow-only predictor combining residence
    # time (age) and delivery rate (flux). Computed as percentile 0–100.
    joint_rank = None
    if age is not None and flux is not None:
        from scipy.stats import rankdata
        r_age = rankdata(age, method="average")
        r_flux = rankdata(flux, method="average")
        joint_rank = 100.0 * rankdata(r_age * r_flux, method="average") / len(xyz)
    _nm = os.path.basename(run_dir.rstrip("/"))
    _mm = re.search(r"(L\d+_)?p32_\d+_s\d+", _nm)
    return dict(name=_mm.group(0) if _mm else _nm,
                age=age, joint_rank=joint_rank,
                xyz=xyz, vol=vol, flux=flux, ca=ca, mg=mg, cat=ca + mg,
                carbon=carbon, ph=ph,
                diss=(anor0 - anor1) if anor0 is not None else None, release=release,
                carb=carb, unphys=unphys,
                n_carb=int(((carb > 0) & ~unphys).sum()),
                carb_per_cell=float(carb[~unphys].sum() / len(xyz))), None



def elevated_panel(ax, xyz, v, *, cmap, thresh_mult=10.0, size=7.0,
                   vmin=None, vmax=None):
    """Grey the background, colour only the cells that are elevated.

    Three quarters of the cells sit at exactly the initial concentration --
    2.7e-8 M for the cation complexes, about 1e-4 M for inorganic carbon --
    because the injected fluid never reached them. On a five-decade colour
    scale those all pin to one end and the panel reads as a single colour,
    hiding the very cells that matter.

    So the background is drawn as small grey points and the colour ramp is
    spent entirely on cells above thresh_mult times the median, which is the
    delivered population. This is the same device already used in the overlay
    panel, and it makes the delivered regions the visual subject rather than
    an afterthought.

    Returns the mappable for the colour bar, and the threshold used.
    """
    from matplotlib.colors import LogNorm
    x, z = xyz[:, 0], xyz[:, 2]
    bg = float(np.median(v[np.isfinite(v)]))
    thr = bg * thresh_mult
    hi = v >= thr
    ax.scatter(x[~hi], z[~hi], s=0.45, c="0.86", linewidths=0, rasterized=True)
    sc = None
    if hi.any():
        lo = vmin if vmin is not None else thr
        top = vmax if vmax is not None else float(np.percentile(v[hi], 99.5))
        sc = ax.scatter(x[hi], z[hi], s=size, c=v[hi], cmap=cmap,
                        norm=LogNorm(vmin=lo, vmax=max(top, lo * 1.5)),
                        linewidths=0, rasterized=True)
    return sc, thr, int(hi.sum())



def upper_tail_panel(ax, xyz, v, *, cmap, pct, size=7.0, vmin=None, vmax=None):
    """Grey the lower part of a field, colour only its upper tail.

    Unlike the delivered concentrations, the age field is defined in every
    cell, so there is no background population in the sense used by
    elevated_panel. Colouring all of it makes the panel dense and gives most of
    the scale to water that is too young to matter. The cells below the
    percentile cut are therefore drawn grey and the ramp is spent on the tail.

    Using the same cut as the co-location overlay means this panel and the grey
    population in that overlay are the same set of cells.
    """
    from matplotlib.colors import LogNorm
    x, z = xyz[:, 0], xyz[:, 2]
    ok = np.isfinite(v) & (v > 0)
    thr = float(np.percentile(v[ok], pct))
    hi = ok & (v >= thr)
    ax.scatter(x[~hi], z[~hi], s=0.45, c="0.86", linewidths=0, rasterized=True)
    sc = None
    if hi.any():
        lo = vmin if vmin is not None else thr
        top = vmax if vmax is not None else float(np.percentile(v[hi], 99.5))
        sc = ax.scatter(x[hi], z[hi], s=size, c=v[hi], cmap=cmap,
                        norm=LogNorm(vmin=lo, vmax=max(top, lo * 1.5)),
                        linewidths=0, rasterized=True)
    return sc, thr, int(hi.sum())


def lognorm(v, lo_pct=5, hi_pct=99.9, vmin=None, vmax=None):
    from matplotlib.colors import LogNorm
    p = v[np.isfinite(v) & (v > 0)]
    lo = vmin if vmin is not None else np.percentile(p, lo_pct)
    hi = vmax if vmax is not None else np.percentile(p, hi_pct)
    return LogNorm(vmin=max(lo, hi * 1e-8), vmax=hi)


def mode_maps(cases, out, thresh_mult=10.0):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 13, "axes.linewidth": 0.5})
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    # shared scales so the two rows are comparable
    fl = np.concatenate([c["flux"][c["flux"] > 0] for c in cases])
    ct = np.concatenate([c["cat"][c["cat"] > 0] for c in cases])
    cb = np.concatenate([c["carbon"][c["carbon"] > 0] for c in cases])
    # Scales are shared BETWEEN THE TWO ROWS so the realisations can be
    # compared, but set independently FOR EACH FIELD from that field's own
    # distribution: flux spans about three decades, the cations five, the
    # carbon three, so one common scale would flatten most of them. Each colour
    # bar therefore carries its own numbers while the RdBu ramp is the same in
    # every panel. The cation and carbon fields are strongly skewed with the
    # informative values in the top few percent, hence the tighter low clip.
    # For the two delivered fields, fix the colour limits across both rows so
    # the plumes are directly comparable: low end at thresh_mult times the
    # pooled background, high end at the 99.5th percentile of the delivered
    # cells. Reported below so the threshold can be checked against the data.
    def delivered_limits(key):
        pooled = np.concatenate([c[key] for c in cases])
        bg = float(np.median(pooled))
        thr = bg * thresh_mult
        hi = pooled[pooled >= thr]
        top = float(np.percentile(hi, 99.5)) if hi.size else thr * 10
        print(f"    {key:<7} background {bg:.3e}, threshold {thr:.3e}, "
              f"{hi.size:,} of {pooled.size:,} cells above "
              f"({100*hi.size/pooled.size:.2f}%)")
        return thr, top
    print("  delivered-field thresholds:")
    cat_lim = delivered_limits("cat")
    cb_lim = delivered_limits("carbon")

    n_fl = lognorm(fl, lo_pct=15, hi_pct=99.5)
    n_ct = lognorm(ct, lo_pct=40, hi_pct=99.9)
    n_cb = lognorm(cb, lo_pct=40, hi_pct=99.9)
    # pH on a shared linear scale: carbonate stability is what matters, so the
    # range is fixed to the full span rather than clipped to percentiles
    global PH_LABEL
    phv = np.concatenate([c["ph"] for c in cases if c["ph"] is not None])
    if PH_MODE == "mean":
        # pivot at the volume-weighted network mean, the definition of a high pH in
        # the co-location mapping; range from the 1st to the 99th percentile
        pivot = float(np.mean([np.average(c["ph"], weights=c["vol"]) for c in cases]))
        ph_lo = float(np.floor(np.percentile(phv, 1) * 10) / 10)
        ph_hi = float(np.ceil(np.percentile(phv, 99) * 10) / 10)
        ph_lo = min(ph_lo, pivot - 0.1); ph_hi = max(ph_hi, pivot + 0.1)
        PH_LABEL = f"pH\n(pivot at the network mean, {pivot:.2f})"
    else:
        pivot = PH_PIVOT
        ph_lo = float(np.floor(phv.min() * 2) / 2)
        ph_hi = float(np.ceil(phv.max() * 2) / 2)
        # TwoSlopeNorm needs the pivot strictly inside the range
        ph_lo = min(ph_lo, PH_PIVOT - 0.5)
        ph_hi = max(ph_hi, PH_PIVOT + 0.5)
        PH_LABEL = f"pH  (pivot {PH_PIVOT:g})"
    print(f"    pH scale {ph_lo:g} to {ph_hi:g}, pivot {pivot:.2f} ({PH_MODE})")

    nrow = len(cases)
    npanel = 6 if all(c.get("joint_rank") is not None for c in cases) else 5
    # shared scale of the dissolution panel (every cell dissolves, so a full-field log scale as for the flux)
    rl = np.concatenate([c["release"][c["release"] > 0] for c in cases]) if all(c.get("release") is not None for c in cases) else np.array([1.0])
    n_rl = lognorm(rl, lo_pct=5, hi_pct=99.5)
    if npanel == 5 and any(c.get("joint_rank") is not None for c in cases):
        print("    joint_rank field missing for at least one case; "
              "falling back to five panels")

    if npanel == 6:
        titles = ("(a) Flow",
                  "(b) pH",
                  "(c) Ca/Mg complexes",
                  "(d) Age × Flux rank",
                  "(e) Carbon complexes",
                  "(f) Co-location")
    else:
        npanel = 6
        titles = ("(a) Flow",
                  "(b) pH",
                  "(c) Dissolution",
                  "(d) Ca/Mg complexes",
                  "(e) Carbon complexes",
                  "(f) Co-location")
    diss_layout = all(c.get("joint_rank") is None for c in cases)

    # Six panels in a 2x3 grid rather than a 1x6 strip: on a journal page the
    # strip forces each panel down to a width at which the fracture traces are
    # not resolvable. GRID_ROWS applies per realisation, so a second
    # realisation adds a further block of rows rather than a second strip.
    # One realisation is easier to read as a 2x3 grid, since a 1x6 strip on a
    # journal page leaves each panel too narrow to resolve the fracture traces.
    # Two realisations invert that: the point is to compare the same field
    # between them, so each becomes a single row of six and the rows sit
    # directly above one another.
    GRID_ROWS, GRID_COLS = (2, 3) if nrow == 1 else (1, npanel)
    nrow_ax = nrow * GRID_ROWS
    fig, axes = plt.subplots(nrow_ax, GRID_COLS,
                             figsize=(4.6 * GRID_COLS if GRID_COLS <= 3 else 2.5 * GRID_COLS,
                                      4.15 * nrow_ax if GRID_COLS <= 3
                                      else 3.2 * nrow_ax),
                             squeeze=False)
    # panel index -> (row, col) within this realisation's block
    def rc(row, k):
        return row * GRID_ROWS + k // GRID_COLS, k % GRID_COLS
    for row, c in enumerate(cases):
        x, z = c["xyz"][:, 0], c["xyz"][:, 2]
        sel = (c["carb"] > 0) & ~c["unphys"]
        # One sequential colour map for the three concentration/magnitude
        # fields so that "darker means more" reads the same in each. pH is the
        # exception: it is not a magnitude but a scale with a chemically
        # meaningful midpoint, so it gets a diverging map centred on the pH at
        # which carbonate becomes stable. Using the same sequential map for pH
        # would imply that low pH is simply "less of something", which is the
        # wrong reading.
        # (a) and (d) are fields every cell has, so every cell is coloured.
        # (b) and (e) are delivered quantities that most cells never receive,
        # so the background is greyed and the ramp is spent on the rest.
        # (d) is the age × flux rank product: a linear 0–100 field like pH,
        # with a midpoint at 50 (a cell at the median of both age and flux).
        from matplotlib.colors import TwoSlopeNorm as _TSN
        full = {0: (c["flux"], CMAP_SEQ, n_fl),
                1: (c["ph"], CMAP_PH,
                    TwoSlopeNorm(vmin=ph_lo, vcenter=pivot, vmax=ph_hi))}
        if npanel == 6 and c.get("joint_rank") is not None:
            full[3] = (c["joint_rank"], CMAP_SEQ,
                       Normalize(vmin=0, vmax=100))
        tail = {}
        if diss_layout:
            full[2] = (np.where(c["release"] > 0, c["release"], np.nan), CMAP_SEQ, n_rl)
            cat_col, carbon_col = 3, 4
        else:
            cat_col, carbon_col = 2, (4 if npanel == 6 else 3)
        elev = {cat_col: c["cat"], carbon_col: c["carbon"]}
        elim = {cat_col: cat_lim, carbon_col: cb_lim}
        for col in range(npanel - 1):
            ax = axes[rc(row, col)]
            if col in full:
                v, cmap, norm = full[col]
                sc = ax.scatter(x, z, s=0.55, c=v, cmap=cmap, norm=norm,
                                linewidths=0, rasterized=True)
                if col == 1:
                    # Redraw the strongly alkaline cells larger. They are few
                    # and they mark where carbonate is stable, so they should
                    # not be lost among the acidic majority.
                    em = np.isfinite(v) & (v >= PH_EMPHASIS)
                    if em.any():
                        ax.scatter(x[em], z[em], s=5.0, c=v[em], cmap=cmap,
                                   norm=norm, linewidths=0, rasterized=True,
                                   zorder=3)
                    # count for the caption (no in-panel note: it overlaps the data)
                    print(f"    {c['name']}: {int(em.sum()):,} cells at pH >= {PH_EMPHASIS:g}")
                    c["n_alkaline"] = int(em.sum())
                n_hi = None
            elif col in tail:
                sc, thr, n_hi = upper_tail_panel(
                    ax, c["xyz"], tail[col], cmap=CMAP_SEQ,
                    pct=OVERLAY_AGE_PCT,
                    vmin=age_lim[0], vmax=age_lim[1])
            else:
                sc, thr, n_hi = elevated_panel(
                    ax, c["xyz"], elev[col], cmap=CMAP_SEQ,
                    thresh_mult=thresh_mult,
                    vmin=elim[col][0], vmax=elim[col][1])
            style(ax, x, z, titles[col] if row == 0 else None)
            if col == 0 and row == 0:
                # x arrow BELOW the panel
                ax.annotate("", xy=(0.26, -0.06), xytext=(0.0, -0.06),
                            xycoords="axes fraction",
                            textcoords="axes fraction", annotation_clip=False,
                            arrowprops=dict(arrowstyle="-|>", lw=0.9,
                                            color="0.15"))
                ax.text(0.285, -0.06, "$x$", transform=ax.transAxes,
                        fontsize=13, va="center", ha="left", color="0.15",
                        clip_on=False)
                # z arrow on the right of panel (e), the last column
                ax_e = axes[rc(row, npanel - 1)]
                ax_e.annotate("", xy=(1.03, 0.30), xytext=(1.03, 0.04),
                            xycoords="axes fraction",
                            textcoords="axes fraction", annotation_clip=False,
                            arrowprops=dict(arrowstyle="-|>", lw=0.9,
                                            color="0.15"))
                ax_e.text(1.03, 0.325, "$z$", transform=ax_e.transAxes,
                        fontsize=13, ha="center", va="bottom", color="0.15",
                        clip_on=False)
            if n_hi is not None:
                # count for the caption (no in-panel note: it overlaps the data)
                print(f"    {c['name']}, {titles[col]}: {n_hi:,} cells above background")
            if row == 0 and sc is not None:
                globals().setdefault("_scs", {})[col] = sc

        # (e) the overlay: grey background + precipitating cells only
        ax = axes[rc(row, npanel - 1)]
        ax.scatter(x, z, s=0.15, c="0.90", linewidths=0, rasterized=True,
                   zorder=1)
        if sel.any():
            n_sel = int(sel.sum())
            if n_sel <= 150:
                ring, lw, halo = 46.0, 1.3, True
            elif n_sel <= 600:
                ring, lw, halo = 20.0, 0.9, True
            else:
                ring, lw, halo = 7.0, 0.5, False
            if halo:
                ax.scatter(x[sel], z[sel], s=ring, facecolors="none",
                           edgecolors="w", linewidths=lw + 1.1, zorder=4)
            ax.scatter(x[sel], z[sel], s=ring, facecolors="none",
                       edgecolors="#d62728", linewidths=lw, zorder=5,
                       label="Precipitating")
        else:
            ax.scatter([], [], s=20, facecolors="none",
                       edgecolors="#d62728", linewidths=0.9,
                       label="Precipitating")
        style(ax, x, z, titles[npanel - 1] if row == 0 else None)
        if row == 0:
            _overlay_handles, _overlay_labels = ax.get_legend_handles_labels()

        # the realisation name goes on the left of the row as a text label,
        # because panel (a)'s y axis now carries the z direction instead
        lab = f"{c['name']}\n({c['tag']})" if nrow > 1 else c["name"]
        a0 = axes[rc(row, 0)]
        a0.text(-0.115, 0.5, lab, transform=a0.transAxes,
                fontsize=13, rotation=90, ha="center", va="center")

    if npanel == 6 and not diss_layout:
        labs = ("Flow rate (relative)",
                PH_LABEL,
                "Ca/Mg complexes (mol/L)",
                "Age × Flux rank (%)",
                "HCO$_3^{-}$ + CO$_3^{2-}$ (mol/L)")
    else:
        labs = ("Flow rate (relative)",
                PH_LABEL,
                DISS_LABEL,
                "Ca/Mg complexes (mol/L)",
                "HCO$_3^{-}$ + CO$_3^{2-}$ (mol/L)")
    # The panel geometry must be final before the colour bar axes are
    # positioned from it: subplots_adjust moves the panels but leaves any
    # axes already added by add_axes where they were.
    fig.subplots_adjust(left=0.045, right=0.99, top=0.965, bottom=0.11,
                        wspace=0.06, hspace=0.02)

    # Shared legend for the co-location overlay, below the bottom overlay panel
    if '_overlay_handles' in dir() and _overlay_handles:
        bb = axes[rc(nrow - 1, npanel - 1)].get_position()
        from matplotlib.lines import Line2D
        clean_handles = []
        for h, l in zip(_overlay_handles, _overlay_labels):
            clean_handles.append(Line2D([0], [0], marker="o", color="w",
                                        markerfacecolor="none",
                                        markeredgecolor="#d62728",
                                        markeredgewidth=1.0,
                                        markersize=6, label=l))
        fig.legend(handles=clean_handles, loc="lower center",
                   bbox_to_anchor=(bb.x0 + bb.width / 2, bb.y0 - 0.065),
                   ncol=1, fontsize=12, frameon=False,
                   handletextpad=0.3, borderpad=0.3)
    for col, lab in enumerate(labs[:npanel - 1]):
        sc = globals().get("_scs", {}).get(col)
        if sc is None:
            continue
        # Each bar sits immediately below the panel it belongs to, which on a
        # grid is no longer the bottom of the figure.
        bb = axes[rc(nrow - 1, col)].get_position()
        bw = bb.width * 0.74
        cax = fig.add_axes([bb.x0 + (bb.width - bw) / 2,
                            bb.y0 - 0.030, bw, 0.010])
        cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cbar.set_label(lab, fontsize=11)
        cbar.ax.tick_params(labelsize=10, length=2, pad=1.5)
        # choose the ticks explicitly: the default locator either labels one
        # value on a narrow log range or collides on a wide one
        from matplotlib.ticker import LogFormatterMathtext, FixedLocator
        norm = sc.norm
        if hasattr(norm, "vmin") and type(norm).__name__ == "LogNorm":
            t = log_ticks(norm.vmin, norm.vmax, want=4)
            if t:
                cbar.ax.xaxis.set_major_locator(FixedLocator(t))
                cbar.ax.xaxis.set_minor_locator(FixedLocator([]))
                import math as _m
                def _lab(v):
                    e = _m.floor(_m.log10(v))
                    man = v / 10.0 ** e
                    # whole decades as 10^e, otherwise m x 10^e; the default
                    # LogFormatterMathtext writes 10^{-2.30} for the latter
                    if abs(man - 1.0) < 1e-6:
                        return rf"$10^{{{e:d}}}$"
                    return rf"${man:.0f}\times10^{{{e:d}}}$"
                cbar.ax.set_xticklabels([_lab(v) for v in t])
        else:
            cbar.ax.locator_params(axis="x", nbins=5)
    # No footnote inside the figure: everything it said belongs in the
    # caption, which is printed at the end of this function so it can be
    # pasted straight into the manuscript.
    save(fig, out)



def log_ticks(vmin, vmax, want=4):
    """Tick values for a logarithmic colour bar, chosen not to collide.

    Left to itself matplotlib puts a single label on a narrow log range and
    several overlapping ones on a wide range, because the decade locator does
    not know how little horizontal space a colour bar has. This picks at most
    `want` values: whole decades when the range covers enough of them,
    otherwise 1-2-5 multiples inside the range.
    """
    import math
    if not (vmin > 0 and vmax > vmin):
        return None
    lo, hi = math.log10(vmin), math.log10(vmax)
    decades = [10.0 ** e for e in range(math.ceil(lo), math.floor(hi) + 1)]
    if len(decades) >= 2:
        if len(decades) > want:                    # thin them evenly
            step = math.ceil(len(decades) / want)
            decades = decades[::step]
        return decades
    cand = [m * 10.0 ** e
            for e in range(math.floor(lo) - 1, math.ceil(hi) + 1)
            for m in (1, 2, 5)]
    inside = [c for c in cand if vmin <= c <= vmax]
    if len(inside) > want:
        step = math.ceil(len(inside) / want)
        inside = inside[::step]
    return inside or None


def style(ax, x, z, title):
    ax.axvspan(x.min(), x.min() + 0.2 * (x.max() - x.min()),
               color="#4a7fb5", alpha=0.06, lw=0, zorder=0)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    if title:
        ax.set_title(title, fontsize=14, pad=4)


def mode_profile(cases, out, nbins=45):
    """Cation and carbon supply against distance downstream.

    Plotted as the 95th percentile in each bin rather than the median: the
    quantities of interest are confined to a small fraction of cells, so a
    median over ~2000 cells per bin reports the background and hides them.
    """
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                                "xtick.labelsize": 7, "ytick.labelsize": 7})
    import matplotlib.pyplot as plt

    COL = {"most carbonate": "#1b5e9c", "least carbonate": "#b8532a"}
    x0 = min(c["xyz"][:, 0].min() for c in cases)
    x1 = max(c["xyz"][:, 0].max() for c in cases)
    edges = np.linspace(x0, x1, nbins + 1)
    ctr = 0.5 * (edges[:-1] + edges[1:]) - x0

    rows = [("dissolved Ca + Mg complexes", "concentration (M)", "cat", True),
            ("dissolved inorganic carbon", "concentration (M)", "carbon", True),
            ("pH", "pH", "ph", False),
            ("carbonate formed", "VF", "carb", True)]
    fig, axes = plt.subplots(len(rows), 1, figsize=(6.4, 1.45 * len(rows)),
                             sharex=True)
    for ax, (title, ylab, key, logy) in zip(axes, rows):
        for c in cases:
            v = c.get(key)
            if v is None:
                continue
            keep = ~c["unphys"]
            idx = np.digitize(c["xyz"][:, 0], edges) - 1
            hi = np.full(nbins, np.nan); md = np.full(nbins, np.nan)
            for b in range(nbins):
                s = (idx == b) & keep
                if s.sum() >= 3:
                    hi[b] = np.percentile(v[s], 95)
                    md[b] = np.median(v[s])
            ax.plot(ctr, hi, lw=1.3, color=COL[c["tag"]],
                    label=f"{c['name']} ({c['tag']})")
            ax.plot(ctr, md, lw=0.7, color=COL[c["tag"]], ls=":", alpha=0.8)
        if logy:
            ax.set_yscale("log")
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.text(0.006, 0.90, title, transform=ax.transAxes, fontsize=8,
                va="top", fontweight="bold")
        ax.grid(alpha=0.25, lw=0.4)
        ax.axvspan(0, 0.2 * (x1 - x0), color="#4a7fb5", alpha=0.06, lw=0)
    axes[0].legend(fontsize=6.4, loc="lower right", framealpha=0.9)
    axes[-1].set_xlabel("distance from the injection face (m)", fontsize=8)
    fig.text(0.5, 0.013,
             "Solid: 95th percentile in each bin. Dotted: median. The "
             "quantities are confined to a small fraction of cells, so the "
             "median reports the background.",
             fontsize=6.3, ha="center", color="0.3")
    fig.subplots_adjust(left=0.125, right=0.985, top=0.99, bottom=0.115,
                        hspace=0.12)
    save(fig, out)


def caption(cases, thresh_mult):
    """The caption text, printed so it can be pasted into the manuscript."""
    names = " and ".join(c["name"] for c in cases)
    n = cases[0]
    sel = (n["carb"] > 0) & ~n["unphys"]
    return (
        f"Figure. Flow, chemical conditions, dissolution and precipitation in realisations {names}. "
        f"(a) Volumetric throughflow of each cell, from the simulated pressure field by Darcy's law, in relative units. "
        f"(b) {PH_LABEL.replace(chr(10), ' ')}, on a diverging scale; the cells at pH {PH_EMPHASIS:g} and above are drawn larger "
        f"({', '.join(str(c.get('n_alkaline', 0)) for c in cases)} cells). "
        f"(c) Ca and Mg released by the dissolution of anorthite, diopside and forsterite, per cubic metre of rock; dissolution is strongest at the injection face and decreases downstream. "
        f"(d) Dissolved Ca and Mg, as the sum of the complexes CaCO3(aq), CaHCO3+, MgCO3(aq) and MgHCO3+. "
        f"(e) Dissolved inorganic carbon, HCO3- + CO3(2-). "
        f"(f) The cells that form carbonate, ringed in red ({', '.join(str(int(((c['carb'] > 0) & ~c['unphys']).sum())) for c in cases)} cells). "
        f"{TIME_TEXT} The projection is onto the x-z plane, with x running from the injection face to the outflow; "
        f"the shaded strip marks the injection region. In panels (d) and (e), grey denotes cells below "
        f"{thresh_mult:g} times the field median, and the colour scale is reserved for the cells above it. "
        f"Red denotes the high end of every field, including the alkaline side of the pH scale. "
        f"Colour limits are shared between the realisations and set independently for each field.")


def save(fig, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=300)
    fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\nwrote {out}")
    print(f"      {out.rsplit('.',1)[0]}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--high", default="C_baseline__p32_150_s117")
    ap.add_argument("--low", default=None,
                    help="optional second realisation, drawn as a second row")
    ap.add_argument("--out", default=None)
    ap.add_argument("--profile", action="store_true",
                    help="downstream profiles instead of maps")
    ap.add_argument("--library", default="../dfn_library",
                    help="DFN mesh directory, for the groundwater age solve")
    ap.add_argument("--no-age", action="store_true",
                    help="omit the age panel and use flux in the overlay")
    ap.add_argument("--thresh", type=float, default=3.0,
                    help="a cell counts as having received the injectate when "
                         "its concentration exceeds this multiple of the "
                         "field median; lower to show more, raise to show less")
    ap.add_argument("--time", type=float, default=None,
                    help="draw the fields at this time (yr) and ring the cells that form "
                         "carbonate between --time and --time-end; default: 50 yr, net carbonate")
    ap.add_argument("--time-end", type=float, default=None)
    ap.add_argument("--ph-pivot", choices=("mean", "calcite"), default="mean",
                    help="pivot of the pH scale: the volume-weighted network mean (default) "
                         f"or the calcite pivot {PH_PIVOT:g}")
    a = ap.parse_args()
    global PH_MODE, TIME_TEXT
    PH_MODE = a.ph_pivot
    if a.time is not None:
        if a.time_end is None or a.time_end <= a.time:
            sys.exit("--time-end must be given and later than --time")
        TIME_TEXT = (f"All fields are at t = {a.time:g} yr; the dissolution and the ringed cells "
                     f"cover {a.time:g} to {a.time_end:g} yr.")
        global DISS_LABEL
        DISS_LABEL = f"Ca + Mg released\n({a.time:g}-{a.time_end:g} yr, mol m$^{{-3}}$)"
    out = a.out or ("figures/fig_cation_profile.pdf" if a.profile
                    else "figures/fig_cations.pdf")

    wanted = [("most carbonate", a.high)]
    if a.low:
        wanted.append(("least carbonate", a.low))
    cases = []
    for tag, name in wanted:
        c, why = gather(os.path.join(a.runs, name),
                        want_age=not a.no_age, library=a.library,
                        t_start=a.time, t_end=a.time_end)
        if c is None:
            sys.exit(f"{name}: {why}")
        c["tag"] = tag
        cases.append(c)
        sel = (c["carb"] > 0) & ~c["unphys"]
        r = ((np.median(c["cat"][sel]) / np.median(c["cat"][~sel]))
             if sel.any() else float("nan"))
        print(f"  {c['name']:<18} {int(sel.sum()):>5} precipitating cells, "
              f"Ca+Mg contrast x{r:,.0f}" if sel.any() else
              f"  {c['name']:<18} no precipitation")

    (mode_profile(cases, out) if a.profile
     else mode_maps(cases, out, thresh_mult=a.thresh))
    if not a.profile:
        print("\n--- caption " + "-" * 58)
        print(caption(cases, a.thresh))
        print("-" * 70)


if __name__ == "__main__":
    main()
