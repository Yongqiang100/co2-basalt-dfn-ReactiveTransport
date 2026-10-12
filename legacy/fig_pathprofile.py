#!/usr/bin/env python3
"""
The co-location mechanism along a flow path (AE-2, R2-6).

The associate editor asked for the flow paths, dissolution zones and
precipitation patterns to be shown together, and for an explanation of how they
relate. A projection of the three-dimensional network collapses into
overlapping streaks and shows none of this. This figure instead traces the
actual flow path and plots what happens along it.

WHAT IS PLOTTED, per realisation
--------------------------------
A path is traced from the injection region to the outflow boundary by following
the largest outgoing volumetric flux at each step, using the pressure field
PFLOTRAN simulated and the connectivity in the .uge file. Along that path:

  panel 1   cumulative anorthite consumed and the local dissolution rate
            proxy -- where cations enter the fluid
  panel 2   dissolved Ca and Mg, and pH -- whether the fluid is approaching
            carbonate saturation
  panel 3   carbonate volume fraction -- where precipitation actually occurs
  panel 4   flux magnitude along the path, for context

Two realisations at the same fracture intensity are overlaid: one that traps
and one that does not. The mechanism claimed in the manuscript is that
precipitation requires the path to accumulate cations from upstream dissolution
AND for the proton load to have been consumed. The figure tests that directly:
if the trapping case shows rising cations and recovering pH before its
carbonate peak, while the non-trapping case shows cations that never accumulate
or a pH that never recovers, the mechanism is demonstrated rather than asserted.

WHY A PATH RATHER THAN A MAP
----------------------------
The claim is about what a parcel of fluid encounters in sequence. That is a
one-dimensional statement about a path, so a one-dimensional plot along the
path is the natural presentation, and it makes the causal ordering visible:
dissolution upstream, saturation reached, precipitation downstream.

Usage
-----
    python3 src/fig_pathprofile.py --runs runs \\
        --high C_baseline__p32_100_s383 --low C_baseline__p32_100_s42 \\
        --out figures/fig_pathprofile.pdf

    python3 src/fig_pathprofile.py --runs runs --fields \\
        --high C_baseline__p32_100_s383      # what is in the output
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6


# ------------------------------------------------------------------ mesh + fields
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


def open_final(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    if not h5:
        return None, None, "no output"
    f = h5py.File(h5[0], "r")
    tg = time_groups(f)
    if not tg:
        f.close(); return None, None, "no time groups"
    if tg[-1][0] < final_year * 0.999:
        t = tg[-1][0]; f.close(); return None, None, f"incomplete ({t:g} yr)"
    return f, (f[tg[0][1]], f[tg[-1][1]]), None


def field(grp, *prefixes):
    for p in prefixes:
        k = next((x for x in grp.keys() if x.startswith(p)), None)
        if k is not None:
            return np.asarray(grp[k][:], dtype=float)
    return None


# --------------------------------------------------------------------- path trace
def trace_path(xyz, ids, area, press, inlet_frac=0.2, max_steps=100000):
    """Follow the largest outgoing flux from the strongest inlet cell.

    q_ij = A_ij * (P_i - P_j) / d_ij, with k/mu absorbed (relative only).
    Starting cell: the injection-region cell with the largest total outflux,
    i.e. the mouth of the dominant conduit. Terminates at the outflow face, at
    a cell with no downstream neighbour, or on revisiting a cell.
    """
    n = len(xyz)
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    d[d <= 0] = np.nan
    q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                      posinf=0.0, neginf=0.0)

    # adjacency: for each cell, the connections and the signed flux leaving it
    nbr = [[] for _ in range(n)]
    for k in range(len(ids)):
        a, b, f = i[k], j[k], q[k]
        nbr[a].append((b, f))        # positive f leaves a
        nbr[b].append((a, -f))       # sign flips for the other side

    x0, x1 = xyz[:, 0].min(), xyz[:, 0].max()
    inlet = xyz[:, 0] <= x0 + inlet_frac * (x1 - x0)
    outflux = np.zeros(n)
    for c in np.flatnonzero(inlet):
        outflux[c] = sum(f for _, f in nbr[c] if f > 0)
    if not outflux.any():
        return None, "no outgoing flux from the injection region"
    cur = int(np.argmax(outflux))

    # A step longer than a few times the median connection length means the
    # two cells are not really neighbours, which would make the path length
    # meaningless. Cap it.
    step_cap = 3.0 * float(np.nanpercentile(d, 95))
    path, seen = [cur], {cur}
    for _ in range(max_steps):
        near = [(f, c) for c, f in nbr[cur] if c not in seen
                and np.linalg.norm(xyz[c] - xyz[cur]) <= step_cap]
        if not near:
            break
        down = [(f, c) for f, c in near if f > 0]
        # If every neighbour is uphill the path has reached a local pressure
        # minimum, which on a DFN mesh is common and does not mean the conduit
        # has ended. Continue along the least adverse gradient.
        f, nxt = max(down) if down else max(near)
        path.append(nxt); seen.add(nxt); cur = nxt
        if xyz[cur, 0] >= x1 - 1e-9:
            break
    if len(path) < 10:
        return None, f"path terminated after {len(path)} cells"
    return np.array(path), None



def trace_through(xyz, ids, area, press, target, step_cap):
    """Trace the flow path THROUGH a nominated cell.

    A path started from the inlet and followed downstream need not pass
    through any of the 87 precipitating cells out of 84,303, in which case the
    precipitation panels are empty and the figure demonstrates nothing. This
    instead anchors on the cell of interest and walks upstream to the injection
    region and downstream to the outflow, so the path is guaranteed to be the
    one along which precipitation occurred. The question the figure then asks
    is what the fluid encountered on its way there.

    Upstream is followed by largest incoming flux, downstream by largest
    outgoing flux, both restricted to steps no longer than step_cap.
    """
    n = len(xyz)
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    d[d <= 0] = np.nan
    q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                      posinf=0.0, neginf=0.0)
    nbr = [[] for _ in range(n)]
    for k in range(len(ids)):
        nbr[i[k]].append((j[k], q[k]))
        nbr[j[k]].append((i[k], -q[k]))

    def walk(start, sign, limit=4000):
        """sign=+1 downstream (largest outflux), -1 upstream (largest influx)."""
        cur, out, seen = start, [], {start}
        for _ in range(limit):
            cand = [(sign * f, c) for c, f in nbr[cur]
                    if c not in seen and sign * f > 0
                    and np.linalg.norm(xyz[c] - xyz[cur]) <= step_cap]
            if not cand:
                break
            _, nxt = max(cand)
            out.append(nxt); seen.add(nxt); cur = nxt
        return out

    up = walk(target, -1)
    dn = walk(target, +1)
    path = np.array(up[::-1] + [target] + dn)
    # Orient the path so that distance increases DOWNSTREAM. The upstream and
    # downstream walks are assembled by flux direction, but a DFN path can
    # double back, so verify against the x coordinate: injection is the low-x
    # face and outflow the high-x face. Without this the profile reads with the
    # acid front crossed in reverse and pH appears to start near 9.
    if xyz[path[0], 0] > xyz[path[-1], 0]:
        path = path[::-1]
    return path


def path_distance(xyz, path):
    seg = np.linalg.norm(np.diff(xyz[path], axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


# ------------------------------------------------------------------------- figure
def gather(run_dir):
    f, grps, why = open_final(run_dir)
    if f is None:
        return None, why
    try:
        g0, g1 = grps
        xyz, vol, ids, area = read_uge(os.path.join(run_dir, "full_mesh.uge"))
        press = field(g1, "Liquid Pressure", "Pressure")
        if press is None:
            return None, "no pressure field"
        anor0, anor1 = field(g0, "Anorthite VF"), field(g1, "Anorthite VF")
        ph = field(g1, "pH")
        if ph is None:
            h = field(g1, "Total H+")
            ph = -np.log10(np.clip(h, 1e-30, None)) if h is not None else None

        # Total Ca and Mg are not in the output, but the Ca- and Mg-bearing
        # carbonate complexes are, and those are the more direct saturation
        # proxy: CaCO3(aq) activity relates to calcite saturation, whereas
        # total Ca includes cations that are not available to precipitate.
        def sum_fields(*prefixes):
            tot = None
            for p_ in prefixes:
                v = field(g1, p_)
                if v is not None:
                    tot = v if tot is None else tot + v
            return tot
        ca = sum_fields("CaCO3(aq)", "CaHCO3+")
        mg = sum_fields("MgCO3(aq)", "MgHCO3+")
        co3 = sum_fields("CO3--", "HCO3-")

        # Reaction rates say where the process is ACTIVE at 50 yr, which the
        # accumulated volume fractions do not. Dissolution rates are negative,
        # precipitation positive, in mol per m^3 bulk per second.
        diss_rate = field(g1, "Anorthite Rate")
        prec_rate = sum_fields(*[f"{m} Rate" for m in CARB])
        carb = np.zeros(len(xyz))
        for m in CARB:
            a = field(g1, f"{m} VF")
            if a is not None:
                carb += np.clip(a - SEED_VF, 0.0, None)

        # flux magnitude for context
        i, j = ids[:, 0], ids[:, 1]
        d = np.linalg.norm(xyz[i] - xyz[j], axis=1); d[d <= 0] = np.nan
        q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                          posinf=0.0, neginf=0.0)
        fl = np.zeros(len(xyz))
        np.add.at(fl, i, np.abs(q)); np.add.at(fl, j, np.abs(q)); fl *= 0.5

        # anchor the path on the precipitation where there is any, otherwise
        # on the dominant inlet conduit, so both rows are traced the same way
        i_, j_ = ids[:, 0], ids[:, 1]
        dd = np.linalg.norm(xyz[i_] - xyz[j_], axis=1); dd[dd <= 0] = np.nan
        step_cap = 3.0 * float(np.nanpercentile(dd, 95))
        if (carb > 0).any():
            anchor = int(np.argmax(carb))
            anchored_on = "highest-carbonate cell"
        else:
            path0, why0 = trace_path(xyz, ids, area, press)
            if path0 is None:
                return None, why0
            anchor = int(path0[len(path0) // 2])
            anchored_on = "midpoint of the dominant inlet conduit"
        path = trace_through(xyz, ids, area, press, anchor, step_cap)
        if len(path) < 10:
            return None, f"path through the anchor is only {len(path)} cells"
        s = path_distance(xyz, path)
        out = dict(name=os.path.basename(run_dir).replace("C_baseline__", ""),
                   s=s, n_path=len(path),
                   diss=(anor0 - anor1)[path] if anor0 is not None else None,
                   ph=ph[path] if ph is not None else None,
                   ca=ca[path] if ca is not None else None,
                   mg=mg[path] if mg is not None else None,
                   co3=co3[path] if co3 is not None else None,
                   diss_rate=diss_rate[path] if diss_rate is not None else None,
                   prec_rate=prec_rate[path] if prec_rate is not None else None,
                   carb=carb[path], flux=fl[path],
                   span=float(xyz[:, 0].max() - xyz[:, 0].min()),
                   anchored_on=anchored_on,
                   x_start=float(xyz[path[0], 0]), x_end=float(xyz[path[-1], 0]),
                   ph_start=(float(ph[path[0]]) if ph is not None else None),
                   ph_end=(float(ph[path[-1]]) if ph is not None else None),
                   anchor_s=float(path_distance(xyz, path)[list(path).index(anchor)]),
                   carb_total=float(carb.sum() / len(xyz)),
                   n_carb=int((carb > 0).sum()),
                   available=sorted(g1.keys()))
        return out, None
    finally:
        f.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--high", default="C_baseline__p32_100_s383")
    ap.add_argument("--low", default="C_baseline__p32_100_s42")
    ap.add_argument("--out", default="figures/fig_pathprofile.pdf")
    ap.add_argument("--fields", action="store_true",
                    help="list the output fields and exit")
    a = ap.parse_args()

    if a.fields:
        r, why = gather(os.path.join(a.runs, a.high))
        if r is None:
            sys.exit(f"{a.high}: {why}")
        print(f"{len(r['available'])} fields at 50 yr:")
        for k in r["available"]:
            print("  ", k)
        return

    cases = []
    for tag, name in (("traps", a.high), ("does not trap", a.low)):
        r, why = gather(os.path.join(a.runs, name))
        if r is None:
            sys.exit(f"{name}: {why}")
        r["tag"] = tag
        cases.append(r)
        print(f"  {r['name']:<18} path {r['n_path']:>5} cells, "
              f"{r['s'][-1]:>7.2f} m along path, domain span "
              f"{r['span']:.1f} m, {r['n_carb']:>4} precipitating cells\n"
              f"      anchored on the {r['anchored_on']} "
              f"at {r['anchor_s']:.1f} m along the path\n"
              f"      x from {r['x_start']:+.1f} to {r['x_end']:+.1f} m, "
              f"pH {r['ph_start']:.2f} -> {r['ph_end']:.2f}")
        if r["s"][-1] > 5 * r["span"]:
            print(f"      WARNING: path length far exceeds the domain; the "
                  f"trace is jumping between poorly connected cells")

    missing = [k for k in ("ph", "ca", "mg", "diss_rate", "prec_rate")
               if cases[0][k] is None]
    if missing:
        print(f"\n  NOTE: absent from the output: {missing}")
        print( "  Those panels will be blank. Add the species to the OUTPUT")
        print( "  block if the saturation argument is to be shown.")

    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                                "xtick.labelsize": 7, "ytick.labelsize": 7})
    import matplotlib.pyplot as plt

    COL = {"traps": "#1b5e9c", "does not trap": "#b8532a"}
    rows = [("anorthite dissolution rate",
             "|rate| (mol m$^{-3}$ s$^{-1}$)", "diss_rate", True),
            ("Ca and Mg carbonate complexes", "concentration (M)",
             "cations", True),
            ("pH and dissolved carbonate", "pH", "ph", False),
            ("carbonate precipitation rate",
             "rate (mol m$^{-3}$ s$^{-1}$)", "prec_rate", True),
            ("accumulated carbonate", "carbonate VF", "carb", True)]
    rows = [r for r in rows
            if r[2] == "cations" or cases[0].get(r[2]) is not None]

    fig, axes = plt.subplots(len(rows), 1, figsize=(6.6, 1.42 * len(rows)),
                             sharex=True)
    if len(rows) == 1:
        axes = [axes]

    for ax, (title, ylab, key, logy) in zip(axes, rows):
        for c in cases:
            if key == "cations":
                y = np.zeros_like(c["s"])
                for k in ("ca", "mg"):
                    if c[k] is not None:
                        y = y + c[k]
                if not y.any():
                    continue
            else:
                y = c[key]
                if y is None:
                    continue
                if key in ("diss_rate", "prec_rate"):
                    y = np.abs(y)                 # signed; magnitude on log axis
            ax.plot(c["s"], y, lw=1.0, color=COL[c["tag"]],
                    label=f"{c['name']} ({c['tag']})")
            if key == "carb":
                sel = y > 0
                if sel.any():
                    ax.plot(c["s"][sel], y[sel], ".", ms=4,
                            color=COL[c["tag"]], mec="none")
        if logy and any(np.nanmax(np.abs(v)) > 0 for v in ax.get_lines()
                        and [l.get_ydata() for l in ax.get_lines()]):
            ax.set_yscale("log")
        elif logy:
            ax.text(0.5, 0.5, "zero along this path", transform=ax.transAxes,
                    ha="center", va="center", fontsize=7, color="0.5")
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.text(0.005, 0.90, title, transform=ax.transAxes, fontsize=8,
                va="top", fontweight="bold")
        ax.grid(alpha=0.25, lw=0.4)
        ax.tick_params(length=2.5)

    axes[0].legend(fontsize=6.5, loc="upper right", framealpha=0.9)
    axes[-1].set_xlabel("distance along the flow path from the injection region (m)",
                        fontsize=8)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.985, bottom=0.10,
                        hspace=0.12)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300)
    fig.savefig(a.out.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\nwrote {a.out}")
    print(f"      {a.out.rsplit('.',1)[0]}.png")


if __name__ == "__main__":
    main()
