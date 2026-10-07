#!/usr/bin/env python3
"""
Lagrangian particle tracking on a discrete fracture network, with the flow
field taken from either source.

Two modes:

  --field pflotran   UNVERIFIED. Reconstructs the face fluxes from the Liquid
                     Pressure and
                     Permeability that PFLOTRAN writes to its HDF5 output,
                     using Darcy's law on the connectivity list in the .uge
                     file. This is the field the reactive transport actually
                     used.

  --field solver     Uses flowfield.solve(), the separate Darcy solve in this
                     directory. This is the field behind the mean groundwater
                     age and the pre-reaction metric.

  --field both       Tracks the same particles through both fields and reports
                     whether the trajectories agree.

At each cell the next cell is drawn at random with probability proportional to
the outgoing flux, so a junction splits particles in the same ratio as it
splits the water. The time to cross a cell is its pore volume divided by its
total outflow; the distance is the separation of the two cell centroids. Each
particle therefore carries a travel time and a path length, and the fields it
passes through are recorded against both.

    python3 src/ptrack.py --dfn A_p32_100_s1181
    python3 src/ptrack.py --dfn A_p32_100_s1181 --field both
    python3 src/ptrack.py --dfn A_p32_100_s1181 --particles 5000 --csv tracks.csv
"""
from __future__ import annotations
import argparse
import os
import re
import sys

import numpy as np
import h5py

CARB = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
CATION = ["CaCO3(aq) [M]", "CaHCO3+ [M]", "MgCO3(aq) [M]", "MgHCO3+ [M]"]
SEED = 1.0e-6
YEAR = 3.15576e7


# ----------------------------------------------------------------------
# mesh and output
# ----------------------------------------------------------------------
def read_uge(path):
    """Cell centroids and volumes, connection pairs and face areas."""
    with open(path) as f:
        h = f.readline().split()
        if h[0].upper() != "CELLS":
            raise ValueError(f"{path}: expected 'CELLS n', got {h[:2]}")
        n = int(h[1])
        cells = np.empty((n, 4))
        for i in range(n):
            p = f.readline().split()
            cells[i] = (float(p[1]), float(p[2]), float(p[3]), float(p[4]))
        h = f.readline().split()
        if h[0].upper() != "CONNECTIONS":
            raise ValueError(f"{path}: expected 'CONNECTIONS m', got {h[:2]}")
        m = int(h[1])
        ids = np.empty((m, 2), dtype=np.int64)
        area = np.empty(m)
        for j in range(m):
            p = f.readline().split()
            ids[j] = (int(p[0]) - 1, int(p[1]) - 1)
            area[j] = float(p[5])
    return cells, ids, area


def time_groups(f):
    out = []
    for k in f.keys():
        m = re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)
        if m:
            out.append((float(m.group(1)), k))
    out.sort()
    return out


def read_h5(run_dir, which="last"):
    """Pressure, permeability, cation content, pH and carbonate per cell."""
    path = os.path.join(run_dir, "pflotran_co2.h5")
    with h5py.File(path, "r") as f:
        tg = time_groups(f)
        if not tg:
            raise KeyError("no time groups")
        key = tg[0][1] if which == "first" else tg[-1][1]
        g = f[key]
        out = {"time_y": tg[0][0] if which == "first" else tg[-1][0]}
        for name, prefix in (("P", "Liquid Pressure"),
                             ("k", "Permeability"),
                             ("pH", "pH")):
            kk = next((x for x in g if x.startswith(prefix)), None)
            out[name] = np.asarray(g[kk][:], float).flatten() if kk else None
        cat = None
        for s in CATION:
            kk = next((x for x in g if x == s or
                       x.startswith(s.split(" [")[0])), None)
            if kk:
                v = np.asarray(g[kk][:], float).flatten()
                cat = v if cat is None else cat + v
        out["cation"] = cat
        carb = None
        for p in CARB:
            kk = next((x for x in g if x.startswith(p + " VF")), None)
            if kk:
                v = np.maximum(np.asarray(g[kk][:], float).flatten() - SEED, 0.0)
                carb = v if carb is None else carb + v
        out["carb"] = carb
    return out


# ----------------------------------------------------------------------
# flow fields
# ----------------------------------------------------------------------
def flux_from_pressure(P, kperm, xyz, ids, area, mu):
    """Darcy flux per connection from PFLOTRAN's pressure and permeability.

    Q_ij = -(k_ij A_ij / mu) (P_j - P_i) / d_ij,  positive means i -> j.
    """
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[j] - xyz[i], axis=1)
    d[d <= 0] = np.nan
    kh = 2.0 * kperm[i] * kperm[j] / np.maximum(kperm[i] + kperm[j], 1e-300)
    return -(kh * area / mu) * (P[j] - P[i]) / d, i, j


def flux_from_solver(run_dir):
    sys.path.insert(0, "src")
    import flowfield
    sol = flowfield.solve(run_dir, verbose=False)
    return (np.asarray(sol["flux"], float),
            np.asarray(sol["up"], int),
            np.asarray(sol["dn"], int),
            sol)


# ----------------------------------------------------------------------
# tracking
# ----------------------------------------------------------------------
def build_outflow(i, j, q, ncell):
    """Per cell, the cells it sends water to and the flux to each."""
    dest = [[] for _ in range(ncell)]
    for e in range(len(q)):
        f = q[e]
        if not np.isfinite(f) or f == 0.0:
            continue
        a, b = int(i[e]), int(j[e])
        if f > 0:
            dest[a].append((b, f))
        else:
            dest[b].append((a, -f))
    return dest


def inlet_cells(xyz, frac=0.02):
    """Cells within frac of the domain length of the minimum x face."""
    x = xyz[:, 0]
    return np.where(x <= x.min() + frac * (x.max() - x.min()))[0]


def track_one(start, dest, xyz, vol, poro, rng, max_steps=50000):
    cell = int(start)
    t = s = 0.0
    cells, times, dists = [cell], [0.0], [0.0]
    seen = {cell}
    for _ in range(max_steps):
        out = dest[cell]
        if not out:
            break
        js = np.fromiter((o[0] for o in out), int, len(out))
        fs = np.fromiter((o[1] for o in out), float, len(out))
        tot = fs.sum()
        if tot <= 0:
            break
        nxt = int(rng.choice(js, p=fs / tot))
        t += poro * vol[cell] / tot
        s += float(np.linalg.norm(xyz[nxt] - xyz[cell]))
        cell = nxt
        cells.append(cell)
        times.append(t)
        dists.append(s)
        if cell in seen:
            break
        seen.add(cell)
    return np.array(cells), np.array(times), np.array(dists)


def run_field(label, q, i, j, xyz, vol, poro, starts, rng, fld, show):
    ncell = len(vol)
    dest = build_outflow(i, j, q, ncell)

    tracks, reached = [], []
    for s in starts:
        c, t, x = track_one(s, dest, xyz, vol, poro, rng)
        if len(c) < 3:
            continue
        hit = fld["carb"][c] > 0
        tracks.append((c, t, x))
        if hit.any():
            reached.append((c, t, x, int(np.argmax(hit))))

    print(f"--- {label} ---")
    print(f"  tracks longer than two cells   : {len(tracks)} of {len(starts)}")
    print(f"  tracks reaching a carbonate cell: {len(reached)}")
    if tracks:
        L = np.array([len(c) for c, _, _ in tracks])
        T = np.array([t[-1] for _, t, _ in tracks]) / YEAR
        X = np.array([x[-1] for _, _, x in tracks])
        print(f"  median track: {int(np.median(L))} cells, "
              f"{np.median(T):.3f} y, {np.median(X):.2f} m")

    if reached:
        frac = np.array([o / max(len(c) - 1, 1) for c, _, _, o in reached])
        print(f"  carbonate first appears at {np.median(frac):.2f} of the "
              f"track")
        cat = fld["cation"]
        st = np.array([np.nanmedian(cat[c[:max(1, len(c) // 4)]])
                       for c, _, _, _ in reached])
        on = np.array([cat[c[o]] for c, _, _, o in reached])
        print(f"  cation content: {np.median(st):.3e} over the first quarter, "
              f"{np.median(on):.3e} at the onset cell, "
              f"ratio {np.median(on) / max(np.median(st), 1e-30):.1f}")

        from scipy.stats import spearmanr
        rs = []
        for c, t, x, _ in reached:
            v = cat[c]
            ok = np.isfinite(v) & np.isfinite(x)
            if ok.sum() > 5:
                r, _ = spearmanr(x[ok], v[ok])
                if np.isfinite(r):
                    rs.append(r)
        if rs:
            rs = np.array(rs)
            print(f"  cation against distance along the track: median rho "
                  f"{np.median(rs):+.3f}, rising on {(rs > 0).mean() * 100:.0f}% "
                  f"of tracks")

    for c, t, x, o in reached[:show]:
        idx = np.linspace(0, len(c) - 1, 8).astype(int)
        print(f"    track of {len(c)} cells, carbonate first at step {o}")
        print("      time (y)  " + " ".join(f"{t[m]/YEAR:8.3f}" for m in idx))
        print("      dist (m)  " + " ".join(f"{x[m]:8.2f}" for m in idx))
        print("      cation    " + " ".join(f"{fld['cation'][c[m]]:8.1e}"
                                            for m in idx))
        if fld["pH"] is not None:
            print("      pH        " + " ".join(f"{fld['pH'][c[m]]:8.2f}"
                                                for m in idx))
    print()
    return tracks, reached


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn", required=True)
    ap.add_argument("--field", choices=("pflotran", "solver", "both"),
                    default="solver",
                    help="solver is verified: its inlet discharge matches the "
                         "prescribed injection rate")
    ap.add_argument("--particles", type=int, default=500)
    ap.add_argument("--porosity", type=float, default=0.50,
                    help="pore fraction, if the .uge volume is bulk")
    ap.add_argument("--mu", type=float, default=5.47e-4,
                    help="dynamic viscosity, Pa s, water at 50 C")
    ap.add_argument("--pressure-at", choices=("first", "last"), default="first",
                    help="timestep for the pressure and permeability fields")
    ap.add_argument("--show", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--csv")
    a = ap.parse_args()

    d = os.path.join("runs", a.dfn)
    fld = read_h5(d, "last")                       # chemistry at 50 years
    pf = read_h5(d, a.pressure_at)                 # flow field
    cells, ids, area = read_uge(os.path.join(d, "full_mesh.uge"))
    xyz, vol = cells[:, :3], cells[:, 3]
    ncell = len(cells)

    print(f"case {a.dfn}: {ncell} cells, {len(ids)} connections")
    print(f"chemistry at {fld['time_y']:.1f} y, flow field at "
          f"{pf['time_y']:.2f} y")
    if fld["cation"] is None or fld["carb"] is None:
        print("cation or carbonate field absent")
        return 1
    print(f"carbonate-bearing cells: {(fld['carb'] > 0).sum()}")
    print()

    rng = np.random.default_rng(a.seed)

    # flowfield.solve() carries the inlet set used by the verified Darcy
    # solve, whose inlet discharge matches the prescribed injection rate.
    # Prefer it over a guess from cell position.
    inl = None
    try:
        sys.path.insert(0, "src")
        import flowfield
        _sol = flowfield.solve(d, verbose=False)
        cand = np.asarray(_sol.get("inlet", []), int).flatten()
        if cand.size:
            inl = cand
            print(f"inlet cells from flowfield.solve(): {len(inl)}")
    except Exception as e:
        print(f"could not read the inlet set from flowfield.solve(): "
              f"{str(e)[:50]}")
    if inl is None:
        inl = inlet_cells(xyz)
        print(f"inlet cells identified by position: {len(inl)}")

    results = {}

    if a.field in ("pflotran", "both"):
        if pf["P"] is None or pf["k"] is None:
            print("Liquid Pressure or Permeability absent; cannot use this "
                  "field")
        else:
            q, i, j = flux_from_pressure(pf["P"], pf["k"], xyz, ids, area, a.mu)
            dest = build_outflow(i, j, q, ncell)
            send = np.array([sum(f for _, f in dest[c]) for c in inl], float)
            if send.sum() <= 0:
                print("no outgoing flux at the inlet in the PFLOTRAN field")
            else:
                starts = rng.choice(inl, size=a.particles,
                                    p=send / send.sum())
                results["pflotran"] = run_field(
                    "flow field from the PFLOTRAN pressure and permeability",
                    q, i, j, xyz, vol, a.porosity, starts, rng, fld, a.show)

    if a.field in ("solver", "both"):
        q2, up, dn, sol = flux_from_solver(d)
        if len(q2) != len(ids):
            print(f"solver returns {len(q2)} connections, mesh has {len(ids)}")
        else:
            dest = build_outflow(up, dn, q2, ncell)
            send = np.array([sum(f for _, f in dest[c]) for c in inl], float)
            if send.sum() <= 0:
                print("no outgoing flux at the inlet in the solver field")
            else:
                starts = rng.choice(inl, size=a.particles,
                                    p=send / send.sum())
                results["solver"] = run_field(
                    "flow field from flowfield.solve()",
                    q2, up, dn, xyz, vol, a.porosity, starts, rng, fld,
                    a.show)

    if a.field == "both" and len(results) == 2:
        print("=== comparison ===")
        for k, (tr, re_) in results.items():
            L = np.array([len(c) for c, _, _ in tr]) if tr else np.array([0])
            print(f"  {k:10s} {len(tr):5d} tracks, "
                  f"{len(re_):4d} reach carbonate, "
                  f"median {int(np.median(L))} cells")
        a_, b_ = results["pflotran"], results["solver"]
        if a_[0] and b_[0]:
            fa = len(a_[1]) / max(len(a_[0]), 1)
            fb = len(b_[1]) / max(len(b_[0]), 1)
            print(f"  fraction reaching carbonate: {fa:.3f} against {fb:.3f}")
            if abs(fa - fb) < 0.05:
                print("  The two fields give the same result.")
            else:
                print("  The two fields give different results. The flow field "
                      "matters,")
                print("  and the PFLOTRAN field is the one the chemistry saw.")

    if a.csv and results:
        key = "pflotran" if "pflotran" in results else "solver"
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["track", "step", "cell", "time_y", "dist_m",
                        "cation", "pH", "carbonate"])
            for n, (c, t, x, _) in enumerate(results[key][1]):
                for m in range(len(c)):
                    w.writerow([n, m, int(c[m]), t[m] / YEAR, x[m],
                                fld["cation"][c[m]],
                                "" if fld["pH"] is None else fld["pH"][c[m]],
                                fld["carb"][c[m]]])
        print(f"wrote {a.csv} from the {key} field")
    return 0


if __name__ == "__main__":
    sys.exit(main())
