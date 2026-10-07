#!/usr/bin/env python3
"""
Traces representative flow routes from the PFLOTRAN pressure and permeability
fields, after checking the reconstructed flux against a known quantity.

The flux on the connection between cells i and j follows Darcy's law,

    Q_ij = -(k_ij * A_ij / mu) * (P_j - P_i) / d_ij,

with k_ij the harmonic mean of the two cell permeabilities, A_ij the face area
from the .uge file, d_ij the distance between cell centroids, and mu the
dynamic viscosity. The route follows the largest flux at each step, not the
steepest pressure gradient: in a fracture network the apertures vary by orders
of magnitude, so the permeability and the face area decide where the water goes
as much as the pressure drop does.

An earlier attempt at this reconstruction was wrong by four orders of
magnitude. This version therefore validates itself before tracing anything: it
sums the reconstructed flux leaving the injection boundary and compares that
with the discharge reported by flowfield.solve(), whose mass balance closes to
one part in 10^13 and whose inlet discharge matches the prescribed injection
rate. If the two disagree by more than the tolerance, the script reports the
ratio and stops rather than producing routes from a wrong field.

What the routes are and are not: each is the sequence of cells obtained by
following the largest outgoing flux, so it traces the dominant route through
the network. Water divides at every junction, so a route is representative
rather than the path of any particular parcel. Velocity is not reconstructed
inside a cell, so a route is a sequence of cells rather than a streamline. For
trajectories in the proper sense, use dfnTrans, which reconstructs the local
velocity field first \\(Makedonska et al., 2015\\).

    python3 src/flowpath.py --dfn A_p32_100_s1181
    python3 src/flowpath.py --dfn A_p32_100_s1181 --time 10.0 --routes 20
    python3 src/flowpath.py --dfn A_p32_100_s1181 --to-carbonate --csv routes.csv
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
# input
# ----------------------------------------------------------------------
def read_uge(path):
    """Cell centroids and volumes, connection pairs, face areas."""
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


def snapshots(path):
    out = []
    with h5py.File(path, "r") as f:
        for k in f.keys():
            m = re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)
            if m:
                out.append((float(m.group(1)), k))
    out.sort()
    return out


def read_snapshot(path, group):
    """Every field this script needs, from one snapshot."""
    want = {"P": "Liquid Pressure", "k": "Permeability", "pH": "pH",
            "rho": "Liquid Density"}
    out = {}
    with h5py.File(path, "r") as f:
        g = f[group]
        for name, prefix in want.items():
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
# flux
# ----------------------------------------------------------------------
def darcy_flux(P, kperm, xyz, ids, area, mu, rho=None, g=9.80665):
    """Volumetric flux per connection, positive from ids[:,0] to ids[:,1].

    PFLOTRAN reports the total liquid pressure, which carries the hydrostatic
    component. Gravity balances that component, so the flux follows the
    gradient of the pressure with the hydrostatic part removed,

        q = -(k A / mu) * [(P_j - P_i) - rho g (z_i - z_j)] / d.

    The hydrostatic gradient here is rho*g, near 9700 Pa/m, whereas the
    gradient that drives the flow is of order 1e-3 Pa across a cell. Leaving
    the term out therefore overestimates the flux by several orders of
    magnitude, which is what an earlier version of this script did.
    """
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[j] - xyz[i], axis=1)
    d[d <= 0] = np.nan
    kh = 2.0 * kperm[i] * kperm[j] / np.maximum(kperm[i] + kperm[j], 1e-300)
    dP = P[j] - P[i]
    if rho is not None:
        # z increases upward, so the hydrostatic difference between j and i
        # is -rho g (z_j - z_i)
        dP = dP + rho * g * (xyz[j, 2] - xyz[i, 2])
    q = -(kh * area / mu) * dP / d
    q[~np.isfinite(q)] = 0.0
    return q


def calibrate(q, ids, inlet, q_ref):
    """Scale factor that matches the reconstructed inlet discharge to q_ref.

    Returns the factor and the discharge before scaling. A factor far from one
    means the reconstruction has a unit or formula error, and the caller should
    stop rather than trace routes through a wrong field.
    """
    inl = set(int(c) for c in inlet)
    tot = 0.0
    for e in range(len(q)):
        a, b = int(ids[e, 0]), int(ids[e, 1])
        f = q[e]
        if a in inl and b not in inl and f > 0:
            tot += f
        elif b in inl and a not in inl and f < 0:
            tot += -f
    if tot <= 0:
        return None, tot
    return q_ref / tot, tot


def outflow_map(q, ids, ncell):
    dest = [[] for _ in range(ncell)]
    for e in range(len(q)):
        f = q[e]
        if f == 0.0 or not np.isfinite(f):
            continue
        a, b = int(ids[e, 0]), int(ids[e, 1])
        if f > 0:
            dest[a].append((b, f))
        else:
            dest[b].append((a, -f))
    return dest


def inflow_map(q, ids, ncell):
    src = [[] for _ in range(ncell)]
    for e in range(len(q)):
        f = q[e]
        if f == 0.0 or not np.isfinite(f):
            continue
        a, b = int(ids[e, 0]), int(ids[e, 1])
        if f > 0:
            src[b].append((a, f))
        else:
            src[a].append((b, -f))
    return src


# ----------------------------------------------------------------------
# routes
# ----------------------------------------------------------------------
def trace(start, table, xyz, vol, poro, max_steps=20000):
    """Follow the largest flux. Returns cells, cumulative time and distance."""
    cell = int(start)
    t = s = 0.0
    cells, times, dists = [cell], [0.0], [0.0]
    seen = {cell}
    for _ in range(max_steps):
        opts = [(c, f) for c, f in table[cell] if c not in seen]
        if not opts:
            break
        nxt, _ = max(opts, key=lambda o: o[1])
        tot = sum(f for _, f in table[cell])
        if tot > 0:
            t += poro * vol[cell] / tot
        s += float(np.linalg.norm(xyz[nxt] - xyz[cell]))
        cell = nxt
        cells.append(cell)
        times.append(t)
        dists.append(s)
        seen.add(cell)
    return np.array(cells), np.array(times), np.array(dists)


def report(cells, times, dists, fld, xyz, label):
    c = cells
    idx = np.linspace(0, len(c) - 1, min(8, len(c))).astype(int)
    net = float(np.linalg.norm(xyz[c[-1]] - xyz[c[0]]))
    print(f"  {label}: {len(c)} cells, {dists[-1]:.2f} m travelled, "
          f"{net:.2f} m net, {times[-1]/YEAR:.4f} y")
    print("    position  " + " ".join(f"{dists[m]:8.2f}" for m in idx))
    if fld["cation"] is not None:
        print("    cation    " + " ".join(f"{fld['cation'][c[m]]:8.1e}"
                                          for m in idx))
    if fld["pH"] is not None:
        print("    pH        " + " ".join(f"{fld['pH'][c[m]]:8.2f}"
                                          for m in idx))
    if fld["carb"] is not None:
        print("    carbonate " + " ".join(f"{fld['carb'][c[m]]:8.1e}"
                                          for m in idx))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn", required=True)
    ap.add_argument("--time", type=float,
                    help="snapshot in years; the last one by default")
    ap.add_argument("--routes", type=int, default=10)
    ap.add_argument("--to-carbonate", action="store_true",
                    help="trace upstream from the carbonate-bearing cells "
                         "instead of downstream from the inlet")
    ap.add_argument("--porosity", type=float, default=0.50)
    ap.add_argument("--mu", type=float, default=5.47e-4,
                    help="dynamic viscosity, Pa s, water at 50 C")
    ap.add_argument("--rho", type=float, default=None,
                    help="liquid density for the hydrostatic correction; "
                         "taken from the output if not given")
    ap.add_argument("--no-gravity", action="store_true",
                    help="omit the hydrostatic correction")
    ap.add_argument("--tol", type=float, default=2.0,
                    help="largest acceptable factor between the reconstructed "
                         "and reference inlet discharge")
    ap.add_argument("--force", action="store_true",
                    help="trace even if the flux check fails")
    ap.add_argument("--csv")
    a = ap.parse_args()

    d = os.path.join("runs", a.dfn)
    h5 = os.path.join(d, "pflotran_co2.h5")
    snaps = snapshots(h5)
    if not snaps:
        print("no snapshots in the output")
        return 1
    if a.time is None:
        t_use, grp = snaps[-1]
    else:
        t_use, grp = min(snaps, key=lambda s: abs(s[0] - a.time))
    fld = read_snapshot(h5, grp)

    cells, ids, area = read_uge(os.path.join(d, "full_mesh.uge"))
    xyz, vol = cells[:, :3], cells[:, 3]
    ncell = len(cells)

    print(f"case {a.dfn}")
    print(f"  {ncell} cells, {len(ids)} connections")
    print(f"  snapshot at {t_use:g} y")
    if fld["P"] is None or fld["k"] is None:
        print("  Liquid Pressure or Permeability absent from the output")
        return 1
    if fld["carb"] is not None:
        print(f"  carbonate-bearing cells: {int((fld['carb'] > 0).sum())}")
    print()

    # the reference field, used for the inlet set and the discharge
    sys.path.insert(0, "src")
    import flowfield
    sol = flowfield.solve(d, verbose=False)
    inlet = np.asarray(sol["inlet"], int).flatten()
    q_ref = float(sol["q_inlet_m3s"])
    print(f"reference: {len(inlet)} inlet cells, "
          f"inlet discharge {q_ref:.4e} m^3/s")

    rho = None if a.no_gravity else a.rho
    if rho is None and not a.no_gravity and fld["rho"] is not None:
        rho = float(np.nanmedian(fld["rho"]))
        print(f"  liquid density from the output: {rho:.1f} kg/m^3")
    q = darcy_flux(fld["P"], fld["k"], xyz, ids, area, a.mu, rho=rho)
    factor, raw = calibrate(q, ids, inlet, q_ref)
    if factor is None:
        print("reconstruction gives no outflow at the inlet; stopping")
        return 1
    print(f"reconstruction: inlet discharge {raw:.4e} m^3/s, "
          f"factor {1.0/factor:.4g} times the reference")
    if not (1.0 / a.tol <= 1.0 / factor <= a.tol):
        print()
        print(f"  The reconstructed flux differs from the reference by more "
              f"than {a.tol:g} times.")
        print(f"  The Darcy reconstruction is therefore wrong, most likely in "
              f"the units of")
        print(f"  the face area or in the treatment of the interface "
              f"permeability. Routes")
        print(f"  traced through this field would not be meaningful.")
        if not a.force:
            print(f"  Stopping. Use --force to trace anyway.")
            return 1
        print(f"  Continuing because --force was given.")
    else:
        print(f"  Within tolerance; the reconstruction reproduces the "
              f"reference discharge.")
    print()

    table = (inflow_map(q, ids, ncell) if a.to_carbonate
             else outflow_map(q, ids, ncell))

    if a.to_carbonate:
        if fld["carb"] is None or (fld["carb"] > 0).sum() == 0:
            print("no carbonate-bearing cells to trace back from")
            return 1
        prec = np.where(fld["carb"] > 0)[0]
        starts = prec[np.argsort(fld["carb"][prec])[::-1]][:a.routes]
        print(f"tracing upstream from the {len(starts)} cells holding the "
              f"most carbonate")
    else:
        send = np.array([sum(f for _, f in outflow_map(q, ids, ncell)[c])
                         for c in inlet], float)
        order = np.argsort(send)[::-1]
        starts = inlet[order[:a.routes]]
        print(f"tracing downstream from the {len(starts)} inlet cells with "
              f"the largest discharge")
    print()

    routes = []
    for n, s in enumerate(starts, 1):
        c, t, x = trace(int(s), table, xyz, vol, a.porosity)
        if len(c) < 3:
            continue
        if a.to_carbonate:                     # read from the inlet end
            c, t, x = c[::-1], t[-1] - t[::-1], x[-1] - x[::-1]
        routes.append((c, t, x))
        report(c, t, x, fld, xyz, f"route {n}")
        print()

    if routes and fld["cation"] is not None:
        from scipy.stats import spearmanr
        rs = []
        for c, _, x in routes:
            v = fld["cation"][c]
            ok = np.isfinite(v) & np.isfinite(x)
            if ok.sum() > 5:
                r, _ = spearmanr(x[ok], v[ok])
                if np.isfinite(r):
                    rs.append(r)
        if rs:
            rs = np.array(rs)
            print(f"cation content against distance along the route: "
                  f"median rho {np.median(rs):+.3f}, "
                  f"rising on {(rs > 0).mean()*100:.0f}% of routes")
        tortuosity = [x[-1] / max(np.linalg.norm(xyz[c[-1]] - xyz[c[0]]), 1e-12)
                      for c, _, x in routes]
        print(f"route length over net displacement: "
              f"median {np.median(tortuosity):.2f}")

    if a.csv and routes:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["route", "step", "cell", "x", "y", "z",
                        "dist_m", "time_y", "cation", "pH", "carbonate"])
            for n, (c, t, x) in enumerate(routes):
                for m in range(len(c)):
                    w.writerow([n, m, int(c[m]),
                                xyz[c[m], 0], xyz[c[m], 1], xyz[c[m], 2],
                                x[m], t[m] / YEAR,
                                "" if fld["cation"] is None else fld["cation"][c[m]],
                                "" if fld["pH"] is None else fld["pH"][c[m]],
                                "" if fld["carb"] is None else fld["carb"][c[m]]])
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
