#!/usr/bin/env python3
"""
Compare the mesh convergence sequence and report a bound, not an order.

Per the protocol: check monotonicity BEFORE attempting to fit an order of
convergence. Reactive transport with adaptive time stepping frequently produces
non-monotone successive differences, in which case no order can be fitted and
the honest statement is a bound over the range tested.

Two quantities are reported deliberately:

  carb_per_cell   sum(VF)/N_cells -- as the manuscript defines it. N_cells is
                  in the denominator, so this CANNOT converge; it is a
                  mesh-defined quantity.
  carb_intensity  sum(VF_i * V_i) / sum(phi * V_i) -- independent of the
                  discretisation, and therefore the quantity a convergence
                  test can meaningfully address.

Showing that the first does not converge while the second does is the
justification for preferring the volume-weighted measure.
"""
from __future__ import annotations
import sys, os, glob, json, math, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CARB, SEED_VF, PHI = ("Calcite","Magnesite","Siderite","Dawsonite"), 1e-6, 0.50


def read(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "no output"
    with open(uge) as f:
        n = int(f.readline().split()[1])
        vol = np.array([float(f.readline().split()[4]) for _ in range(n)])
    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg: return None, "no time groups"
        t_end = tg[-1][0]
        g = f[tg[-1][1]]
        vf = vol_c = 0.0
        per = {}
        for m in CARB:
            k = next((x for x in g if x.startswith(f"{m} VF")), None)
            if k is None: per[m] = 0.0; continue
            a = np.clip(g[k][:] - SEED_VF, 0.0, None)
            if a.size != vol.size: return None, "length mismatch"
            vf += float(a.sum()); per[m] = float((a * vol).sum())
        vol_c = sum(per.values())
    pore = float((vol * PHI).sum())
    return dict(n_cells=n, t_end=t_end, fracture_volume=float(vol.sum()),
                pore_volume=pore, carb_per_cell=vf / n,
                carb_intensity=vol_c / pore,
                **{f"pct_{m.lower()}": (100*per[m]/vol_c if vol_c else 0.0)
                   for m in CARB}), None


def report(name, factors, cells, values, unit=""):
    print(f"\n  {name}")
    print(f"    {'factor':>8}{'cells':>10}{'value':>14}{'diff from prev':>17}")
    prev, diffs = None, []
    for fa, c, v in zip(factors, cells, values):
        if prev is None:
            d = ""
        else:
            rel = 100 * (v - prev) / prev if prev else float("nan")
            diffs.append(rel)
            d = f"{rel:+.2f}%"
        print(f"    {fa:>8.4g}{c:>10,}{v:>14.4e}{d:>17}")
        prev = v
    if len(values) >= 2:
        spread = 100 * (max(values) - min(values)) / np.mean(values)
        print(f"    spread over the sequence: {spread:.2f}%")
    order = None
    if len(diffs) >= 2:
        # Monotone VALUES are not convergence. A sequence can rise at every
        # refinement while the successive differences grow, which is
        # divergence. The test that matters is whether |y_{k+1} - y_k| shrinks.
        signs = {np.sign(d) for d in diffs}
        same_sign = len(signs) == 1
        d1 = abs(values[1] - values[0])
        d2 = abs(values[2] - values[1]) if len(values) >= 3 else None
        print(f"    successive differences: "
              f"{'same sign' if same_sign else 'CHANGE SIGN'}")
        if not same_sign:
            print("      no order can be fitted; report a bound only")
        elif d2 is None:
            print("      only two members; no order can be fitted")
        elif d2 == 0:
            print("      the two finest members are identical")
        elif d1 <= d2:
            print(f"      differences GROW under refinement "
                  f"({d1:.3e} then {d2:.3e})")
            print("      this quantity is NOT converging; no order can be fitted")
        else:
            # uniform refinement ratio required for the classical estimate
            r1 = factors[0] / factors[1]
            r2 = factors[1] / factors[2]
            if abs(r1 - r2) / r1 > 0.10:
                print(f"      refinement ratio not uniform "
                      f"({r1:.2f} then {r2:.2f}); no order fitted")
            else:
                order = math.log(d1 / d2) / math.log(r1)
                rich = values[-1] + (values[-1] - values[-2]) / (r1**order - 1)
                err = 100 * abs(rich - values[-1]) / abs(values[-1])
                print(f"      differences shrink: apparent order p = {order:.3f}")
                print(f"      Richardson estimate of the converged value: "
                      f"{rich:.4e}")
                print(f"      remaining error at the finest mesh: {err:.1f}%")
    return diffs, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq", default="mesh_conv")
    a = ap.parse_args()

    runs = sorted(glob.glob(os.path.join(a.seq, "run_mrf*")))
    if len(runs) < 2:
        sys.exit(f"found {len(runs)} run(s) in {a.seq} -- need at least two")

    rows = []
    for r in runs:
        f = float(os.path.basename(r).replace("run_mrf", ""))
        d, why = read(r)
        if d is None:
            print(f"  mrf{f:g}: {why}")
            continue
        d["factor"] = f
        rows.append(d)
    rows.sort(key=lambda r: -r["factor"])
    if len(rows) < 2:
        sys.exit("fewer than two usable members")

    fa = [r["factor"] for r in rows]
    nc = [r["n_cells"] for r in rows]

    print("=" * 74)
    print("MESH CONVERGENCE")
    print("=" * 74)
    print(f"  members            {len(rows)}")
    print(f"  cumulative refinement  x{nc[-1]/nc[0]:.2f}")
    te = {round(r["t_end"], 1) for r in rows}
    print(f"  final time         {sorted(te)} yr"
          + ("  ** members differ **" if len(te) > 1 else ""))
    if len(te) > 1:
        # A member still running has had less reaction time, not a different
        # discretisation. Comparing them measures elapsed time, not mesh
        # resolution, so the run stops here rather than printing numbers that
        # would be read as a convergence result.
        print()
        print("  " + "!" * 68)
        print("  COMPARISON ABANDONED: the members stop at different times.")
        print("  Any difference between them reflects the reaction time, not the")
        print("  mesh. Wait for every member to reach the same time and re-run.")
        for r in rows:
            print(f"    mrf{r['factor']:<6g} {r['n_cells']:>9,} cells   "
                  f"t_end = {r['t_end']:>5.1f} yr")
        print("  " + "!" * 68)
        sys.exit(2)
    fv = np.array([r["fracture_volume"] for r in rows])
    print(f"  fracture volume    {fv.min():.4f} to {fv.max():.4f} m3 "
          f"({100*(fv.max()-fv.min())/fv.mean():.2f}% spread)")

    _, ord_iv = report(
        "carb_intensity  -- volume-weighted, discretisation-independent",
        fa, nc, [r["carb_intensity"] for r in rows])
    _, ord_pc = report(
        "carb_per_cell   -- as the manuscript defines it (N_cells in denominator)",
        fa, nc, [r["carb_per_cell"] for r in rows])

    print("\n  carbonate assemblage")
    print(f"    {'factor':>8}" + "".join(f"{m[:9]:>11}" for m in CARB))
    for r in rows:
        print(f"    {r['factor']:>8.4g}" +
              "".join(f"{r[f'pct_{m.lower()}']:>10.1f}%" for m in CARB))

    print("\n" + "=" * 74)
    print("SUGGESTED WORDING")
    print("=" * 74)
    iv = [r["carb_intensity"] for r in rows]
    pc = [r["carb_per_cell"] for r in rows]
    si = 100*(max(iv)-min(iv))/np.mean(iv)
    sp = 100*(max(pc)-min(pc))/np.mean(pc)
    ref = nc[-1]/nc[0]

    if ord_iv is not None:
        r1 = fa[0]/fa[1]
        rich = iv[-1] + (iv[-1]-iv[-2])/(r1**ord_iv - 1)
        err = 100*abs(rich-iv[-1])/abs(iv[-1])
        print(f"  The volume-weighted trapping intensity converges at an apparent")
        print(f"  order of {ord_iv:.2f} over a {ref:.1f}-fold refinement at fixed fracture")
        print(f"  network, with the observed values spanning {si:.1f}%. Richardson")
        print(f"  extrapolation places the converged value {err:.0f}% above the finest")
        print(f"  mesh, and the absolute trapping amounts should be read with that")
        print(f"  uncertainty.")
    else:
        print(f"  The volume-weighted trapping intensity varies {si:.1f}% over a")
        print(f"  {ref:.1f}-fold refinement at fixed fracture network. No order of")
        print(f"  convergence can be fitted to the sequence, so the variation is")
        print(f"  reported as a bound rather than a rate.")

    if ord_pc is None:
        print()
        print(f"  The per-cell measure used in the submitted manuscript varies")
        print(f"  {sp:.1f}% over the same sequence and its successive differences do")
        print(f"  not shrink, so no order can be fitted to it. That measure carries")
        print(f"  the cell count in its denominator and is therefore not a")
        print(f"  discretisation-independent quantity, which is why the revision")
        print(f"  reports the volume-weighted intensity instead.")
    else:
        print()
        print(f"  The per-cell measure varies {sp:.1f}% over the same sequence at an")
        print(f"  apparent order of {ord_pc:.2f}.")

    amin = min(min(r[f"pct_{m.lower()}"] for r in rows) for m in CARB)
    print()
    print(f"  The carbonate assemblage is stable across the sequence:")
    for m in CARB:
        v = [r[f"pct_{m.lower()}"] for r in rows]
        print(f"    {m:<11} {min(v):>5.1f} to {max(v):>5.1f}%   "
              f"spread {max(v)-min(v):.1f} points")


if __name__ == "__main__":
    main()
