#!/usr/bin/env python3
"""
Health check for a rerun block, run once its jobs have finished.

For every network: finished to 50 years, gravity off in its deck, remaining
circulation (gross/net flow across x = -4 m, 1 = none), runaway cells (share
of the carbonate in cells with volume fraction above 1), carbonate per cell,
and the change from the old run of the same network.

    cd $MYSCRATCH/co2-basalt/revision
    python3 src/check_block.py --prefix A_
    python3 src/check_block.py --prefix C_ --old runs_dirichlet_20260924/runs_C050
"""
import argparse, glob, os, sys
import numpy as np, h5py
from scipy.stats import spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import compare_hydrotest as c


def last_time(h5):
    with h5py.File(h5) as h:
        t = [float(k.split("Time")[1].split()[0]) for k in h if "Time" in k]
    return max(t) if t else None


def metrics(d):
    r = c.load(d)
    x, q, i, j = r["xyz"][:, 0], r["q"], r["i"], r["j"]
    lo, hi = np.minimum(x[i], x[j]), np.maximum(x[i], x[j]); cr = (lo < -4) & (hi >= -4)
    net = (q[cr] * np.sign(x[j][cr] - x[i][cr])).sum()
    t = r["total"]; pos = t[t > 0]; tot = pos.sum()
    return dict(gv=r["gv"], gn=np.abs(q[cr]).sum() / abs(net) if net else float("inf"),
                maxvf=float(pos.max()) if pos.size else 0.0,
                runaway=100 * pos[pos > 1].sum() / tot if tot else 0.0,
                carb=float(np.clip(t, 0, None).mean()), cells=int(pos.size))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="A_")
    ap.add_argument("--new", default="runs_gravityoff")
    ap.add_argument("--old", default="runs_dirichlet_20260924/runs",
                    help="superseded runs; the flat runs_dirichlet_20260924/ layout is also searched")
    a = ap.parse_args()
    new_root = a.new if os.path.isabs(a.new) else os.path.join(REV, a.new)
    old_root = a.old if os.path.isabs(a.old) else os.path.join(REV, a.old)
    dirs = sorted(d for d in glob.glob(os.path.join(new_root, a.prefix + "p32_*")) if os.path.isdir(d))
    if not dirs:
        sys.exit(f"no {a.prefix}p32_* directories in {new_root}")

    W = max(22, max(len(os.path.basename(d)) for d in dirs) + 2)
    print(f"{'network':<{W}}{'status':<14}{'gross/net':>10}{'max VF':>9}{'runaway%':>10}{'carb/cell':>12}{'new/old':>9}")
    rows, notdone, gravity_on = [], [], []
    for d in dirs:
        name = os.path.basename(d); h5 = os.path.join(d, "pflotran_co2.h5")
        t = last_time(h5) if os.path.isfile(h5) else None
        if t is None or t < 49.99:
            notdone.append(name)
            print(f"{name:<{W}}{('no output' if t is None else f'at {t:.3g} y'):<14}"); continue
        m = metrics(d)
        if np.any(m["gv"] != 0):
            gravity_on.append(name)
        old = next((o for o in (os.path.join(old_root, name),
                                os.path.join(REV, "runs_dirichlet_20260924", name))
                    if os.path.isfile(os.path.join(o, "pflotran_co2.h5"))), os.path.join(old_root, name))
        ratio = None
        if os.path.isfile(os.path.join(old, "pflotran_co2.h5")):
            oc = float(np.clip(c.load(old)["total"], 0, None).mean())
            ratio = m["carb"] / oc if oc > 0 else None
            m["old"] = oc
        rows.append((name, m, ratio))
        print(f"{name:<{W}}{'done':<14}{m['gn']:>10.2f}{m['maxvf']:>9.3g}{m['runaway']:>10.1f}{m['carb']:>12.3e}"
              f"{(f'{ratio:9.2f}' if ratio else '        -')}")

    print(f"\n{len(rows)} of {len(dirs)} finished" + (f"; {len(notdone)} not finished (running, queued or stopped)" if notdone else ""))
    if not rows:
        return
    gn = np.array([m["gn"] for _, m, _ in rows]); run = np.array([m["runaway"] for _, m, _ in rows])
    print(f"gravity off in every finished deck: {'yes' if not gravity_on else 'NO: ' + ', '.join(gravity_on)}")
    print(f"circulation (gross/net at x = -4 m): median {np.median(gn):.2f}, max {gn.max():.2f}   (1 = none)")
    print(f"runaway cells (VF > 1): in {int((run > 0).sum())} networks; above 5% of the carbonate in "
          f"{int((run > 5).sum())}; median share {np.median(run):.1f}%")
    carb = np.array([m["carb"] for _, m, _ in rows])
    print(f"carbonate per cell: mean {carb.mean():.3e}, median {np.median(carb):.3e}")
    paired = [(m["old"], m["carb"]) for _, m, _ in rows if "old" in m]
    if paired:
        o, n = np.array(paired).T
        ok = o > 0
        print(f"change from the old runs ({len(paired)} pairs): median ratio {np.median(n[ok] / o[ok]):.2f}, "
              f"Spearman(old, new) = {spearmanr(o, n).correlation:.3f}")


if __name__ == "__main__":
    main()
