#!/usr/bin/env python3
"""
Lowest porosity each coupled run reached: porosity = 1 - (sum of all mineral
volume fractions), at every saved snapshot. Shows whether the porosity floor
(MINIMUM_POROSITY) was, or would have been, active.

    python3 src/min_porosity.py lists/A_feedback.txt [--floor 0.01]

Per run: snapshots read, last time, lowest porosity and when, cells that fell
below the floor at any snapshot, and whether the deck carries the floor.
A run WITHOUT the floor whose porosity never fell below it gives the same
result as a floored run. Only saved snapshots are seen: a dip between two
snapshots is missed, so treat "never below" as "not below at any snapshot".
"""
import glob, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE


def tval(k):
    return float(k.split("Time")[1].split()[0])


def scan(d, floor):
    h5 = os.path.join(d, "pflotran_co2.h5")
    if not os.path.isfile(h5):
        return None
    lo, t_lo, below, times = np.inf, None, None, []
    with h5py.File(h5, "r") as f:
        for k in sorted((k for k in f if "Time" in k), key=tval):
            g = f[k]; vf = [n for n in g if re.search(r"\bVF\b", n)]
            if not vf:
                continue
            phi = 1.0 - sum(np.asarray(g[n][:], float).ravel() for n in vf)
            times.append(tval(k))
            below = (phi < floor) if below is None else (below | (phi < floor))
            if phi.min() < lo:
                lo, t_lo = float(phi.min()), tval(k)
    if not times:
        return None
    return dict(n=len(times), last=times[-1], lo=lo, t_lo=t_lo, below=int(below.sum()))


def main():
    a = sys.argv[1:]
    floor = float(a[a.index("--floor") + 1]) if "--floor" in a else 0.01
    lists = [x for x in a if x.endswith(".txt")]
    if not lists:
        sys.exit(__doc__)
    rows = []
    for lst in lists:
        for rid in (l.strip() for l in open(os.path.join(REV, lst) if not os.path.isabs(lst) else lst) if l.strip()):
            d = os.path.join(REV, rid) if "/" in rid else os.path.join(REV, "runs", rid)
            deck = next(iter(sorted(glob.glob(os.path.join(d, "*.in")))), None)
            fl = bool(deck) and re.search(r"^\s*MINIMUM_POROSITY", open(deck).read(), re.M | re.I) is not None
            rows.append((os.path.basename(d), fl, scan(d, floor)))
    W = max(22, max(len(r[0]) for r in rows) + 2)
    print(f"{'run':<{W}}{'floor':>6}{'snaps':>7}{'last y':>8}{'min phi':>10}{'at y':>7}{'cells<' + str(floor):>12}")
    for name, fl, s in rows:
        if s is None:
            print(f"{name:<{W}}{'yes' if fl else 'no':>6}   no output yet"); continue
        print(f"{name:<{W}}{'yes' if fl else 'no':>6}{s['n']:>7}{s['last']:>8.3g}{s['lo']:>10.4f}{s['t_lo']:>7.3g}{s['below']:>12}")
    done = [(n, fl, s) for n, fl, s in rows if s]
    nf = [(n, s) for n, fl, s in done if not fl]
    hit = [n for n, s in nf if s["lo"] < floor]
    print(f"\n{len(done)} runs with output; {len(nf)} without the floor, of which {len(hit)} went below {floor} "
          f"(the floor would have acted there){': ' + ', '.join(hit) if hit else ''}")
    fl_hit = [n for n, fl, s in done if fl and s["lo"] <= floor + 1e-9]
    if fl_hit:
        print(f"floored runs that reached the floor (it was active): {len(fl_hit)}")


if __name__ == "__main__":
    main()
