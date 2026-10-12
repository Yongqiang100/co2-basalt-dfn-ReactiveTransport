"""pH against fracture intensity, and acid-front arrival (Section 3.1).

    python3 src/ph_statistics.py [--root runs_gravityoff]

For each intensity-ensemble network (A_feedback__*): the volume-weighted mean pH at 0.1, 1, 10
and 50 years, tested against P32 with Spearman and Kruskal-Wallis; and the first output time at
which 90% of the fracture volume has pH below 5 (acid-front arrival), tested against P32.
"""
import argparse, glob, os, re, sys
import numpy as np, h5py
from scipy.stats import spearmanr, kruskal

sys.path.insert(0, os.path.dirname(__file__))
import fig_cations as fc
from domain_mean_ph import cell_volumes

ap = argparse.ArgumentParser(); ap.add_argument("--root", default="runs_gravityoff"); a = ap.parse_args()

def ph_series(d):
    out = {}
    for t, h5, g in fc._snaps(d):
        with h5py.File(h5, "r") as f:
            k = next((k for k in f[g].keys() if re.match(r"\s*pH\b", k)), None)
            if k is not None: out[round(t, 3)] = np.asarray(f[g][k][:], float)
    return out

def volumes(d, n):
    for args in ([d, n], [d]):
        try:
            v = np.asarray(cell_volumes(*args), float)
            if v.size == n: return v
        except Exception:
            pass
    raise SystemExit(f"cell volumes not found for {d}")

rows = []
for d in sorted(glob.glob(os.path.join(a.root, "A_feedback__p32_*"))):
    lvl = float(re.search(r"p32_(\d+)_", d).group(1)) / 100
    P = ph_series(d)
    if not P: continue
    v = volumes(d, len(next(iter(P.values())))); v = v / v.sum()
    mean = {t: float((p * v).sum()) for t, p in P.items()}
    front = next((t for t in sorted(P) if v[P[t] < 5].sum() >= 0.9), np.nan)
    rows.append((lvl, mean, front))
print(f"networks: {len(rows)}")
for t in (0.1, 1, 10, 50):
    x = [r[0] for r in rows if t in r[1]]; y = [r[1][t] for r in rows if t in r[1]]
    g = {}; [g.setdefault(p, []).append(q) for p, q in zip(x, y)]
    rho, p = spearmanr(x, y); _, pk = kruskal(*g.values())
    print(f"t={t:>4} y: mean pH {min(y):.2f}-{max(y):.2f} | level medians " + ", ".join(f"x{k:.2f} {np.median(v):.2f}" for k, v in sorted(g.items()))
          + f" | Spearman rho={rho:+.2f} p={p:.2g} | Kruskal-Wallis p={pk:.2g} | n={len(y)}")
x = [r[0] for r in rows if not np.isnan(r[2])]; y = [r[2] for r in rows if not np.isnan(r[2])]
g = {}; [g.setdefault(p, []).append(q) for p, q in zip(x, y)]
rho, p = spearmanr(x, y)
print(f"acid-front arrival (pH<5 in 90% of the fracture volume): {min(y)}-{max(y)} y | level medians "
      + ", ".join(f"x{k:.2f} {np.median(v)}" for k, v in sorted(g.items())) + f" | Spearman rho={rho:+.2f} p={p:.2g} | n={len(y)} of {len(rows)}")
