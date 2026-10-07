"""Shut-in gain against fracture intensity (Section 3.5).

    python3 src/shutin_by_intensity.py

Ratio of the carbonate at 50 years with a shut-in at 10 years (E_feedback) to that under
continuous injection (A_feedback), for the 50 paired networks, tested against P32.
Reads results/A_coupled_networks.csv and results/E_coupled_networks.csv (summarize_block.py).
"""
import csv, collections, os, sys
from scipy.stats import spearmanr, kruskal

def read(path):
    if not os.path.exists(path): sys.exit(f"{path} missing: run summarize_block.py first")
    return {r["network"].split("__")[-1]: r for r in csv.DictReader(open(path))}

A = read("results/A_coupled_networks.csv"); E = read("results/E_coupled_networks.csv")
x, y, g = [], [], collections.defaultdict(list)
for k, r in A.items():
    if k in E:
        lv = float(r["p32_factor"]); q = float(E[k]["co2_kg"]) / float(r["co2_kg"])
        x.append(lv); y.append(q); g[lv].append(q)
rho, p = spearmanr(x, y); _, pk = kruskal(*g.values())
print(f"pairs {len(y)} | shut-in/continuous ratio by level (median): "
      + ", ".join(f"x{k:.2f} {sorted(v)[len(v) // 2]:.2f}" for k, v in sorted(g.items())))
print(f"Spearman rho = {rho:+.2f}, p = {p:.2g} | Kruskal-Wallis p = {pk:.2g}")
