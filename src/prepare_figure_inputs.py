#!/usr/bin/env python3
"""
Inputs for ensemble_statistics_figures.py from the corrected coupled runs.

ensemble_statistics_figures.py finds its columns by name, so its code and layout stay unchanged.
This writes the two files it reads, with the column names it accepts:

  results/figure_input_blockA.csv    for fig_carbonate_dissolution_variability.pdf: run_id, p32_mult,
                                  carb_per_cell, forsterite_dissolved_pct, diopside_dissolved_pct,
                                  anorthite_dissolved_pct
      carb_per_cell from results/flow_index_test.csv (Block A)
      p32 and dissolution from results/A_coupled_networks.csv
  results/figure_input_sensitivity.csv      for fig_parameter_sensitivity.pdf: variant, ratio
      from results/sensitivity_ratios.csv (carbonate per cell as in blockc2.py),
      without vf_consistent, identical to the baseline since the correction

    python3 src/prepare_figure_inputs.py
    python3 src/ensemble_statistics_figures.py --runs runs_gravityoff --out figures_R2 \\
        --blocka results/figure_input_blockA.csv --sens results/figure_input_sensitivity.csv
"""
import csv, os, sys
REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda f: os.path.join(REV, "results", f)


def rows(f):
    p = R(f)
    if not os.path.isfile(p):
        sys.exit(f"missing {p}")
    return list(csv.DictReader(open(p)))


def main():
    cpc = {r["network"]: r["carb_per_cell"] for r in rows("flow_index_test.csv") if r.get("block") == "A"}
    out = []
    for r in rows("A_coupled_networks.csv"):
        n = r["network"]
        if n not in cpc:
            print(f"  no carbonate per cell for {n}, skipped"); continue
        out.append(dict(run_id=n, p32_mult=float(r["p32_factor"]), carb_per_cell=cpc[n],
                        forsterite_dissolved_pct=r["forsterite_dissolved_pct"], diopside_dissolved_pct=r["diopside_dissolved_pct"],
                        anorthite_dissolved_pct=r["anorthite_dissolved_pct"]))
    with open(R("figure_input_blockA.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"  wrote results/figure_input_blockA.csv ({len(out)} networks)")
    sens = [dict(variant=r["variant"], ratio=r["ratio"]) for r in rows("sensitivity_ratios.csv") if r["variant"] != "vf_consistent"]
    with open(R("figure_input_sensitivity.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["variant", "ratio"]); w.writeheader(); w.writerows(sens)
    print(f"  wrote results/figure_input_sensitivity.csv ({len(sens)} rows, {len({r['variant'] for r in sens})} variants)")


if __name__ == "__main__":
    main()
