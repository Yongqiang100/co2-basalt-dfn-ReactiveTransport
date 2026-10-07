#!/usr/bin/env python3
"""
Percolation / outflow-boundary audit -- runs on ARCHIVED METADATA ONLY.

Why this matters
----------------
run_pflotran.py:226-234 attaches the outflow boundary only when
boundary_right_e.ex exists AND has more than one line. Otherwise the deck gets
no outflow coupler at all and the realisation runs as a CLOSED SYSTEM under
continuous mass injection: no throughflow, rising pressure, batch-like
chemistry. That is a different physical problem, not a poorly-connected one.

prepare_dfn.py:211-216 already stored the per-file .ex line counts in
dfn_summary.json under "boundary_cells", so every realisation can be classified
without re-running anything.

If the low-trapping realisations turn out to be the closed-system ones, that is
a confound with the paper's central mechanism and has to become a stated scope
condition rather than a clause in section 2.3.2.

Usage
-----
    python3 src/audit.py                     # table + confound test
    python3 src/audit.py --csv audit.csv     # also write the SI table
"""
from __future__ import annotations
import os, sys, json, glob, argparse, csv

HERE = os.path.dirname(os.path.abspath(__file__))
DEPOSIT = os.environ.get("DEPOSIT_ROOT", os.path.abspath(os.path.join(HERE, "..", "..")))
OUTFLOW_EX = "boundary_right_e.ex"


def load():
    rows = []
    pat = os.path.join(DEPOSIT, "dfn_library", "*", "dfn_summary.json")
    for f in sorted(glob.glob(pat)):
        try:
            d = json.load(open(f))
        except Exception as e:
            print(f"  WARN unreadable {f}: {e}", file=sys.stderr)
            continue
        bc = d.get("boundary_cells") or {}
        n_out = bc.get(OUTFLOW_EX, 0)
        dom = d.get("domain") or []
        rows.append({
            "name": d.get("name") or os.path.basename(os.path.dirname(f)),
            "seed": d.get("seed"),
            "p32_mult": d.get("p32_mult"),
            "num_fractures": d.get("num_fractures"),
            "n_nodes": d.get("n_nodes"),
            "n_elements": d.get("n_elements"),
            "n_connections": d.get("n_connections"),
            "conn_per_frac": d.get("conn_per_frac"),
            "domain_L": dom[0] if dom else None,
            "h": d.get("h"),
            "ex_files": d.get("ex_files"),
            "outflow_ex_lines": n_out,
            # run_pflotran.py's own condition, reproduced exactly
            "has_outflow": bool(n_out and n_out > 1),
            "dfnworks_version_recorded": any(
                k for k in d if "version" in k.lower() or "commit" in k.lower()),
        })
    return rows


def join_outcomes(rows):
    """Attach published trapping outcomes if betweenness_results.csv is present."""
    p = os.path.join(DEPOSIT, "betweenness_results.csv")
    if not os.path.exists(p):
        return False
    by = {}
    with open(p) as f:
        for r in csv.DictReader(f):
            by[r.get("name", "").strip()] = r
    n = 0
    for row in rows:
        r = by.get(row["name"])
        if not r:
            continue
        n += 1
        for k_src, k_dst in [("carb_per_cell_final", "carb_final"),
                             ("carb_per_cell_peak", "carb_peak"),
                             ("ncells", "ncells_csv")]:
            v = r.get(k_src)
            try:
                row[k_dst] = float(v) if v not in (None, "") else None
            except ValueError:
                row[k_dst] = None
    return n > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv")
    a = ap.parse_args()

    rows = load()
    if not rows:
        sys.exit(f"No dfn_summary.json found under {DEPOSIT}/dfn_library/\n"
                 f"  set DEPOSIT_ROOT if the deposit is elsewhere")
    have_out = join_outcomes(rows)

    print(f"{len(rows)} realisation(s) from {DEPOSIT}/dfn_library/\n")
    hdr = f"{'name':<18}{'p32':>5}{'frac':>6}{'cells':>8}{'conn':>6}{'outEX':>7}{'outflow':>9}"
    if have_out:
        hdr += f"{'carb/cell':>12}"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: (x["p32_mult"] or 0, str(x["seed"]))):
        line = (f"{r['name']:<18}{r['p32_mult'] or 0:>5.2f}{r['num_fractures'] or 0:>6}"
                f"{r['n_nodes'] or 0:>8}{r['n_connections'] or 0:>6}"
                f"{r['outflow_ex_lines']:>7}{'yes' if r['has_outflow'] else 'CLOSED':>9}")
        if have_out:
            c = r.get("carb_final")
            line += f"{(f'{c:.3e}' if c is not None else 'n/a'):>12}"
        print(line)

    closed = [r for r in rows if not r["has_outflow"]]
    print(f"\n{'='*72}")
    print(f"CLOSED-SYSTEM REALISATIONS: {len(closed)} of {len(rows)}")
    if closed:
        for r in closed:
            print(f"  {r['name']}  (P32x{r['p32_mult']}, {OUTFLOW_EX} = "
                  f"{r['outflow_ex_lines']} line(s))")
        print("\n  These ran with NO outflow coupler -> closed domain under")
        print("  continuous injection. Report them as a separate flow regime.")
    else:
        print("  None. Every realisation had a usable outflow boundary, so the")
        print("  closed-system confound does not apply to this ensemble.")

    # ---- the confound test -------------------------------------------------
    if have_out and closed:
        oc = [r["carb_final"] for r in rows if r["has_outflow"] and r.get("carb_final") is not None]
        cc = [r["carb_final"] for r in closed if r.get("carb_final") is not None]
        if oc and cc:
            print(f"\n  CONFOUND TEST")
            print(f"    open   n={len(oc):2d}  median carb/cell = {sorted(oc)[len(oc)//2]:.3e}")
            print(f"    closed n={len(cc):2d}  median carb/cell = {sorted(cc)[len(cc)//2]:.3e}")
            zc = sum(1 for v in cc if v == 0)
            zo = sum(1 for v in oc if v == 0)
            print(f"    zero-trapping: {zc}/{len(cc)} closed vs {zo}/{len(oc)} open")
            if zc and zc / len(cc) > zo / max(len(oc), 1):
                print("    -> zero trapping is CONCENTRATED in closed-system cases.")
                print("       Treat them separately; do not attribute this to")
                print("       flow-reaction overlap.")

    if not any(r["dfnworks_version_recorded"] for r in rows):
        print(f"\n{'='*72}")
        print("NOTE: no dfnWorks/pydfnworks version or commit is recorded in any")
        print("dfn_summary.json. The archive cannot attest which build produced")
        print("the published ensemble -- the manuscript text is the only record.")
        print("provenance.py now captures this for every new realisation.")

    if a.csv:
        keys = sorted({k for r in rows for k in r})
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {a.csv} ({len(rows)} rows) -- basis for the SI table.")


if __name__ == "__main__":
    main()
