#!/usr/bin/env python3
"""
Build a Block C variant deck: emit the baseline via the published code, apply
declarative edits, save the deck + its diff + provenance, and register the run.

  python src/build_variant.py --list
  python src/build_variant.py --variant anor_as30 --dry
  python src/build_variant.py --variant anor_as30 --index 0

Always read the emitted .diff once per variant before queueing the runs. The
engine guarantees an op MATCHED; only you can confirm it matched the right thing.
"""
import argparse, os, sys, json, shutil, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "config")]
from deckmod import Deck, apply_variant, DeckError
import provenance, registry
import variants as V

REPO = os.environ.get("DEPOSIT_ROOT", os.path.join(HERE, "..", ".."))
OUT  = os.environ.get("REVISION_RUNS", os.path.join(HERE, "..", "runs"))

def baseline_deck(case, workdir):
    """Emit the unmodified deck using the PUBLISHED run_pflotran.py.

    `--write_only` (run_pflotran.py:882) writes the deck and stages full_mesh.uge,
    full_mesh.inp, the thermodynamic databases and every *.ex boundary file into
    pflotran_results/<case>/, then returns without launching MPI. We copy that
    staged directory to the variant directory and edit the deck there, so the
    baseline is produced entirely by the published code path.
    """
    import glob as _glob
    dfn = os.path.join(REPO, "dfn_library", case)
    if not os.path.isdir(dfn):
        sys.exit(f"FATAL: {dfn} not found (is dfn_library/ unpacked from Zenodo?)")

    # run_pflotran.py misdetects any Slurm job as Setonix (see slurm/env.sh),
    # so these must be exported or it resolves non-existent Pawsey paths.
    for v in ("PFLOTRAN_EXE", "PFLOTRAN_DB", "CO2_DB"):
        if not os.environ.get(v) or not os.path.isfile(os.environ[v]):
            sys.exit(f"FATAL: {v} unset or not a file.\n"
                     f"  Run inside a Slurm job that sources slurm/env.sh, or:\n"
                     f"    source revision/slurm/env.sh\n"
                     f"  before calling build_variant.py")
    staged = os.path.join(REPO, "pflotran_results", case)
    r = subprocess.run([sys.executable, "run_pflotran.py",
                        "--dfn", case, "--write_only"],
                       cwd=REPO, capture_output=True, text=True)
    decks = _glob.glob(os.path.join(staged, "*.in"))
    if not decks:
        sys.exit("FATAL: run_pflotran.py --write_only produced no *.in in\n"
                 f"  {staged}\nstdout:\n{r.stdout[-800:]}\nstderr:\n{r.stderr[-1200:]}")

    # Copy ONLY what a run needs. A blanket copy would clone any *.h5 output
    # from a previous published run into all 29 variant directories.
    os.makedirs(workdir, exist_ok=True)
    KEEP_EXT = {".in", ".uge", ".inp", ".ex", ".dat"}
    for f in os.listdir(staged):
        src = os.path.join(staged, f)
        if os.path.isfile(src) and os.path.splitext(f)[1].lower() in KEEP_EXT:
            shutil.copy2(src, workdir)
    return os.path.join(workdir, os.path.basename(decks[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant"); ap.add_argument("--index", type=int)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list or not a.variant:
        tot = 0
        for k, v in V.VARIANTS.items():
            print(f"  {k:<18s} {len(v['cases'])} case(s)  {v['comment'][:58]}")
            tot += len(v["cases"])
        print(f"\n  {len(V.VARIANTS)} variants, {tot} runs")
        print("\n  not yet buildable (need more than a deck edit):")
        for k, why in V.NEEDS_MORE_WORK.items():
            print(f"    {k:<20s} {why[:66]}")
        return

    spec = V.VARIANTS.get(a.variant) or sys.exit(f"unknown variant {a.variant!r}")
    cases = spec["cases"]
    if a.index is not None:
        if not 0 <= a.index < len(cases):
            sys.exit(f"--index out of range 0..{len(cases)-1}")
        cases = [cases[a.index]]

    for case in cases:
        rid = f"C_{a.variant}__{case}"
        wd  = os.path.join(OUT, rid)
        print(f"\n=== {rid} ===")
        if a.dry:
            print(f"  would apply {len(spec['ops'])} op(s) to {case}")
            for o in spec["ops"]: print(f"    {o[0]}{o[1]}")
            continue
        baseline_path = baseline_deck(case, wd)
        deck = Deck.from_file(baseline_path)
        deck.name = rid
        try:
            apply_variant(deck, spec["ops"])
        except DeckError as e:
            registry.upsert(rid, block="C", variant=a.variant, case=case,
                            state="invalid", note=str(e))
            sys.exit(f"  DECK EDIT FAILED -- run marked invalid:\n  {e}")
        deck.write(baseline_path)  # overwrite the staged copy in place
        with open(os.path.join(wd, "variant.diff"), "w") as f:
            f.write(deck.diff())
        provenance.write(os.path.join(wd, "provenance.json"),
                         extra={"run_id": rid, "block": "C",
                                "variant": a.variant, "case": case,
                                "comment": spec["comment"],
                                "deck_ops": deck.provenance()})
        registry.upsert(rid, block="C", variant=a.variant, case=case,
                        run_dir=wd, state="planned", n_ops=len(deck.ops))
        nd = len(deck.diff().splitlines())
        print(f"  {len(deck.ops)} op(s) applied | diff {nd} lines -> {wd}/variant.diff")
        print(f"  READ THE DIFF before queueing.")

if __name__ == "__main__":
    main()
