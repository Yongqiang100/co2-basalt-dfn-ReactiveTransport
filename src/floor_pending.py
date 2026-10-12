#!/usr/bin/env python3
"""
Add the porosity floor (MINIMUM_POROSITY 1.d-2) to every coupled deck whose
job has not started yet, then optionally release the jobs.

    python3 src/floor_pending.py 49948902:lists/A_feedback.txt ...                     # plan only
    python3 src/floor_pending.py 49948902:lists/A_feedback.txt ... --apply             # edit
    python3 src/floor_pending.py 49948902:lists/A_feedback.txt ... --apply --release   # edit, release
    python3 src/floor_pending.py --dirs runs/A_feedback__p32_100_s727 --apply          # named directories

Only PENDING tasks are touched: a running or finished run keeps the deck it
ran with, so deck and result always agree. A pending task that is not held
is held before its deck is edited, so it cannot start mid-edit. Fixed-porosity
decks are left alone. Each edited run is appended to lists/minpor_applied.txt.
"""
import datetime, os, re, subprocess, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import add_coupling as ac

RECORD = os.path.join(REV, "lists", "minpor_applied.txt")


def sh(args):
    return subprocess.run(args, capture_output=True, text=True)


def pending(job):
    """{array index: reason} for the job's pending tasks."""
    out = {}
    for line in sh(["squeue", "-r", "-h", "-j", job, "-t", "PD", "-o", "%K %r"]).stdout.splitlines():
        t = line.split()
        if len(t) >= 2 and t[0].isdigit():
            out[int(t[0])] = t[1]
    return out


def resolve(rid):
    return os.path.join(REV, rid) if "/" in rid else os.path.join(REV, "runs", rid)


def deck_of(d):
    ins = sorted(f for f in os.listdir(d) if f.endswith(".in")) if os.path.isdir(d) else []
    return os.path.join(d, ins[0]) if ins else None


def has(deck, key):
    return re.search(rf"^\s*{key}\b", open(deck).read(), re.M | re.I) is not None


def main():
    a = sys.argv[1:]
    apply, release = "--apply" in a, "--release" in a
    dirs = []
    if "--dirs" in a:
        dirs = [x for x in a[a.index("--dirs") + 1].split(",") if x]
    pairs = [x for x in a if ":" in x and not x.startswith("--")]
    if not pairs and not dirs:
        sys.exit(__doc__)

    targets, jobs = [], []
    for pair in pairs:
        job, lst = pair.split(":", 1)
        ids = [l.strip() for l in open(os.path.join(REV, lst) if not os.path.isabs(lst) else lst) if l.strip()]
        pend = pending(job); jobs.append(job)
        n = dict(coupled=0, floored=0, fixed=0, held=0)
        for idx, reason in sorted(pend.items()):
            if idx >= len(ids):
                print(f"  WARN {job}_{idx}: no entry in {lst}"); continue
            deck = deck_of(resolve(ids[idx]))
            if deck is None:
                print(f"  WARN {job}_{idx}: no deck for {ids[idx]}"); continue
            if not has(deck, "UPDATE_POROSITY"):
                n["fixed"] += 1; continue
            if has(deck, "MINIMUM_POROSITY"):
                n["floored"] += 1; continue
            n["coupled"] += 1; n["held"] += reason == "JobHeldUser"
            targets.append((job, idx, reason, ids[idx], deck))
        print(f"  {job} ({lst}): {len(pend)} pending | {n['coupled']} coupled to floor "
              f"({n['held']} already held) | {n['floored']} already floored | {n['fixed']} fixed, left alone")
    for d in dirs:
        deck = deck_of(resolve(d))
        if deck is None:
            print(f"  WARN {d}: no deck"); continue
        if has(deck, "UPDATE_POROSITY") and not has(deck, "MINIMUM_POROSITY"):
            targets.append((None, None, "named", d, deck))
        else:
            print(f"  {d}: {'already floored' if has(deck, 'MINIMUM_POROSITY') else 'not coupled'}; left alone")

    print(f"decks to receive the floor: {len(targets)}")
    if not apply:
        print("plan only: nothing changed. Re-run with --apply (and --release to free the jobs)."); return

    done, stamp = [], datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    for job, idx, reason, rid, deck in targets:
        if job and reason != "JobHeldUser":
            if sh(["scontrol", "hold", f"{job}_{idx}"]).returncode:
                print(f"  SKIPPED {rid}: could not hold {job}_{idx} (it may have just started)"); continue
        msg = ac.edit(deck)
        if not has(deck, "MINIMUM_POROSITY"):
            print(f"  FAILED {rid}: {msg}"); continue
        done.append(f"{rid}\t{stamp}\t{job + '_' + str(idx) if job else 'named'}")
    os.makedirs(os.path.dirname(RECORD), exist_ok=True)
    with open(RECORD, "a") as f:
        f.writelines(line + "\n" for line in done)
    print(f"floor added to {len(done)} of {len(targets)} decks; recorded in {os.path.relpath(RECORD, REV)}")
    if release:
        for job in jobs:
            r = sh(["scontrol", "release", job])
            print(f"  released {job}" if r.returncode == 0 else f"  release {job}: {r.stderr.strip() or 'no held tasks'}")


if __name__ == "__main__":
    main()
