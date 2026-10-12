#!/usr/bin/env python3
"""
Floor and resubmit every outstanding coupled failure that ran without the floor.

    python3 src/collect_failures.py            # refresh lists/failures.tsv first
    python3 src/resubmit_coupled.py            # plan
    python3 src/resubmit_coupled.py --apply    # floor, resubmit, record
    python3 src/resubmit_coupled.py --apply --time 16:00:00 --stop 14   # longer limit, with checkpoints

For each failure in lists/failures.tsv with no later completion, either
  config "coupled" (no floor at run time: a clogging failure the floor addresses), or
  a resource limit (TIMEOUT, OUT_OF_MEMORY, NODE_FAIL) in any coupled run, floored or not.
A floored run that failed numerically is left alone: it needs analysis, not a retry.
For each one: adds MINIMUM_POROSITY to its deck, finds the
list it belongs to, and resubmits the failed runs as one array per list, sized
for the block (D40: 128 ranks, 112G, 24 h; D30: 64, 56G, 16 h; others: 32, 56G,
8 h; --time overrides the limit, and --stop H adds checkpointing with a
graceful stop after H hours, via add_checkpoint.py). New job IDs are appended to the set's line in lists/jobs.txt, and each
run to lists/resubmitted.txt, so a run is never resubmitted twice.
"""
import datetime, glob, os, re, subprocess, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import add_coupling as ac
import add_checkpoint as ack

LISTS = os.path.join(REV, "lists")
RECORDS = {"failures.tsv", "jobs.txt", "minpor_applied.txt", "resubmitted.txt"}
SIZE = {"D40": ("128", "112G", "24:00:00"), "D30": ("64", "56G", "16:00:00")}


def read_tsv(p):
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    head = rows[0]
    return [dict(zip(head, r + [""] * (len(head) - len(r)))) for r in rows[1:]]


def list_of(run):
    """(list file, index) of the run; lists with 'feedback' in the name first."""
    cands = sorted((p for p in glob.glob(os.path.join(LISTS, "*.txt")) if os.path.basename(p) not in RECORDS),
                   key=lambda p: ("feedback" not in os.path.basename(p), p))
    for p in cands:
        ids = [l.strip() for l in open(p) if l.strip()]
        for i, rid in enumerate(ids):
            if rid == run or os.path.basename(rid) == run:
                return p, i
    return None, None


def main():
    apply = "--apply" in sys.argv
    a = sys.argv
    t_over = a[a.index("--time") + 1] if "--time" in a else None
    stop = float(a[a.index("--stop") + 1]) if "--stop" in a else None
    tsv = os.path.join(LISTS, "failures.tsv")
    if not os.path.isfile(tsv):
        sys.exit("no lists/failures.tsv: run src/collect_failures.py first")
    done_before = set()
    rec = os.path.join(LISTS, "resubmitted.txt")
    if os.path.isfile(rec):
        done_before = {l.split("\t")[0] for l in open(rec) if l.strip()}
    # STOPPED_EARLY runs are continued with make_restart.py, never restarted from scratch
    RES = ("TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL")
    rows = [r for r in read_tsv(tsv) if not r["status"] and r["state"] != "STOPPED_EARLY"]
    todo = [r for r in rows if r["config"] == "coupled" or (r["config"].startswith("coupled") and r["state"] in RES)]
    held = [r["run"] for r in rows if r["config"] == "coupled+floor" and r["state"] not in RES]
    if held:
        print(f"  floored runs that failed numerically, left for analysis: {', '.join(held)}")
    groups, skipped = {}, []
    for r in todo:
        run = r["run"]
        if run in done_before:
            skipped.append(run); continue
        lst, idx = list_of(run)
        if lst is None:
            print(f"  WARN {run}: in no list under lists/; resubmit by hand"); continue
        groups.setdefault((lst, r["set"]), []).append((idx, run))
    if skipped:
        print(f"  already resubmitted earlier (lists/resubmitted.txt), left alone: {', '.join(skipped)}")
    if not groups:
        print("nothing to resubmit"); return
    for (lst, sset), items in sorted(groups.items()):
        key = next((k for k in SIZE if os.path.basename(lst).startswith(k)), None)
        n, mem, t = SIZE.get(key, ("32", "56G", "08:00:00"))
        t = t_over or t
        arr = ",".join(str(i) for i, _ in sorted(items))
        print(f"  {sset}: {len(items)} run(s) from {os.path.relpath(lst, REV)}, indices {arr}, {n} ranks, {mem}, {t}"
              + (f", checkpoints with a stop at {stop:g} h" if stop else ""))
        for _, run in sorted(items):
            print(f"      {run}")
        if not apply:
            continue
        stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        for _, run in items:
            deck = os.path.join(REV, "runs", run, "pflotran_co2.in")
            had_floor = bool(re.search(r"^\s*MINIMUM_POROSITY", open(deck).read(), re.M | re.I))
            msg = ac.edit(deck)
            if not had_floor:                     # recorded, so collect_failures.py knows it came after the failure
                with open(os.path.join(LISTS, "minpor_applied.txt"), "a") as f:
                    f.write(f"{run}\t{stamp}\tresubmit_coupled\n")
            if not re.search(r"^\s*MINIMUM_POROSITY", open(deck).read(), re.M | re.I):
                sys.exit(f"  FATAL {run}: floor not added ({msg}); nothing submitted for this group")
            if stop:
                m2 = ack.edit(deck, stop, 5.0)
                if m2.startswith("NOT WRITTEN"):
                    sys.exit(f"  FATAL {run}: checkpointing not added ({m2}); nothing submitted for this group")
        r = subprocess.run(["sbatch", "--parsable", f"--array={arr}", f"--ntasks={n}", f"--mem={mem}", f"--time={t}",
                            f"--export=ALL,LIST={os.path.relpath(lst, REV)}", "slurm/run_dirs.sh"],
                           cwd=REV, capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"  FATAL sbatch: {r.stderr.strip()}")
        job = r.stdout.strip().split(";")[0]
        print(f"      submitted as job {job}")
        jt = os.path.join(LISTS, "jobs.txt")
        lines = open(jt).read().split("\n")
        hit = False
        for k, l in enumerate(lines):
            t_ = l.split()
            if len(t_) >= 2 and t_[0] == sset:
                lines[k] = l.rstrip() + "," + job; hit = True
        if not hit:
            lines.append(f"{sset:<14} {job}")
        open(jt, "w").write("\n".join(lines))
        with open(rec, "a") as f:
            for _, run in items:
                f.write(f"{run}\t{job}\t{stamp}\n")
    if not apply:
        print("plan only: nothing changed. Re-run with --apply.")


if __name__ == "__main__":
    main()
