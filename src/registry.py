"""Run registry: one JSON file tracking every run's state. Needed because the
run set is ~120 runs across four blocks and Slurm's own accounting is purged."""
from __future__ import annotations
import json, os, fcntl, datetime, glob

REG = os.environ.get("REVISION_REGISTRY", "runs/registry.json")
STATES = ("planned","queued","running","done","failed","invalid")

def _load():
    if not os.path.exists(REG): return {}
    with open(REG) as f: return json.load(f)

def _lock():
    """Lock the REGISTRY, not the temp file. Slurm array tasks call upsert()
    concurrently; locking the temp file protected nothing and lost records."""
    os.makedirs(os.path.dirname(REG) or ".", exist_ok=True)
    fh = open(REG + ".lock", "a")
    fcntl.flock(fh, fcntl.LOCK_EX)
    return fh

def _save(d, fh=None):
    tmp = f"{REG}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(d, f, indent=2, sort_keys=True)
    os.replace(tmp, REG)

def upsert(run_id, **fields):
    fh = _lock()          # held until this function returns
    d = _load()
    r = d.setdefault(run_id, {"run_id": run_id, "state": "planned",
                              "history": []})
    if "state" in fields and fields["state"] != r.get("state"):
        r["history"].append({
            "at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "state": fields["state"]})
    r.update(fields)
    _save(d)
    fh.close()
    return r

def summary():
    d = _load()
    by = {}
    for r in d.values():
        by.setdefault(r.get("block","?"), {}).setdefault(r.get("state","?"), 0)
        by[r.get("block","?")][r.get("state","?")] += 1
    return by

def verify(results_root="pflotran_results", final_year=50.0):
    """Mark runs done/failed by inspecting output, not by trusting Slurm exit."""
    import h5py
    d = _load(); changed = 0
    for rid, r in d.items():
        rd = r.get("run_dir") or os.path.join(results_root, rid)
        h5 = sorted(glob.glob(os.path.join(rd, "*.h5")))
        if not h5:
            if r["state"] in ("running","queued"): r["state"]="failed"; changed+=1
            continue
        try:
            with h5py.File(h5[-1],"r") as f:
                ts = [k for k in f.keys() if k.startswith("Time")]
                yrs = max(float(k.split()[1]) for k in ts) if ts else 0.0
            r["last_output_year"] = yrs
            r["state"] = "done" if yrs >= final_year*0.999 else "failed"
            if yrs < final_year*0.999:
                r["note"] = f"stalled at t={yrs:g} yr (expected {final_year:g})"
            changed += 1
        except Exception as e:
            r["state"], r["note"] = "failed", f"unreadable: {e}"; changed += 1
    _save(d); return changed

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1 and sys.argv[1]=="verify":
        print(f"updated {verify()} record(s)")
    print(f"{'BLOCK':<10}", " ".join(f"{s:>9}" for s in STATES))
    for b, st in sorted(summary().items()):
        print(f"{b:<10}", " ".join(f"{st.get(s,0):>9}" for s in STATES))
