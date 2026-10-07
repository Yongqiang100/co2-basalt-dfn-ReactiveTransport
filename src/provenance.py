"""Capture what produced a run. Given the dfnWorks version findings, this is
the difference between an auditable ensemble and an unauditable one."""
from __future__ import annotations
import os, sys, json, subprocess, platform, datetime, hashlib

def _git(path, *args):
    try:
        return subprocess.run(["git","-C",str(path),*args], capture_output=True,
                              text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None

def _sha256(path, limit=None):
    h = hashlib.sha256()
    try:
        with open(path,"rb") as f:
            for chunk in iter(lambda: f.read(1<<20), b""):
                h.update(chunk)
                if limit and f.tell() > limit: break
        return h.hexdigest()[:16]
    except Exception:
        return None

DFNWORKS_SRC = os.environ.get("DFNWORKS_SRC", "/opt/sw/dfnworks/dfnWorks")
# Upstream commits (2026-08-01) that changed seed -> NumPy/random mapping
SEED_COMMITS = ["210a910", "fb5f0c2", "fce27ef"]

def capture(extra: dict | None = None) -> dict:
    p = {
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "host": platform.node(),
        "python": sys.version.split()[0],
        "slurm_job": os.environ.get("SLURM_JOB_ID"),
        "slurm_task": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "cwd": os.getcwd(),
        "revision_pkg_build": (lambda f: open(f).read().strip()
            if os.path.exists(f) else "unknown")(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "VERSION")),
    }
    p["deposit_commit"] = _git(".", "rev-parse", "HEAD")
    # Distinguish MODIFIED TRACKED files (which would break the reproducibility
    # guarantee) from untracked additions like revision/ and dfn_library/.
    # Reporting a bare "dirty: true" in an archived record implies the published
    # code was edited.
    _st = _git(".", "status", "--porcelain") or ""
    _mod = [l[3:] for l in _st.splitlines() if l[:2].strip() and not l.startswith("??")]
    p["deposit_modified_tracked"] = _mod          # MUST be empty
    p["deposit_untracked_count"] = sum(1 for l in _st.splitlines() if l.startswith("??"))
    p["deposit_clean"] = (len(_mod) == 0)

    d = {"src": DFNWORKS_SRC,
         "commit": _git(DFNWORKS_SRC, "rev-parse", "HEAD"),
         "describe": _git(DFNWORKS_SRC, "describe", "--tags", "--always", "--dirty"),
         "date": _git(DFNWORKS_SRC, "log", "-1", "--format=%ad", "--date=short")}
    # Does this build contain the 2026-08-01 seeding change?
    if d["commit"]:
        anc = {}
        for c in SEED_COMMITS:
            # Distinguish "not an ancestor" from "unknown to this clone".
            # A shallow clone reports neither, and conflating them would
            # silently record has_seeding_change=False for a build that has it.
            known = subprocess.run(["git","-C",DFNWORKS_SRC,"cat-file","-e",c+"^{commit}"],
                                   capture_output=True).returncode == 0
            if not known:
                anc[c] = "unknown_to_clone"
                continue
            r = subprocess.run(["git","-C",DFNWORKS_SRC,"merge-base",
                                "--is-ancestor",c,"HEAD"], capture_output=True)
            anc[c] = (r.returncode == 0)
        d["has_2026_08_01_seeding_change"] = (
            True if any(v is True for v in anc.values())
            else ("indeterminate" if "unknown_to_clone" in anc.values() else False))
        d["seed_commit_ancestry"] = anc
    # /opt/sw/dfnworks/dfnWorks may not be a git checkout (the setup guide's
    # clone can be discarded after building). Detect the 2026-08-01 seeding
    # change from the INSTALLED SOURCE instead, which is authoritative anyway:
    # commit fce27ef added DFN.set_seed(), called before hydraulic properties
    # are assigned.
    if d.get("has_2026_08_01_seeding_change") is None:
        d["has_2026_08_01_seeding_change"] = "unknown"
        try:
            import pydfnworks as _pd
            root = os.path.dirname(_pd.__file__)
            found = False
            for dp, _, fns in os.walk(root):
                for fn in fns:
                    if not fn.endswith(".py"):
                        continue
                    try:
                        txt = open(os.path.join(dp, fn), errors="ignore").read()
                    except Exception:
                        continue
                    if "def set_seed" in txt:
                        found = True
                        d["set_seed_defined_in"] = os.path.relpath(
                            os.path.join(dp, fn), root)
                        break
                if found:
                    break
            d["has_2026_08_01_seeding_change"] = found
            d["detected_via"] = "installed source (def set_seed)"
        except Exception as e:
            d["seed_detect_error"] = str(e)
    p["dfnworks"] = d

    try:
        import pydfnworks
        p["pydfnworks"] = {"version": getattr(pydfnworks,"__version__",None),
                           "path": os.path.dirname(pydfnworks.__file__)}
    except Exception:
        p["pydfnworks"] = None

    for mod in ("numpy","scipy","h5py","networkx"):
        try:
            p.setdefault("libs",{})[mod] = __import__(mod).__version__
        except Exception:
            p.setdefault("libs",{})[mod] = None

    exe = os.environ.get("PFLOTRAN_EXE")
    if exe and os.path.exists(exe):
        p["pflotran"] = {"exe": exe, "sha256_16": _sha256(exe, limit=1<<26)}
    p["modules"] = os.environ.get("LOADEDMODULES")
    if extra: p.update(extra)
    return p

def write(path, extra=None):
    rec = capture(extra)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path,"w") as f: json.dump(rec, f, indent=2, sort_keys=True)
    return rec

if __name__ == "__main__":
    print(json.dumps(capture(), indent=2, sort_keys=True))
