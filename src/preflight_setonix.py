#!/usr/bin/env python3
"""
Preflight for the option-B (gravity off) reruns on Setonix. Changes nothing.

    cd $MYSCRATCH/co2-basalt/revision
    python3 src/preflight_setonix.py --skip A:43     # A:43 is p32_200_s941, left out on purpose
    python3 src/preflight_setonix.py --blocks A --skip A:43   # only what block A needs

Each check reports PASS / WARN / FAIL; exit status 1 if anything FAILs.
  1. launchers   run_rt.sh, run_pulse.sh, run_dirs.sh: hook once, before PFLOTRAN,
                 Setonix Slurm settings (work, pawsey1284, <= 24 h, no pmix)
  2. databases   the copy the hook would supply equals the copy run_rt.sh stages
  3. staging     the run_rt.sh correction chain on a template deck, in a temporary
                 directory, diffed against the old production deck
  4. reuse       run directories run_rt.sh would reuse (double VF rescale), and inputs
  5. pre-staged  every directory in lists/*.txt (or runs/F_* if no lists): corrected
                 once, Richards mode, one GRAVITY line, databases reachable
"""
import argparse, difflib, glob, os, re, shutil, subprocess, sys, tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
ROOT = os.path.dirname(REV)
SRC = os.path.join(REV, "src")
PRIM = ("Anorthite", "Albite", "Diopside", "Forsterite", "Fayalite", "Enstatite")
GRAV = re.compile(r"^[ \t]*GRAVITY[ \t]+(\S+)[ \t]+(\S+)[ \t]+(\S+)", re.M | re.I)
OUT = []
sys.path.insert(0, SRC)


def detect_machine():
    host = os.uname().nodename.lower()
    if os.environ.get("PAWSEY_CLUSTER", "").lower() == "setonix" or "setonix" in host or host.startswith("nid"):
        return "setonix"
    return "hpc01"


def rep(status, what, detail=""):
    OUT.append(status)
    print(f"  {status:4s}  {what}" + (f"\n        {detail}" if detail else ""))


def num(t):
    return float(t.replace("d", "e").replace("D", "E"))


def facts(text):
    vf = 0.0
    for line in text.splitlines():
        t = line.split()
        # plain decimals (0.176471) as well as Fortran exponents (0.30d0, 1.d-1)
        if len(t) >= 2 and t[0] in PRIM and re.fullmatch(r"[0-9]*\.?[0-9]+([dDeE][+-]?[0-9]+)?", t[1]):
            vf += num(t[1])
    return dict(richards=bool(re.search(r"^[ \t]*MODE[ \t]+RICHARDS\b", text, re.M | re.I)),
                gravity=[tuple(num(v) for v in g) for g in GRAV.findall(text)], vf=vf,
                scaled="SCALED_MASS_RATE" in text, coupled=bool(re.search(r"UPDATE_POROSITY", text, re.I)))


def live(lines):
    return [(k, l) for k, l in enumerate(lines) if l.strip() and not l.lstrip().startswith("#")]


# ------------------------------------------------------------------ 1
def check_launchers(machine):
    print(f"\n1. launchers ({machine})")
    if machine == "setonix" and os.path.isfile(os.path.join(REV, "slurm", "run_F.sh")):
        rep("WARN", "slurm/run_F.sh is an hpc01 one-off (partition normal, pmix, one network); "
                    "run the F directories through run_dirs.sh instead")
    names = (("slurm/run_rt.sh", "slurm/run_pulse.sh", "slurm/run_dirs.sh") if machine == "setonix"
             else ("slurm/run_rt.sh", "slurm/run_dirs_hpc01.sh"))
    for name in names:
        p = os.path.join(REV, name); tag = f"{name}:"
        if not os.path.isfile(p):
            rep("FAIL", f"{tag} not found"); continue
        if subprocess.run(["bash", "-n", p], capture_output=True).returncode:
            rep("FAIL", f"{tag} bash -n reports a syntax error"); continue
        lines = open(p).read().split("\n"); L = live(lines); text = "\n".join(lines)
        hook = [k for k, l in L if "fix_gravity_off.py" in l]
        srun = [k for k, l in L if re.match(r"\s*srun\b", l) and "pflotran" in l]
        corr = [k for k, l in L if "apply_corrections.py" in l]
        if any("fix_boundary.py" in l for _, l in L):
            rep("FAIL", f"{tag} still calls fix_boundary.py (option A)")
        if len(hook) != 1 or len(srun) != 1:
            rep("FAIL", f"{tag} gravity hook x{len(hook)}, pflotran srun x{len(srun)} (expected 1 and 1)")
        else:
            h = hook[0]; m = re.search(r'"([^"]*fix_gravity_off\.py)"', lines[h]); path = m.group(1) if m else ""
            if path.startswith("/") and (not os.path.isfile(path) or not path.startswith(REV)):
                rep("FAIL", f"{tag} hook points to a file that is missing or outside this tree", path)
            ok = h < srun[0] and (not corr or corr[0] < h)
            rep("PASS" if ok else "FAIL", f"{tag} hook runs {'after apply_corrections and ' if corr else ''}before PFLOTRAN")
            if not any("unset SLURM_CPUS_PER_TASK" in lines[k] for k in range(max(0, srun[0] - 6), srun[0])):
                rep("WARN", f"{tag} no 'unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK' before srun; "
                            "with --hint=nomultithread on 2-thread cores srun can refuse to start")
        pmix = bool(srun) and "--mpi=pmix" in " ".join(lines[k] for k in srun)
        if machine == "setonix":
            if pmix:
                rep("FAIL", f"{tag} srun uses --mpi=pmix, the hpc01 flag")
            if re.search(r"#SBATCH\s+--partition=normal", text):
                rep("FAIL", f"{tag} partition normal is hpc01's; Setonix uses work")
            if not re.search(r"#SBATCH\s+--account=", text):
                rep("FAIL", f"{tag} no #SBATCH --account")
            m = re.search(r"#SBATCH\s+--time=(\d+):(\d+):(\d+)", text)
            if m and int(m.group(1)) + int(m.group(2)) / 60 > 24:
                rep("FAIL", f"{tag} --time={m.group(1)}:{m.group(2)}:{m.group(3)} exceeds the 24 h work limit")
        else:
            if srun and not pmix:
                rep("FAIL", f"{tag} srun lacks --mpi=pmix, which hpc01 needs")
            if re.search(r"#SBATCH\s+--partition=work", text):
                rep("FAIL", f"{tag} partition work is Setonix's; hpc01 uses normal")
            if re.search(r"#SBATCH\s+--account=pawsey", text):
                rep("FAIL", f"{tag} carries the Setonix account, which hpc01 does not have")
    rt = os.path.join(REV, "slurm", "run_rt.sh")
    if os.path.isfile(rt) and "prepare_dfn.py" in open(rt).read() and not any(
            os.path.isfile(os.path.join(d, "prepare_dfn.py")) for d in (REV, ROOT)):
        rep("WARN", "run_rt.sh finds the project root through prepare_dfn.py, which is not here; "
                    f"add DEPOSIT_ROOT={ROOT} to every --export")
    if os.path.isfile(rt):
        skip = [l.strip() for l in open(rt) if re.search(r"skip|already|complete", l, re.I) and not l.lstrip().startswith("#")]
        if skip:
            rep("WARN", "run_rt.sh has lines that may skip finished cases; confirm they cannot skip reruns", " | ".join(skip[:4]))


# ------------------------------------------------------------------ 2
def check_databases(case):
    print("\n2. databases")
    tpl_dir = os.path.join(ROOT, "pflotran_results", case)
    tpl = os.path.join(tpl_dir, "pflotran_co2.in")
    if not os.path.isfile(tpl):
        rep("FAIL", f"template deck missing: {tpl}"); return
    import fix_gravity_off as fg
    env_sh = os.path.join(REV, "slurm", "env.sh")
    pfv6 = re.search(r"PFV6=(\S+)", open(env_sh).read()).group(1) if os.path.isfile(env_sh) and \
        re.search(r"PFV6=(\S+)", open(env_sh).read()) else None
    for raw in fg.DB.findall(open(tpl).read()):
        name = raw.strip("'\"")
        # a copy run_rt.sh actually staged: a new run first, then an old one
        staged = next((p for p in sorted(glob.glob(os.path.join(REV, "runs", "A_*", name))) +
                       sorted(glob.glob(os.path.join(REV, "runs_dirichlet_20260924", "runs", "*", name))) +
                       sorted(glob.glob(os.path.join(REV, "runs_dirichlet_20260924", "*", name)))
                       if os.path.isfile(p)), None)
        if staged is None:
            rep("WARN", f"{name}: no staged run yet to compare the hook's copy against"); continue
        hook = fg.db_sources(name)
        if hook is None and pfv6:
            cand = os.path.join(pfv6, "pflotran", "database", name)
            hook = cand if os.path.isfile(cand) else None
        if hook is None:
            rep("FAIL", f"{name}: the hook finds no reference copy to supply to directories lacking it")
            continue
        same = open(staged, "rb").read() == open(hook, "rb").read()
        rep("PASS" if same else "FAIL",
            f"{name}: the hook's copy {'is identical to' if same else 'DIFFERS from'} the copy run_rt.sh stages",
            f"hook: {hook}\n        staged: {staged}")


# ------------------------------------------------------------------ 3
def check_staging(case):
    print(f"\n3. staging chain on {case} (temporary copy)")
    tpl = os.path.join(ROOT, "pflotran_results", case, "pflotran_co2.in")
    lib = os.path.join(ROOT, "dfn_library", case)
    if not os.path.isfile(tpl) or not os.path.isdir(lib):
        rep("FAIL", f"template or library missing for {case}"); return
    f0 = facts(open(tpl).read())
    rep("PASS" if abs(f0["vf"] - 0.85) < 1e-3 else "WARN", f"template is uncorrected (VF sum {f0['vf']:.4f}, expect 0.8500)")
    with tempfile.TemporaryDirectory() as tmp:
        for f in glob.glob(os.path.join(lib, "*")):
            os.symlink(os.path.realpath(f), os.path.join(tmp, os.path.basename(f)))
        for f in glob.glob(os.path.join(os.path.dirname(tpl), "*.dat")):
            shutil.copy(f, tmp)
        shutil.copy(tpl, os.path.join(tmp, "pflotran_co2.in"))
        for cmd, env in (([sys.executable, os.path.join(SRC, "apply_corrections.py"), "pflotran_co2.in"], dict(os.environ, RESCALE_VF="1")),
                         ([sys.executable, os.path.join(SRC, "fix_gravity_off.py"), "pflotran_co2.in"], dict(os.environ))):
            r = subprocess.run(cmd, cwd=tmp, env=env, capture_output=True, text=True)
            if r.returncode:
                rep("FAIL", f"{os.path.basename(cmd[1])} failed on the template", (r.stderr or r.stdout)[-300:]); return
        new = open(os.path.join(tmp, "pflotran_co2.in")).read()
    f1 = facts(new)
    rep("PASS" if f1["richards"] else "FAIL", "flow mode is RICHARDS")
    rep("PASS" if f1["gravity"] == [(0.0, 0.0, 0.0)] else "FAIL", f"one GRAVITY line, all zero ({f1['gravity']})")
    rep("PASS" if abs(f1["vf"] - 0.50) < 1e-3 else "FAIL", f"primary VF rescaled once (sum {f1['vf']:.4f}, expect 0.5000)")
    rep("PASS" if f1["scaled"] else "FAIL", "injection uses SCALED_MASS_RATE")
    old = next((p for p in (os.path.join(REV, "runs_dirichlet_20260924", "runs", f"A_{case}", "pflotran_co2.in"),
                            os.path.join(REV, "runs_dirichlet_20260924", f"A_{case}", "pflotran_co2.in"),
                            os.path.join(REV, "runs", f"A_{case}", "pflotran_co2.in")) if os.path.isfile(p)), None)
    if not old:
        rep("WARN", f"no old production deck for A_{case} to diff against"); return
    norm = lambda s: [" ".join(l.split()) for l in s.splitlines() if l.strip()]
    diff = [d for d in difflib.unified_diff(norm(open(old).read()), norm(new), lineterm="", n=0)
            if d[:1] in "+-" and not d.startswith(("+++", "---")) and "GRAVITY" not in d.upper()]
    rep("PASS" if not diff else "WARN",
        "staged deck differs from the old production deck only in GRAVITY" if not diff else
        f"staged deck differs from {old} in {len(diff)} other line(s); review them", "\n        ".join(diff[:12]))


# ------------------------------------------------------------------ 4
def check_blocks(skip, blocks):
    print("\n4. reuse hazard and inputs, per run_rt.sh block")
    for block in blocks:
        r = subprocess.run([sys.executable, os.path.join(SRC, "blocks.py"), "--block", block, "--count"],
                           cwd=REV, capture_output=True, text=True)
        try:
            n = int(r.stdout.strip().split()[-1])
        except (ValueError, IndexError):
            rep("FAIL", f"block {block}: blocks.py --count gave '{r.stdout.strip()}'"); continue
        reuse, missing = [], []
        skipped = sorted(i for b, i in skip if b == block)
        for i in range(n):
            if i in skipped:
                continue
            s = subprocess.run([sys.executable, os.path.join(SRC, "blocks.py"), "--block", block, "--index", str(i)],
                               cwd=REV, capture_output=True, text=True).stdout.strip().split("\t")
            if len(s) < 3:
                missing.append(f"index {i}: no spec"); continue
            rd = s[2] if os.path.isabs(s[2]) else os.path.join(REV, s[2])
            if os.path.isfile(os.path.join(rd, "pflotran_co2.in")):
                reuse.append(os.path.basename(rd))
            m = re.search(r"(L\d+_)?p32_\d+_s\d+", "\t".join(s))
            if not m:
                missing.append(f"index {i}: no network name"); continue
            for need in (os.path.join(ROOT, "pflotran_results", m.group(0), "pflotran_co2.in"),
                         os.path.join(ROOT, "dfn_library", m.group(0), "full_mesh.uge")):
                if not os.path.isfile(need):
                    missing.append(os.path.relpath(need, ROOT))
        rep("PASS" if not reuse else "FAIL", f"block {block} ({n} cases): " +
            ("no run directory would be reused" if not reuse else f"{len(reuse)} directories would be reused"),
            ", ".join(reuse[:6]) if reuse else "")
        if missing:
            rep("FAIL", f"block {block}: {len(missing)} inputs missing", "; ".join(missing[:6]))
        if skipped:
            rep("WARN", f"block {block}: indices {', '.join(map(str, skipped))} skipped by request; "
                        f"leave them out of the sbatch array")


# ------------------------------------------------------------------ 5
def check_prestaged():
    print("\n5. pre-staged directories")
    import fix_gravity_off as fg
    lists = sorted(glob.glob(os.path.join(REV, "lists", "*.txt")))
    groups = {}
    for L in lists:
        ids = [l.strip() for l in open(L) if l.strip()]
        groups[os.path.relpath(L, REV)] = [os.path.join(REV, i) if "/" in i else os.path.join(REV, "runs", i) for i in ids]
    if not groups:
        groups["runs/F_* (no lists/ yet)"] = sorted(d for d in glob.glob(os.path.join(REV, "runs", "F_*")) if os.path.isdir(d))
    for label, dirs in groups.items():
        if not dirs:
            rep("WARN", f"{label}: empty"); continue
        bad, nodb, coupled, outs, nofloor = [], [], 0, 0, 0
        for d in dirs:
            p = os.path.join(d, "pflotran_co2.in")
            if not os.path.isfile(p):
                bad.append(f"{os.path.relpath(d, REV)}: no deck"); continue
            text = open(p).read(); f = facts(text); why = []
            if not f["richards"]: why.append("not RICHARDS")
            if len(f["gravity"]) != 1: why.append(f"{len(f['gravity'])} GRAVITY lines")
            if abs(f["vf"] - 0.50) > 1e-3: why.append(f"VF sum {f['vf']:.4f}")
            if not f["scaled"]: why.append("no SCALED_MASS_RATE")
            if not os.path.isfile(os.path.join(d, "full_mesh.uge")): why.append("no mesh")
            coupled += f["coupled"]
            nofloor += f["coupled"] and not re.search(r"^[ \t]*MINIMUM_POROSITY", text, re.M | re.I)
            outs += bool(glob.glob(os.path.join(d, "pflotran_co2*.h5")))
            for raw in fg.DB.findall(text):
                name = raw.strip("'\"")
                if not os.path.isabs(name) and not os.path.isfile(os.path.join(d, name)) and fg.db_sources(name) is None:
                    nodb.append(f"{os.path.relpath(d, REV)}: {name}")
            if why:
                bad.append(f"{os.path.relpath(d, REV)}: {', '.join(why)}")
        rep("PASS" if not bad else "FAIL", f"{label}: {len(dirs) - len(bad)} of {len(dirs)} ready "
            f"({coupled} coupled, {len(dirs) - coupled} fixed porosity)", "; ".join(bad[:6]))
        if nodb:
            rep("FAIL", f"{label}: {len(nodb)} databases missing with no reference copy", "; ".join(nodb[:4]))
        if outs:
            rep("WARN", f"{label}: {outs} directories still hold old outputs; archive them first (run_dirs.sh deletes them)")
        if nofloor:
            rep("WARN", f"{label}: {nofloor} coupled decks have no MINIMUM_POROSITY (expected only for runs that "
                        "started before the floor was introduced)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="p32_100_s1181")
    ap.add_argument("--skip", default="", help="indices left out on purpose, e.g. A:43,C:5")
    ap.add_argument("--machine", choices=("setonix", "hpc01"), default=None,
                    help="default: detected from PAWSEY_CLUSTER / the hostname")
    ap.add_argument("--blocks", default="all",
                    help="only these, e.g. A or A,B; add 'pre' for the pre-staged directories")
    a = ap.parse_args()
    skip = [(t.split(":")[0], int(t.split(":")[1])) for t in a.skip.split(",") if ":" in t]
    machine = a.machine or detect_machine()
    print(f"revision tree: {REV}   machine: {machine}")
    for f in ("fix_gravity_off.py", "apply_corrections.py", "blocks.py"):
        if not os.path.isfile(os.path.join(SRC, f)):
            print(f"FAIL  src/{f} missing"); sys.exit(1)
    check_launchers(machine); check_databases(a.case); check_staging(a.case)
    sel = ["A", "B", "C", "D30", "D40", "pre"] if a.blocks == "all" else [b.strip() for b in a.blocks.split(",")]
    check_blocks(skip, [b for b in ("A", "B", "C", "D30", "D40") if b in sel])
    if "pre" in sel:
        check_prestaged()
    n = {s: OUT.count(s) for s in ("PASS", "WARN", "FAIL")}
    print(f"\nsummary: {n['PASS']} pass, {n['WARN']} warn, {n['FAIL']} fail")
    sys.exit(1 if n["FAIL"] else 0)


if __name__ == "__main__":
    main()
