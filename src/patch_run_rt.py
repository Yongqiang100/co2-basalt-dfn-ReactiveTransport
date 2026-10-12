#!/usr/bin/env python3
"""
Hook src/fix_gravity_off.py (option B) into every PFLOTRAN launcher, so no run
starts with gravity on against a uniform-pressure boundary.

  * a launcher that calls apply_corrections.py (run_rt.sh) gets the step
    directly after that call, cloned from it;
  * a launcher that only runs a pre-staged deck (run_F.sh, run_pulse.sh) gets
    it directly before the srun that starts PFLOTRAN, on the deck srun reads;
  * a launcher already carrying the earlier fix_boundary.py hook is converted
    in place, so none ends up with both.

fix_gravity_off.py leaves a gravity-off deck unchanged, so repeats are safe.

    cd <revision>
    python3 src/patch_run_rt.py slurm/run_rt.sh slurm/run_F.sh slurm/run_pulse.sh          # dry run
    python3 src/patch_run_rt.py slurm/run_rt.sh slurm/run_F.sh slurm/run_pulse.sh --apply  # write
With no script named it patches slurm/run_rt.sh. Each file keeps a .bak_dirichlet.
"""
import os, re, shutil, subprocess, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
NAME = "fix_gravity_off.py"
FIX = os.path.join(REV, "src", NAME)
NOTE = ["# Gravity off (24 Sep 2026): in Richards mode gravity drives no flow, and with it",
        "# on the uniform-pressure outflow face drives circulation. Every run, every block."]
OLD_NOTE = ["# Hydrostatic boundary (24 Sep 2026): one uniform Dirichlet pressure under",
            "# gravity drives circulation through the outflow face. Every run, every block."]
MSG, OLD_MSG = "gravity correction failed", "boundary correction failed"


def live(l):
    return not l.lstrip().startswith("#")


def statement(lines, a):
    """Last line of the shell statement starting at line a."""
    b, text = a, lines[a]
    while text.rstrip().endswith("\\") or text.count("{") > text.count("}"):
        b += 1
        text += "\n" + lines[b]
    return b


def convert_old(lines):
    new = []
    for l in lines:
        for o, n in zip(OLD_NOTE, NOTE):
            if l.strip() == o:
                l = l.replace(o, n)
        new.append(l.replace("fix_boundary.py", NAME).replace(OLD_MSG, MSG))
    changed = [n for o, n in zip(lines, new) if o != n]
    return new, ("in place of the earlier fix_boundary.py hook",
                 [o for o, n in zip(lines, new) if o != n], changed)


def plan(path):
    lines = open(path).read().split("\n")
    if any(NAME in l for l in lines):
        return None, f"already calls {NAME}; nothing to do"
    if any("fix_boundary.py" in l for l in lines):
        return convert_old(lines)
    ac = [k for k, l in enumerate(lines) if "apply_corrections.py" in l and live(l)]
    if len(ac) == 1:
        a = ac[0]; b = statement(lines, a)
        ind = lines[a][:len(lines[a]) - len(lines[a].lstrip())]
        clone = [l.replace("apply_corrections.py", NAME).replace("deck corrections failed", MSG)
                 for l in lines[a:b + 1]]
        block = [ind + n for n in NOTE] + clone
        return lines[:b + 1] + block + lines[b + 1:], ("after the apply_corrections.py call", lines[a:b + 1], block)
    if len(ac) > 1:
        return None, f"{len(ac)} apply_corrections.py calls; patch by hand"
    sr = [k for k, l in enumerate(lines) if live(l) and re.match(r"\s*srun\b", l) and "pflotran" in l]
    if len(sr) != 1:
        return None, f"no apply_corrections.py call and {len(sr)} pflotran srun lines; patch by hand"
    k = sr[0]
    m = re.search(r"-input_prefix\s+(\S+)", lines[k])
    if not m or not re.fullmatch(r"[\w./-]+", m.group(1)):
        return None, "cannot read the deck name from the srun line; patch by hand"
    deck = m.group(1) + ".in"
    ind = lines[k][:len(lines[k]) - len(lines[k].lstrip())]
    block = [ind + n for n in NOTE] + \
            [f'{ind}python3 "{FIX}" {deck} || {{ echo "FATAL: {MSG}"; exit 1; }}']
    return lines[:k] + block + lines[k:], ("before the srun that starts PFLOTRAN", [lines[k]], block)


def main():
    apply = "--apply" in sys.argv
    scripts = [a for a in sys.argv[1:] if a != "--apply"] or [os.path.join(REV, "slurm", "run_rt.sh")]
    if not os.path.isfile(FIX):
        sys.exit(f"{FIX} is missing; copy {NAME} into src/ first")
    todo = []
    for p in scripts:
        if not os.path.isfile(p):
            sys.exit(f"not found: {p}")
        new, info = plan(p)
        print(f"\n=== {p}")
        if new is None:
            print(f"  {info}"); continue
        where, anchor, block = info
        print(f"  change {where}:")
        print("\n".join("      " + l for l in anchor))
        print("\n".join("    + " + l for l in block))
        todo.append((p, new))
    if not apply:
        print("\nDry run only. Re-run with --apply to write."); return
    for p, new in todo:
        orig, tmp = p + ".bak_dirichlet", p + ".bak_tmp"
        keep = os.path.exists(orig)            # never overwrite the first backup: it is the original
        if not keep:
            shutil.copy(p, orig)
        shutil.copy(p, tmp)                    # restore point for this write only
        open(p, "w").write("\n".join(new))
        if subprocess.run(["bash", "-n", p]).returncode:
            shutil.copy(tmp, p); os.remove(tmp)
            sys.exit(f"bash -n failed on {p}; restored to its state before this run")
        os.remove(tmp)
        print(f"wrote {p}, syntax OK  (original kept in {os.path.basename(orig)}"
              f"{', from an earlier run' if keep else ''})")


if __name__ == "__main__":
    main()
