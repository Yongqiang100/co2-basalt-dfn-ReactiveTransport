#!/usr/bin/env python3
"""Mesh convergence sequence: vary the discretisation, hold the network fixed.

h is an input to the GENERATOR, not only to the mesher: FRAM discards fractures
whose intersections are shorter than h (computationalGeometry.cpp:661). A
sequence at h = 0.4, 0.2, 0.1 would vary discretisation and topology together.
This holds h and the seed fixed and varies max_resolution_factor, which is
applied at meshing only. dfn_gen() calls mesh_network() with no arguments, so the
stages are called separately here to pass the parameter.

VERIFICATION -- all four must hold or the sequence is invalid:
  1. md5(radii_Final.dat) identical    -- accepted fracture set unchanged
  2. cell counts differ                -- the discretisation actually changed
  3. total fracture VOLUME agrees ~1%  -- same surfaces x same apertures must
                                          give the same volume however subdivided
  4. mean aperture identical           -- the seeded draw is reproducible
"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, sys
import numpy as np

DEPOSIT = os.environ.get("DEPOSIT_ROOT", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..")))
PHI = 0.50

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""): h.update(blk)
    return h.hexdigest()

def uge_stats(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        vol = np.array([float(f.readline().split()[4]) for _ in range(n)])
    return n, float(vol.sum()), float(vol.min()), float(vol.max())

def generate(case, factor, outroot, seed, p32_mult):
    sys.path.insert(0, DEPOSIT)
    import prepare_dfn as P
    from pydfnworks import DFNWORKS
    # ABSOLUTE path: create_network() does os.chdir(self.jobname), and
    # check_input() has already changed the working directory by that point,
    # so a relative jobname no longer resolves. dfn_gen() does not hit this
    # because it runs the stages in one uninterrupted sequence.
    jobname = os.path.abspath(os.path.join(outroot, f"mrf{factor:g}"))
    cfg = P.make_config(p32_mult)
    DFN = DFNWORKS(jobname=jobname, ncpu=1)
    DFN.params["domainSize"]["value"] = cfg["domain"]
    DFN.params["h"]["value"] = cfg["h"]
    DFN.params["orientationOption"]["value"] = None
    DFN.params["keepOnlyLargestCluster"]["value"] = True
    DFN.params["seed"]["value"] = seed
    for k in ("family1", "family2", "family3"):
        fam = cfg[k]
        DFN.add_fracture_family(
            shape=fam["shape"], distribution=fam["distribution"],
            kappa=fam["kappa"], orientation_distribution="fisher",
            theta=fam["theta"], phi=fam["phi"], log_mean=fam["log_mean"],
            log_std=fam["log_std"], min_radius=fam["min_radius"],
            max_radius=fam["max_radius"], number_of_points=8,
            p32=fam["p32"], aspect=fam["aspect"], beta_distribution=1, beta=0.0,
            hy_variable="aperture", hy_function="log-normal",
            hy_params={"mu": fam["aperture_mu"], "sigma": fam["aperture_sigma"]})
    here = os.getcwd()
    try:
        DFN.make_working_directory()
        DFN.check_input()
        DFN.create_network()
        DFN.mesh_network(max_resolution_factor=factor)
    finally:
        os.chdir(here)
    rng = np.random.default_rng(seed)
    fam_id = np.asarray(DFN.families, dtype=int)
    mu = np.array([cfg[f"family{f}"]["aperture_mu"] for f in fam_id])
    sd = np.array([cfg[f"family{f}"]["aperture_sigma"] for f in fam_id])
    b = np.exp(rng.normal(mu, sd))
    DFN.aperture, DFN.perm = b, b**2/12.0
    DFN.transmissivity = DFN.perm * b
    DFN.cell_based_aperture = False
    DFN.inp_file, DFN.uge_file, DFN.flow_solver = "full_mesh.inp", "full_mesh.uge", "PFLOTRAN"
    cwd = os.getcwd(); os.chdir(jobname)
    try:
        DFN.lagrit2pflotran()
        if not os.path.exists("full_mesh_vol_area.uge"):
            raise RuntimeError("correct_uge_file() produced no output")
        shutil.move("full_mesh.uge", "full_mesh_area_uncorrected.uge")
        shutil.copy("full_mesh_vol_area.uge", "full_mesh.uge")
        json.dump(dict(case=case, seed=seed, p32_mult=p32_mult, h=cfg["h"],
                       max_resolution_factor=factor, num_fractures=int(DFN.num_frac),
                       aperture_mean_m=float(b.mean()), uge_corrected=True),
                  open("mesh_conv.json", "w"), indent=2)
    finally:
        os.chdir(cwd)
    return jobname

def verify(outroot):
    dirs = sorted(d for d in os.listdir(outroot)
                  if d.startswith("mrf") and os.path.isdir(os.path.join(outroot, d)))
    if len(dirs) < 2:
        print(f"  only {len(dirs)} member(s)"); return False
    rows = []
    for d in dirs:
        p = os.path.join(outroot, d)
        radii = os.path.join(p, "dfnGen_output", "radii_Final.dat")
        uge = os.path.join(p, "full_mesh.uge")
        if not (os.path.isfile(radii) and os.path.isfile(uge)):
            print(f"  {d}: incomplete"); return False
        meta = json.load(open(os.path.join(p, "mesh_conv.json")))
        n, tot, vmin, vmax = uge_stats(uge)
        rows.append(dict(dir=d, factor=meta["max_resolution_factor"],
                         radii_md5=md5(radii), n_cells=n, volume=tot,
                         n_frac=meta["num_fractures"], aperture=meta["aperture_mean_m"]))
    rows.sort(key=lambda r: -r["factor"])
    print(f"\n{'dir':<8}{'factor':>8}{'fractures':>11}{'cells':>10}"
          f"{'volume m3':>12}{'aperture mm':>13}{'radii md5':>12}")
    for r in rows:
        print(f"{r['dir']:<8}{r['factor']:>8.4g}{r['n_frac']:>11}{r['n_cells']:>10,}"
              f"{r['volume']:>12.4f}{1000*r['aperture']:>13.4f}{r['radii_md5'][:10]:>12}")
    ok = True
    print("\n  CHECK 1  identical accepted fracture set")
    if len({r["radii_md5"] for r in rows}) == 1:
        print("    PASS  radii_Final.dat identical")
    else:
        print("    FAIL  radii_Final.dat differs -- topology varied, sequence invalid"); ok = False
    print("  CHECK 2  discretisation actually changed")
    if len({r["n_cells"] for r in rows}) == len(rows):
        print(f"    PASS  cell counts differ; refinement x{rows[-1]['n_cells']/rows[0]['n_cells']:.2f}")
    else:
        print("    FAIL  cell counts repeat -- parameter had no effect"); ok = False
    print("  CHECK 3  fracture volume invariant under remeshing")
    v = np.array([r["volume"] for r in rows])
    sp = 100*(v.max()-v.min())/v.mean()
    if sp < 2.0: print(f"    PASS  volume varies {sp:.2f}%")
    else:
        print(f"    FAIL  volume varies {sp:.2f}% -- check correct_uge_file() at each resolution"); ok = False
    print("  CHECK 4  apertures unchanged")
    a = np.array([r["aperture"] for r in rows])
    if np.allclose(a, a[0], rtol=1e-9): print(f"    PASS  mean aperture {1000*a[0]:.4f} mm")
    else: print("    FAIL  apertures differ"); ok = False
    print(f"\n  sequence is {'USABLE' if ok else 'NOT usable'}")
    if ok:
        print("  next: run at a FIXED rank count. Changing ranks between members")
        print("        alters the decomposition and introduces a second variable.")
    return ok

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="p32_150_s117")
    ap.add_argument("--factors", nargs="+", type=float, default=[10, 5, 2.5])
    ap.add_argument("--out", default="mesh_conv")
    ap.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if not a.verify_only:
        m = a.case.split("_")
        p32_mult, seed = int(m[1])/100, int(m[2][1:])
        print(f"case {a.case}: P32 x{p32_mult}, seed {seed}, factors {a.factors}\n")
        for f in a.factors:
            print(f"--- max_resolution_factor = {f:g}")
            generate(a.case, f, a.out, seed, p32_mult)
    verify(a.out)

if __name__ == "__main__":
    main()
