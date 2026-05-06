# -*- coding: utf-8 -*-
"""
run_pflotran.py -- Run PFLOTRAN reactive transport on DFN meshes

CarbFix1-style carbonated brine injection into fractured basalt.
CO2 is dissolved in water before injection (no free gas phase).

Kinetics from Palandri & Kharaka (2004) USGS OFR 2004-1068.
  - Rate constants in mol/m^2-sec (SI)
  - Activation energies in kJ/mol

CarbFix1 reservoir conditions:
  T = 50 C, P = 5 MPa (Gislason et al., 2010; Aradottir et al., 2009)

Supports both Setonix (srun) and local WSL2/Linux (mpirun).

Usage:
    python run_pflotran.py --dfn p32_100_s42 --nprocs 4
    python run_pflotran.py --dfn all --nprocs 4
    python run_pflotran.py --dfn all --post_only
    python run_pflotran.py --list
    python run_pflotran.py --dfn all --write_only
    python run_pflotran.py --dfn all --write_only --total_years 100
"""

import os
import sys
import json
import glob
import shutil
import argparse
import numpy as np

# =============================================================
# AUTO-DETECT ENVIRONMENT: HPC (Setonix) vs Local (WSL2/Linux)
# =============================================================
def _detect_environment():
    on_setonix = (
        os.path.exists("/software/projects/pawsey1284") or
        "PAWSEY_PROJECT" in os.environ or
        "SLURM_JOB_ID" in os.environ or
        os.environ.get("HOSTNAME", "").startswith("setonix") or
        os.environ.get("HOSTNAME", "").startswith("nid")
    )
    if on_setonix:
        _soft = os.environ.get("PFLOTRAN_SOFT",
            "/software/projects/pawsey1284/ychen6/pflotran")
        return {
            "platform": "setonix",
            "launcher": "srun",
            "pflotran_exe": os.environ.get("PFLOTRAN_EXE",
                f"{_soft}/pflotran/src/pflotran/pflotran"),
            "pflotran_db": os.environ.get("PFLOTRAN_DB",
                f"{_soft}/pflotran/database/hanford.dat"),
            "co2_db": os.environ.get("CO2_DB",
                f"{_soft}/pflotran/database/co2_sw.dat"),
        }
    else:
        home = os.path.expanduser("~")
        search_paths = [
            f"{home}/pflotran/pflotran-src",
            f"{home}/pflotran/pflotran",
            f"{home}/pflotran",
            "/usr/local",
            "/opt/pflotran",
        ]
        exe = os.environ.get("PFLOTRAN_EXE", "")
        if not exe or not os.path.isfile(exe):
            for base in search_paths:
                candidate = os.path.join(base, "src", "pflotran", "pflotran")
                if os.path.isfile(candidate):
                    exe = candidate
                    break
            if not exe:
                exe = shutil.which("pflotran") or "pflotran"
        db = os.environ.get("PFLOTRAN_DB", "")
        if not db or not os.path.isfile(db):
            for base in search_paths:
                candidate = os.path.join(base, "database", "hanford.dat")
                if os.path.isfile(candidate):
                    db = candidate
                    break
        co2 = os.environ.get("CO2_DB", "")
        if not co2 or not os.path.isfile(co2):
            for base in search_paths:
                candidate = os.path.join(base, "database", "co2_sw.dat")
                if os.path.isfile(candidate):
                    co2 = candidate
                    break
        return {
            "platform": "local",
            "launcher": "mpirun",
            "pflotran_exe": exe,
            "pflotran_db": db,
            "co2_db": co2,
        }


ENV = _detect_environment()
PFLOTRAN_EXE = ENV["pflotran_exe"]
PFLOTRAN_DB = ENV["pflotran_db"]
CO2_DB = ENV["co2_db"]

DFN_ROOT = os.path.join(os.getcwd(), "dfn_library")
RESULTS_ROOT = os.path.join(os.getcwd(), "pflotran_results")
SEC_PER_YEAR = 3.156e7


# =============================================================
# HELPER: Parse PFLOTRAN HDF5 time groups
# =============================================================
def parse_time_groups(h5f):
    """Find and sort time groups in PFLOTRAN HDF5.
    Handles format: '   0 Time  0.00000E+00 y'
    Returns list of (group_key, time_in_years)."""
    results = []
    for k in h5f.keys():
        if "Time" not in k:
            continue
        if "failure" in k.lower() or "cut" in k.lower():
            continue
        parts = k.strip().split()
        try:
            tidx = parts.index("Time")
            t_val = float(parts[tidx + 1])
            if len(parts) > tidx + 2 and parts[tidx + 2] == "y":
                t_yr = t_val
            else:
                t_yr = t_val / SEC_PER_YEAR
            results.append((k, t_yr))
        except (ValueError, IndexError):
            continue
    results.sort(key=lambda x: x[1])
    return results


def find_h5_var(grp, short_name):
    """Find HDF5 variable key starting with short_name.
    e.g. find_h5_var(grp, 'Calcite VF') matches
    'Calcite VF [m^3 mnrl_m^3 bulk]'."""
    for k in grp.keys():
        if k.startswith(short_name):
            return k
    return None


# =============================================================
# MESH UTILITIES
# =============================================================
def get_mesh_bounds(uge_path):
    xs, ys, zs = [], [], []
    with open(uge_path) as f:
        header = f.readline().strip().split()
        n_cells = int(header[1])
        for _ in range(n_cells):
            parts = f.readline().strip().split()
            xs.append(float(parts[1]))
            ys.append(float(parts[2]))
            zs.append(float(parts[3]))
    return {
        "x": [min(xs), max(xs)],
        "y": [min(ys), max(ys)],
        "z": [min(zs), max(zs)],
        "n_cells": n_cells,
    }


def estimate_domain_volume(bounds):
    dx = bounds["x"][1] - bounds["x"][0]
    dy = bounds["y"][1] - bounds["y"][0]
    dz = bounds["z"][1] - bounds["z"][0]
    return dx * dy * dz


def scale_injection_rate(base_rate, bounds, ref_volume=8000.0):
    vol = estimate_domain_volume(bounds)
    scale = vol / ref_volume
    return base_rate * max(scale, 0.01)


# =============================================================
# WRITE PFLOTRAN INPUT
# =============================================================
def write_input(dfn_dir, run_dir, params):
    os.makedirs(run_dir, exist_ok=True)

    for src_name in ["full_mesh.uge", "full_mesh.inp"]:
        src = os.path.join(dfn_dir, src_name)
        dst = os.path.join(run_dir, src_name)
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)

    for db in [PFLOTRAN_DB, CO2_DB]:
        dst = os.path.join(run_dir, os.path.basename(db))
        if not os.path.exists(dst):
            shutil.copy2(db, dst)

    for ex_file in glob.glob(os.path.join(dfn_dir, "*.ex")):
        dst = os.path.join(run_dir, os.path.basename(ex_file))
        if not os.path.exists(dst):
            shutil.copy2(ex_file, dst)

    bounds = get_mesh_bounds(os.path.join(dfn_dir, "full_mesh.uge"))
    margin = 5.0
    xmin = bounds["x"][0] - margin
    xmax = bounds["x"][1] + margin
    ymin = bounds["y"][0] - margin
    ymax = bounds["y"][1] + margin
    zmin = bounds["z"][0] - margin
    zmax = bounds["z"][1] + margin
    x_inj = xmin + 0.2 * (xmax - xmin)

    T = params["temperature"]
    P = params["pressure"]
    water_rate_base = params["water_rate_kg_s"]
    inject_years = params["inject_years"]
    total_years = params["total_years"]
    co2_molal = params["co2_molal"]

    water_rate = scale_injection_rate(water_rate_base, bounds)
    vol = estimate_domain_volume(bounds)
    print(f"  Domain volume: {vol:.0f} m^3")
    print(f"  Water rate: {water_rate:.3f} kg/s (scaled from {water_rate_base:.3f})")
    print(f"  CO2(aq) in injectate: {co2_molal} mol/kgw")

    has_right_ex = os.path.exists(os.path.join(run_dir, "boundary_right_e.ex"))
    right_size = 0
    if has_right_ex:
        with open(os.path.join(run_dir, "boundary_right_e.ex")) as f:
            right_size = len(f.readlines())

    bc_region = ""
    bc_coupler = ""
    if has_right_ex and right_size > 1:
        bc_region = """
REGION outflow_east
  FILE boundary_right_e.ex
END
"""
        bc_coupler = """
BOUNDARY_CONDITION outflow
  FLOW_CONDITION initial
  TRANSPORT_CONDITION initial_brine
  REGION outflow_east
END
"""
        print(f"  Boundary: outflow on right_e ({right_size} cells)")
    else:
        print(f"  Boundary: closed system")

    out_times = []
    for t in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0,
              10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0]:
        if t <= total_years:
            out_times.append(f"{t}")
    out_times_str = " ".join(out_times)

    water_rate_10pct = water_rate * 0.1

    input_text = f"""
#=============================================================
# CO2 Reactive Transport in Fractured Basalt
# T={T}C | P={P/1e6:.1f} MPa
# Carbonated brine: {water_rate:.3f} kg/s, CO2(aq)={co2_molal} mol/kgw
# Injection for {inject_years} yr, total {total_years} yr
# Mesh: {bounds['n_cells']} cells
#=============================================================

SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE RICHARDS
    /
    SUBSURFACE_TRANSPORT transport
      MODE GIRT
    /
  /
END

SUBSURFACE

NUMERICAL_METHODS FLOW
  NEWTON_SOLVER
    ATOL 1.d-6
    RTOL 1.d-6
    STOL 1.d-30
    MAXIMUM_NUMBER_OF_ITERATIONS 25
  /
  LINEAR_SOLVER
    SOLVER GMRES
    PRECONDITIONER ASM
    ATOL 1.d-12
    RTOL 1.d-6
  /
  TIMESTEPPER
    TS_ACCELERATION 8
    TIMESTEP_MAXIMUM_GROWTH_FACTOR 2.0
    MAXIMUM_CONSECUTIVE_TS_CUTS 32
  /
END

NUMERICAL_METHODS TRANSPORT
  NEWTON_SOLVER
    ATOL 1.d-12
    RTOL 1.d-6
    STOL 1.d-30
    MAXIMUM_NUMBER_OF_ITERATIONS 250
  /
  LINEAR_SOLVER
    SOLVER GMRES
    PRECONDITIONER ASM
    ATOL 1.d-12
    RTOL 1.d-5
  /
  TIMESTEPPER
    TS_ACCELERATION 16
    TIMESTEP_MAXIMUM_GROWTH_FACTOR 2.0
    MAXIMUM_CONSECUTIVE_TS_CUTS 40
  /
END

GRID
  TYPE unstructured_explicit full_mesh.uge
  GRAVITY 0.d0 0.d0 -9.8068d0
END

CHEMISTRY
  PRIMARY_SPECIES
    Al+++
    Ca++
    Fe++
    K+
    Mg++
    Na+
    SiO2(aq)
    H+
    CO2(aq)
    Cl-
    O2(aq)
  /
  SECONDARY_SPECIES
    HCO3-
    OH-
    CO3--
    AlOH++
    Al(OH)2+
    AlO2-
    HAlO2(aq)
    CaCO3(aq)
    CaHCO3+
    MgCO3(aq)
    MgHCO3+
    HSiO3-
  /
  MINERALS
    Anorthite
    Albite
    Diopside
    Forsterite
    Fayalite
    Enstatite
    Calcite
    Magnesite
    Siderite
    Dawsonite
    Kaolinite
    Chalcedony
  /
  MINERAL_KINETICS
    Anorthite
      PREFACTOR
        RATE_CONSTANT 3.16e-4 mol/m^2-sec
        ACTIVATION_ENERGY 16.6d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 1.411d0
        /
      /
      PREFACTOR
        RATE_CONSTANT 7.59e-10 mol/m^2-sec
        ACTIVATION_ENERGY 17.8d0 kJ/mol
      /
    /
    Albite
      PREFACTOR
        RATE_CONSTANT 6.92e-11 mol/m^2-sec
        ACTIVATION_ENERGY 65.0d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 0.457d0
        /
      /
      PREFACTOR
        RATE_CONSTANT 2.75e-13 mol/m^2-sec
        ACTIVATION_ENERGY 69.8d0 kJ/mol
      /
    /
    Diopside
      PREFACTOR
        RATE_CONSTANT 4.37e-7 mol/m^2-sec
        ACTIVATION_ENERGY 96.1d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 0.71d0
        /
      /
      PREFACTOR
        RATE_CONSTANT 7.76e-12 mol/m^2-sec
        ACTIVATION_ENERGY 40.6d0 kJ/mol
      /
    /
    Forsterite
      PREFACTOR
        RATE_CONSTANT 1.41e-7 mol/m^2-sec
        ACTIVATION_ENERGY 67.2d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 0.47d0
        /
      /
      PREFACTOR
        RATE_CONSTANT 2.29e-11 mol/m^2-sec
        ACTIVATION_ENERGY 79.0d0 kJ/mol
      /
    /
    Fayalite
      PREFACTOR
        RATE_CONSTANT 1.58e-13 mol/m^2-sec
        ACTIVATION_ENERGY 94.4d0 kJ/mol
      /
    /
    Enstatite
      PREFACTOR
        RATE_CONSTANT 9.55e-14 mol/m^2-sec
        ACTIVATION_ENERGY 80.0d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 0.6d0
        /
      /
      PREFACTOR
        RATE_CONSTANT 1.91e-13 mol/m^2-sec
        ACTIVATION_ENERGY 80.0d0 kJ/mol
      /
    /
    Calcite
      PREFACTOR
        RATE_CONSTANT 1.55e-6 mol/m^2-sec
        ACTIVATION_ENERGY 23.5d0 kJ/mol
      /
      PREFACTOR
        RATE_CONSTANT 5.01e-1 mol/m^2-sec
        ACTIVATION_ENERGY 14.4d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 1.0d0
        /
      /
    /
    Magnesite
      PREFACTOR
        RATE_CONSTANT 4.57e-10 mol/m^2-sec
        ACTIVATION_ENERGY 23.5d0 kJ/mol
      /
      PREFACTOR
        RATE_CONSTANT 4.17e-7 mol/m^2-sec
        ACTIVATION_ENERGY 14.4d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 1.0d0
        /
      /
    /
    Siderite
      PREFACTOR
        RATE_CONSTANT 1.26e-9 mol/m^2-sec
        ACTIVATION_ENERGY 21.0d0 kJ/mol
      /
    /
    Dawsonite
      PREFACTOR
        RATE_CONSTANT 1.0e-7 mol/m^2-sec
        ACTIVATION_ENERGY 62.8d0 kJ/mol
      /
    /
    Kaolinite
      PREFACTOR
        RATE_CONSTANT 6.61e-14 mol/m^2-sec
        ACTIVATION_ENERGY 22.2d0 kJ/mol
      /
      PREFACTOR
        RATE_CONSTANT 4.90e-12 mol/m^2-sec
        ACTIVATION_ENERGY 65.9d0 kJ/mol
        PREFACTOR_SPECIES H+
          ALPHA 0.777d0
        /
      /
    /
    Chalcedony
      PREFACTOR
        RATE_CONSTANT 5.89e-13 mol/m^2-sec
        ACTIVATION_ENERGY 74.5d0 kJ/mol
      /
    /
  /
  DATABASE hanford.dat
  LOG_FORMULATION
  ACTIVITY_COEFFICIENTS TIMESTEP
  OUTPUT
    SECONDARY_SPECIES
    FREE_ION
    PH
    TOTAL
    MINERALS
  /
END

FLUID_PROPERTY
  DIFFUSION_COEFFICIENT 2.d-9 m^2/s
END

MATERIAL_PROPERTY fracture
  ID 1
  CHARACTERISTIC_CURVES cc_frac
  POROSITY 0.50d0
  TORTUOSITY 1.d0
  ROCK_DENSITY 2900.d0 kg/m^3
  THERMAL_CONDUCTIVITY_DRY 1.5d0 W/m-C
  THERMAL_CONDUCTIVITY_WET 2.0d0 W/m-C
  HEAT_CAPACITY 880.d0 J/kg-C
  PERMEABILITY
    PERM_ISO 8.33d-8
  /
END

CHARACTERISTIC_CURVES cc_frac
  DEFAULT
END

REGION all
  COORDINATES
    {xmin:.1f}d0 {ymin:.1f}d0 {zmin:.1f}d0
    {xmax:.1f}d0 {ymax:.1f}d0 {zmax:.1f}d0
  /
END

REGION injection_zone
  COORDINATES
    {xmin:.1f}d0 {ymin:.1f}d0 {zmin:.1f}d0
    {x_inj:.1f}d0 {ymax:.1f}d0 {zmax:.1f}d0
  /
END

{bc_region}

FLOW_CONDITION initial
  TYPE
    LIQUID_PRESSURE DIRICHLET
  /
  LIQUID_PRESSURE {P}
END

FLOW_CONDITION carbonated_water_injection
  TYPE
    RATE MASS_RATE
  /
  RATE LIST
    TIME_UNITS y
    DATA_UNITS kg/s
    0.d0       0.d0
    1.d-4      {water_rate_10pct:.6e}
    1.d-2      {water_rate:.6e}
  /
END

TRANSPORT_CONDITION initial_brine
  TYPE ZERO_GRADIENT
  CONSTRAINT_LIST
    0.d0 basalt_brine
  /
END

TRANSPORT_CONDITION carbonated_water
  TYPE DIRICHLET
  CONSTRAINT_LIST
    0.d0 co2_rich_water
  /
END

CONSTRAINT basalt_brine
  CONCENTRATIONS
    Al+++       1.d-9      T
    Ca++        2.7d-5     T
    Fe++        1.72d-5    T
    K+          6.01d-5    T
    Mg++        4.53d-6    T
    Na+         4.03d-3    T
    SiO2(aq)    1.28d-4    T
    H+          7.5d0      P
    CO2(aq)     1.d-4      T
    Cl-         5.d-1      Z
    O2(aq)      1.d-64     T
  /
  MINERALS
    Anorthite     0.30d0    10.d0  cm^2/cm^3
    Albite        0.20d0    100.d0 cm^2/cm^3
    Diopside      0.25d0    100.d0 cm^2/cm^3
    Forsterite    0.03d0    100.d0 cm^2/cm^3
    Fayalite      0.02d0    100.d0 cm^2/cm^3
    Enstatite     0.05d0    100.d0 cm^2/cm^3
    Calcite       1.d-6     1.d0   cm^2/cm^3
    Magnesite     1.d-6     1.d0   cm^2/cm^3
    Siderite      1.d-6     1.d0   cm^2/cm^3
    Dawsonite     1.d-6     1.d0   cm^2/cm^3
    Kaolinite     1.d-6     1.d0   cm^2/cm^3
    Chalcedony    1.d-6     1.d0   cm^2/cm^3
  /
END

CONSTRAINT co2_rich_water
  CONCENTRATIONS
    Al+++       1.d-12     T
    Ca++        1.d-6      T
    Fe++        1.d-9      T
    K+          1.d-6      T
    Mg++        1.d-6      T
    Na+         1.d-3      T
    SiO2(aq)    1.d-6      T
    H+          3.4d0      P
    CO2(aq)     {co2_molal}    T
    Cl-         1.d-3      Z
    O2(aq)      1.d-64     T
  /
END

INITIAL_CONDITION
  FLOW_CONDITION initial
  TRANSPORT_CONDITION initial_brine
  REGION all
END

SOURCE_SINK carbonated_water_injector
  FLOW_CONDITION carbonated_water_injection
  TRANSPORT_CONDITION carbonated_water
  REGION injection_zone
END

{bc_coupler}

STRATA
  REGION all
  MATERIAL fracture
END

OUTPUT
  SNAPSHOT_FILE
    TIMES y {out_times_str}
    FORMAT HDF5
    PRINT_COLUMN_IDS
  /
  MASS_BALANCE_FILE
    PERIODIC TIME 1.d-1 y
    TOTAL_MASS_REGIONS
      all
    /
  /
  VARIABLES
    LIQUID_PRESSURE
    LIQUID_SATURATION
    LIQUID_DENSITY
    PERMEABILITY
    MATERIAL_ID
  /
END

TIME
  FINAL_TIME {total_years:.1f}d0 y
  INITIAL_TIMESTEP_SIZE 1.d-7 y
  MAXIMUM_TIMESTEP_SIZE 1.d-1 y
END

END_SUBSURFACE
"""

    input_path = os.path.join(run_dir, "pflotran_co2.in")
    with open(input_path, "w") as f:
        f.write(input_text)

    params_out = dict(params)
    params_out["water_rate_kg_s_scaled"] = water_rate
    params_out["domain_volume_m3"] = vol
    with open(os.path.join(run_dir, "simulation_params.json"), "w") as f:
        json.dump(params_out, f, indent=2)

    print(f"  Input: {input_path}")
    print(f"  Cells: {bounds['n_cells']}")


# =============================================================
# EXECUTE PFLOTRAN
# =============================================================
def execute_pflotran(run_dir, nprocs):
    """Run PFLOTRAN using auto-detected launcher."""
    launcher = ENV["launcher"]
    cmd = (f"cd {run_dir} && "
           f"{launcher} -n {nprocs} {PFLOTRAN_EXE} "
           f"-pflotranin pflotran_co2.in")

    print(f"  PFLOTRAN running ({nprocs} procs via {launcher})...")
    ret = os.system(cmd)

    if ret == 0:
        print("  PFLOTRAN completed successfully.")
    else:
        print(f"  WARNING: PFLOTRAN exited with code {ret}")

    h5 = glob.glob(os.path.join(run_dir, "pflotran_co2*.h5"))
    mas = glob.glob(os.path.join(run_dir, "*-mas.dat"))
    print(f"  HDF5: {len(h5)}, Mass balance: {len(mas)}")
    return ret == 0


# =============================================================
# POSTPROCESS
# =============================================================
def postprocess(run_dir):
    """Generate summary plots from PFLOTRAN HDF5 output."""
    try:
        import h5py
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams['text.usetex'] = False
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"  Missing: {e}")
        return

    # Only read pflotran_co2*.h5, skip pflotran_clean.h5
    h5_files = sorted(glob.glob(os.path.join(run_dir, "pflotran_co2*.h5")))
    if not h5_files:
        print("  No HDF5 output")
        return

    print(f"  Post-processing: {os.path.basename(h5_files[0])}")

    with h5py.File(h5_files[0], "r") as f:
        time_groups = parse_time_groups(f)
        if not time_groups:
            print("  No time groups found")
            return

        first_grp = f[time_groups[0][0]]
        first_vars = sorted(first_grp.keys())
        print(f"  Variables ({len(first_vars)}): {first_vars[:8]}...")
        print(f"  Time steps: {len(time_groups)}")

        times_yr = []
        avg_pH = []
        sum_calc, sum_mag, sum_sid, sum_daw = [], [], [], []
        sum_anor, sum_alb, sum_diop = [], [], []
        sum_forst, sum_fay, sum_enst = [], [], []
        sum_kaol, sum_chal = [], []

        for tg_key, t_yr in time_groups:
            grp = f[tg_key]
            times_yr.append(t_yr)

            # pH
            ph_key = find_h5_var(grp, "pH")
            if ph_key:
                avg_pH.append(float(np.nanmean(grp[ph_key][:])))
            else:
                avg_pH.append(np.nan)

            # Mineral volume fractions
            mineral_map = {
                "Calcite VF": sum_calc, "Magnesite VF": sum_mag,
                "Siderite VF": sum_sid, "Dawsonite VF": sum_daw,
                "Anorthite VF": sum_anor, "Albite VF": sum_alb,
                "Diopside VF": sum_diop, "Forsterite VF": sum_forst,
                "Fayalite VF": sum_fay, "Enstatite VF": sum_enst,
                "Kaolinite VF": sum_kaol, "Chalcedony VF": sum_chal,
            }
            for short_name, store in mineral_map.items():
                var_key = find_h5_var(grp, short_name)
                if var_key:
                    store.append(float(np.nansum(grp[var_key][:])))
                else:
                    store.append(0.0)

    if not times_yr:
        return

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    axes[0, 0].plot(times_yr, avg_pH, "k-", lw=2)
    axes[0, 0].set_xlabel("Time [yr]"); axes[0, 0].set_ylabel("pH")
    axes[0, 0].set_title("pH evolution"); axes[0, 0].grid(alpha=0.3)

    axes[0, 1].set_visible(False)

    for data, lbl, clr in [
        (sum_calc, "Calcite", "blue"), (sum_mag, "Magnesite", "green"),
        (sum_sid, "Siderite", "red"), (sum_daw, "Dawsonite", "orange"),
    ]:
        if any(v > 1e-10 for v in data):
            axes[0, 2].plot(times_yr, data, color=clr, lw=2, label=lbl)
    axes[0, 2].set_xlabel("Time [yr]"); axes[0, 2].set_ylabel("Total VF")
    axes[0, 2].set_title("Carbonate precipitation")
    axes[0, 2].legend(fontsize=8); axes[0, 2].grid(alpha=0.3)

    for data, lbl, clr in [
        (sum_anor, "Anorthite", "darkblue"),
        (sum_alb, "Albite", "cornflowerblue"),
        (sum_diop, "Diopside", "teal"),
        (sum_forst, "Forsterite", "darkgreen"),
        (sum_fay, "Fayalite", "brown"),
        (sum_enst, "Enstatite", "navy"),
    ]:
        if data:
            axes[1, 0].plot(times_yr, data, color=clr, lw=2, label=lbl)
    axes[1, 0].set_xlabel("Time [yr]"); axes[1, 0].set_ylabel("Total VF")
    axes[1, 0].set_title("Primary mineral dissolution")
    axes[1, 0].legend(fontsize=7); axes[1, 0].grid(alpha=0.3)

    for data, lbl, clr in [
        (sum_kaol, "Kaolinite", "sienna"),
        (sum_chal, "Chalcedony", "purple"),
    ]:
        if any(v > 1e-10 for v in data):
            axes[1, 1].plot(times_yr, data, color=clr, lw=2, label=lbl)
    axes[1, 1].set_xlabel("Time [yr]"); axes[1, 1].set_ylabel("Total VF")
    axes[1, 1].set_title("Clay and silica")
    axes[1, 1].legend(fontsize=8); axes[1, 1].grid(alpha=0.3)

    c = np.array(sum_calc); m = np.array(sum_mag)
    s = np.array(sum_sid); d = np.array(sum_daw)
    total = c + m + s + d
    if np.any(total > 1e-10):
        axes[1, 2].fill_between(times_yr, 0, c,
            alpha=0.4, label="Calcite", color="blue")
        axes[1, 2].fill_between(times_yr, c, c+m,
            alpha=0.4, label="Magnesite", color="green")
        axes[1, 2].fill_between(times_yr, c+m, c+m+s,
            alpha=0.4, label="Siderite", color="red")
        axes[1, 2].fill_between(times_yr, c+m+s, total,
            alpha=0.4, label="Dawsonite", color="orange")
        axes[1, 2].legend(fontsize=8)
    axes[1, 2].set_xlabel("Time [yr]")
    axes[1, 2].set_ylabel("Cumulative VF")
    axes[1, 2].set_title("CO2 mineral trapping")
    axes[1, 2].grid(alpha=0.3)

    run_name = os.path.basename(run_dir)
    fig.suptitle(f"CarbFix1 Carbonated Brine Injection: {run_name}",
                 fontsize=14)
    fig.tight_layout()
    path = os.path.join(run_dir, "results_summary.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")

    if len(times_yr) >= 2:
        print(f"\n  Summary (t={times_yr[0]:.2f} -> t={times_yr[-1]:.1f} yr):")
        print(f"    pH:         {avg_pH[0]:.2f} -> {avg_pH[-1]:.2f}")
        print(f"    Calcite:    {sum_calc[0]:.2e} -> {sum_calc[-1]:.2e}")
        print(f"    Magnesite:  {sum_mag[0]:.2e} -> {sum_mag[-1]:.2e}")
        print(f"    Siderite:   {sum_sid[0]:.2e} -> {sum_sid[-1]:.2e}")
        print(f"    Anorthite:  {sum_anor[0]:.4f} -> {sum_anor[-1]:.4f}")
        print(f"    Forsterite: {sum_forst[0]:.4f} -> {sum_forst[-1]:.4f}")


# =============================================================
# MAIN
# =============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Run PFLOTRAN CO2 reactive transport on DFN meshes")
    parser.add_argument("--dfn", type=str, default="p32_100_s42")
    parser.add_argument("--nprocs", type=int, default=4)
    parser.add_argument("--temperature", type=float, default=50.0)
    parser.add_argument("--pressure", type=float, default=5.0e6)
    parser.add_argument("--water_rate", type=float, default=0.01)
    parser.add_argument("--co2_molal", type=float, default=0.82)
    parser.add_argument("--inject_years", type=float, default=1.0)
    parser.add_argument("--total_years", type=float, default=50.0)
    parser.add_argument("--post_only", action="store_true")
    parser.add_argument("--write_only", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    if args.list:
        print(f"\nDFNs in {DFN_ROOT}:")
        if os.path.exists(DFN_ROOT):
            for d in sorted(os.listdir(DFN_ROOT)):
                uge = os.path.join(DFN_ROOT, d, "full_mesh.uge")
                if os.path.exists(uge):
                    size = os.path.getsize(uge) / 1e6
                    n_ex = len(glob.glob(os.path.join(DFN_ROOT, d, "*.ex")))
                    print(f"  {d:<25s} {size:>6.1f} MB  {n_ex} .ex")
        return

    print(f"\n  Platform:    {ENV['platform']}")
    print(f"  Launcher:    {ENV['launcher']}")
    print(f"  PFLOTRAN:    {PFLOTRAN_EXE}")
    print(f"  Database:    {PFLOTRAN_DB}")
    print(f"  CO2 DB:      {CO2_DB}")

    missing = []
    if not os.path.isfile(PFLOTRAN_EXE):
        missing.append(f"  PFLOTRAN_EXE: {PFLOTRAN_EXE}")
    if not os.path.isfile(PFLOTRAN_DB):
        missing.append(f"  PFLOTRAN_DB:  {PFLOTRAN_DB}")
    if not os.path.isfile(CO2_DB):
        missing.append(f"  CO2_DB:       {CO2_DB}")
    if missing and not args.post_only:
        print("\n  ERROR: Missing files:")
        for m in missing:
            print(f"    {m}")
        sys.exit(1)

    params = {
        "temperature": args.temperature,
        "pressure": args.pressure,
        "water_rate_kg_s": args.water_rate,
        "co2_molal": args.co2_molal,
        "inject_years": args.inject_years,
        "total_years": args.total_years,
    }

    os.makedirs(RESULTS_ROOT, exist_ok=True)

    if args.dfn == "all":
        if not os.path.exists(DFN_ROOT):
            print(f"ERROR: DFN library not found: {DFN_ROOT}")
            sys.exit(1)
        dfn_names = sorted([
            d for d in os.listdir(DFN_ROOT)
            if os.path.exists(os.path.join(DFN_ROOT, d, "full_mesh.uge"))
        ])
    else:
        dfn_names = [args.dfn]

    print(f"  Found {len(dfn_names)} DFN(s) to process")

    for dfn_name in dfn_names:
        dfn_dir = os.path.join(DFN_ROOT, dfn_name)
        if not os.path.exists(os.path.join(dfn_dir, "full_mesh.uge")):
            print(f"ERROR: {dfn_dir}/full_mesh.uge not found")
            continue

        run_dir = os.path.join(RESULTS_ROOT, dfn_name)

        print(f"\n{'='*60}")
        print(f"DFN: {dfn_name}")
        print(f"{'='*60}")

        if args.post_only:
            postprocess(run_dir)
        elif args.write_only:
            write_input(dfn_dir, run_dir, params)
        else:
            write_input(dfn_dir, run_dir, params)
            execute_pflotran(run_dir, args.nprocs)
            postprocess(run_dir)

    print(f"\nDone. Results at: {RESULTS_ROOT}/")


if __name__ == "__main__":
    main()