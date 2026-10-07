# CO₂ Basalt DFN Reactive Transport

Source code accompanying:
> Chen, Y., Xie, Q., Kang, Q., & Regenauer-Lieb, K. (2026). Fracture network and geochemical controls on CO₂ mineral trapping in basalt: A stochastic reactive transport study. *Water Resources Research*. DOI: [pending]

This repository contains the code to reproduce the simulations and every result of the article. The archive of the simulation outputs, DFN meshes and a copy of this code is deposited on Zenodo: [10.5281/zenodo.20047873](https://doi.org/10.5281/zenodo.20047873).

## Branches

| Branch | Content |
|---|---|
| `Revision1` (default) | Code of the revised manuscript (revision R1), release `v2.0-R1` |
| `original` | Code of the original submission |

## Overview

The pipeline couples stochastic discrete fracture network (DFN) generation with reactive transport modeling to investigate the controls on CO₂ mineral trapping in fractured basalt. Seventy three-dimensional DFN realizations are simulated across seven fracture intensity levels under CarbFix-like conditions (50 °C, 5 MPa, continuous injection of carbonated water at pH 3.4 for 50 years), with a full basalt mineralogy, kinetic rate laws and evolving fracture apertures. The simulations comprise 431 runs: the intensity ensemble (50 networks), an independent ensemble (20 networks), domain-size tests, a shut-in ensemble, an injection-duration series and geochemical sensitivity tests.

The principal findings are:

- Fracture intensity increases the amount of carbonate through the pore volume, but mineralization efficiency shows no trend with fracture intensity.
- Carbonate forms where Ca/Mg complexes, carbonate species and a higher local pH coincide. During injection, 3 to 6% of the fracture volume is carbonate-supersaturated and contains more than 90% of the new carbonate.
- The injection scheme controls the mineralization efficiency. A shut-in after 10 years increases carbonate formation by a median factor of 6.2.

## Repository contents

| Path | Purpose |
|---|---|
| `prepare_dfn.py` | DFN generation and meshing via dfnWorks (Steps B1 and B2) |
| `run_pflotran.py` | Template PFLOTRAN decks per network (Step B4); also used by `src/deckmod.py` and `src/build_variant.py` |
| `verify_matrix.py` | Quality control of the DFN library (Step B1) |
| `slurm/` | Job scripts: DFN generation, deck staging, PFLOTRAN runs, job status |
| `config/variants.py` | Definitions of the geochemical sensitivity variants |
| `src/` | Deck building and corrections, run checks, analysis, statistics and figures |
| `legacy/` | Analyses of the original submission not used in the revised article |
| `legacy/original/` | Scripts of the original pipeline replaced by the revision (betweenness, XDMF visualization, original figure and run scripts) |
| `requirements.txt` | Python packages for the analysis |

This repository hosts source code only. Simulation inputs, reduced HDF5 outputs and DFN meshes (15 GB in `.tar.zst` archives) are deposited on Zenodo.

## Archived simulation results

The fastest path to verify the analysis is to run the post-processing on the outputs in the Zenodo deposit. PFLOTRAN is not required for this step.

```bash
git clone -b Revision1 https://github.com/Yongqiang100/co2-basalt-dfn-ReactiveTransport.git revision
cd revision
pip install -r requirements.txt
# download the archives from Zenodo into this folder, then unpack
for f in runs_gravityoff runs_gravityoff_setonix_overlap runs_verify_shutin runs_gravityoff_prefloor archive_hpc01 results; do
  zstd -d --long=31 -c $f.tar.zst | tar -xf -
done
```

Then run the steps of [Post-processing](#post-processing-from-outputs-to-the-article).

## Reproducing the simulations

Reproducing the simulations rather than the post-processing requires the following external software:

- **dfnWorks** (v2.7, with LaGriT): <https://github.com/lanl/dfnWorks>
- **PFLOTRAN** (v6) built with PETSc, with the databases `hanford.dat` and `co2_sw.dat`: <https://www.pflotran.org>
- An MPI runtime and Slurm

dfnWorks and PFLOTRAN are not included in `requirements.txt`. Install or compile them from their official websites following their user guides, which are updated continuously.

### Layout

```
<root>/
├── prepare_dfn.py           DFN generation and meshing (dfnWorks)
├── dfn_library/             networks and meshes, written by step B1
├── pflotran_results/        template decks per network, written by step B1
└── revision/
    ├── src/                 deck building, checks, analysis and figures
    ├── slurm/               job scripts
    ├── config/              sensitivity variants (variants.py)
    ├── lists/               run lists for the job arrays
    └── runs_gravityoff/     one folder per run
```

Clone this branch into `<root>/revision` and place `prepare_dfn.py` in `<root>`. Set the root before submitting jobs:

```bash
export DEPOSIT_ROOT=<root>
cd $DEPOSIT_ROOT/revision
```

All commands below run from `<root>/revision`.

### Cluster settings

The job scripts were written for Pawsey Setonix and the hpc01/hpc02 clusters. Adapt the account, partition and module lines in `slurm/env.sh` and `slurm/env_pflotran.sh`, and the `#SBATCH` lines of each job script, to your system. `slurm/run_dirs.sh` is the Setonix version of `slurm/run_dirs_hpc01.sh`; use the one that matches your cluster.

### Pipeline

The pipeline has three stages: build the fracture networks and their meshes (Steps B1 to B3), prepare and correct the PFLOTRAN decks (Steps B4 and B5), and run each simulation set (Steps B6 to B12). Step B13 checks the runs.

#### Step B1. Original DFN library (25 networks)

**Article results:** the networks of the original submission, which the intensity ensemble extends (Section 2.2).

```bash
cd $DEPOSIT_ROOT
python3 prepare_dfn.py matrix        # 5 intensity levels x 5 seeds: dfnWorks generation + LaGriT meshing
python3 verify_matrix.py             # quality control of the DFN library
```

For each network, `prepare_dfn.py` generates the fractures with dfnWorks from the fracture-family parameters (Table S4), meshes them with LaGriT, and writes `dfn_library/p32_<level>_s<seed>/`:

| File | Content |
|---|---|
| `full_mesh.uge` | unstructured mesh: cell volumes and connections (read by PFLOTRAN) |
| `full_mesh.inp` | the same mesh in AVS format (visualization) |
| `boundary_*.ex` | cells on the domain faces (inflow and outflow boundaries) |
| `dfn_properties.h5` | fracture properties: apertures, permeabilities |

#### Step B2. Networks of the revision (Blocks A, B and D)

**Article results:** the intensity ensemble (10 networks per level), the independent ensemble and the domain-size networks; Figure 2; Tables S4 and S5.

```bash
cd $DEPOSIT_ROOT/revision
sbatch --array=0-49 --export=ALL,BLOCK=A slurm/gen_dfn.sh   # intensity ensemble: 5 levels x 10 seeds
sbatch --array=0-19 --export=ALL,BLOCK=B slurm/gen_dfn.sh   # independent ensemble: P32 x0.90 and x1.75, 10 seeds each
sbatch --array=0-15 --export=ALL,BLOCK=D slurm/gen_dfn.sh   # domain size: 30 and 40 m, 8 seeds each
```

`gen_dfn.sh` loads dfnWorks and runs `src/extend_matrix.py`, which wraps `prepare_dfn.py` without changing its matrix logic, so the 25 original networks stay reproducible. Block D needs the one-line patch `01_matrix_ext.patch` to `prepare_dfn.py` applied first. Add `%4` to the array (for example `--array=0-49%4`) if LaGriT runs short of memory. To list or test a block without generating it:

```bash
python3 src/extend_matrix.py --block A --dry
```

#### Step B3. Percolation check of the networks

**Article results:** Table S2.

```bash
python3 src/check_percolation.py --csv percolation.csv
```

The check reads the archived metadata only. A network without an outflow boundary (`boundary_right_e.ex` with fewer than two lines) runs as a closed system. The network `p32_200_s941` does not percolate and is replaced by `p32_200_s117`: array index 43 is skipped and `--extra p32_200_s117` is added in the build steps below.

#### Step B4. Template decks

**Article results:** the model setup of Section 2.3 (Table 1, Table S1).

```bash
cd $DEPOSIT_ROOT
python3 run_pflotran.py --dfn all --write_only      # one PFLOTRAN deck per network in pflotran_results/
```

Each template deck `pflotran_results/p32_<level>_s<seed>/pflotran_co2.in` sets the Richards flow and GIRT reactive transport, the mineral assemblage and kinetics (Table 1), the formation water and the injectate (`CONSTRAINT basalt_brine` and `co2_rich_water`, Table S1), 50 °C and 5 MPa, the injection region (left 20% of the domain) and the uniform-pressure outflow boundary.

#### Step B5. Deck corrections (applied to every run)

**Article results:** the corrected model of the revision (Section 2.3); every result of the article.

The launchers apply these corrections to each deck before PFLOTRAN starts. They can also be run by hand on a staged deck:

| Correction | Script | What it changes |
|---|---|---|
| Injection rate | `src/apply_corrections.py` | `RATE MASS_RATE` becomes `RATE SCALED_MASS_RATE VOLUME`, so the stated rate is the total over the injection region rather than a rate per cell |
| Porosity and mineral volume fractions | `src/apply_corrections.py` with `RESCALE_VF=1` | the primary-mineral volume fractions are rescaled once so that they sum to 1 − φ = 0.50 for a fracture porosity of 0.50 (previously 0.85) |
| Gravity | `src/fix_gravity_off.py` | `GRAVITY 0.d0 0.d0 0.d0`, so the uniform-pressure outflow boundary drives no circulation, and the database paths |
| Porosity-permeability coupling (evolving aperture only) | `src/add_coupling.py` | `UPDATE_POROSITY`, `UPDATE_PERMEABILITY`, permeability power 3, critical porosity 0.01, minimum scale factor 1e-6, and a porosity floor `MINIMUM_POROSITY 0.01` |
| Injection stop (shut-in and duration runs) | `src/add_shutin.py` | the injection rate drops to zero after the given time |
| Evolving surface areas (Step B12 only) | `src/add_surface_area.py` | `SURFACE_AREA_FUNCTION` for the primary and secondary minerals |

Check a deck after the corrections:

```bash
python3 src/check_porosity_and_vf.py      # porosity setting and volume-fraction sum of the decks
grep -i "GRAVITY\|SCALED_MASS_RATE\|UPDATE_POROSITY\|MINIMUM_POROSITY" <run folder>/pflotran_co2.in
```

Two launch routes are used:

- **Fixed porosity:** `slurm/run_rt.sh` stages each run from `pflotran_results/`, applies the corrections and starts PFLOTRAN.
- **Evolving aperture:** `slurm/build_dirs.sh` builds each run folder with the corrections and the coupling (default `--variant feedback`), and `slurm/run_dirs.sh` starts the folders listed in a file.

#### Step B6. Intensity ensemble (50 networks)

**Article results:** Sections 3.1 to 3.4; Figures 3 to 8; Tables 2, S6, S7 and S8; Figures S1 to S4; Key Points 1 and 2.

```bash
# fixed porosity (comparison, Table S8)
sbatch --array=0-49%6 --export=ALL,BLOCK=A,RESCALE_VF=1 slurm/run_rt.sh
# evolving aperture (primary model)
bash slurm/build_dirs.sh --skip 43 --extra p32_200_s117 --apply
ls -d runs_gravityoff/A_feedback__* | sed 's#^runs_gravityoff/##' > lists/A_feedback.txt
sbatch --array=0-49 --export=ALL,LIST=lists/A_feedback.txt slurm/run_dirs.sh
```

#### Step B7. Independent ensemble (20 networks)

**Article results:** the independent test of the intensity trends (Section 3.4).

```bash
sbatch --array=0-19 --export=ALL,BLOCK=B,RESCALE_VF=1 slurm/run_rt.sh
bash slurm/build_dirs.sh --block B --apply
ls -d runs_gravityoff/B_feedback__* | sed 's#^runs_gravityoff/##' > lists/B_feedback.txt
sbatch --array=0-19 --export=ALL,LIST=lists/B_feedback.txt slurm/run_dirs.sh
```

#### Step B8. Domain size (16 networks)

**Article results:** the domain-size test (Section 3.4).

```bash
sbatch --array=0-7 --ntasks=64  --mem=112G --export=ALL,BLOCK=D30,RESCALE_VF=1 slurm/run_rt.sh
sbatch --array=0-7 --ntasks=128 --mem=224G --export=ALL,BLOCK=D40,RESCALE_VF=1 slurm/run_rt.sh
bash slurm/build_dirs.sh --block D30 --apply
bash slurm/build_dirs.sh --block D40 --apply
sbatch --array=0-7 --ntasks=64  --mem=112G --export=ALL,LIST=lists/D30_feedback.txt slurm/run_dirs.sh
sbatch --array=0-7 --ntasks=128 --mem=224G --export=ALL,LIST=lists/D40_feedback.txt slurm/run_dirs.sh
```

#### Step B9. Shut-in after 10 years (50 paired networks)

**Article results:** Section 3.5; Key Point 3; the shut-in results of the Abstract and Conclusions.

```bash
bash slurm/build_dirs.sh --name E --shutin 10 --skip 43 --extra p32_200_s117 --apply                   # evolving aperture
bash slurm/build_dirs.sh --name E --shutin 10 --variant fixed --skip 43 --extra p32_200_s117 --apply   # fixed porosity
ls -d runs_gravityoff/E_* | sed 's#^runs_gravityoff/##' > lists/E.txt
sbatch --array=0-$(($(wc -l < lists/E.txt)-1)) --export=ALL,LIST=lists/E.txt slurm/run_dirs.sh
```

**Verification of the dissolved carbon at shut-in (3 runs):**

```bash
python3 src/prepare_verify_shutin.py        # networks p32_075_s1597, p32_200_s1289, p32_100_s1063
ls -d runs_verify_shutin/* > lists/verify.txt
sbatch --array=0-2 --export=ALL,LIST=lists/verify.txt slurm/run_dirs.sh
```

#### Step B10. Injection duration (1 network, 6 durations)

**Article results:** Table 3; Sections 3.5 and 4.2.

Each duration is built with and without `--variant fixed`. The shut-in time is in years: 10 days = 0.0273785, 30 days = 0.0821355, 45 days = 0.1232033, 2 years = 2, 5 years = 5. The 1-day injection ends inside the start-up ramp and uses the older offset of 1e-4.

```bash
for d in "F_10d 0.0273785" "F_30d 0.0821355" "F_45d 0.1232033" "F_2yr 2" "F_5yr 5"; do
  set -- $d
  bash slurm/build_dirs.sh --name $1 --shutin $2 --networks p32_100_s1181 --apply
  bash slurm/build_dirs.sh --name $1 --shutin $2 --variant fixed --networks p32_100_s1181 --apply
done
bash slurm/build_dirs.sh --name F_1d --shutin 2.74e-3 --offset 1e-4 --networks p32_100_s1181 --apply
bash slurm/build_dirs.sh --name F_1d --shutin 2.74e-3 --offset 1e-4 --variant fixed --networks p32_100_s1181 --apply
ls -d runs_gravityoff/F_*_feedback__* | sed 's#^runs_gravityoff/##' > lists/F_all_feedback.txt
ls -d runs_gravityoff/F_*__* | grep -v _feedback__ | sed 's#^runs_gravityoff/##' > lists/F_all_fixed.txt
sbatch --array=0-5 --export=ALL,LIST=lists/F_all_feedback.txt slurm/run_dirs.sh
sbatch --array=0-5 --export=ALL,LIST=lists/F_all_fixed.txt slurm/run_dirs.sh
```

The 10-year and 50-year durations are the E and A runs of steps B5 and B2.

#### Step B11. Geochemical sensitivity (16 variants, 25 reference networks)

**Article results:** Section 3.6 (surface areas, anorthite surface area, secondary minerals, seed volume fraction, dawsonite rate, dawsonite removal, analcime).

`python3 src/build_variant.py --list` lists the variants (defined in `config/variants.py`). Build and run the fixed-porosity decks, then make each evolving-aperture twin from its fixed deck:

```bash
python3 src/build_variant.py --variant <variant> --index <i>     # for each variant and network
sbatch --array=0-68 --export=ALL,BLOCK=C,RESCALE_VF=1 slurm/run_rt.sh
ls -d runs_gravityoff/C_* | grep -v _feedback__ > lists/C_sensitivity.txt
: > lists/C_feedback.txt
for d in $(grep -v "C_feedback__" lists/C_sensitivity.txt); do
  n=${d##*__}; v=$(basename "${d%__*}"); t="runs_gravityoff/${v}_feedback__$n"
  mkdir -p "$t" && cp "$d"/pflotran_co2.in "$d"/*.uge "$d"/*.ex "$d"/*.dat "$t"/ 2>/dev/null
  python3 src/add_coupling.py "$t/pflotran_co2.in" > /dev/null && echo "$t" >> lists/C_feedback.txt
done
sbatch --array=0-$(($(wc -l < lists/C_feedback.txt)-1)) --export=ALL,LIST=lists/C_feedback.txt slurm/run_dirs.sh
```

The 4 cases of the variant `feedback` are already coupled and get no twin, giving 65 coupled decks.

#### Step B12. Evolving surface areas (10 networks)

**Article results:** the evolving-surface-area test (Section 3.6).

`S0` keeps the surface areas constant (controls); `S1` lets them evolve with the mineral volumes. `$N` is the comma-separated list of the 10 intensity-ensemble networks:

```bash
bash slurm/build_dirs.sh --block A --name S0 --networks $N --apply
bash slurm/build_dirs.sh --block A --name S1 --networks $N --apply
python3 src/add_surface_area.py runs_gravityoff/S1_feedback__*/pflotran_co2.in
ls -d runs_gravityoff/S[01]_feedback__* | sed 's#^runs_gravityoff/##' > lists/S.txt
sbatch --array=0-19 --export=ALL,LIST=lists/S.txt slurm/run_dirs.sh
```

#### Step B13. Check the runs

```bash
bash slurm/status.sh                                          # job states
python3 src/collect_failures.py                               # failed tasks
python3 src/validate_run.py --final-year 50 <run folder>      # run reached 50 years
python3 src/check_block.py --prefix A_feedback__              # finished, gravity off, carbonate per cell
python3 src/min_porosity.py lists/A_feedback.txt              # lowest porosity reached
python3 src/make_restart.py <run folder>                      # continue a run that stopped early
```

## Post-processing: from outputs to the article

Run from `<root>/revision`, with the run folders in `runs_gravityoff/` (from the simulations above or from the Zenodo archives). Each step writes CSV files to `results/`. Run the steps in this order, since later steps read the earlier CSV files.

### Step C1. Block summaries

**Article results:** Table 2; Tables S6, S7 and S8; the efficiencies of the Abstract and Conclusions.

```bash
python3 src/summarize_block.py --prefix A_feedback__ --name A_coupled --root runs_gravityoff
python3 src/summarize_block.py --prefix A_ --name A_fixed --root runs_gravityoff
python3 src/summarize_block.py --prefix B_feedback__ --name B_coupled --root runs_gravityoff
python3 src/summarize_block.py --prefix E_feedback__ --name E_coupled --root runs_gravityoff
python3 src/time_series.py --prefix A_feedback__ --name A_coupled --root runs_gravityoff
python3 src/time_series.py --prefix A_ --name A_fixed --root runs_gravityoff
```

### Step C2. pH and acid front

**Article results:** Section 3.1.

```bash
python3 src/ph_statistics.py
```

### Step C3. Dissolution

**Article results:** Section 3.2; Figure 4.

```bash
python3 src/dissolution_capture_analysis.py --root runs_gravityoff
```

### Step C4. Carbonate formation and localization

**Article results:** Section 3.3; Figure 6; Key Point 2; Table S3.

```bash
python3 src/colocation_mapping.py                       # co-location, supersaturated volume
python3 src/colocation_mapping.py --stages              # full history
python3 src/colocation_mapping.py --cations freeion     # free-ion proxies
python3 src/local_cation_balance.py
python3 src/local_cation_balance.py --stages
python3 src/check_redissolution.py --all                # Table S3
```

### Step C5. Fracture intensity, domain size and aperture

**Article results:** Section 3.4; Figure 7; Key Point 1.

```bash
python3 src/intensity_trend_test.py
python3 src/revision_statistics.py --only volume_effect domain_size_variability clogging_cells porosity_limit_pairs
```

### Step C6. Injection scheme

**Article results:** Section 3.5; Table 3; Key Point 3.

```bash
python3 src/revision_statistics.py --only injection_duration shutin_porewater_carbon
python3 src/shutin_carbon_balance.py --root runs_gravityoff --prefix E_feedback__
python3 src/verify_shutin_dic.py
python3 src/shutin_by_intensity.py
```

### Step C7. Geochemical sensitivity

**Article results:** Section 3.6.

```bash
python3 src/revision_statistics.py --only sensitivity_ratios anorthite_surface_area shrinking_surface_ratios silicate_regrowth
```

### Step C8. Platform checks

**Article results:** the reproducibility statement of Section 2.

```bash
python3 src/revision_statistics.py --only platform_comparison
```

This needs the folders `runs_gravityoff_setonix_overlap/`, `runs_gravityoff_prefloor/` and `archive_hpc01/results/` from the Zenodo record.

### Step C9. Figures

**Article results:** Figures 1 to 8; Figures S1 to S4.

```bash
python3 src/fig_study_design.py --out figures/fig_study_design.pdf                                   # Figure 1
python3 src/generate_figures.py --results-dir runs_gravityoff --dfn-dir ../dfn_library --output-dir figures   # Figures 2 to 8, S1 to S4
python3 src/prepare_figure_inputs.py
python3 src/ensemble_statistics_figures.py --out figures                                             # figures of the response letter
```

`generate_figures.py --only <n>` regenerates a single figure. Figure 2 needs the DFN library of Steps B1 and B2.

## Map of the article results

| Article item | Simulation step | Post-processing step | Output |
|---|---|---|---|
| Figure 1 | — | C9 `fig_study_design.py` | `figures/fig_study_design.pdf` |
| Figure 2, Tables S4, S5 | B1, B2 | C9 `generate_figures.py` | `figures/` |
| Table S2 | B3 | `check_percolation.py` | `runs_gravityoff/registry.json` |
| Section 3.1 | B6 | C2 `ph_statistics.py` | printed |
| Section 3.2, Figure 4 | B6 | C1, C3 | `A_coupled_by_intensity.csv` |
| Section 3.3, Figure 6, Key Point 2 | B6 | C4 | `colocation_mapping_*.csv`, `local_cation_balance*.csv` |
| Table S3 | B6 | C4 `check_redissolution.py` | printed |
| Table 2 | B6 | C1 | `A_coupled_by_intensity.csv`, `figure_input_blockA.csv` |
| Tables S6, S7 | B6 | C1 | `A_coupled_networks.csv` |
| Table S8 | B6 | C1 | `A_fixed_*.csv`, `A_coupled_*.csv` |
| Section 3.4, Figure 7, Key Point 1 | B6, B7, B8 | C1, C5 | `intensity_trend_tests.csv`, `volume_effect.csv`, `domain_size_variability.csv` |
| Section 3.5, Table 3, Key Point 3 | B6, B9, B10 | C1, C6 | `injection_duration.csv`, `shutin_*.csv`, `E_coupled_*.csv` |
| Section 3.6 | B11, B12 | C7 | `sensitivity_ratios.csv`, `anorthite_surface_area.csv`, `shrinking_surface_ratios.csv`, `silicate_regrowth.csv` |
| Figure 3, Figure 5, Figure 8, Figures S1 to S4 | B6 | C9 `generate_figures.py` | `figures/` |
| Platform checks (Section 2) | repeated runs | C8 | `platform_comparison.csv` |

## Folder `legacy/`

Scripts of analyses from the original submission that the revised article does not use (betweenness, connectivity, particle tracking, flow-path profiles, earlier co-location figures). `legacy/original/` holds the scripts of the original pipeline that the revision replaced. They are kept for reference and expect the folder layout of the original runs. The complete original code is on the `original` branch.

## Citation

Use of this code or the associated data should cite both the paper and the Zenodo archive:

```bibtex
@article{Chen2026_Paper,
  author  = {Chen, Yongqiang and Xie, Quan and Kang, Qinjun
             and Regenauer-Lieb, Klaus},
  title   = {Fracture network and geochemical controls on {CO_2} mineral
             trapping in basalt: A stochastic reactive transport study},
  journal = {Water Resources Research},
  year    = {2026}
}

@dataset{Chen2026_Code,
  author    = {Chen, Yongqiang and Xie, Quan and Kang, Qinjun
               and Regenauer-Lieb, Klaus},
  title     = {Code and data for: Fracture network and geochemical controls
               on {CO_2} mineral trapping in basalt},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20047873}
}
```

## License

Apache License, Version 2.0.

## External software dependencies

| Software | License | Reference |
|---|---|---|
| [dfnWorks](https://github.com/lanl/dfnWorks) | BSD 3-Clause | Hyman et al. (2015) |
| [PFLOTRAN](https://www.pflotran.org) | BSD 3-Clause | Lichtner et al. (2015) |
| [NetworkX](https://networkx.org/) | BSD 3-Clause | Hagberg et al. (2008) |

## Funding

This work was supported by the Australian Research Council Centre of Excellence for Carbon Science and Innovation (CE230100032) and the Australian Research Council Discovery Early Career Researcher Award (DE250100674).

## Acknowledgements

Computational resources were provided by the Pawsey Supercomputing Research Centre's Setonix supercomputer (<https://doi.org/10.48569/18sb-8s43>). Claude Opus 4.6 was used to assist code debugging and figure formatting; the authors reviewed and edited the code as needed and take full responsibility for it.

## Contact

Yongqiang Chen — <yongqiang.chen@curtin.edu.au>, Curtin University, Perth, WA, Australia
