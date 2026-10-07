# CO₂ Mineral Trapping in Fractured Basalt: Discrete Fracture Network and Reactive Transport Code

Source code for:

> Chen, Y., Xie, Q., Kang, Q., & Regenauer-Lieb, K. (2026). Fracture network and geochemical controls on CO₂ mineral trapping in basalt: A stochastic reactive transport study. *Water Resources Research*. DOI: [pending]

This repository contains the source code used to generate the simulations and analyses presented in the article. The simulation outputs, DFN meshes, and a copy of the source code are archived on Zenodo: [10.5281/zenodo.20047873](https://doi.org/10.5281/zenodo.20047873).

## Repository Branches: Revised and Original Code

| Branch | Content |
|---|---|
| `Revision1` (default) | Code of the revised manuscript (revision R1), release `v2.0-R1` |
| `original` | Code of the original submission |

## Study Overview: Simulation Design and Principal Findings

The pipeline couples stochastic discrete fracture network (DFN) generation with reactive transport modeling to investigate CO₂ mineral trapping in fractured basalt. The simulations use three-dimensional DFN realizations under CarbFix-like conditions: 50 °C, 5 MPa, and continuous injection of carbonated water at pH 3.4 for 50 years. The model includes basalt mineral kinetics and evolving fracture apertures.

The simulation set contains 431 runs. It includes the intensity ensemble, independent ensemble, domain-size tests, shut-in simulations, injection-duration tests, and geochemical sensitivity tests.

The principal findings are:

- Fracture intensity increases total carbonate formation through its effect on pore volume, but mineralization efficiency shows no trend with fracture intensity.
- Carbonate forms where Ca/Mg complexes, carbonate species, and elevated local pH occur together. During injection, 3–6% of the fracture volume is carbonate-supersaturated and contains more than 90% of the new carbonate.
- Injection strategy controls mineralization efficiency. A 10-year injection followed by shut-in increases carbonate formation by a median factor of 6.2.

## Repository Contents: Scripts, Job Files, and Folders

| Path | Purpose |
|---|---|
| `prepare_dfn.py` | DFN generation and meshing with dfnWorks (Steps B1 and B2) |
| `run_pflotran.py` | Template PFLOTRAN decks for each network (Step B4); also used by `src/deckmod.py` and `src/build_variant.py` |
| `verify_matrix.py` | Quality control of the DFN library (Step B1) |
| `slurm/` | Slurm job scripts for DFN generation, deck staging, PFLOTRAN runs, and job status |
| `config/variants.py` | Definitions of geochemical sensitivity variants |
| `src/` | Scripts for deck construction, model corrections, run checks, analysis, statistics, and figures |
| `legacy/` | Analyses from the original submission that are not used in the revised article |
| `legacy/original/` | Scripts from the original pipeline that were replaced during revision |
| `requirements.txt` | Python packages required for the analysis |

This repository contains source code. Simulation inputs, reduced HDF5 outputs, and DFN meshes are archived on Zenodo in `.tar.zst` archives with a total size of approximately 15 GB.

## Quick Start: Post-Processing the Archived Zenodo Outputs

The post-processing workflow uses the archived simulation outputs. PFLOTRAN is not required for this step.

```bash
git clone -b Revision1 https://github.com/Yongqiang100/co2-basalt-dfn-ReactiveTransport.git revision
cd revision
pip install -r requirements.txt
# download the archives from Zenodo into this folder, then unpack
for f in runs_gravityoff runs_gravityoff_setonix_overlap runs_verify_shutin runs_gravityoff_prefloor archive_hpc01 results; do
  zstd -d --long=31 -c $f.tar.zst | tar -xf -
done
```

Then run the steps in [Post-Processing](#post-processing-simulation-outputs-to-article-figures-tables-and-results).

## Full Reproduction: Running the Simulations

Full reproduction of the simulations requires the following external software:

- **dfnWorks** (v2.7, with LaGriT): <https://github.com/lanl/dfnWorks>
- **PFLOTRAN** (v6) built with PETSc, with the databases `hanford.dat` and `co2_sw.dat`: <https://www.pflotran.org>
- An MPI runtime and Slurm

dfnWorks and PFLOTRAN are not included in `requirements.txt`. Install or compile these packages using their official documentation.

### Expected Folder Layout

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

Clone this branch into `<root>/revision` and place `prepare_dfn.py` in `<root>`. Set the root directory before submitting jobs:

```bash
export DEPOSIT_ROOT=<root>
cd $DEPOSIT_ROOT/revision
```

All commands below are executed from `<root>/revision`.

### Cluster Settings for the Slurm Job Scripts

The job scripts were developed for Pawsey Setonix and the hpc01/hpc02 clusters. Modify the account, partition, module, and `#SBATCH` settings in `slurm/env.sh`, `slurm/env_pflotran.sh`, and the relevant job scripts for the target cluster. `slurm/run_dirs.sh` is configured for Setonix. Use `slurm/run_dirs_hpc01.sh` for hpc01.

## Simulation Pipeline: DFN Generation, Model Preparation, and PFLOTRAN Runs

The simulation pipeline has three stages:

1. DFN generation and meshing (Steps B1–B3)
2. PFLOTRAN deck preparation and correction (Steps B4–B5)
3. Simulation execution and run checks (Steps B6–B13)

### Step B1. Generate and Mesh the Original DFN Library

**Article results:** Networks from the original submission, which form the basis of the revised intensity ensemble (Section 2.2).

```bash
cd $DEPOSIT_ROOT
python3 prepare_dfn.py matrix        # 5 intensity levels x 5 seeds: dfnWorks generation + LaGriT meshing
python3 verify_matrix.py             # quality control of the DFN library
```

For each network, `prepare_dfn.py` generates the fracture network with dfnWorks using the fracture-family parameters in Table S4, meshes the network with LaGriT, and writes the following files to `dfn_library/p32_<level>_s<seed>/`:

| File | Content |
|---|---|
| `full_mesh.uge` | unstructured mesh: cell volumes and connections (read by PFLOTRAN) |
| `full_mesh.inp` | the same mesh in AVS format (visualization) |
| `boundary_*.ex` | cells on the domain faces (inflow and outflow boundaries) |
| `dfn_properties.h5` | fracture properties: apertures, permeabilities |

### Step B2. Generate and Mesh the Revision Networks

**Article results:** Intensity ensemble, independent ensemble, domain-size networks, Figure 2, and Tables S4 and S5.

```bash
cd $DEPOSIT_ROOT/revision
sbatch --array=0-49 --export=ALL,BLOCK=A slurm/gen_dfn.sh   # intensity ensemble: 5 levels x 10 seeds
sbatch --array=0-19 --export=ALL,BLOCK=B slurm/gen_dfn.sh   # independent ensemble: P32 x0.90 and x1.75, 10 seeds each
sbatch --array=0-15 --export=ALL,BLOCK=D slurm/gen_dfn.sh   # domain size: 30 and 40 m, 8 seeds each
```

`gen_dfn.sh` loads dfnWorks and runs `src/extend_matrix.py`, which calls `prepare_dfn.py` without changing its matrix logic. Block D requires the one-line patch `01_matrix_ext.patch` before `prepare_dfn.py` is run. The Slurm array limit can be reduced, for example `--array=0-49%4`, if LaGriT exceeds the available memory.

To list or test a block without generating the networks:

```bash
python3 src/extend_matrix.py --block A --dry
```

### Step B3. Check Percolation and Outflow Boundaries

**Article results:** Table S2.

```bash
python3 src/check_percolation.py --csv percolation.csv
```

The script reads the archived network metadata. A network without an outflow boundary (`boundary_right_e.ex` with fewer than two lines) is treated as a closed system. Network `p32_200_s941` does not percolate and is replaced by `p32_200_s117`. Array index 43 is skipped, and `--extra p32_200_s117` is used in the subsequent build steps.

### Step B4. Generate the Template PFLOTRAN Decks

**Article results:** Model setup in Section 2.3, Table 1, and Table S1.

```bash
cd $DEPOSIT_ROOT
python3 run_pflotran.py --dfn all --write_only      # one PFLOTRAN deck per network in pflotran_results/
```

Each template deck `pflotran_results/p32_<level>_s<seed>/pflotran_co2.in` defines Richards flow, GIRT reactive transport, the mineral assemblage and kinetic parameters in Table 1, formation water and injectate compositions in Table S1, temperature and pressure, the injection region, and the uniform-pressure outflow boundary.

### Step B5. Apply Model Corrections

**Article results:** Corrected model used for the revised manuscript (Section 2.3).

The launch scripts apply the following corrections before PFLOTRAN is started. The corrections can also be applied manually to a staged deck.

| Correction | Script | What it changes |
|---|---|---|
| Injection rate | `src/apply_corrections.py` | `RATE MASS_RATE` becomes `RATE SCALED_MASS_RATE VOLUME`, so the stated rate is the total over the injection region rather than a rate per cell |
| Porosity and mineral volume fractions | `src/apply_corrections.py` with `RESCALE_VF=1` | the primary-mineral volume fractions are rescaled once so that they sum to 1 − φ = 0.50 for a fracture porosity of 0.50 (previously 0.85) |
| Gravity | `src/fix_gravity_off.py` | `GRAVITY 0.d0 0.d0 0.d0`, so the uniform-pressure outflow boundary drives no circulation, and the database paths |
| Porosity-permeability coupling (evolving aperture only) | `src/add_coupling.py` | `UPDATE_POROSITY`, `UPDATE_PERMEABILITY`, permeability power 3, critical porosity 0.01, minimum scale factor 1e-6, and a porosity floor `MINIMUM_POROSITY 0.01` |
| Injection stop (shut-in and duration runs) | `src/add_shutin.py` | the injection rate drops to zero after the given time |
| Evolving surface areas (Step B12 only) | `src/add_surface_area.py` | `SURFACE_AREA_FUNCTION` for the primary and secondary minerals |

Check a corrected deck with:

```bash
python3 src/check_porosity_and_vf.py      # porosity setting and volume-fraction sum of the decks
grep -i "GRAVITY\|SCALED_MASS_RATE\|UPDATE_POROSITY\|MINIMUM_POROSITY" <run folder>/pflotran_co2.in
```

Two execution routes are used:

- **Fixed porosity:** `slurm/run_rt.sh` stages each run from `pflotran_results/`, applies the corrections, and starts PFLOTRAN.
- **Evolving aperture:** `slurm/build_dirs.sh` creates the run directory with the required corrections and porosity-permeability coupling. `slurm/run_dirs.sh` then starts the runs listed in the corresponding run list.

### Step B6. Run the Intensity Ensemble

**Article results:** Sections 3.1–3.4, Figures 3–8, Tables 2 and S6–S8, Figures S1–S4, and Key Points 1 and 2.

```bash
# fixed porosity (comparison, Table S8)
sbatch --array=0-49%6 --export=ALL,BLOCK=A,RESCALE_VF=1 slurm/run_rt.sh
# evolving aperture (primary model)
bash slurm/build_dirs.sh --skip 43 --extra p32_200_s117 --apply
ls -d runs_gravityoff/A_feedback__* | sed 's#^runs_gravityoff/##' > lists/A_feedback.txt
sbatch --array=0-49 --export=ALL,LIST=lists/A_feedback.txt slurm/run_dirs.sh
```

### Step B7. Run the Independent Ensemble

**Article results:** Independent test of the fracture-intensity results (Section 3.4).

```bash
sbatch --array=0-19 --export=ALL,BLOCK=B,RESCALE_VF=1 slurm/run_rt.sh
bash slurm/build_dirs.sh --block B --apply
ls -d runs_gravityoff/B_feedback__* | sed 's#^runs_gravityoff/##' > lists/B_feedback.txt
sbatch --array=0-19 --export=ALL,LIST=lists/B_feedback.txt slurm/run_dirs.sh
```

### Step B8. Run the Domain-Size Test

**Article results:** Domain-size analysis (Section 3.4).

```bash
sbatch --array=0-7 --ntasks=64  --mem=112G --export=ALL,BLOCK=D30,RESCALE_VF=1 slurm/run_rt.sh
sbatch --array=0-7 --ntasks=128 --mem=224G --export=ALL,BLOCK=D40,RESCALE_VF=1 slurm/run_rt.sh
bash slurm/build_dirs.sh --block D30 --apply
bash slurm/build_dirs.sh --block D40 --apply
sbatch --array=0-7 --ntasks=64  --mem=112G --export=ALL,LIST=lists/D30_feedback.txt slurm/run_dirs.sh
sbatch --array=0-7 --ntasks=128 --mem=224G --export=ALL,LIST=lists/D40_feedback.txt slurm/run_dirs.sh
```

### Step B9. Run the Shut-In Ensemble

**Article results:** Section 3.5, Key Point 3, and the shut-in results reported in the Abstract and Conclusions.

```bash
bash slurm/build_dirs.sh --name E --shutin 10 --skip 43 --extra p32_200_s117 --apply                   # evolving aperture
bash slurm/build_dirs.sh --name E --shutin 10 --variant fixed --skip 43 --extra p32_200_s117 --apply   # fixed porosity
ls -d runs_gravityoff/E_* | sed 's#^runs_gravityoff/##' > lists/E.txt
sbatch --array=0-$(($(wc -l < lists/E.txt)-1)) --export=ALL,LIST=lists/E.txt slurm/run_dirs.sh
```

**Verification of dissolved carbon at shut-in:**

```bash
python3 src/prepare_verify_shutin.py        # networks p32_075_s1597, p32_200_s1289, p32_100_s1063
ls -d runs_verify_shutin/* > lists/verify.txt
sbatch --array=0-2 --export=ALL,LIST=lists/verify.txt slurm/run_dirs.sh
```

### Step B10. Run the Injection-Duration Series

**Article results:** Table 3 and Sections 3.5 and 4.2.

The injection-duration series uses one network. Each duration is simulated with fixed and evolving apertures. The shut-in times are given in years: 10 days = 0.0273785, 30 days = 0.0821355, 45 days = 0.1232033, 2 years = 2, and 5 years = 5. The 1-day simulation ends during the start-up ramp and uses the time offset of 1e-4.

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

The 10-year and 50-year durations correspond to the E and A runs of Steps B9 and B6.

### Step B11. Run the Geochemical Sensitivity Tests

**Article results:** Section 3.6, including surface-area, anorthite surface-area, secondary-mineral, seed-volume-fraction, dawsonite-rate, dawsonite-removal, and analcime tests.

`python3 src/build_variant.py --list` lists the variants defined in `config/variants.py`. The fixed-porosity decks are built and run first. Each evolving-aperture deck is then created from the corresponding fixed-porosity deck:

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

The 4 cases of the `feedback` variant already include aperture coupling and do not require a second coupled deck. This produces 65 coupled sensitivity decks.

### Step B12. Run the Evolving Surface-Area Test

**Article results:** Evolving surface-area analysis in Section 3.6.

`S0` uses constant surface areas as the control case. `S1` allows surface areas to evolve with mineral volume. `$N` is the comma-separated list of the 10 intensity-ensemble networks.

```bash
bash slurm/build_dirs.sh --block A --name S0 --networks $N --apply
bash slurm/build_dirs.sh --block A --name S1 --networks $N --apply
python3 src/add_surface_area.py runs_gravityoff/S1_feedback__*/pflotran_co2.in
ls -d runs_gravityoff/S[01]_feedback__* | sed 's#^runs_gravityoff/##' > lists/S.txt
sbatch --array=0-19 --export=ALL,LIST=lists/S.txt slurm/run_dirs.sh
```

### Step B13. Check Simulation Completion and Model Conditions

```bash
bash slurm/status.sh                                          # job states
python3 src/collect_failures.py                               # failed tasks
python3 src/validate_run.py --final-year 50 <run folder>      # run reached 50 years
python3 src/check_block.py --prefix A_feedback__              # finished, gravity off, carbonate per cell
python3 src/min_porosity.py lists/A_feedback.txt              # lowest porosity reached
python3 src/make_restart.py <run folder>                      # continue a run that stopped early
```

## Post-Processing: Simulation Outputs to Article Figures, Tables, and Results

Run the following commands from `<root>/revision` with the simulation folders in `runs_gravityoff/`. The folders can be generated by the simulation steps above or extracted from the Zenodo archive. Each step writes CSV files to `results/`. Run the steps in the order given because later steps use CSV files generated by earlier steps.

### Step C1. Summarize Simulation Blocks

**Article results:** Table 2, Tables S6–S8, and the mineralization efficiencies reported in the Abstract and Conclusions.

```bash
python3 src/summarize_block.py --prefix A_feedback__ --name A_coupled --root runs_gravityoff
python3 src/summarize_block.py --prefix A_ --name A_fixed --root runs_gravityoff
python3 src/summarize_block.py --prefix B_feedback__ --name B_coupled --root runs_gravityoff
python3 src/summarize_block.py --prefix E_feedback__ --name E_coupled --root runs_gravityoff
python3 src/time_series.py --prefix A_feedback__ --name A_coupled --root runs_gravityoff
python3 src/time_series.py --prefix A_ --name A_fixed --root runs_gravityoff
```

### Step C2. Analyze pH and Acid-Front Arrival

**Article results:** Section 3.1.

```bash
python3 src/ph_statistics.py
```

### Step C3. Quantify Primary-Mineral Dissolution

**Article results:** Section 3.2 and Figure 4.

```bash
python3 src/dissolution_capture_analysis.py --root runs_gravityoff
```

### Step C4. Analyze Carbonate Formation and Chemical Co-Location

**Article results:** Section 3.3, Figure 6, Key Point 2, and Table S3.

```bash
python3 src/colocation_mapping.py                       # co-location, supersaturated volume
python3 src/colocation_mapping.py --stages              # full history
python3 src/colocation_mapping.py --cations freeion     # free-ion proxies
python3 src/local_cation_balance.py
python3 src/local_cation_balance.py --stages
python3 src/check_redissolution.py --all                # Table S3
```

### Step C5. Analyze Fracture Intensity, Domain Size, and Aperture Evolution

**Article results:** Section 3.4, Figure 7, and Key Point 1.

```bash
python3 src/intensity_trend_test.py
python3 src/revision_statistics.py --only volume_effect domain_size_variability clogging_cells porosity_limit_pairs
```

### Step C6. Analyze Continuous Injection, Shut-In, and Injection Duration

**Article results:** Section 3.5, Table 3, and Key Point 3.

```bash
python3 src/revision_statistics.py --only injection_duration shutin_porewater_carbon
python3 src/shutin_carbon_balance.py --root runs_gravityoff --prefix E_feedback__
python3 src/verify_shutin_dic.py
python3 src/shutin_by_intensity.py
```

### Step C7. Calculate Geochemical Sensitivity Ratios

**Article results:** Section 3.6.

```bash
python3 src/revision_statistics.py --only sensitivity_ratios anorthite_surface_area shrinking_surface_ratios silicate_regrowth
```

### Step C8. Compare Results Across Computing Systems

**Article results:** Reproducibility analysis in Section 2.

```bash
python3 src/revision_statistics.py --only platform_comparison
```

This analysis requires `runs_gravityoff_setonix_overlap/`, `runs_gravityoff_prefloor/`, and `archive_hpc01/results/` from the Zenodo archive.

### Step C9. Generate Article and Supporting Information Figures

**Article results:** Figures 1–8 and Figures S1–S4.

All figures of the article and the response letter are generated with one command, run from `<root>/revision`:

```bash
bash src/make_figures.sh               # all figures, written to figures_R2/
bash src/make_figures.sh --skip-3d     # skips the slower 3D figures
```

`make_figures.sh` reads the runs in `runs_gravityoff/` and the DFN library in `../dfn_library`. It selects the two networks of Figure 6, the evolving-aperture networks at P32 × 1.00 with the highest and lowest carbonate per cell at 50 years, and calls `generate_figures.py`, `fig_study_design.py`, `fig_cations.py`, `prepare_figure_inputs.py` and `ensemble_statistics_figures.py`. The figures are written to `figures_R2/` for review before they are copied into `figures/`.

Individual figures are generated with the following commands:

```bash
python3 src/fig_study_design.py --out figures/fig_study_design.pdf                                   # Figure 1
python3 src/generate_figures.py --results-dir runs_gravityoff --dfn-dir ../dfn_library --output-dir figures   # Figures 2 to 8, S1 to S4
python3 src/prepare_figure_inputs.py
python3 src/ensemble_statistics_figures.py --out figures                                             # figures of the response letter
```

`generate_figures.py --only <n>` regenerates an individual figure. Figure 2 requires the DFN library generated in Steps B1 and B2.

## Map of Article Results to Simulation and Post-Processing Steps

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

## Legacy Scripts from the Original Submission

The `legacy/` directory contains analyses from the original submission that are not used in the revised article. These include betweenness, connectivity, particle tracking, flow-path profiles, and earlier co-location analyses. The `legacy/original/` directory contains scripts from the original simulation pipeline that were replaced during revision. These scripts are retained for reference and use the folder structure of the original simulation archive. The complete original source code is available on the `original` branch.

## Citation of the Article and Zenodo Archive

Use of the source code or associated data should cite both the article and the Zenodo archive:

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

## External Software Dependencies and Their Licenses

| Software | License | Reference |
|---|---|---|
| [dfnWorks](https://github.com/lanl/dfnWorks) | BSD 3-Clause | Hyman et al. (2015) |
| [PFLOTRAN](https://www.pflotran.org) | BSD 3-Clause | Lichtner et al. (2015) |
| [NetworkX](https://networkx.org/) | BSD 3-Clause | Hagberg et al. (2008) |

## Funding

This work was supported by the Australian Research Council Centre of Excellence for Carbon Science and Innovation (CE230100032) and the Australian Research Council Discovery Early Career Researcher Award (DE250100674).

## Computing Resources and AI Assistance

Computational resources were provided by the Pawsey Supercomputing Research Centre's Setonix supercomputer (<https://doi.org/10.48569/18sb-8s43>).

Claude Opus 4.6 was used for code debugging and figure formatting. The authors reviewed and edited the resulting code and are responsible for the final implementation.

## Contact

Yongqiang Chen — <yongqiang.chen@curtin.edu.au>, Curtin University, Perth, WA, Australia
