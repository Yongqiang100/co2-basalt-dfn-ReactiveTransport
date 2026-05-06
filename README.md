# CO₂ Basalt DFN Reactive Transport

Source code accompanying:

> Chen, Y., Xie, Q., Kang, Q., & Regenauer-Lieb, K. (2026). Fracture network connectivity controls on CO₂ mineral trapping efficiency in basalt: A stochastic reactive transport study. *Water Resources Research*. DOI: [pending]

This repository contains the active development version of the analysis pipeline. The frozen archive corresponding to the published version of the paper, including simulation outputs and DFN meshes, is deposited on Zenodo: [10.5281/zenodo.20047873](https://doi.org/10.5281/zenodo.20047873).

## Overview

The pipeline couples stochastic discrete fracture network (DFN) generation with reactive transport modeling to investigate how fracture network topology controls CO₂ mineral trapping efficiency in basalt. Twenty-five three-dimensional DFN realizations are simulated across five fracture intensity levels under representative reservoir conditions (50°C, 5 MPa, continuous injection of carbonated brine for 50 years), using a full basalt mineralogy and experimentally derived kinetic rate laws. The principal finding is that at fixed fracture intensity, carbonate precipitation varies by four orders of magnitude, and that a flow-weighted topological metric — the flow-reactive co-location index — predicts the rank order of trapping outcomes (Spearman ρ = 0.865).

## Repository contents

| File | Purpose |
|---|---|
| `prepare_dfn.py` | DFN library generation via dfnWorks |
| `verify_matrix.py` | DFN library quality control |
| `run_pflotran.py` | PFLOTRAN input generation and simulation execution |
| `run_local.sh` | Sequential and parallel `mpirun` driver for the 25 PFLOTRAN runs |
| `compute_betweenness.py` | Computation of the flow-reactive co-location index (Table 8) |
| `stress_test_analysis.py` | Computation of backbone, dead-end, and finite-size metrics (Table 8) |
| `generate_figures.py` | Single entry point producing all nine manuscript figures |
| `generate_xdmf.py` | Optional export of PFLOTRAN HDF5 to ParaView XDMF |

This repository hosts source code only. Simulation inputs, HDF5 outputs, and DFN meshes (approximately 6–13 GB) are deposited on Zenodo.

## Reproducing the published figures

The fastest path to verify the analysis is to regenerate the figures from the HDF5 outputs supplied in the Zenodo deposit. PFLOTRAN itself is not required for this step.

```bash
git clone https://github.com/Yongqiang100/co2-basalt-dfn-ReactiveTransport.git
cd co2-basalt-dfn-ReactiveTransport
pip install -r requirements.txt

# Link to the simulation outputs in the Zenodo archive
ln -s /path/to/zenodo/archive/dfn_library dfn_library
mkdir pflotran_results
for f in /path/to/zenodo/archive/pflotran_outputs/*_pflotran_co2.h5; do
    name=$(basename "$f" _pflotran_co2.h5)
    mkdir -p "pflotran_results/$name"
    ln -s "$f" "pflotran_results/$name/pflotran_co2.h5"
done

python generate_figures.py
```

Output figures are written to `paper_figures/` in both PNG and PDF format.

## Reproducing the simulations

Reproducing the underlying simulations rather than the post-processing requires the following external software:

- **dfnWorks** (v2.10): https://github.com/lanl/dfnWorks
- **PFLOTRAN** built with PETSc: https://www.pflotran.org
- An MPI runtime (OpenMPI or MPICH)

The complete pipeline is then:

```bash
python prepare_dfn.py matrix          # construct 25 DFN meshes (~30 min)
python verify_matrix.py               # quality control
python run_pflotran.py --dfn all --write_only   # generate 25 PFLOTRAN inputs

# PFLOTRAN simulations via mpirun (~5 hr per realization at 16 ranks)
bash run_local.sh --nprocs 16                   # sequential, 16 MPI ranks per run
# or, if memory permits, run several realizations concurrently:
bash run_local.sh --parallel --nprocs 16 --max-parallel 4

python compute_betweenness.py         # Table 8 source data
python stress_test_analysis.py        # additional Table 8 metrics
python generate_figures.py            # all nine manuscript figures
```

`run_local.sh` invokes `mpirun` and auto-detects the executable; the `MPIRUN` environment variable or `--mpirun` flag can override this. The 25-realization ensemble requires approximately 125 CPU-hours in serial; concurrent execution via `--parallel` reduces wall-clock time at the cost of higher memory usage.

## Figure generation

`generate_figures.py` is the single entry point for all nine manuscript figures. It supersedes the earlier multi-script workflow.

```bash
python generate_figures.py                      # all nine figures
python generate_figures.py --only 2 4 8         # selected figures
python generate_figures.py --skip-3d            # omit Figs 1 and 6 (slow 3D)
python generate_figures.py --output-dir myfigs  # custom output directory
python generate_figures.py --help               # full options
```

## Citation

Use of this code or the associated data should cite both the paper and the Zenodo archive:

```bibtex
@article{Chen2026_Paper,
  author  = {Chen, Yongqiang and Xie, Quan and Kang, Qinjun
             and Regenauer-Lieb, Klaus},
  title   = {Fracture network connectivity controls on {CO_2} mineral
             trapping efficiency in basalt: A stochastic reactive
             transport study},
  journal = {Water Resources Research},
  year    = {2026}
}

@dataset{Chen2026_Code,
  author    = {Chen, Yongqiang and Xie, Quan and Kang, Qinjun
               and Regenauer-Lieb, Klaus},
  title     = {Code and data for: Fracture network connectivity controls
               on {CO_2} mineral trapping efficiency in basalt},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20047873}
}
```

## License

Apache License, Version 2.0 — see [LICENSE](LICENSE).

## External software dependencies

| Software | License | Reference |
|---|---|---|
| [dfnWorks](https://github.com/lanl/dfnWorks) | BSD 3-Clause | Hyman et al. (2015) |
| [PFLOTRAN](https://www.pflotran.org) | BSD 3-Clause | Lichtner et al. (2015) |
| [NetworkX](https://networkx.org/) | BSD 3-Clause | Hagberg et al. (2008) |

## Funding

This work was supported by the Australian Research Council Centre of Excellence for Carbon Science and Innovation (CE230100032), the Australian Research Council Discovery Early Career Researcher Award (DE250100674), and the Pawsey Supercomputing Research Centre (Setonix).

## Contact

Yongqiang Chen — yongqiang.chen@curtin.edu.au
Curtin University, Perth, WA, Australia
