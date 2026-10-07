# Revision run set driver. Run from the deposit root's revision/ directory.
SHELL := /bin/bash
.PHONY: help build check version bench audit validate dfn-a dfn-b dfn-d variants status verify test clean

help:
	@echo "Revision run set targets"
	@echo "  make build       - print the package build stamp"
	@echo "  make test        - self-tests for the deck override engine"
	@echo "  make version     - dfnWorks version + seed-reproducibility audit"
	@echo "  make check       - provenance snapshot of this environment"
	@echo "  make bench       - throughput calibration (4/8/16/24 ranks)"
	@echo "  make dfn-a       - generate Block A DFNs (expansion, 50)"
	@echo "  make dfn-b       - generate Block B DFNs (held-out, 20)"
	@echo "  make dfn-d       - generate Block D DFNs (domain size, 16)"
	@echo "  make audit       - closed-system / percolation audit (no compute)"
	@echo "  make validate H5=<file> [CASE=<name>] - check a run finished + reproduces"
	@echo "  make variants    - list Block C variants"
	@echo "  make status      - registry summary"
	@echo "  make verify      - mark runs done/failed by INSPECTING OUTPUT"

build:     ; @cat VERSION
test:      ; python3 tests/test_deckmod.py
version:   ; sbatch slurm/version_check.sh
check:     ; python3 src/provenance.py
bench:     ; for n in 4 8 16 24; do sbatch --ntasks=$$n --export=ALL,CASE=p32_100_s383,RANKS=$$n slurm/benchmark.sh; done
dfn-a:     ; sbatch --array=0-49%4 --export=ALL,BLOCK=A slurm/gen_dfn.sh
dfn-b:     ; sbatch --array=0-19%4 --export=ALL,BLOCK=B slurm/gen_dfn.sh
dfn-d:     ; sbatch --array=0-15%2 --export=ALL,BLOCK=D slurm/gen_dfn.sh
variants:  ; python3 src/build_variant.py --list
audit:     ; python3 src/audit.py --csv audit_boundaries.csv
validate:  ; @test -n "$(H5)" || { echo "usage: make validate H5=path/to/x.h5 [CASE=name]"; exit 1; }; \
	     python3 src/validate_run.py "$(H5)" $(if $(CASE),--case $(CASE),)
status:    ; python3 src/registry.py
verify:    ; python3 src/registry.py verify
clean:     ; rm -rf runs/*/carbfix.in runs/*/*.h5
