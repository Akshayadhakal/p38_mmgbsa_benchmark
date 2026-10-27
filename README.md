# p38 MM/GBSA Benchmark

A computational workflow for evaluating molecular docking scores and MM/GBSA binding-energy estimates against experimental binding-affinity data for p38 protein–ligand complexes.

## Project contents

- `scripts/` — Python and shell scripts for ligand preparation, docking, molecular dynamics, MM/GBSA calculations, statistical analysis, and figure generation.
- `config/` — GROMACS molecular-dynamics parameters and MM/GBSA input files.
- `data/` — Ligand information, docking scores, pose/RMSD measurements, MM/GBSA results, analysis tables, and figures.
- `docs/` — Analysis protocol and project documentation.
- `environment_*.yml` — Exported Conda environment specifications.

## Methods overview

The workflow includes the following stages, as applicable to each ligand and simulation:

1. Prepare protein and ligand structures.
2. Perform molecular docking and collect docking scores.
3. Compare docking scores with experimental affinity data.
4. Evaluate docked-pose and molecular-dynamics stability.
5. Calculate MM/GBSA binding-energy estimates using the provided input files.
6. Compare estimates across ligands and simulation replicas.
7. Generate summary tables and figures.

Consult the scripts and configuration files for the actual settings used. Do not assume that every script must be run for every analysis.

## Requirements

The workflow uses a combination of Python, Conda, AutoDock Vina, GROMACS, and an MM/GBSA analysis tool. Exact dependencies and compatible versions depend on the workflow stage.

Environment exports are provided in `environment_bench.yml`, `environment_md.yml`, and `environment_mmpbsa.yml`. Review them before creating environments, and verify that any external executables and required force-field files are installed separately.

## Data and provenance

The files in `data/` include derived analysis outputs. Before reusing or redistributing any underlying experimental measurements, protein structures, benchmark datasets, or third-party files, consult their original sources and applicable licenses.

**Before publication:** complete the data-provenance record with the original dataset/benchmark citation, structure identifiers, experimental-affinity source, any preprocessing performed, and the licenses or terms applicable to redistributed files.

## Reproducibility

The scripts, configuration files, and environment exports are provided to document the computational workflow. Reproducing the results may require original input structures, ligand files, external force fields, software executables, and additional simulation files that are not included in this repository.

Check all paths and software versions before running scripts. Some scripts may require adaptation to your local installation.

## License

The MIT License in `LICENSE` applies to original code in this repository, subject to the notices in that file. It does not automatically apply to third-party data, structures, software, force fields, or other materials. See `DATA_LICENSES.md` for data-specific terms.
