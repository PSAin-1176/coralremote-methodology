# Installation

This repository is code only. Data, local outputs, and machine-local configuration stay outside Git.

## Tested Environment

- MATLAB `R2025a`
- Windows

## Required MATLAB Products

- MATLAB
- Mapping Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox

## Likely Needed Depending on Workflow

- Image Processing Toolbox

## Expected Repository Setup

1. Clone the repository.
2. Keep all PlanetScope RAW data and processed local datasets outside the repository.
3. Create `config/paths_local.m` from `config/paths_example.m`.
4. Point `paths_local.m` to your local data roots and output roots.

## Required Local Configuration

Create:

- `config/paths_local.m`

Based on:

- `config/paths_example.m`

Set at least:

- `REPO_ROOT`
- `DATA_ROOT`
- `RAW_DIR`
- `ACA_DIR`
- `OUT_ROOT`
- `PREP_SOURCE_ROOT`
- `MASTERLOG_DIR`

## Required External Inputs

The production pipeline expects local raster/vector support files outside Git, including data referenced through `ACA_DIR` and related configured paths.

Typical required support inputs include:

- reef polygon GeoJSON
- deep-water mask GeoJSON
- sand mask GeoJSON
- wave/deep mask GeoJSON
- road mask raster
- bathymetry raster

## MATLAB Notes

Main runner:

- `code/RunAll_LogRatio.m`

Shared runner:

- `code/RunAll_LogRatio_Mode.m`

The main runner comment header is the authoritative short runtime note for this repository.

## Optional AROSICS Dependency

If you use:

- `code/01_preprocessing/Coreg_AROSICS_runner.m`

you may need a local Python/Conda environment with AROSICS available.

This is an optional preprocessing dependency, not a Git dependency.

## First Run

Minimal first-run order:

1. edit `config/paths_local.m`
2. review the comment header in `code/RunAll_LogRatio.m`
3. run `code/RunAll_LogRatio.m`

## Output Behavior

The pipeline writes outputs to the configured local output root.
Generated outputs should remain outside Git.
