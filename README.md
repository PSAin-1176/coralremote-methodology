# CoralRemote Methodology

This repository contains the code-only production methodology pipeline for coral reef bleaching analysis using a log-ratio workflow.

## What is included

- production MATLAB pipeline code
- shared utility functions
- lightweight repository documentation
- a 4-part code layout: preprocessing, baseline construction, classification, evaluation

## What is intentionally excluded

- PlanetScope RAW data
- processed datasets
- model files
- generated results
- local archives
- manuscript-output scripts
- machine-local configuration

## Entry Point

- main runner: `code/RunAll_LogRatio.m`
- shared runner: `code/RunAll_LogRatio_Mode.m`
- preprocessing entry point: `code/01_preprocessing/Preprocess_Run.m`

All required runtime notes are maintained in the comment header of `code/RunAll_LogRatio.m`.

## Quick Start

1. create `config/paths_local.m` from `config/paths_example.m`
2. point it to your local data and output roots
3. review `code/RunAll_LogRatio.m`
4. run `code/RunAll_LogRatio.m`

Detailed setup:

- `docs/INSTALLATION.md`

## Requirements

- MATLAB
- Mapping Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox
- local configuration file based on `config/paths_example.m`
- local external data stored outside this repository
- local geospatial support needed by the MATLAB pipeline
- Image Processing Toolbox is likely needed depending on workflow details
- Python/Conda environment for `code/01_preprocessing/Coreg_AROSICS_runner.m` if AROSICS-based coregistration is used

## Tested Environment

- tested on `MATLAB R2025a`
- tested on `Windows`
- full methodology pipeline execution-tested on the local workstation used to prepare this repository

## Repository Structure

The main production path is organized and locally configurable:

- `code/01_preprocessing`
- `code/02_baseline_construction`
- `code/03_classification`
- `code/04_evaluation`
- `config/paths_example.m`
- `config/paths_local.m` (local only, ignored by Git)

The production runner and the main methodology sections no longer depend on private machine-root defaults.
The current public candidate has no remaining machine-root references in code or docs.

## Local Setup

The production code should run against local data that is stored outside the repository.

The repository uses:

- a public config template in `config/paths_example.m`
- a machine-local config file in `config/paths_local.m`

Additional setup detail:

- `docs/LOCAL_SETUP.md`
- `docs/INSTALLATION.md`

## Publishing Rule

Do not add RAW imagery, processed `.mat` datasets, generated outputs, or local archives to this repository.

## Sample Folder

- sample/ contains lightweight public examples only and must not contain real imagery, real .mat products, or private run outputs.

## Processed Examples

- sample/processed_examples/ may contain a small number of curated public processed outputs, such as tiny baseline or classification examples.
- Full processed datasets and full run outputs are still excluded.

