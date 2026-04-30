# CoralRemote Production Candidate

This repository candidate contains the current active methodology pipeline first.

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

## Current state

This folder is a GitHub publication candidate, not the final published repository yet.

The main production path is now reorganized and locally configurable:

- `code/01_preprocessing`
- `code/02_baseline_construction`
- `code/03_classification`
- `code/04_evaluation`
- `config/paths_example.m`
- `config/paths_local.m` (local only, ignored by Git)

The production runner and the main methodology sections no longer depend on private machine-root defaults.
The current public candidate has no remaining machine-root references in code or docs.

## Expected local setup

The production code should run against local data that is stored outside the repository.

The repository uses:

- a public config template in `config/paths_example.m`
- a machine-local config file in `config/paths_local.m`

## Publishing rule

Do not add RAW imagery, processed `.mat` datasets, generated outputs, or local archives to this repository.




## Sample Folder

- sample/ contains lightweight public examples only and must not contain real imagery, real .mat products, or private run outputs.


## Processed Examples

- sample/processed_examples/ may contain a small number of curated public processed outputs, such as tiny baseline or classification examples.
- Full processed datasets and full run outputs are still excluded.

