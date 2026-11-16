# Local Setup

The repository uses a split config model.

## Public file

- `config/paths_example.m`

## Local file

- `config/paths_local.m`

Create `config/paths_local.m` by copying `config/paths_example.m` and editing the machine-specific paths.

Required fields:

- `RAW_DIR`
- `ACA_DIR`
- `OUT_ROOT`
- `PREP_SOURCE_ROOT`
- `MASTERLOG_DIR`

## Runner

Main entry point:

- `code/RunAll_LogRatio.m`

The runner loads local paths from `config/paths_local.m` and writes `pipeline_logratio_config.mat` into the selected output root.
