# Sample Folder

This folder is for lightweight public sample materials only.

## Allowed contents

- small dummy config examples
- tiny mock CSV manifests
- minimal fake metadata for folder structure examples
- toy rasters created only for testing repository wiring
- redacted examples of expected file naming conventions

## Do not place here

- PlanetScope RAW imagery
- real reef analysis rasters
- processed `.mat` datasets
- trained models
- generated manuscript figures
- private polygons or private run outputs

## Recommended contents for first release

- `sample_manifest.csv`
  - example rows showing date, product name, and expected local path pattern
- `sample_dates.txt`
  - example date list for a tiny dry-run workflow
- `sample_folder_tree.txt`
  - expected local data layout outside Git
- `sample_config_notes.md`
  - short explanation of how sample names map to `config/paths_local.m`

## Purpose

The sample folder should help a new user understand input structure without publishing real data.
