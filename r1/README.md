# R1: code behind the revised (round 1) tables

This folder holds the code that produced the numbers in the R1 manuscript and its calibration supplement. It uses the `coralremote` package of this repository (corrected input scaling, `fixed` mode). The submitted (pre-revision) tables came from the earlier configuration, which `compat` mode reproduces for provenance (see the top-level README).

| File | Purpose |
|---|---|
| `revision_runner.py` | Calibration, Stage 1, consensus mask, Stage 2, date-omission and leave-one-polygon-out tests |
| `legacy_replay.py` | Replays the archived (submitted) Stage 2 models in `compat` mode; not the revised analysis |
| `table3_retest.py` | Table 3 feature-family ablation on the Table 4A held-out design |
| `table3_spatial_cv.py` | Table 3 spatial-CV column and the earlier exploratory test column (S2) |
| `make_tables.py` | Summary tables and context statistics from a completed run (Tables 2, 4A, 4B, 5; S1, S3, S4, S5) |
| `calibration_sensitivity.py` | S1 sensitivities: largest polygon, Lyzenga estimator, deep-water subtraction |
| `scenes.csv`, `scenes.py` | S7 scene identifiers and acquisition times; `scenes.py` rebuilds them from licensed metadata |
| `revision_config.json`, `revision_config_seed1.json`, `revision_config_seed2.json` | Input identity (SHA-256 and size) and fixed parameters for seeds 0, 1 and 2 |

No imagery, reflectance stacks or derived arrays are included. PlanetScope data are licensed and cannot be redistributed. The configurations identify every input by SHA-256 and size, so a reader with authorised copies can check that they hold the same files.

## Settings

Paths in the configurations are placeholders. Set these environment variables before running:

| Variable | Meaning |
|---|---|
| `CORALREMOTE_ROOT` | Data root with `Source/` and `Asset/` (required) |
| `CORALREMOTE_A3_ROOT` | Folder of the A3 pipeline outputs (`Phase2_Masking`, `Phase3_LogRatio`) |
| `REVISION_UDM` | Folder of the coregistered UDM2 flag files `udm2_<date>.npz` (default `<root>/coralremote-python-cache/udm2`) |
| `CORALREMOTE_STABLESAND` | The `StableSand.geojson` file |
| `REVISION_PLAN` | Optional: a plan document whose hash is recorded in `manifest/design.json` |
| `REVISION_RUN` | `table3_retest.py`, `make_tables.py`: output folder of the completed revision run |
| `REVISION_OLD_RUN`, `TABLE3_OUT` | `table3_retest.py`: output folder of `table3_spatial_cv.py` (`TABLE3_CV_OUT`), and its own output folder |
| `CORALREMOTE_ACA_BATHY` | `make_tables.py`: Allen Coral Atlas bathymetry GeoTIFF (mask depth summary) |
| `DEEP_MASK_BM`, `DEEP_MASK_MDEEP` | `calibration_sensitivity.py`: the Phase 2 deep-water masks (`bm_deep_mask.tif`, `mdeep_mask.tif`) |
| `PLANET_SCENES` | `scenes.py`: folder of the Planet STAC item records `<scene id>.json` |
| `TABLES_OUT`, `TABLE3_CV_OUT`, `CALIB_OUT`, `SCENES_OUT` | Optional output folders (default under `r1/output/`) |

The runner stops with a clear error if `CORALREMOTE_ROOT` is missing. Before any output is written, `primary` checks every input hash in the configuration and the code identity (`vendor/coralremote/*` entries are compared with this repository's `coralremote/*`, and the two entry-point scripts with their released versions; code hashes use LF line endings, so Windows and Unix checkouts give the same result). Existing inconsistent caches are rejected.

## Commands

Run from the repository root. The full pipeline takes hours; use new output and cache folders.

```bash
# Primary run for one seed (calibration, Stage 1 mask and accuracy models, mask, Stage 2)
python r1/revision_runner.py primary --config r1/revision_config.json --out <OUT> --cache <OUT>/cache
python r1/revision_runner.py primary --config r1/revision_config_seed1.json --out <OUT_SEED1> --cache <OUT_SEED1>/cache
python r1/revision_runner.py primary --config r1/revision_config_seed2.json --out <OUT_SEED2> --cache <OUT_SEED2>/cache

# Staged commands for the full revision run (all policies, seeds and sensitivity tests)
# in this order, each with --out <OUT>: freeze, calibrate, prepare, stage1, select, stage2, loo, omissions
python r1/revision_runner.py calibrate --out <OUT>

# Table 3: spatial-CV column first, then the test columns (needs the completed run in REVISION_RUN)
python r1/table3_spatial_cv.py
python r1/table3_retest.py

# Summary tables, S1 sensitivities and the scene list
python r1/make_tables.py
python r1/calibration_sensitivity.py
python r1/scenes.py

# Archived (submitted) Stage 2 replay in compat mode
python r1/legacy_replay.py --json legacy_replay.json
```

Software used for the R1 runs: Python 3.13.9, numpy 2.4.4, scipy 1.16.3, scikit-learn 1.7.2, rasterio 1.5.0, h5py 3.15.1, matplotlib 3.10.6 and, for `legacy_replay.py`, mat-io 1.0.2. Runs used one thread per worker.

## Which command produces which manuscript item

Runner outputs are relative to `--out`; summary tables are written by `make_tables.py` (to `TABLES_OUT`). "Verified" means the release script was run on 21 Nov 2026 against the stored outputs of the R1 runs and reproduced them (numbers to 1e-9, model fits bit for bit).

| Manuscript item | Script | Output | Check |
|---|---|---|---|
| Pooled StableSand pair 0.381 / 0.374 (printed 0.38 / 0.37), 7,487 observations, R2, bootstrap intervals, single-date fits and eligibility, leave-one-date-out pairs (S1) | `revision_runner.py primary` or `calibrate` | `calibration/calibration.json`, `calibration/date_coefficients.csv` | Verified |
| S1 common-pixel table (2,112 locations) | `make_tables.py` | `calm_date_common_support.csv` | Verified |
| S1 largest polygon (51.6% of pixels; 0.203 / 0.209 without it), single-polygon range, Lyzenga estimator 0.393 / 0.386, deep-water subtraction table, 1,006 of 7,487 lost, deep-water red 59 to 85% of sand red, reef pixels below deep water 0.2 to 21% | `calibration_sensitivity.py` | `results.csv`, `per_date.csv`, `deepwater.csv`, `paired_bootstrap.csv`, `reef_below_deepwater.csv` | Verified |
| 20,460-pixel mask (seed 0); 20,304 and 20,326 (seeds 1, 2); 10.4% of the reef, median atlas depth 0.8 m (Section 3.5.6) | `revision_runner.py primary`; `make_tables.py` | `stage1/pooled_mask_seed<s>_masks.npz`; `paired_masks.csv`, `context_statistics.json` | Verified |
| Table 2 (gate 6,000 = 3,000 / 1,064 / 1,936; stress 3,101 = 1,098 / 2,003; OOB errors; polygons HPY 374/293, PLE 43/29, BLD 71/48) | `make_tables.py` | `candidate_Table2.csv`, `OOB_and_drawn_training_counts.csv`, `reference_polygons.json` | Verified |
| Table 3 test columns (OA 57.9, 82.0, 76.0, 81.5; CORAL recall) | `table3_retest.py` | `results.csv`, `results_per_seed.csv` | Verified |
| Table 3 spatial-CV macro-F1 (0.741, 0.867, 0.860, 0.871); S2 earlier exploratory test column (70.0, 96.9, 95.0, 96.7%) | `table3_spatial_cv.py` | `report/ablation_summary.csv`, `stage1_ablation/*.json` | Verified |
| Table 4A (905-pixel test: OA 81.5%, kappa 0.701, CORAL recall 0.325) | `make_tables.py` | `candidate_Table4A.csv`, `stage1_accuracy.csv` | Verified |
| Leave-one-polygon-out mean recall 0.928 (0.522 to 1.000) | `revision_runner.py loo`; `make_tables.py` | `loo/*.json`; `LOO_polygons.csv` | Verified |
| Stage 1 transfer (0.90 on 20240307, 0.035 on 20240317) | `make_tables.py` | `Stage1_transfer_metrics.csv` | Verified |
| Table 4B (8,929 pixels: core OA 97.1%, 96.3% including UNK); anchor checks (0.810 on 1,262; 1.000 on 418) | `make_tables.py` | `candidate_Table4B.csv`, `context_statistics.json` | Verified |
| Table 5 (14 dates, n = 20,460) and S3 seed ranges | `make_tables.py` | `candidate_Table5.csv` | Verified |
| S3 date-specific policy (20 June PLE 48.54 to 49.94, BLD 43.12 to 45.21), mask IoU and area change | `make_tables.py` | `paired_stage2_dates.csv`, `paired_masks.csv` | Verified |
| S4 June date omission | `revision_runner.py omissions`; `make_tables.py` | `omissions/*.json`; `June_omission_results.csv` | Verified |
| S4 common (8,854) and location-disjoint (8,320) subsets, PLE F1 0.776 / 0.783 | `make_tables.py` | `stage2_evaluation_and_training_counts.csv` | Verified |
| S4 uncurated pale candidates; NIR-threshold shares; haze on 20240616 (90.1% of the mask) | `make_tables.py` | `context_statistics.json` | Verified |
| S4 surface-state thresholds (median + 1.5 MAD per band) | `scripts/run_r1_task2.py` (wave diagnostic; does not depend on the calibration) | values in `docs/R1_LOO_DATE.md` | Verified against the documented run |
| S5 May 2025 outside the axis range (53.91 to 54.13%); BLD-DHW r 0.8151 (0.8146 to 0.8155), without June 0.156 to 0.479 | `make_tables.py` | `context_statistics.json`, `DHW_correlations.csv` | Verified |
| S7 scene identifiers and times | `scenes.py` | `scenes.csv` (committed) | Verified against the Planet item records |

Not rebuilt by this release: the S2 ten-configuration table of the earlier binary analysis. It reports archived results of the submitted (MATLAB) workflow and is labelled as such in the supplement.


## Notes

- Seed 0 is the primary run; seeds 1 and 2 report fitting variability, not confidence intervals.
- The date-specific calibration policy is a sensitivity analysis. The pooled pair is the operational choice (S1).
- Every output is written to a folder given by `--out` (default `r1/output/`, which is git-ignored). Originals are never modified.
