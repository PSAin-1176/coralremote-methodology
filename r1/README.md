# R1: code behind the revised (round 1) tables

This folder holds the code that produced the numbers in the R1 manuscript and its calibration supplement. It uses the `coralremote` package of this repository (corrected input scaling, `fixed` mode). The submitted (pre-revision) tables came from the earlier configuration, which `compat` mode reproduces for provenance (see the top-level README).

| File | Purpose |
|---|---|
| `revision_runner.py` | Calibration, Stage 1, consensus mask, Stage 2, date-omission and leave-one-polygon-out tests |
| `legacy_replay.py` | Replays the archived (submitted) Stage 2 models in `compat` mode; not the revised analysis |
| `table3_retest.py` | Table 3 feature-family ablation on the Table 4A held-out design |
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
| `REVISION_RUN`, `REVISION_OLD_RUN`, `TABLE3_OUT` | `table3_retest.py` only: the completed revision run, the earlier exploratory ablation run (read for the spatial-CV column) and the output folder |

The runner stops with a clear error if `CORALREMOTE_ROOT` is missing. Before any output is written, `primary` checks every input hash in the configuration and the code identity (`vendor/coralremote/*` entries are compared with this repository's `coralremote/*`; the recorded hashes of the two entry-point scripts predate the path edits and are reported, not enforced). Existing inconsistent caches are rejected.

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

# Table 3 (needs the completed run in REVISION_RUN)
python r1/table3_retest.py

# Archived (submitted) Stage 2 replay in compat mode
python r1/legacy_replay.py --json legacy_replay.json
```

Software used for the R1 runs: Python 3.13.9, numpy 2.4.4, scipy 1.16.3, scikit-learn 1.7.2, rasterio 1.5.0, h5py 3.15.1, matplotlib 3.10.6 and, for `legacy_replay.py`, mat-io 1.0.2. Runs used one thread per worker.

## Which command produces which manuscript item

Output names are relative to `--out`. Entries under `report/` are summary tables written by post-processing scripts of the revision workspace that are not part of this release; the same values can be read from the runner outputs listed next to them. "Verified" means the manuscript value was compared with the stored run output on 21 Nov 2026.

| Manuscript item | Command | Output file | Check |
|---|---|---|---|
| Pooled StableSand OLS pair 0.381 / 0.374 (printed 0.38 / 0.37), 7,487 observations, R2 0.799 / 0.802, bootstrap intervals, per-date fits and eligibility (3 of 8 pass), common-support slopes, leave-one-calibration-date-out pairs (Supplement S1) | `primary` or `calibrate` | `calibration/calibration.json`, `calibration/date_coefficients.csv`, `report/calm_date_common_support.csv` | Verified |
| 20,460-pixel baseline mask (seed 0); 20,304 and 20,326 (seeds 1, 2) | `primary` with each seed configuration | `stage1/pooled_mask_seed<s>_masks.npz` (`maskFinal`) | Verified |
| Table 2 (training composition: gate 6,000 = 3,000 / 1,064 / 1,936; stress 3,101 = 1,098 / 2,003; OOB errors) | `primary` | `stage2/pooled_seed0.json`, `report/candidate_Table2.csv`, `report/OOB_and_drawn_training_counts.csv` | Counts and OOB verified; polygon counts (HPY 374/293, PLE 43/29, BLD 71/48) not checked |
| Table 3 (feature ablation: test OA 57.9, 82.0, 76.0, 81.5; CORAL recall; OA ranges over three seeds) | `table3_retest.py` | `results.csv`, `results_per_seed.csv`, `sanity_gate.json` in `TABLE3_OUT` | Test OA and recall verified. The spatial-CV column (0.741, 0.867, 0.860, 0.871) is carried from the earlier exploratory run (`REVISION_OLD_RUN`), which is not in this release: unverified here |
| Table 4A (Stage 1, 905-pixel test: OA 81.5%, kappa 0.701, CORAL recall 0.325) | `primary` | `stage1/pooled_accuracy_seed0.json`, `report/stage1_accuracy.csv`, `report/candidate_Table4A.csv` | Verified |
| Leave-one-polygon-out mean recall 0.928 (36 polygons, 3 seeds) | `loo` | `loo/polygon<id>_seed<s>.json`, `report/LOO_polygons.csv` | Verified (mean 0.9282; range 0.522 to 1.000) |
| Stage 1 transfer to other dates (0.90 on 20240307, 0.035 on 20240317) | `primary` (accuracy model), summary by post-processing | `report/Stage1_transfer_metrics.csv` | Verified |
| Table 4B (Stage 2, 8,929 test pixels: core OA 97.1%, including UNK 96.3%) | `primary` | `stage2/pooled_seed0.json`, `report/candidate_Table4B.csv` | Verified |
| Table 5 (14 dates, n = 20,460 each) | `primary` | `stage2/pooled_seed0.json` (`table`), `report/candidate_Table5.csv` | Verified |
| Supplement S3: seed ranges for the primary trajectory | `primary` for seeds 0, 1, 2 | `stage2/pooled_seed<s>.json`, replicate columns of `report/candidate_Table5.csv` | Verified |
| Supplement S3: date-specific policy with pooled fallbacks (20 June order reverses; spatial-CV macro-F1 0.866 pooled, 0.863 date-specific) | `stage1`, `select`, `stage2` | `stage2/per_date_seed<s>.json`, `stage1/per_date_mask_seed<s>_masks.npz`, `manifest/selection.json` | Verified (20 June PLE 48.54 to 49.94, BLD 43.12 to 45.21) |
| Supplement S4: June date omission (held-date BLD, PLE and stressed fractions, 3 seeds) | `omissions` | `omissions/pooled_omit<date>_seed<s>.json`, `report/June_omission_results.csv` | Verified |
| Supplement S5: BLD-DHW correlation 0.8151 (ten dates) | post-processing | `report/DHW_correlations.csv` | Verified |

Items in the supplement that could not be tied to a stored run output are marked **unverified**: the polygon-level leave-one-out values of the pooled pair and the 51.6% polygon share (S1), the Lyzenga variance-covariance estimator comparison (S1), the deep-water subtraction sensitivity table (S1), the ten-configuration table of the earlier binary analysis (S2), the surface-state thresholds and NIR-threshold shares (S4; they come from the wave diagnostic and the `scripts/run_r1_*.py` analyses of the earlier configuration, not from this runner), the uncurated pale-candidate distributions (S4) and the scene identifiers (S7).

## Notes

- Seed 0 is the primary run; seeds 1 and 2 report fitting variability, not confidence intervals.
- The date-specific calibration policy is a sensitivity analysis. The pooled pair is the operational choice (S1).
- Every output is written to a folder given by `--out` (default `r1/output/`, which is git-ignored). Originals are never modified.
