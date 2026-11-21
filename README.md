# CoralRemote

Python pipeline for mapping coral bleaching from PlanetScope SuperDove imagery: a log-ratio substrate screen (Stage 1), a consensus baseline coral mask, a hierarchical healthy / pale / bleached classifier (Stage 2) and an image-only wave diagnostic. Study site: Koh Man Nai, Rayong, Thailand, 14 dates from 2023 to 2025.

The revised (round 1, R1) AJSTR paper (Sa-ngadsup, Haq and Cheewinsiriwat, *Multitemporal PlanetScope-based assessment of spectral signatures of coral bleaching and recovery at Koh Man Nai, Gulf of Thailand*) reports tables rebuilt in Python from the code in this repository. The submitted (pre-revision) tables were produced by the author's original MATLAB workflow, whose surviving scripts and archived outputs are identified in [docs/PROVENANCE.md](docs/PROVENANCE.md). Those scripts are not part of this repository.

**History.** Earlier versions of this repository (up to commit 2c0d556) contained a two-class MATLAB prototype of the workflow. It predates the hierarchical method in the AJSTR paper and was not used to produce its results, so it was removed from the current version to avoid confusion. It remains available in the repository history at that commit (tag `matlab-prototype`).

The MATLAB licence was no longer available when the analyses requested in peer review were run, so the workflow was reimplemented in Python. The package has two unit modes:
- **`compat` mode:** loaded with the archived MATLAB-trained forests, it reproduces the SUBMITTED (pre-revision) Tables 4 and 5 exactly. It is kept for provenance and parity tests, not as the R1 analysis.
- **`fixed` mode:** corrects the input scaling described in Section 3.5.6 of the paper (see [Units](#units-and-known-issue)). The R1 analysis uses it.

## Results of the R1 paper: where they come from

| Manuscript item | Source |
|---|---|
| Tables 2, 4 and 5; the 20,460-pixel baseline mask; the pooled StableSand OLS coefficients (0.381 and 0.374, printed 0.38 and 0.37); Supplement S1 (calibration), S3 (seeds 0 to 2 and the date-specific policy) and S4 (June date omission) | [`r1/revision_runner.py`](r1/revision_runner.py) with the configurations in `r1/` |
| Table 3 (feature ablation): test columns and spatial-CV column | [`r1/table3_retest.py`](r1/table3_retest.py), [`r1/table3_spatial_cv.py`](r1/table3_spatial_cv.py) |
| Summary tables and context statistics (Tables 2, 4A, 4B, 5; S1, S3 to S5) | [`r1/make_tables.py`](r1/make_tables.py) from a completed run |
| S1 calibration sensitivities (largest polygon, Lyzenga estimator, deep-water subtraction) | [`r1/calibration_sensitivity.py`](r1/calibration_sensitivity.py) |
| S7 scene identifiers | [`r1/scenes.csv`](r1/scenes.csv) |
| Submitted (pre-revision) Tables 4 and 5, and the 20,330-pixel mask | `compat` mode ([`scripts/run_frozen_stage2.py`](scripts/run_frozen_stage2.py), [`r1/legacy_replay.py`](r1/legacy_replay.py)) |

[`r1/README.md`](r1/README.md) lists the settings, the commands, and the command and output file behind each table and supplement section, each checked by rerunning the release scripts against the stored R1 outputs. Only the S2 table of the earlier binary analysis is not rebuilt: it reports archived results of the submitted workflow. The earlier review analyses on the submitted configuration (leave-one-date-out, surface-state values, mask depth) are in [docs/R1_LOO_DATE.md](docs/R1_LOO_DATE.md). They are superseded where the R1 paper reports a rebuilt value.

## What's in it
| Module | Content |
|---|---|
| `units`, `io`, `paths` | Unit contract (`assert_reflectance`, `compat` / `fixed` modes), loaders for the source stacks, mask, labels and saved MATLAB runs |
| `preprocess` | Rebuild of the Phase 3–5 Stage 1 features (glint reflectance, UDM2 + wave masking, B/G1, G1/G2, B/R, DII) |
| `stage1` | Stage 1 substrate classifier (CORAL / SND / RBL), evaluation and cross-date transfer |
| `mask` | Baseline coral mask: Stage 1 votes, then the morphology chain |
| `features`, `stage2` | 21-feature hybrid set, PLE curation on the healthy→bleached axis, hierarchical gate/stress prediction with the 0.55 UNK rule |
| `cart`, `forest`, `train` | TreeBagger-equivalent random forests (MATLAB split rules); a reader for saved MATLAB TreeBagger models |
| `wave` | NDWI-tile FFT diagnostic: CALM / TURBULENT / WAVE |
| `evaluate`, `regression` | Metrics, report parsing, P8 regression tile |

## Verification against the submitted paper's archived MATLAB runs (`compat` mode)
These checks show that the Python port reproduces the submitted (pre-revision) results. They are parity checks, not the R1 tables (see the table above).

| Check | Result | Report |
|---|---|---|
| Paper's Stage 2 forests on Python features | Table 5 (56 cells) and the test confusion matrix identical (OA 0.974841, κ 0.934418) | `scripts/run_frozen_stage2.py` |
| Paper's Stage 1 forest on rebuilt features | OOB 0.1004, test confusion matrix and the 8-date transfer table identical | [PY4_TRAINING_PARITY.md](docs/PY4_TRAINING_PARITY.md) |
| Retraining (Stage 1, Stage 2) | Within seed-to-seed variation of the paper | [PY4_TRAINING_PARITY.md](docs/PY4_TRAINING_PARITY.md) |
| Baseline mask chain | Every saved intermediate identical (0 px), 20,330 px | [PY5_MASK_WAVE.md](docs/PY5_MASK_WAVE.md) |
| Wave diagnostic | All 14 verdicts and confidence tiers identical | [PY5_MASK_WAVE.md](docs/PY5_MASK_WAVE.md) |
| Units fix re-run | See below | [PY6_FIX_AND_RERUN.md](docs/PY6_FIX_AND_RERUN.md) |

Where each MATLAB script is reproduced: [docs/PROVENANCE.md](docs/PROVENANCE.md).

## Units and known issue
All arrays are surface reflectance on a 0–1 scale. Reflectance features use an explicit floor of 1e-4 reflectance (`units.FIXED`). Details and Methods wording: [docs/UNITS.md](docs/UNITS.md).

The source stacks behind the AJSTR paper (built 5 Apr 2026) divided reflectance by 10,000 twice. The 1e-6 floor in Stage 2 then acted as a 0.01 reflectance floor and clipped NIR on 23–84% of coral-mask pixels on five dates. `units.COMPAT` reproduces that behaviour exactly, for parity only and to regenerate the submitted tables. The R1 paper discloses this (Section 3.5.6) and reports tables rebuilt with `units.FIXED` by `r1/revision_runner.py`; the submitted values are not carried over.

Paired comparison: five Python retrainings in each mode, same procedure; the corrected floor is 1e-4 and the PLE curation is recomputed (81 vs 88 polygons).
- Mean OA moves from 97.3% to 97.1%; the archived MATLAB model scores 97.5%.
- Every date's dominant class is unchanged.
- Mean class fractions move by up to 9 pp, most on 20240514 (PLE 25.9% → 35.0%).

See [PY6_FIX_AND_RERUN.md](docs/PY6_FIX_AND_RERUN.md). The units fix itself does not affect Stage 1, the baseline mask, the wave diagnostic or the training labels. The R1 mask (20,460 pixels, against 20,330 submitted) differs because R1 re-estimates the sand-regression coefficients (pooled 0.381 and 0.374 in place of the earlier fixed 0.82 and 0.80), not because of the units.

## Install
```bash
conda env create -f environment.yml
conda activate coralremote
pip install -e .
```

## Data
No imagery, reflectance stacks, feature stacks or model files are in this repository: PlanetScope imagery and the reflectance-derived stacks are licensed and can't be redistributed. Inputs are identified by SHA-256 and size (see [docs/INPUTS.md](docs/INPUTS.md) and the `r1/` configurations), so a reader with authorised copies can confirm they hold the same files. Point the package to local data by:
- copying `paths.example.json` to `paths.local.json` (git-ignored), or
- setting `CORALREMOTE_ROOT`, `CORALREMOTE_A3_ROOT`, `CORALREMOTE_RAW_PLANET` and `CORALREMOTE_CACHE`.

The R1 commands also take placeholders for the input folders; see [r1/README.md](r1/README.md). The full list of inputs (files, variables, units, grid, original run IDs, SHA-256 checksums, and how each can be obtained or rebuilt) is in [docs/INPUTS.md](docs/INPUTS.md). Most commands need archived products that only the author can supply; that file says which.

Every write by the package and scripts goes through `paths.assert_writable`, which refuses paths inside the protected source folders (`paths.READ_ONLY_DIRS`). This covers the cache folder, forest exports and the scripts' output arguments (`--json` and similar). The `r1/` scripts also refuse output inside any configured input location (A3 outputs, UDM2 flags, StableSand, deep-water masks, bathymetry, scene records). Derived arrays go to the cache folder. The cache is keyed on the source file's path, size and modification time, so a replaced input is re-read rather than served from a stale cache.

## Tests
```bash
pytest -m "not slow"          # about 30 s; tests marked `data` skip when the data is absent
pytest                        # adds Stage 2 retraining parity (several minutes)
```
The regression tile (`tests/test_regression_tile.py`, plan P8) is 400 baseline-mask pixels on all 14 dates. It is cut from the local data into the cache. Only aggregate golden numbers are committed (`tests/golden/`); regenerate them with `scripts/make_regression_golden.py` when a change is intended.

## Licence and citation
MIT (see `LICENSE`). Cite the release DOI once it's archived (see [CHANGELOG](CHANGELOG.md)).
