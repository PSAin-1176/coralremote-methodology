# CoralRemote

Python pipeline for mapping coral bleaching from PlanetScope SuperDove imagery: a log-ratio substrate screen (Stage 1), a consensus baseline coral mask, a hierarchical healthy / pale / bleached classifier (Stage 2) and an image-only wave diagnostic. Study site: Koh Man Nai, Rayong, Thailand, 14 dates from 2023 to 2025.

Version 1.0 replaces the earlier MATLAB code, which is kept unchanged in [`legacy/matlab/`](legacy/matlab/README.md). The Python code was first checked against the MATLAB runs behind the AJSTR paper, then a reflectance-units error was fixed (see [Units](#units-and-known-issue)).

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

## Verification against the paper of record
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

The source stacks behind the AJSTR paper (built 5 Apr 2026) divided reflectance by 10,000 twice. The 1e-6 floor in Stage 2 then acted as a 0.01 reflectance floor and clipped NIR on 23–84% of coral-mask pixels on five dates. `units.COMPAT` reproduces that behaviour exactly, for parity only. Re-running with correct units:
- leaves OA (−0.4 pp), F1 (≤ 0.015) and every date's dominant class unchanged;
- moves some Table 5 cells by more than 2 pp, for example 20240620 BLD +4.7 pp and 20240514 PLE +9 pp / UNK −9 pp.

See [PY6_FIX_AND_RERUN.md](docs/PY6_FIX_AND_RERUN.md). Stage 1, the baseline mask, the wave diagnostic and the training labels aren't affected.

## Install
```bash
conda env create -f environment.yml
conda activate coralremote
pip install -e .
```

## Data
No imagery, reflectance stacks, feature stacks or model files are in this repository: PlanetScope data can't be redistributed. Point the package to local data by:
- copying `paths.example.json` to `paths.local.json` (git-ignored), or
- setting `CORALREMOTE_ROOT`, `CORALREMOTE_A3_ROOT`, `CORALREMOTE_RAW_PLANET` and `CORALREMOTE_CACHE`.

The package never writes inside the source folders (`paths.READ_ONLY_DIRS`). Derived arrays go to the cache folder.

## Tests
```bash
pytest -m "not slow"          # about 30 s; tests marked `data` skip when the data is absent
pytest                        # adds Stage 2 retraining parity (several minutes)
```
The regression tile (`tests/test_regression_tile.py`, plan P8) is 400 baseline-mask pixels on all 14 dates. It is cut from the local data into the cache. Only aggregate golden numbers are committed (`tests/golden/`); regenerate them with `scripts/make_regression_golden.py` when a change is intended.

## Licence and citation
MIT (see `LICENSE`). Cite the release DOI once it's archived (see [CHANGELOG](CHANGELOG.md)).
