# Changelog

## 1.0.0 (prepared 17 Nov 2026; date set when tagged)

### Peer-review analyses (AJSTR R1)
- Leave-one-June-date-out Stage 2 test, surface-state values, mask depth and Stage 1 LOO provenance: `docs/R1_LOO_DATE.md`, `scripts/run_r1_*.py`.
- Machine-local Allen Coral Atlas bathymetry path (`aca_bathymetry` / `CORALREMOTE_ACA_BATHY`); no hard-coded local paths remain.


### Python becomes the maintained implementation
- New package `coralremote/`: Stage 1, the baseline coral mask, Stage 2 (hierarchical, curated PLE), the wave FFT diagnostic, and a TreeBagger-compatible random forest (MATLAB split rules) with a reader for saved MATLAB models.
- Verified against the AJSTR paper-of-record runs before any fix:
  - Table 5 and the Stage 2 confusion matrix are exact;
  - Stage 1 OOB, test and transfer results are exact;
  - the baseline mask is exact to 0 px;
  - all 14 wave verdicts match.
- The previous MATLAB pipeline moved unchanged to `legacy/matlab/`. It is frozen.

### Fixed: reflectance double-scaling in the Stage 2 inputs
- The source stacks used for the AJSTR Stage 2 analysis (built 5 Apr 2026) divided PlanetScope SR by 10,000 twice. The bare 1e-6 feature floor therefore acted as a 0.01-reflectance floor and clipped NIR on 23–84% of coral-mask pixels on five dates (most on 20240307).
- New default: true reflectance with an explicit 1e-4 floor (`units.FIXED`). `units.COMPAT` reproduces the paper-of-record behaviour for comparison only.
- Effect on the paper's Stage 2 results (5 seeds, same procedure):
  - paired mean OA 97.1% vs 97.3% (archived MATLAB model 97.5%);
  - F1 within 0.015;
  - dominant class unchanged on every date;
  - some class fractions move by more than 2 pp: 20240620 BLD +4.7 pp, and on 20240514 PLE +9 pp / UNK −9 pp.
- Stage 1, the baseline mask, the wave diagnostic and the label polygons aren't affected.
- The earlier public MATLAB code (now `legacy/matlab/`) never had this bug: it divides once, on GeoTIFF input.

### Added
- Unit contract: `assert_reflectance` on load; scale-invariance and NIR-identity tests.
- Regression tile (P8): 400 pixels × 14 dates, built locally; aggregate golden values in `tests/golden/`.
- Machine-local paths via the git-ignored `paths.local.json` (template `paths.example.json`) or environment variables.
- Py8 seed-sensitivity diagnosis (`docs/PY8_SEED_SENSITIVITY.md`, `scripts/run_py8_*.py`).

### Fixed: endless loop in tree growing
- When a node's two lowest values were neighbouring single-precision numbers, their midpoint could round onto the lower one. That cut sent no rows left, and the node split the same way forever. Seen once, in a 1,000-tree fixed-units stress forest. Such a node is now a leaf. The guard only acts where growth never finished before, so no earlier result changes.

### Known limitations
- Retraining isn't bit-identical to MATLAB, because MATLAB's random stream can't be reproduced. On 20240514 and 20250518 the PLE/UNK split varies 17–30 pp between seeds even within MATLAB-rule training. Py8 traces 20240514 to bootstrap noise in the stress model near the 0.55 cut (more trees or seed averaging narrow it). 20250518 depends on the gate's training-row draw and is outside the training range on `axis_t` (more trees don't help). Report 20240514 as a seed mean with its range and 20250518 as a range.
- Stage 1 trees grown in Python are about 9% larger than MATLAB's. This doesn't change any metric beyond seed noise.
- Not ported:
  - co-registration (AROSICS);
  - the autolabel scripts (frozen labels are used);
  - Stage 1 leave-one-polygon-out and spatial-block checks;
  - k-means diagnostics;
  - figure scripts.
