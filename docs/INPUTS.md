# Inputs: what the scripts need and how to get it

None of the files below is in this repository. Point the package to them with `paths.local.json` or environment variables (see the README, Data section).

**Path settings:**
- `<ROOT>` = `CORALREMOTE_ROOT`;
- `<A3>` = `CORALREMOTE_A3_ROOT`;
- `<RAW>` = `CORALREMOTE_RAW_PLANET`;
- `<ACA>` = the folder of `CORALREMOTE_ACA_BATHY`.

**Analysis grid:** 1,366 rows × 1,078 columns, 3 m cells, WGS 84 / UTM zone 47N (EPSG:32647). All stacks, masks and class maps below are on this grid unless stated otherwise.

**Dates:** 14 acquisition dates, 20230404 to 20250518.

**Checksums:** SHA-256 values were computed on 18 Nov 2026 from the author's archived copies, the ones the verification reports in `docs/` were run on.

## 1. Licensed imagery (PlanetScope): not redistributable

| Item | Files | Units / scale | Needed by | How to obtain |
|---|---|---|---|---|
| PlanetScope SuperDove 8-band surface reflectance, 14 dates | `<RAW>/<scene>_3B_AnalyticMS_SR_8b_harmonized_clip.tif` | Reflectance × 10,000 (int16); native UTM 47N | Stage 1 rebuild (`preprocess`, `tests/test_stage1_parity.py`) | Order the same scenes from Planet under your own licence. Scene IDs and acquisition metadata are listed in the author's `RAW_PlanetScope_Metadata_Summary.txt`, available on request. |
| UDM2 usable-data masks | `<RAW>/<scene>_3B_udm2_clip.tif` | Bit flags | Stage 1 rebuild | Delivered with each Planet order |
| Scene metadata | `<RAW>/<scene>_3B_AnalyticMS_8b_metadata_clip.xml`, `<scene>.json` | — | Sun and view angles in the reports | Delivered with each Planet order |

No checksums are listed for the licensed scenes: they come from Planet, not from the author.

## 2. Derived stacks (built by the author from the imagery)

These contain reflectance derived from Planet data. They are shared only under Planet's terms, so ask the author.

| File | Variable(s) | Units / scale | Built by | Size (bytes) | SHA-256 |
|---|---|---|---|---|---|
| `<ROOT>/Source/SourceStack_HedleyGlintReflectance.mat` (v7.3 / HDF5) | `stack_glint` (stored N × 8 × W × H = 14 × 8 × 1078 × 1366, float32), `dates_col`, `rawBandNames`, `R0` (MapCellsReference), `sourceGlintFiles` | **Reflectance ÷ 10,000** (divided twice; declared `on_disk_scale = 10,000`) | Source stack builder, 5 Apr 2026 (script lost) | 581,581,649 | `002e47e5bc025df43c52767c3bfe153f66ecc324eaf53fe756177deefe842a04` |
| `<ROOT>/Source/SourceStack_CoregRawReflectance.mat` (v7.3) | `stack_raw` (same layout), `dates_col`, `rawBandNames`, `R0`, `sourceRawFiles` | Reflectance 0–1 (`on_disk_scale = 1`) | Same builder (co-registered raw reflectance) | 560,941,970 | `53bc867d6e0375d3f90082fc87ac7be6fc03617d5bed525b972c5378659da02f` |
| `<A3>/Phase3_LogRatio/run_20260401_102508/GLINT_<date>.tif` | 8 bands, float32 | Reflectance 0–1, glint-corrected | A3 pipeline, Phase 3 | see below | see below |

Surviving Phase 3 GLINT GeoTIFFs (8 of the original 14):

| File | Size (bytes) | SHA-256 |
|---|---|---|
| `GLINT_20230404.tif` | 47,264,924 | `e58fd2f9d1609d885bde6d2a056d926193a5ee56eb9414e9b7e5b73fd5bfae81` |
| `GLINT_20230509.tif` | 47,130,040 | `c28db91aa9c919ccf273e1b540d5eae230123df836713059e4698338aad42ad5` |
| `GLINT_20240307.tif` | 47,117,150 | `8cdaf059befb9effb25555d10c1c3f95abc5135200a2749d89c2a4a1db45c9b3` |
| `GLINT_20240514.tif` | 46,410,304 | `cd8572d4b740dc7c6e2dcf73fd9a46421a63b17324e17ebee5cc55ac6b09a287` |
| `GLINT_20240613.tif` | 46,521,746 | `d2478d490021f1c7b9dbd4f629c219040ab986c2b0b069f43dab6ff3878fc45f` |
| `GLINT_20240620.tif` | 46,344,782 | `47349d4d27859e081877217e6871a8da7f79800983cf7bc42855e3305cab96a4` |
| `GLINT_20240926.tif` | 47,338,076 | `783102766120abd3ec91524a4f6f263cfdb2609c3469510ef3ccac8d8be82cfa` |
| `GLINT_20241111.tif` | 47,355,208 | `d1ce2849d10411ee6986080e2feeb1b227d4a7774f631637292dc2710de5f356` |
| `logratio_report.txt` | 1,940 | `312c6f5c890831a88b275850c9b87db716588bdfe416fab7885a3700de0cf238` |

**Rebuilding:** the Python `preprocess` module rebuilds the Stage 1 features from the raw stack plus the UDM2 masks. The 5 Apr Source stack builder itself is lost, so the two Source stacks can't be regenerated bit for bit from Planet scenes with this repository alone. They need Phase 1 co-registration (AROSICS, not ported) and the Hedley glint step on the exact scenes.

## 3. Masks and class maps (no reflectance)

| File | Variable(s) / content | Units | Original run | Size (bytes) | SHA-256 |
|---|---|---|---|---|---|
| `<ROOT>/Asset/baseline_coral_mask.mat` | `baseline_coral_mask`, `H`, `W`, `reef_mask`, `vote_map`, `tier_map`, `mask0`, `maskClean`, `maskPeel`, `maskFinal` | Logical / uint8 | `Baseline_Coral/run_20260401_103612` + 6 Apr clean-up | 252,461 | `e1108ff43bd81e90bab023d8ef14c83e9005e0e88c5cb98b2953ddb3c611b57d` |
| `<A3>/Phase2_Masking/run_20260401_102458/reef_mask.tif` | Reef polygon raster | uint8 | Phase 2 | 31,818 | `44320c4de1cb6024bbe2d851efa0c58b24cdb4ac27fadad70075a118879ddace` |
| `<A3>/Phase9_Classify/run_20260401_102734/CLASS_20230404.tif` | Stage 1 class map (0–4) | uint8 | Phase 9 (from the lost Stage 1 model `run_20260401_102559`) | 121,524 | `d04f18e76abb17db89cb3f40d1f21efd06c4c188d1c76f84c6b5584aaebee136` |
| `…/CLASS_20230509.tif` | as above | | | 114,174 | `ce190334eb94ebc536c405b2674cd1adcf8a1b4a935f4b37aa6eceb68686ad0c` |
| `…/CLASS_20240307.tif` | as above | | | 131,218 | `988123212c3909d3abcef92907f55621898de96bd2d838c2c553a1e80a5f18c4` |
| `…/classification_report.txt` | Phase 9 report | | | 2,728 | `e3896331a843c5faf27c3c8874c3d5293b392b9ede9e7e68ce6154491d714090` |
| `<ACA>/bathymetry_0.tif` | Allen Coral Atlas bathymetry | Depth in cm, positive down, int16, nodata 0; EPSG:4326, about 10 m | Downloaded from the Allen Coral Atlas | 221,822 | `25da665fe83f223bb0820d20ce44d54f9f6d1417ec17e6df2a1b1d35ae5c19ea` |

**Obtaining these:**
- **Masks and class maps:** from the author.
- **Bathymetry:** download it from allencoralatlas.org for the Koh Man Nai area. A new download may differ from this checksum.

## 4. Label GeoJSONs (automatic labels; no reflectance)

All are EPSG:4326 polygons with `Class`, `date`, `spatial_block` and `train_eligible` properties. They come from the author's `Autolabel/` scripts (see [R1_LOO_DATE.md](R1_LOO_DATE.md), BLD label source), which are not ported. Obtain them from the author.

| File | Use | Size (bytes) | SHA-256 |
|---|---|---|---|
| `<ROOT>/Asset/stricklabel_Training_SpatialDefault.geojson` | Stage 2 BLD training | 7,270,602 | `d6135d530b35d9c30abeb8b6e0a02d222234a0f3c8cd06401fe99ef293b8b76d` |
| `<ROOT>/Asset/stricklabel_Testing_SpatialDefault.geojson` | Stage 2 BLD test, bleached anchor | 8,005,605 | `d26b023e1b322589eb1cf8899e87185260c71e862d0a956282cecd96fa1c7e4e` |
| `<ROOT>/Asset/stricklabel_Holdout_20240514.geojson` | 20240514 hold-out | 3,104,337 | `306ef02d8566e41782de5253e75e0d2c0ed21f17673894a2cd60a97b38958224` |
| `<ROOT>/Asset/Unified_Training_SpatialDefault.geojson` | Stage 2 HPY training and PLE candidates | 7,271,859 | `f2b638ed2eae71fec60ff96a69c82db6eea83c4ccf82af6b659e13819e1c8ee6` |
| `<ROOT>/Asset/Unified_Testing_SpatialDefault.geojson` | Stage 2 HPY test, healthy anchor | 8,006,905 | `e4be71bf71055c734dfd18a5d068b92154ba1e5f4f42d61f7514371cdaa1e3ee` |
| `<ROOT>/Asset/Unified_Holdout_20240514.geojson` | 20240514 hold-out | 3,104,767 | `024de2155a1dc709feabd60722df884ad2d5cb411cb7278dbbc4cc6af1743810` |
| `<ROOT>/Asset/TrainingStage1.geojson` | Stage 1 training (HC / SND / RBL; undated) | 29,861 | `6ed77ff5ecee3f521052ee3a011a512e08d6f1dc6bca24e2d96d2436a9089df5` |
| `<ROOT>/Asset/TestingStage1.geojson` | Stage 1 test | 9,489 | `a3b195dc910c79ac712a6b9b983e303c3ed36047d47f717f80773370604c2ae8` |

## 5. Frozen models (MATLAB TreeBagger objects)

Read with `mat-io` 1.0.2. They contain no imagery, but the saved training matrices (`PrivateX`) hold derived features. Obtain them from the author.

| File | Variable(s) | Original run | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| `<ROOT>/Results/stage2_hierarchical_curated_holdout_quicktest/run_20260408_192522/hierarchical_curated_holdout_output.mat` | `rfGate`, `rfStress`, `muHealthy`, `muBleached`, … | Stage 2 paper run | 3,536,791 | `4f42eb868756cbed0b7ccdfb647ee13ccb687ec37644f67babce94a500d222b6` |
| `<A3>/Phase_Stage1/run_20260401_170150/rf_stage1_model.mat` | `rf_model` | Stage 1 paper run | 10,830,404 | `5e856750054e4a0c4e1e72811f705d73f8d222b3e36d69f955fd94b80818a66a` |

## 6. Aggregate results (archived reports; used as references)

| File | Content | Size (bytes) | SHA-256 |
|---|---|---|---|
| `…/run_20260408_192522/hierarchical_curated_holdout_date_summary.csv` | Table 5 (per-date class fractions) | 796 | `a21a660622968196ac8d206aa894f87ad9cc00fc553adf5125b4ced97d6af6a0` |
| `…/run_20260408_192522/hierarchical_curated_holdout_report.txt` | Stage 2 metrics | 1,869 | `f330ab4df02c11d32332d8edb746b92e5020bbe7945ffb920baa59ee56ac366b` |
| `<A3>/Phase_Stage1/run_20260401_170150/stage1_report.txt` | Stage 1 OOB, test, LOO and transfer | 3,386 | `fd92406a45c6e0435b97fe0858745870f793c6659ba391eb4d3f16b840f2ded5` |

These hold only aggregate numbers. The values they contain are in the paper and in `docs/PY*.md`.

## Which commands need what

| Command | Needs |
|---|---|
| `pytest -m "not slow and not data"` (`tests/test_prepush_safeguards.py` and other data-free tests) | Nothing |
| `scripts/run_frozen_stage2.py`, `tests/test_parity_data.py`, `tests/test_regression_tile.py` | Glint Source stack, baseline mask, label GeoJSONs, Stage 2 output `.mat` (all author-supplied) |
| `scripts/run_py4_stage1.py`, `scripts/run_frozen_stage1.py`, `tests/test_stage1_parity.py` | Raw Source stack, Planet UDM2 scenes (`<RAW>`), Stage 1 model and report, Stage 1 GeoJSONs |
| `scripts/run_py5.py`, `tests/test_mask_wave.py` | Baseline mask, Phase 9 class maps, raw Source stack |
| `scripts/run_py6_*.py`, `run_py8_*.py`, `run_r1_*.py`, `tests/test_fixed_units.py` | Both Source stacks, baseline mask, label GeoJSONs, Stage 2 output. `run_r1_task2.py` also needs the ACA bathymetry and the Stage 1 report. |

**Every command except the data-free tests needs archived products that only the author can supply.** Planet scenes alone are not enough: the Source stacks, labels and frozen models were derived from them by steps whose original scripts or intermediate files are partly lost.

## Lost originals

These no longer exist; the port works around them as described in [PROVENANCE.md](PROVENANCE.md).
- **Phase 1 co-registered GeoTIFFs:** all dates.
- **Phase 5 feature files:** `FEATURES_<date>.tif` for all dates.
- **Phase 3 GLINT files:** 6 of the 14 `GLINT_<date>.tif` (the 8 that survive are listed above).
- **Stage 1 model `run_20260401_102559`:** it produced the Phase 9 class maps behind the baseline mask.
- **The 5 Apr 2026 Source stack builder script.**
