# Provenance: MATLAB scripts of the AJSTR paper → Python

The paper's runs used MATLAB scripts that were never in this repository. They live in the author's working folders: `Code/`, `Autolabel/`, `done/baseline_mask/` and the A3 pipeline. The table below shows which Python module reproduces each script and how closely.

| Paper step (MATLAB) | Saved run checked against | Python | Agreement |
|---|---|---|---|
| `Phase3_LogRatioNormalize.m`, `Phase4_WaveHandling.m`, `Phase5_FeatureStack_LogRatio.m`, `buildFeatureStack_LogRatio.m` (A3) | 17,372 training rows saved in the Stage 1 forest | `preprocess` | Pixel counts and labels identical; features within 1 float32 ulp |
| `Train_Stage1_LogRatio.m` (A3 `Phase_Stage1/run_20260401_170150`) | `rf_stage1_model.mat`, `stage1_report.txt` | `stage1`, `forest`, `cart` | Frozen forest exact; retraining within seed noise |
| `Baseline_Coral_LogRatio.m` (A3 `Baseline_Coral/run_20260401_103612`) | `Asset/baseline_coral_mask.mat` (`vote_map`, `tier_map`, `mask0`) | `mask.consensus` | 0 px |
| `Clean_CoralMask_Conservative.m`, `Peel_CoralMask_OnePixel.m`, `Finalize_CoralMask_MinorFixes.m` (6 Apr) | `maskClean`, `maskPeel`, `maskFinal` | `mask.clean_conservative`, `peel`, `finalize` | 0 px at every step |
| `Run_Stage2_Hierarchical_CuratedHoldout_QuickTest.m` (`run_20260408_192522`) | output `.mat`, report, date-summary CSV | `features`, `stage2`, `train` | Frozen forests: Table 5 and confusion matrix exact; retraining within seed noise |
| `Diag_WaveFFT_ThreeClass.m` | `Wave_Diagnostic_Methodology_and_Results_Note.txt` | `wave` | 14/14 verdicts |
| `Run_Autolabel_*`, `Build_HighConf_*`, `Build_HPU/BPU_Polygons.m`, `Build_StrickLabel_Dataset.m`, `Build_Unified_Polygon_Dataset.m` | label GeoJSONs in `Asset/` | not ported; the frozen GeoJSONs are used | Audited: every threshold is relative, so the labels don't depend on scale ([PY6](PY6_FIX_AND_RERUN.md#p5-labels-no-relabelling-needed)) |
| 5 Apr Source stack builder (lost) | `SourceStack_manifest.txt` | `units.FIXED` (load-time correction) | Glint NIR equals raw NIR within 1 ulp |
| Phase 1 co-registration (AROSICS) | not re-run | not ported (needs `arosics`) | — |

**Inputs that no longer exist.** Most of the original `C:\SAT` workspace is gone:
- the Phase 1 co-registered GeoTIFFs;
- the Phase 5 `FEATURES_<date>.tif` files;
- 6 of the 14 Phase 3 `GLINT_<date>.tif` files;
- the Stage 1 model `run_20260401_102559` that produced the Phase 9 class maps behind the baseline mask.

The Python port rebuilds what it needs from the Source stacks and the surviving copies. The reports in `docs/PY*.md` describe each case.
