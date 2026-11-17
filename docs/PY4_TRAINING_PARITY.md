# Py4: Stage 2 training port, parity with the paper of record

Status 16 Nov 2026, compat mode (double-scaled stack, 1e-6 floor), same training pixels as the paper.
Reference = the paper's frozen forests evaluated in Python, which reproduce Table 5 and the confusion matrix exactly (Py3).

## Set-up
- **Training rows.** The gate uses 3,000 HPY + 3,000 stressed (PLE 1,173 + BLD 1,827); the stress model uses all 1,290 PLE + 2,028 BLD. MATLAB's `randperm` stream can't be reproduced in numpy, so compat mode recovers the exact rows, in order, from the training matrix saved inside the paper forests (`PrivateX`). All 9,318 rows are matched.
- **TreeBagger settings**, read from the saved `ModelParams`: Gini, MinParent 10, MinLeaf 5, all 21 predictors, no pruning, no surrogates, bootstrap of size n, equal weights, symmetric cost.
- **Two engines:**
  - `matlab_cart` (default, `coralremote/cart.py`): exhaustive CART with MATLAB's tie rule, so equal-gain splits go to the lowest predictor index. Cuts are single-precision midpoints.
  - `sklearn`: `RandomForestClassifier(250, gini, min_samples_split=10, min_samples_leaf=5, max_features=None)`.
- **Seeds:** five (0–4) per engine. MATLAB's own run was not reproducible either: it trained with `UseParallel` and no fixed stream.

## Results vs the plan's done-when criteria
| Criterion (plan §7.4) | matlab_cart (5 seeds) | sklearn (5 seeds) |
|---|---|---|
| OA within 1 pp | **Met**: 0.9739–0.9750 vs 0.9748 | **Met**: 0.9742–0.9752 |
| Pixel agreement with the paper maps ≥ 97% (mean over 14 dates) | **Met on average**: 96.6, 98.1, 96.9, 97.9, 98.0% | Not met: 93.3–94.3% |
| Every Table 5 cell within 1 pp | Not met: 50–51/56 per seed; 53/56 for the seed mean | Not met: 46–50/56 |
| Tree size (nodes, gate / stress); MATLAB 73.3 / 127.4 | 73.0–73.6 / 126.6–128.3 | not comparable |

Cells more than 1 pp off (range over seeds, Python minus paper, in pp):

| Date | matlab_cart | sklearn |
|---|---|---|
| 20240514 | PLE −9.3…+8.1, UNK −7.6…+8.7 | PLE +4.7…+12.6, UNK −12.9…−4.6 |
| 20250518 | HPY −4.4…−0.1, PLE −6.3…+23.5, UNK −20.5…+7.9 | **HPY −35.2…−30.3**, PLE +33.1…+50.2 |
| 20240616 / 20240620 | ≤ 2 pp | ≤ 4.5 pp |

## Interpretation
1. **The scikit-learn engine has a systematic bias on 20250518.** Every seed moves about 30 pp of the mask from HPY to PLE. The cause is split tie-breaking: scikit-learn breaks ties at random, MATLAB takes the first predictor. The difference shows only on out-of-distribution pixels. With MATLAB's tie rule, HPY on 20250518 is 40.0–44.3% (paper 44.4%).
2. **20250518 is out of distribution.** For 77.6% of its in-mask pixels, the healthy→bleached axis position lies outside every training value. The paper's own gate scores on that date fall mostly between 0.3 and 0.7, near the decision boundary.
3. **The PLE/UNK split on 20240514 and 20250518 depends on the bootstrap.** Across five seeds of the *same* MATLAB-rule algorithm, the seed-to-seed range is 17 pp (20240514) and 30 pp (20250518). It is ≤ 2.3 pp on every other date. Because those pixels' stress scores sit near the 0.55 UNK threshold, a different bootstrap moves them between PLE and UNK. A MATLAB re-run without a fixed stream would vary the same way. A 1 pp tolerance can't be met on these two dates by any retraining, MATLAB included.
4. **Everything else is stable.** HPY on every date, all cells on the other 12 dates, OA, kappa, per-class F1 (within 0.003) and the anchor recalls are all reproduced within tolerance.

**Implication for the paper's numbers (for the author):** the 20240514 and 20250518 PLE and UNK percentages in Table 5 are one draw from a spread of roughly ±10–15 pp. HPY on those dates, and every other cell, is robust. This holds independently of the scaling bug. P6/P7 (the fixed-units re-run) should report these dates as means over several seeds.

## Reproduce
```bash
cd coralremote-python
PYTHONPATH=. python scripts/run_py4_train.py --seeds 0 1 2 3 4                  # matlab_cart, ~75 s/seed
PYTHONPATH=. python scripts/run_py4_train.py --seeds 0 1 2 3 4 --engine sklearn  # ~7 s/seed
```

---

# Stage 1 (substrate CORAL / SND / RBL): Py2–Py4

Paper run: `ReefCode/.../A3/Phase_Stage1/run_20260401_170150` (archived copy; the original working folder no longer exists). Script: `Code/Train_Stage1_LogRatio.m`. Stage 1 read Phase 5 features built from **true** reflectance on 1 Apr, so the Stage 2 double-scaling bug doesn't affect it.

## Inputs rebuilt (Py2): `coralremote/preprocess.py`
The paper's `FEATURES_<date>.tif` files were lost, so they are rebuilt from code and the A3 config (`pipeline_logratio_config.mat`):
- **Glint-corrected reflectance:** the surviving `GLINT_<date>.tif` files, otherwise the Source glint stack × 10,000.
- **NaN masking**, in this order:
  - any raw band ≤ 0;
  - UDM2 `clear` = 0, read from the raw scenes (local archive);
  - wave pixels, where raw NIR > 0.07 inside the reef;
  - everything outside the Phase 2 reef mask.
- **Features:** B/G1, G1/G2, B/R, DII_BR = log B − 0.82 log R, DII_GR = log G1 − 0.80 log R, with EPS 1e-6.
- **Label rasteriser:** the Stage 1 script has its own rasteriser, with a half-pixel offset and clamping that differ from Stage 2's. It is ported separately as `rule="stage1"`.

**Check against the 17,372 training rows saved in the paper forest:** pixel counts and labels are identical (3,751 / 5,974 / 7,647). Ratio features are bit-identical where the GLINT tif survives and within 2.5e-7 otherwise; DII is within 7.2e-7, one float32 ulp of the log terms. The UDM2 mask matters on one date: 20230629, where it removes exactly the 120 light-haze test pixels.

## Frozen model (Py3): exact
| Item | Result |
|---|---|
| OOB error (from the saved per-tree OOB masks) | 0.1004 = report |
| Secondary test (CORAL 20240408, SND/RBL 20230519; n = 905) | Confusion matrix identical; OA 0.7414, κ 0.5740 |
| Cross-date transfer table (8 dates: n, OA, κ, CORAL P/R/F1) | All identical |

**Prediction rule.** With the asymmetric cost `[0 2 2; 1 0 1; 1 1 0]`, MATLAB grows the trees under cost-adjusted priors, so the stored leaf probabilities already include the cost. Applying the cost again at prediction gives OOB 0.1041 and the wrong confusion matrix. `Forest.predict` maps the scores back to the original prior and then minimises expected cost. For these cost matrices that equals the argmax of the stored scores.

## Retraining (Py4): MATLAB-rule CART with cost-adjusted class weights (3 seeds)
| Metric | Paper | Python (seeds 0–2) |
|---|---|---|
| OOB error | 0.1004 | 0.1002–0.1016 |
| Secondary OA / κ | 0.7414 / 0.574 | 0.7425–0.7448 / 0.576–0.580 |
| Per-pixel agreement with the paper forest (reef, 8 dates) | — | 97.7–97.8% (worst date 96.9%) |
| Transfer-table OA | — | within 0.01 of the paper on 6 dates; up to 0.03 on 20240317 and 20241111 |
| Mean tree size (nodes) | 1,540 | 1,683–1,685 (+9%, open) |

**Not ported (all retrain or are random):** the LOO CORAL-polygon check, the spatial-block holdout (100 trees each) and the k-means diagnostic. The tree-size gap is open. Stage 2 tree sizes match, so the likely cause is a detail of how MATLAB handles weighted (cost-adjusted) nodes, possibly its stopping rule. It doesn't change any reported metric beyond seed noise.

```bash
PYTHONPATH=. python scripts/run_frozen_stage1.py         # exact parity, ~30 s
PYTHONPATH=. python scripts/run_py4_stage1.py --seeds 0 1 2   # ~3 min/seed
```
