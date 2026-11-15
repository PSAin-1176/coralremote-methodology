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
cd E:/CoralRemote/coralremote-python
PYTHONPATH=. python scripts/run_py4_train.py --seeds 0 1 2 3 4                  # matlab_cart, ~75 s/seed
PYTHONPATH=. python scripts/run_py4_train.py --seeds 0 1 2 3 4 --engine sklearn  # ~7 s/seed
```
