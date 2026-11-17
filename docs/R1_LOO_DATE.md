# R1 (Reviewer D): leave-one-June-date-out test and supporting numbers

Branch `r1-loo-date` (off `python-port`).

**Scripts:**
- `scripts/run_r1_loo_date.py`: Task 1 runs.
- `scripts/run_r1_report.py`: Task 1 tables.
- `scripts/run_r1_task2.py`: Task 2 (a), (c), (d) and (e).
- `scripts/run_r1_uncurated_ple.py`: add-on (g).

**Outputs:** raw JSON outputs stay in the local cache folder. This file holds aggregates only, with no imagery and no spectra.

## Task 1: leave-one-June-date-out (Stage 2)

### Setup

**Recipe:** the paper's, for every run.
- **Engine:** MATLAB-rule CART, 250 trees, MinLeaf 5, MinParent 10, all 21 predictors.
- **Gate:** 3,000 HPY + up to 3,000 stressed rows (numpy draw).
- **Stress model:** all PLE and BLD (cap 2,500 each).
- **UNK:** if the max score is < 0.55 in either model.
- **Seeds:** 0–4.
- **Pixels:** classification covers all 20,330 mask pixels; every June date has 20,330 valid pixels.
- **Baseline:** the same recipe and seeds with nothing removed. It is the Py6 retraining; its seed 0 is identical to Py6 seed 0. Effects (Δ) are differences from this baseline.

**Settings, for each held-out date D:**
- **"Rows of D removed" (`rows_D`):** every training pixel from D is removed, in every class and in both models. The features stay as in the paper: the healthy and bleached class means, the axis, the spectral-angle references and the PLE curation band were all derived with D included.
- **"D removed including axis and PLE curation" (`strict_D`):** D's training polygons are removed before anything is computed. The healthy and bleached means, the axis (so `axis_t`, the residual and both spectral angles), the curation thresholds and the PLE curation are all re-derived without D. **This is the version with no leakage through the features.**
  - The healthy mean does not change, because no HPY training pixel is from June.
  - The bleached mean moves by up to 0.009.
- **The paper's PLE curation band** (0.736 < t < 0.915, residual ≤ 0.0956) **was derived using all three June dates.** Its bleached-side threshold comes from the June BLD pixels. In `strict_D` it is re-derived and the curated PLE set changes; for example, `strict_20240620` keeps 121 PLE polygons, including some from 13 and 16 June.

**"All three June dates removed":**
- **No BLD left:** every one of the 2,028 BLD training pixels is from June (692 / 684 / 652), so no BLD training pixels remain and the stress model cannot be trained.
- **What runs instead:** the gate alone, HPY against the 538 remaining PLE pixels (all from 20230509). Its HPY / stressed / UNK split is compared with the baseline gate.
- **Labels:** the June dates are also excluded from the training labels. They still enter feature construction through the paper's axis, which cannot be built without BLD pixels.

**Modes:** compat (the paper's units) is the primary result; fixed units are the check.

**Training pixels by date** are shown in the last table of each mode below.

### Key result (compat, mean over 5 seeds; ranges in the tables)

| Date | Frozen paper | Baseline BLD | BLD with D's rows removed | BLD with D removed incl. features | Stressed (PLE + BLD) with D removed | Stressed with all June removed (gate) |
|---|---|---|---|---|---|---|
| 20240613 | 94.5 | 94.0 | 89.3 | 89.9 | 96.4–96.5 | 99.5 |
| 20240616 | 93.7 | 92.9 | 37.0 | 41.4 | 88.9–90.4 | 99.6 |
| 20240620 | 37.2 | 36.7 | 82.4 | 10.2 | 95.8–97.5 | 98.6 |

**Wave controls.** Two wave-affected non-bleaching dates come out healthy in every setting:
- **20230519 (WAVE: 8/8 roughness flags, 7/8 AI flags):** 97–99% HPY in the baseline, the rows variants and the all-June gate. It falls to 85–90% HPY in two `strict` variants.
- **20230629 (TURBULENT):** 100% HPY.

### Tables

#### compat: held-out date classified (mean [min–max] over 5 seeds, % of mask pixels)

| Held-out date | Class | Frozen paper model | Baseline (nothing removed) | Rows of D removed | Δ | D removed incl. axis + PLE curation | Δ |
|---|---|---|---|---|---|---|---|
| 20240613 | HPY | 0.3 | 0.4 [0.2–0.9] | 0.3 [0.3–0.4] | -0.1 | 0.4 [0.3–0.8] | -0.0 |
|  | PLE | 3.9 | 4.1 [4.0–4.2] | 7.3 [6.8–7.5] | +3.2 | 6.5 [6.1–6.8] | +2.4 |
|  | BLD | 94.5 | 94.0 [93.5–94.3] | 89.3 [88.9–89.7] | -4.7 | 89.9 [89.6–90.3] | -4.1 |
|  | UNK | 1.3 | 1.5 [1.2–1.6] | 3.1 [2.9–3.2] | +1.6 | 3.2 [2.9–3.6] | +1.7 |
| 20240616 | HPY | 0.4 | 0.4 [0.2–0.7] | 0.4 [0.2–0.5] | -0.0 | 0.4 [0.4–0.5] | +0.1 |
|  | PLE | 2.3 | 2.6 [2.3–3.0] | 51.9 [50.4–53.6] | +49.3 | 49.0 [47.9–51.1] | +46.4 |
|  | BLD | 93.7 | 92.9 [91.6–94.0] | 37.0 [34.9–38.1] | -55.9 | 41.4 [39.7–43.0] | -51.5 |
|  | UNK | 3.5 | 4.1 [3.4–5.1] | 10.8 [8.4–11.7] | +6.6 | 9.1 [8.7–9.6] | +5.0 |
| 20240620 | HPY | 1.2 | 1.0 [0.7–1.4] | 1.1 [0.8–1.5] | +0.1 | 0.9 [0.7–1.1] | -0.1 |
|  | PLE | 55.6 | 56.3 [56.0–56.6] | 15.1 [15.0–15.3] | -41.2 | 85.6 [85.3–86.0] | +29.3 |
|  | BLD | 37.2 | 36.7 [36.2–37.1] | 82.4 [80.6–83.6] | +45.7 | 10.2 [9.9–10.4] | -26.5 |
|  | UNK | 6.0 | 5.9 [5.8–6.2] | 1.4 [0.2–3.0] | -4.5 | 3.2 [2.9–3.5] | -2.7 |

#### compat: BLD on every June date when one date is removed (mean over seeds)

| Setting | BLD 20240613 | BLD 20240616 | BLD 20240620 | PLE 20240514 | UNK 20240514 |
|---|---|---|---|---|---|
| baseline | 94.0 | 92.9 | 36.7 | 25.9 | 15.4 |
| rows_20240613 | 89.3 | 93.0 | 38.5 | 36.8 | 4.3 |
| rows_20240616 | 89.0 | 37.0 | 34.1 | 21.5 | 19.3 |
| rows_20240620 | 98.0 | 99.5 | 82.4 | 23.8 | 2.3 |
| strict_20240613 | 89.9 | 91.3 | 42.6 | 38.1 | 16.7 |
| strict_20240616 | 76.1 | 41.4 | 25.8 | 58.4 | 5.1 |
| strict_20240620 | 35.7 | 26.6 | 10.2 | 79.9 | 0.1 |

#### compat: test BLD polygons on the held-out date (prediction shares, mean over seeds)

| Date | n test BLD px | Frozen paper BLD recall | Baseline BLD / PLE / UNK | Rows of D removed BLD / PLE / UNK | Strict BLD / PLE / UNK |
|---|---|---|---|---|---|
| 20240613 | 434 | 100.0 | 100.0 / 0.0 / 0.0 | 100.0 / 0.0 / 0.0 | 100.0 / 0.0 / 0.0 |
| 20240616 | 434 | 100.0 | 100.0 / 0.0 / 0.0 | 96.9 / 1.5 / 1.6 | 98.2 / 1.4 / 0.4 |
| 20240620 | 411 | 94.9 | 95.3 / 3.5 / 1.3 | 100.0 / 0.0 / 0.0 | 69.9 / 26.8 / 3.3 |

#### compat: test-set OA (mean [min–max])

| Setting | OA | F1 HPY / PLE / BLD |
|---|---|---|
| baseline | 0.9728 [0.9702–0.9759] | 0.983 / 0.845 / 0.976 |
| rows_20240613 | 0.9717 [0.9687–0.9752] | 0.983 / 0.842 / 0.973 |
| rows_20240616 | 0.9723 [0.9697–0.9756] | 0.982 / 0.848 / 0.971 |
| rows_20240620 | 0.9512 [0.9473–0.9537] | 0.982 / 0.700 / 0.912 |
| strict_20240613 | 0.9696 [0.9678–0.9717] | 0.981 / 0.809 / 0.978 |
| strict_20240616 | 0.9522 [0.9502–0.9534] | 0.969 / 0.649 / 0.979 |
| strict_20240620 | 0.9558 [0.9523–0.9586] | 0.979 / 0.812 / 0.917 |

#### compat: all three June dates removed (gate only; % HPY / stressed / UNK, mean [min–max])

| Date | Baseline gate | All June removed | Δ stressed |
|---|---|---|---|
| 20240613 | 0.4 [0.2–0.9] / 99.4 [98.8–99.7] / 0.2 [0.1–0.3] | 0.3 [0.2–0.4] / 99.5 [99.2–99.7] / 0.2 [0.1–0.4] | +0.1 |
| 20240616 | 0.4 [0.2–0.7] / 99.5 [99.1–99.7] / 0.2 [0.1–0.4] | 0.3 [0.1–0.4] / 99.6 [99.3–99.9] / 0.1 [0.0–0.3] | +0.1 |
| 20240620 | 1.0 [0.7–1.4] / 98.6 [98.2–98.9] / 0.4 [0.2–0.7] | 1.0 [0.6–1.2] / 98.6 [98.4–99.0] / 0.4 [0.3–0.6] | +0.0 |
| 20240514 | 58.7 [54.2–60.6] / 39.5 [38.1–42.7] / 1.8 [0.8–3.1] | 56.5 [55.2–58.2] / 42.1 [40.8–43.4] / 1.5 [1.0–2.0] | +2.6 |

Test BLD polygons under the all-June gate (HPY / stressed / UNK %): {'20240613': [0.0, 100.0, 0.0], '20240616': [0.0, 100.0, 0.0], '20240620': [0.0, 100.0, 0.0]}

#### compat: training pixels by date

| Setting | BLD pool by date | PLE pool by date | Gate stressed rows drawn (seed 0): BLD by date; PLE by date | Stress model rows PLE / BLD |
|---|---|---|---|---|
| baseline | {'20240613': 692, '20240616': 684, '20240620': 652} (n=2028) | {'20230509': 538, '20240620': 752} (n=1290) | {'20240613': 625, '20240616': 619, '20240620': 588}; {'20230509': 493, '20240620': 675} (HPY 3000, stressed 3000) | 1290 / 2028 |
| rows_20240613 | {'20240616': 684, '20240620': 652} (n=1336) | {'20230509': 538, '20240620': 752} (n=1290) | {'20240616': 684, '20240620': 652}; {'20230509': 538, '20240620': 752} (HPY 3000, stressed 2626) | 1290 / 1336 |
| rows_20240616 | {'20240613': 692, '20240620': 652} (n=1344) | {'20230509': 538, '20240620': 752} (n=1290) | {'20240613': 692, '20240620': 652}; {'20230509': 538, '20240620': 752} (HPY 3000, stressed 2634) | 1290 / 1344 |
| rows_20240620 | {'20240613': 692, '20240616': 684} (n=1376) | {'20230509': 538} (n=538) | {'20240613': 692, '20240616': 684}; {'20230509': 538} (HPY 3000, stressed 1914) | 538 / 1376 |
| strict_20240613 | {'20240616': 684, '20240620': 652} (n=1336) | {'20230509': 488, '20240620': 474} (n=962) | {'20240616': 684, '20240620': 652}; {'20230509': 488, '20240620': 474} (HPY 3000, stressed 2298) | 962 / 1336 |
| strict_20240616 | {'20240613': 692, '20240620': 652} (n=1344) | {'20230404': 19, '20230509': 446, '20240620': 483} (n=948) | {'20240613': 692, '20240620': 652}; {'20230404': 19, '20230509': 446, '20240620': 483} (HPY 3000, stressed 2292) | 948 / 1344 |
| strict_20240620 | {'20240613': 692, '20240616': 684} (n=1376) | {'20230509': 764, '20240613': 324, '20240616': 643} (n=1731) | {'20240613': 670, '20240616': 665}; {'20230509': 735, '20240613': 313, '20240616': 617} (HPY 3000, stressed 3000) | 1731 / 1376 |
| alljune | none (n=0) | {'20230509': 538} (n=538) | none; {'20230509': 538} (HPY 3000, stressed 538) | 538 / 0 |

compat strict_20240613: axis means max shift mu_h 0.00e+00, mu_b 6.58e-03; PLE polygons kept 69

compat strict_20240616: axis means max shift mu_h 0.00e+00, mu_b 8.56e-03; PLE polygons kept 67

compat strict_20240620: axis means max shift mu_h 0.00e+00, mu_b 6.98e-03; PLE polygons kept 121

#### fixed: held-out date classified (mean [min–max] over 5 seeds, % of mask pixels)

| Held-out date | Class | Frozen paper model | Baseline (nothing removed) | Rows of D removed | Δ | D removed incl. axis + PLE curation | Δ |
|---|---|---|---|---|---|---|---|
| 20240613 | HPY | – | 0.4 [0.3–0.5] | 0.4 [0.3–0.4] | -0.0 | 0.4 [0.3–1.0] | +0.0 |
|  | PLE | – | 4.3 [4.1–4.7] | 6.2 [6.0–6.4] | +1.9 | 6.5 [6.1–7.0] | +2.2 |
|  | BLD | – | 94.0 [93.8–94.2] | 91.1 [90.8–91.5] | -2.9 | 90.4 [88.1–91.3] | -3.6 |
|  | UNK | – | 1.4 [1.2–1.5] | 2.4 [2.1–2.7] | +1.0 | 2.7 [2.1–4.3] | +1.3 |
| 20240616 | HPY | – | 0.4 [0.3–0.7] | 0.4 [0.2–0.6] | -0.0 | 0.4 [0.4–0.4] | -0.0 |
|  | PLE | – | 3.4 [2.1–4.0] | 50.9 [50.0–51.7] | +47.5 | 48.4 [46.3–50.5] | +45.0 |
|  | BLD | – | 90.4 [89.8–92.4] | 40.7 [39.6–42.4] | -49.7 | 46.1 [44.6–47.2] | -44.3 |
|  | UNK | – | 5.8 [5.1–6.2] | 8.0 [6.9–8.6] | +2.2 | 5.1 [4.5–6.1] | -0.7 |
| 20240620 | HPY | – | 1.2 [0.9–1.5] | 1.2 [1.0–1.4] | +0.0 | 1.1 [0.9–1.5] | -0.1 |
|  | PLE | – | 52.6 [51.7–53.0] | 15.0 [14.9–15.1] | -37.6 | 85.6 [84.9–86.2] | +33.0 |
|  | BLD | – | 41.4 [41.0–41.8] | 82.3 [80.7–83.4] | +41.0 | 10.3 [10.1–10.5] | -31.1 |
|  | UNK | – | 4.9 [4.6–5.1] | 1.5 [0.4–3.1] | -3.4 | 3.0 [2.7–3.5] | -1.8 |

#### fixed: BLD on every June date when one date is removed (mean over seeds)

| Setting | BLD 20240613 | BLD 20240616 | BLD 20240620 | PLE 20240514 | UNK 20240514 |
|---|---|---|---|---|---|
| baseline | 94.0 | 90.4 | 41.4 | 35.0 | 6.6 |
| rows_20240613 | 91.1 | 90.5 | 42.6 | 33.8 | 8.1 |
| rows_20240616 | 90.5 | 40.7 | 37.5 | 39.5 | 2.9 |
| rows_20240620 | 98.0 | 99.4 | 82.3 | 24.2 | 2.4 |
| strict_20240613 | 90.4 | 98.6 | 51.5 | 35.6 | 15.7 |
| strict_20240616 | 70.3 | 46.1 | 26.8 | 60.3 | 4.1 |
| strict_20240620 | 35.7 | 26.6 | 10.3 | 79.5 | 0.5 |

#### fixed: test BLD polygons on the held-out date (prediction shares, mean over seeds)

| Date | n test BLD px | Frozen paper BLD recall | Baseline BLD / PLE / UNK | Rows of D removed BLD / PLE / UNK | Strict BLD / PLE / UNK |
|---|---|---|---|---|---|
| 20240613 | 434 | – | 100.0 / 0.0 / 0.0 | 100.0 / 0.0 / 0.0 | 100.0 / 0.0 / 0.0 |
| 20240616 | 434 | – | 100.0 / 0.0 / 0.0 | 97.8 / 1.7 / 0.5 | 98.4 / 1.4 / 0.2 |
| 20240620 | 411 | – | 97.2 / 2.4 / 0.4 | 100.0 / 0.0 / 0.0 | 69.7 / 26.8 / 3.5 |

#### fixed: test-set OA (mean [min–max])

| Setting | OA | F1 HPY / PLE / BLD |
|---|---|---|
| baseline | 0.9709 [0.9670–0.9730] | 0.981 / 0.835 / 0.976 |
| rows_20240613 | 0.9703 [0.9667–0.9732] | 0.981 / 0.832 / 0.973 |
| rows_20240616 | 0.9719 [0.9686–0.9745] | 0.981 / 0.839 / 0.977 |
| rows_20240620 | 0.9544 [0.9515–0.9578] | 0.981 / 0.713 / 0.923 |
| strict_20240613 | 0.9688 [0.9670–0.9701] | 0.980 / 0.763 / 0.977 |
| strict_20240616 | 0.9617 [0.9604–0.9634] | 0.976 / 0.674 / 0.983 |
| strict_20240620 | 0.9541 [0.9500–0.9577] | 0.976 / 0.809 / 0.916 |

#### fixed: all three June dates removed (gate only; % HPY / stressed / UNK, mean [min–max])

| Date | Baseline gate | All June removed | Δ stressed |
|---|---|---|---|
| 20240613 | 0.4 [0.3–0.5] / 99.5 [99.3–99.6] / 0.1 [0.0–0.2] | 0.4 [0.2–0.4] / 99.4 [99.2–99.6] / 0.2 [0.1–0.4] | -0.1 |
| 20240616 | 0.4 [0.3–0.7] / 99.4 [98.9–99.6] / 0.1 [0.0–0.4] | 0.4 [0.2–0.5] / 99.3 [99.0–99.5] / 0.2 [0.1–0.6] | -0.1 |
| 20240620 | 1.2 [0.9–1.5] / 98.5 [98.2–98.7] / 0.4 [0.3–0.5] | 1.2 [0.9–1.3] / 98.3 [98.1–98.6] / 0.5 [0.5–0.6] | -0.1 |
| 20240514 | 58.4 [56.6–59.5] / 39.5 [38.3–40.5] / 2.1 [0.7–2.9] | 54.9 [54.7–55.2] / 43.6 [42.1–44.5] / 1.5 [0.6–2.7] | +4.1 |

Test BLD polygons under the all-June gate (HPY / stressed / UNK %): {'20240613': [0.0, 100.0, 0.0], '20240616': [0.0, 100.0, 0.0], '20240620': [0.0, 100.0, 0.0]}

#### fixed: training pixels by date

| Setting | BLD pool by date | PLE pool by date | Gate stressed rows drawn (seed 0): BLD by date; PLE by date | Stress model rows PLE / BLD |
|---|---|---|---|---|
| baseline | {'20240613': 692, '20240616': 684, '20240620': 652} (n=2028) | {'20230509': 570, '20240620': 593} (n=1163) | {'20240613': 646, '20240616': 643, '20240620': 606}; {'20230509': 544, '20240620': 561} (HPY 3000, stressed 3000) | 1163 / 2028 |
| rows_20240613 | {'20240616': 684, '20240620': 652} (n=1336) | {'20230509': 570, '20240620': 593} (n=1163) | {'20240616': 684, '20240620': 652}; {'20230509': 570, '20240620': 593} (HPY 3000, stressed 2499) | 1163 / 1336 |
| rows_20240616 | {'20240613': 692, '20240620': 652} (n=1344) | {'20230509': 570, '20240620': 593} (n=1163) | {'20240613': 692, '20240620': 652}; {'20230509': 570, '20240620': 593} (HPY 3000, stressed 2507) | 1163 / 1344 |
| rows_20240620 | {'20240613': 692, '20240616': 684} (n=1376) | {'20230509': 570} (n=570) | {'20240613': 692, '20240616': 684}; {'20230509': 570} (HPY 3000, stressed 1946) | 570 / 1376 |
| strict_20240613 | {'20240616': 684, '20240620': 652} (n=1336) | {'20230509': 493, '20240620': 338} (n=831) | {'20240616': 684, '20240620': 652}; {'20230509': 493, '20240620': 338} (HPY 3000, stressed 2167) | 831 / 1336 |
| strict_20240616 | {'20240613': 692, '20240620': 652} (n=1344) | {'20230404': 19, '20230509': 431, '20240620': 410} (n=860) | {'20240613': 692, '20240620': 652}; {'20230404': 19, '20230509': 431, '20240620': 410} (HPY 3000, stressed 2204) | 860 / 1344 |
| strict_20240620 | {'20240613': 692, '20240616': 684} (n=1376) | {'20230509': 825, '20240613': 397, '20240616': 600} (n=1822) | {'20240613': 651, '20240616': 645}; {'20230509': 781, '20240613': 371, '20240616': 552} (HPY 3000, stressed 3000) | 1822 / 1376 |
| alljune | none (n=0) | {'20230509': 570} (n=570) | none; {'20230509': 570} (HPY 3000, stressed 570) | 570 / 0 |

fixed strict_20240613: axis means max shift mu_h 0.00e+00, mu_b 6.58e-03; PLE polygons kept 59

fixed strict_20240616: axis means max shift mu_h 0.00e+00, mu_b 8.56e-03; PLE polygons kept 62

fixed strict_20240620: axis means max shift mu_h 0.00e+00, mu_b 6.98e-03; PLE polygons kept 126

#### Frozen paper model (Table 5 reference), % HPY / PLE / BLD / UNK

- 20240613: 0.33 / 3.91 / 94.46 / 1.30
- 20240616: 0.40 / 2.34 / 93.74 / 3.53
- 20240620: 1.16 / 55.62 / 37.23 / 5.99
- 20240514: 58.98 / 22.79 / 0.01 / 18.22
frozen test OA check: 0.9748406943528894

### Reading

**Do the June BLD fractions hold up when their own date is excluded?** Partly. The answer has two layers.

1. **"Stressed" (PLE + BLD, the gate's call) holds up on every June date, in every setting.**
   - **One date removed:** 88–99% of each June date is still stressed, whether only the rows are removed or the axis and curation are re-derived too.
   - **All three removed:** a gate that has never seen a June training pixel calls 98.6–99.6% of each June date stressed (99.3–99.6% on 13 and 16 June). It also calls 100% of the June test BLD polygons stressed.
   - **Waves don't explain it:** the two wave-affected non-bleaching dates stay healthy (20230519 85–99% HPY, 20230629 100% HPY).
   - **Conclusion:** the departure from healthy on the June dates is not an artefact of training on those dates, and wave state alone does not produce it.
2. **The BLD / PLE split within the stressed reef holds up for 13 June only.**
   - **20240613:** BLD 94.0% → 89.3% (rows removed) or 89.9% (features re-derived). That is −4 to −5 pp, with the remainder going to PLE and UNK. Its held-out test BLD polygons stay 100% BLD.
   - **20240616: does not hold up.** BLD 92.9% → 37.0% (rows) or 41.4% (re-derived), and about half the reef becomes PLE. Its test BLD polygons are still 97–98% BLD. So the brightest, label-defining patches stay BLD, but the bulk of the reef on that date sits between the PLE and BLD training clusters. Which side it falls on depends on whether its own pixels define the BLD cluster.
   - **20240620: not identifiable.** BLD 36.7% → 82.4% when only its rows are removed, and → 10.2% when the curation is re-derived.
     - **Rows removed:** this also removes 752 of the 1,290 PLE pixels (all PLE except 20230509), leaving PLE defined only by 20230509.
     - **Re-derived:** the new PLE set includes curated pale polygons from 13 and 16 June.
     - **Range:** the date's BLD fraction spans 10–82% depending on which pale examples exist. Its own test polygons go from 95% to 70% BLD.
   - **Removing one date also moves the others.** `strict_20240620` drops 13 June BLD to 36% and 16 June BLD to 27%, because the PLE class then includes June pale pixels. The PLE / BLD boundary is set by which pale and bleached examples are in training, not by the image alone.
3. **Fixed units give the same picture** (13 June 90–91%, 16 June 41–46%, 20 June 10–82%, all-June gate 98–99%).

**For the response letter:**
- **Robust claim:** "Coral on the three June dates departs strongly from its healthy state (≥ 98% stressed). This holds even when every June pixel is excluded from training. Wave-affected dates without bleaching are classified healthy."
- **Not robust:** the specific ~94% BLD on 16 June (and the PLE / BLD split on 20 June). It depends on the June-derived BLD and PLE training examples. Those examples are themselves automatic labels from a prior that pools all three June dates (see the next section), so the test cannot separate the date from its labels.
- **13 June:** the BLD fraction survives the test (about 90%), and the anchor BLD recall of 1.000 survives too (100% when 13 June is held out).
- **Suggested wording:** present June as "stressed, ≥ 98%", with the BLD / PLE split as label-dependent: a range, or 13 June only.


## Where the BLD training labels come from

**The BLD polygons are automatic labels, not field observations.**
- **Scripts:** `Autolabel/` (local copies of the original autolabel scripts; not in this repository). The driver is `Build_StrickLabel_Dataset.m`, which produced the "stricklabel" polygons used for Stage 2 BLD.
- **Step 1: one prior over all three June dates** (`Build_HighConf_BleachedMask.m`).
  - **Inputs:** for each pixel inside the interior of the coral mask, the script takes the median over 20240613, 20240616 and 20240620 of six spectral quantities: brightness, Y/B, B/R, B/G1, DII_BR, and 2·Blue/(Red + Green1). Each is robustly scaled to 0–1.
  - **Score:** 0.28 brightness + 0.20 Y/B + 0.18 B/R + 0.10 B/G1 + 0.14 DII_BR + 0.04 blue-pale + 0.06 temporal stability. Stability is 1 minus the scaled across-date MADs.
  - **Prior:** the top 16% of the score forms the "strict" prior; the top 32% forms the "broad" prior. Both are cleaned of small objects and holes.
- **Step 2: superpixel labelling on each June date** (`Run_Autolabel_BleachedPerDate.m`, stricklabel settings).
  - **Segments:** about 2,600 superpixels (compactness 8) per date.
  - **Bleached** if the segment has ≥ 12 pixels and either:
    - ≥ 70% of it lies in the strict prior and its spectral distance from the prior's centre is ≤ 3.5; or
    - ≥ 55% of it lies in the (here also strict) prior, the centre distance is ≤ 4.5 and its brightness is ≥ 1.0 robust z above the prior's.
  - **Centre distance:** a weighted robust-z distance of the segment medians of B/G1, Y/B, B/R, DII_BR and blue-pale from the prior's medians.
  - **Other classes:** some remaining segments become Pale; small boundary segments become Uncertain or EdgeIgnore.
- **Step 3: polygons and splits** (`Build_BPU_Polygons.m`, then `Build_Unified_Polygon_Dataset.m`). Segments become polygons with a 4 × 4 spatial-block id. Fixed blocks are set aside for testing. Every class except EdgeIgnore is train-eligible.

**What this means for the leave-one-date-out test:**
- **No date-specific step:** the BLD labels on any one June date are drawn where the three-date median score is in the top 16%, and refined by that date's own spectra.
- **What the test can and can't remove:** leaving out date D removes D's pixels from training, but the labels on the other two dates were still placed using a prior that includes D. The test removes D from model fitting and from feature construction. It cannot remove D from label construction.
- **The rule itself:** it defines "bleached" as the brightest, least brown coral pixels of the June period. No sea-state check is involved.

## Task 2: existing numbers

### (a) Surface-state (wave) diagnostic, all 14 dates

**Source:** the Py5 port of `Diag_WaveFFT_ThreeClass.m`, run on the co-registered raw reflectance. All 14 verdicts match the paper.

**How the numbers are made:**
- **Tiles:** 128 × 128 pixels, at least 80% open water (NDWI > 0.05).
- **AI** = 1 − (min ÷ max) of the FFT power over 12 direction sectors.
- **WBP** = mean power in the 6–96 m wavelength annulus.
- **Per date:** each band's value is the median over tiles. The table shows the median across the 8 bands.
- **Flags:** a band is flagged when its value exceeds the cross-date median + 1.5 × MAD for that band.
- **Verdict:** WAVE if ≥ 4 WBP (roughness) flags and ≥ 6 AI flags; TURBULENT if ≥ 4 roughness flags; otherwise CALM.
- **Confidence:** set by the number of tiles: HIGH ≥ 60, MEDIUM ≥ 40, LOW below that.

| Date | Tiles | AI (median of 8 bands) | WBP (median of 8 bands) | Roughness flags | AI flags | Class | Confidence |
|---|---|---|---|---|---|---|---|
| 20230404 | 76 | 0.280 | 7.07e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20230509 | 76 | 0.331 | 7.41e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20230519 | 75 | 0.693 | 2.51e-2 | 8/8 | 7/8 | WAVE | HIGH |
| 20230629 | 38 | 0.525 | 2.07e-2 | 7/8 | 3/8 | TURBULENT | LOW |
| 20240307 | 76 | 0.293 | 5.33e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20240317 | 75 | 0.291 | 5.32e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20240408 | 60 | 0.254 | 1.26e-2 | 2/8 | 0/8 | CALM | HIGH |
| 20240514 | 76 | 0.564 | 1.12e-2 | 0/8 | 6/8 | CALM | HIGH |
| 20240613 | 76 | 0.829 | 2.67e-2 | 6/8 | 8/8 | WAVE | HIGH |
| 20240616 | 76 | 0.796 | 6.55e-2 | 8/8 | 8/8 | WAVE | HIGH |
| 20240620 | 76 | 0.900 | 8.13e-2 | 8/8 | 8/8 | WAVE | HIGH |
| 20240926 | 76 | 0.326 | 2.46e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20241111 | 76 | 0.309 | 3.58e-3 | 0/8 | 0/8 | CALM | HIGH |
| 20250518 | 76 | 0.402 | 5.87e-3 | 0/8 | 0/8 | CALM | HIGH |

MAD-based thresholds per band. A band is flagged when its value is above the threshold (median + 1.5 × MAD across the 14 dates).

| Band | AI median | AI MAD | AI threshold | WBP median | WBP MAD | WBP threshold |
|---|---|---|---|---|---|---|
| Coastal | 0.3142 | 0.0427 | 0.3783 | 2.463e-2 | 9.057e-3 | 3.822e-2 |
| Blue | 0.4001 | 0.1132 | 0.5699 | 7.122e-3 | 4.177e-3 | 1.339e-2 |
| Green1 | 0.4293 | 0.0970 | 0.5748 | 8.379e-3 | 5.134e-3 | 1.608e-2 |
| Green2 | 0.4589 | 0.1172 | 0.6347 | 9.216e-3 | 5.137e-3 | 1.692e-2 |
| Yellow | 0.3665 | 0.1087 | 0.5296 | 7.703e-3 | 3.934e-3 | 1.360e-2 |
| Red | 0.3658 | 0.1174 | 0.5420 | 5.076e-3 | 2.722e-3 | 9.159e-3 |
| RedEdge | 0.3288 | 0.0841 | 0.4550 | 1.493e-2 | 1.028e-2 | 3.035e-2 |
| NIR | 0.3180 | 0.0541 | 0.3991 | 2.083e-2 | 8.946e-3 | 3.425e-2 |

**20240514 is CALM but has 6/8 AI flags with no roughness.** The fringes are directional but not rough, which fits glint (see Py8).

**Wording mismatch with Table 1.** The Table 1 wording for 20240514 and 20230629 differs from the diagnostic output (noted in Py5).

### (b) Random seeds in the original MATLAB runs

**Stage 2** (`Code/Run_Stage2_Hierarchical_CuratedHoldout_QuickTest.m`): seeds are fixed for the downsampling.
- **Calls:**
  - line 115 `downsampleClass(Xgh, Ygh, opts.DownsampleGate, 11)` (gate HPY, 3,000 of 9,883);
  - line 116 `... 17)` (gate stressed, 3,000 of 3,318);
  - lines 128–129 `... 23)` and `... 29)` (stress model PLE and BLD).
- **Inside `downsampleClass`** (lines 478–485): `rng(seed); idx = randperm(size(X,1), nKeep);`. The function returns early, without calling `rng`, when the class already has no more rows than the cap. The stress caps (2,500) exceed both pools (1,290 PLE, 2,028 BLD), so seeds 23 and 29 are never applied.
- **Bootstrap:** the `TreeBagger` calls (lines 157 and 166) have no `rng` of their own. They continue the stream left by `rng(17)` and the `randperm` after it.
- **Result:** given the same MATLAB version, the Stage 2 run is reproducible. The bootstrap stream is not reproducible outside MATLAB (Py4).

**Stage 1** (`Code/Train_Stage1_LogRatio.m`): no seed is fixed.
- **No `rng` call:** the script contains none. Its `randperm` calls (line 576 for the plot sample, lines 695–696 for the SND/RBL subsample in each LOO fold) and all `TreeBagger` calls use whatever global stream state the MATLAB session had.
- **Driver:** the pipeline driver (`RunAll_LogRatio_A3.m`) does not set a seed either.
- **Rest of the pipeline:** the only fixed seeds elsewhere are in diagnostics: `QC_Alignment.m` `rng(42)`, and `rng(0)` / `rng(42)` in several `Diag_*` scripts.

### (c) Stage 2 sample counts (confirmed)

These are recovered from the training matrices saved inside the paper's forests.

| Set | HPY | PLE | BLD | Total |
|---|---|---|---|---|
| Gate training | 3,000 | 1,173 | 1,827 | 6,000 |
| Stress-split training | – | 1,290 | 2,028 | 3,318 |
| Pools before downsampling | 9,883 | 1,290 | 2,028 | |
| Test (spatial hold-out) | 7,068 | 831 | 1,279 | 9,178 |
| Anchors | 1,262 (20230509) | – | 434 (20240613) | |

**The figures in the R1 docs are correct.** The gate's 3,000 stressed rows break down as 1,173 PLE + 1,827 BLD, and the stress model uses every PLE and BLD row.

**Where the training pixels come from, by date:**
- **BLD:** all from the three June dates: 20240613 692, 20240616 684, 20240620 652.
- **PLE:** 20230509 538 and 20240620 752.
- **HPY:** six dates (20230404, 20230509, 20230629, 20240307, 20240317, 20240408), and none from June.

### (d) The 20,330-pixel mask compared with the reef polygon

**Depth source:** Allen Coral Atlas bathymetry (`bathymetry_0.tif`, centimetres, about 10 m pixels), resampled by nearest neighbour onto the 3 m grid. Both areas are 100% covered.

**Shoreline:** there is no shoreline file on disk, so land is taken from the imagery.
- **Rule:** median NDWI (Green1 vs NIR) over the 9 CALM dates < 0, keeping connected areas of at least 1,000 pixels.
- **Result:** one island, 30,474 pixels. It is the only land in the scene and lies south of the mask (half of the mask pixels are 1.47–1.73 km from it). No mask pixel is on land.

| Area | Pixels | Depth, median [IQR] (m) | Distance from shore, median [IQR] (m) |
|---|---|---|---|
| Coral mask | 20,330 | 0.73 [0.32–1.28] | 1,657 [1,469–1,734] |
| Reef polygon | 197,268 | 1.75 [0.86–2.79] | 1,413 [1,077–1,755] |
| Reef polygon outside the mask | 176,938 | 1.92 [1.00–2.92] | 1,376 [1,065–1,772] |

**Reading:**
- **Coverage:** all 20,330 mask pixels lie inside the reef polygon, which is 10.3% of the polygon.
- **Depth:** the mask is the shallow part of the reef: median 0.7 m against 1.9 m for the rest of the polygon.
- **Distance from shore:** similar for both. The mask sits on the outer, shallower part of the shoal rather than near the island.

**For the author:** please confirm this matches your knowledge of the site.
- **ACA depths:** these are relative to the ACA datum and carry about 1 m uncertainty here (a sibling ICESat-2 analysis found RMSE 1.07 m), so read them as relative depths.
- **Distance:** "distance from shore" means distance from that island.

### (e) 20240514 stressed pixels in the Stage 2 feature space

**Setup:** frozen paper forests, compat units. `axis_t` is the position along the healthy (0) to bleached (1) axis; the residual is the distance off that axis. Values are medians [IQR].

| Group | Pixels | axis_t | Residual |
|---|---|---|---|
| 20240514, gate-stressed (all) | 8,163 | 0.283 [0.242–0.325] | 0.0621 [0.0529–0.0685] |
| 20240514, classified PLE | 4,633 | 0.276 [0.237–0.317] | 0.0668 [0.0615–0.0713] |
| 20240514, classified UNK | 3,704 | 0.289 [0.246–0.332] | 0.0527 [0.0448–0.0609] |
| 20240514, classified HPY | 11,990 | 0.170 [0.102–0.225] | 0.0586 [0.0520–0.0669] |
| June BLD training (all three dates) | 2,028 | 1.006 [0.950–1.055] | 0.0233 [0.0182–0.0291] |
| BLD training 20240613 | 692 | 1.052 [1.019–1.085] | 0.0211 [0.0175–0.0251] |
| BLD training 20240616 | 684 | 1.010 [0.972–1.052] | 0.0249 [0.0190–0.0314] |
| BLD training 20240620 | 652 | 0.928 [0.888–0.974] | 0.0253 [0.0187–0.0330] |
| PLE training 20230509 | 538 | 0.768 [0.714–0.810] | 0.0705 [0.0615–0.0789] |
| PLE training 20240620 | 752 | 0.904 [0.870–0.937] | 0.0372 [0.0318–0.0422] |
| HPY training (all dates) | 9,883 | −0.233 [−0.670–0.669] | 0.0657 [0.0503–0.0862] |

**Reading:**
- **Position on the axis:** 20240514's stressed pixels sit low, at a median `axis_t` of 0.28. That is far from June BLD (1.01) and well below the 20230509 PLE pixels (0.77), but above its own HPY pixels (0.17). On axis position, they are only slightly paler than healthy.
- **Residual:** about 0.06, like 20230509 PLE (0.07) and HPY (0.07), and unlike June BLD (0.02).
- **What this means:** the gate sends them to the stress branch on features other than axis position. They are nowhere near the BLD cluster, which is why the stress model calls them PLE or UNK and almost never BLD (0.01%).

## Add-on (f): how "Stage 1 LOO recall = 1.000 across 36 polygons" was computed

**Source:** `Code/Train_Stage1_LogRatio.m`, section "LOO Validation (CORAL polygon holdout)" (lines 633–723), and A3 `Phase_Stage1/run_20260401_170150/stage1_report.txt`.

- **Unit: one polygon per fold.** The 36 folds are the 36 HC (coral) polygons in `TrainingStage1.geojson` (ids 1–22 and 52–65). The recall for a fold is the share of that polygon's pixels predicted CORAL: `loo_recall(fold)=mean(pred_loo==1);` (line 715). The reported 1.000 is the unweighted mean over the 36 polygons (`mean_recall=mean(loo_recall(valid));`, line 721), and every polygon scored 1.000.
- **Retrained for every held-out polygon: yes**, but not with the production model's settings.
  - **Each fold** grows a new forest: `rf_loo=TreeBagger(50,Xl,double(Yl),...,'Cost',cost_loo);` (line 704). That is 50 trees, not the production 200.
  - **Cost matrix:** `cost_loo=[0 3 3; 1 0 1; 1 1 0]` (line 675), against `[0 2 2; ...]` in production.
  - **Spatial hold-out:** coral training pixels come from all four CORAL dates with the held polygon removed: `train_coral_mask=coral_mask & ~held.mask;` (line 680), sampled per date over `CORAL_DATES` (lines 681–687).
  - **SND/RBL:** taken from the baseline stack and randomly subsampled to 3 × the coral count (`idx_s=randperm(...)`, lines 695–696, with no seed).
- **The held-out polygon is tested on the baseline composite, not on a single date.** `[Xheld,~]=sampleFeatures(Abase,held.mask,...)` (line 709). `Abase` is `BASELINE_FEATURES.tif`, the pixel-wise temporal median of 20230404, 20230509 and 20240307 (`Baseline_Build.m`, Phase 6).
- **Within or across dates?** The HC polygons carry no date (no `date` property), so the same footprint is used on every date.
  - **Training:** the held polygon's footprint is excluded on all four CORAL dates (20230404, 20230509, 20240307, 20241111).
  - **Testing:** done on the median of three of those dates.
  - **So:** the test is spatially independent, meaning those pixels are never in training on any date. It is not independent in time, because the test image is a composite of training dates.
  - **All 36 polygons are baseline-period coral.** None comes from a bleaching or wave date.
- **What it shows and doesn't show:** spatially held-out coral on the calm baseline composite is never confused with sand or rubble. It does not test generalisation to other dates or sea states. The report's separate cross-date test (CORAL from 20240408, SND/RBL from 20230519) covers that.

## Add-on (g): uncurated PLE candidates

**The candidate set is on disk.**
- **What it is:** the 923 train-eligible PLE polygons in the unified training labels, which `select_middle_pale` scores.
- **Curation:** 88 were kept, all of them in block 2 (train) or block 7 (test). 834 were not used, and 922 have mask pixels.
- **Classifier:** each unused polygon is classified on its own date with the frozen paper Stage 2 (compat units).
- **Trajectory position:** the centroid `axis_t` that the curation used. It keeps polygons inside the band 0.736 < t < 0.915 with residual ≤ 0.0956.
- **Overlap:** none of these pixels is also a Stage 2 training pixel.

| Group | Polygons | Pixels | HPY % | PLE % | BLD % | UNK % | Polygons by majority class (HPY / PLE / BLD / UNK) |
|---|---|---|---|---|---|---|---|
| **All unused candidates** | 834 | 20,280 | 56.1 | 8.9 | 33.9 | 1.1 | 448 / 83 / 303 / 0 |
| Used (curated; reference) | 88 | 2,121 | 2.6 | 92.1 | 2.6 | 2.7 | 2 / 83 / 2 / 1 |
| Rejected: t below band (healthy side) | 391 | 9,882 | 93.2 | 6.3 | 0.0 | 0.6 | 361 / 30 / 0 / 0 |
| Rejected: t above band (bleached side) | 439 | 10,292 | 20.1 | 11.5 | 66.8 | 1.7 | 83 / 53 / 303 / 0 |
| Rejected: residual above cut | 4 | 106 | 97.2 | 0.9 | 0.0 | 1.9 | 4 / 0 / 0 / 0 |
| Unused, t < 0 | 334 | 8,413 | 100.0 | 0.0 | 0.0 | 0.0 | 334 / 0 / 0 / 0 |
| Unused, 0 ≤ t < 0.25 | 17 | 488 | 100.0 | 0.0 | 0.0 | 0.0 | 17 / 0 / 0 / 0 |
| Unused, 0.25 ≤ t < 0.5 | 2 | 73 | 39.7 | 56.2 | 0.0 | 4.1 | 1 / 1 / 0 / 0 |
| Unused, 0.5 ≤ t < 0.75 | 39 | 940 | 33.1 | 61.4 | 0.0 | 5.5 | 10 / 29 / 0 / 0 |
| Unused, 0.75 ≤ t < 1 | 106 | 2,419 | 9.0 | 33.3 | 54.4 | 3.3 | 12 / 35 / 59 / 0 |
| Unused, t ≥ 1 | 336 | 7,947 | 24.2 | 4.7 | 69.9 | 1.2 | 74 / 18 / 244 / 0 |
| Unused, 20230404 | 84 | 2,076 | 100.0 | 0.0 | 0.0 | 0.0 | 84 / 0 / 0 / 0 |
| Unused, 20230509 | 44 | 1,099 | 34.2 | 60.6 | 0.0 | 5.2 | 12 / 32 / 0 / 0 |
| Unused, 20230629 | 81 | 1,974 | 100.0 | 0.0 | 0.0 | 0.0 | 81 / 0 / 0 / 0 |
| Unused, 20240307 | 94 | 2,431 | 100.0 | 0.0 | 0.0 | 0.0 | 94 / 0 / 0 / 0 |
| Unused, 20240317 | 83 | 2,167 | 100.0 | 0.0 | 0.0 | 0.0 | 83 / 0 / 0 / 0 |
| Unused, 20240408 | 94 | 2,355 | 100.0 | 0.0 | 0.0 | 0.0 | 94 / 0 / 0 / 0 |
| Unused, 20240613 | 134 | 3,066 | 0.0 | 0.2 | 99.7 | 0.1 | 0 / 0 / 134 / 0 |
| Unused, 20240616 | 141 | 3,280 | 0.0 | 0.8 | 98.0 | 1.2 | 0 / 0 / 141 / 0 |
| Unused, 20240620 | 79 | 1,832 | 0.0 | 60.0 | 32.7 | 7.3 | 0 / 51 / 28 / 0 |

**Reading: the frozen classifier follows the trajectory position, not the "PLE" label.**
- **Healthy side:** candidates on the pre-bleaching dates (2023 and March–April 2024), and those with t < 0.25, come out 100% HPY. These autolabels were not pale.
- **Bleached side:** candidates on 13 and 16 June come out ≥ 98% BLD.
- **In between:** candidates in the middle of the trajectory (0.25–0.75, mostly 20230509) come out mostly PLE. On 20240620 the split is 60% PLE / 33% BLD, close to that date's Table 5 split (55.6 / 37.2).
- **What this means:** the curated 88 are a band of the trajectory, and the classifier treats the uncurated rest consistently with where each polygon sits on it.
