# Py8: why 20240514 and 20250518 swing between seeds

This is a diagnosis only. The paper's method is not changed.

## Setup

- **Runs:** both units (compat = paper units, fixed = corrected units) and retraining seeds 0–4. Each seed trains the Stage 2 gate and stress forests with 1,000 trees.
- **Scores kept:** for every in-mask pixel on all 14 dates, the class-1 score and the per-tree spread are stored after 250 trees and after 1,000.
- **Check:** in both modes, the first 250 trees give exactly the Py6 250-tree maps on every date and seed. So the 250-tree numbers below are the Py6 retraining results.
- **Bootstrap-only runs:** a second set holds the training-row draw fixed at seed 0 and varies only the bootstrap (seeds 0–4, fixed units, 250 trees).
- **Scripts:**
  - `scripts/run_py8_collect.py` trains the forests and caches the scores.
  - `scripts/run_py8_analyze.py` produces the tables below.
  - `scripts/run_py8_scene.py` produces the scene statistics.
- **Outputs:** the per-pixel scores stay in the local cache. Only aggregates are reported here.

### How a seed enters training

The seed changes two things in the gate:

- **Training rows:** the gate draws 3,000 HPY + 3,000 stressed rows from the labelled pool using `rng(seed)`.
- **Bootstrap:** the order in which trees resample those rows.

The stress model always uses every PLE and BLD row, so the seed changes only its bootstrap.

### Bug found and fixed along the way

The fixed-units seed 0 stress forest never finished.

- **Where:** tree 402, a node of 35 rows, feature `log_C`.
- **Cause:** the two lowest values at that node are neighbouring single-precision numbers. Their single-precision midpoint rounded onto the lower one, so `x < cut` sent no rows left. The node then split the same way forever.
- **Fix:** `coralremote/cart.py` now makes such a node a leaf.
- **Effect on earlier results:** none. The guard only acts in a case that previously never finished. All 48 tests pass, and the 250-tree maps are unchanged.

## Results

Ranges are the maximum minus the minimum of the class fraction across the five seeds, in percentage points (pp) of in-mask pixels.

| Date | Mode | PLE / UNK range, 250 trees | PLE / UNK range, 1,000 trees | Bootstrap-only PLE / UNK range | UNK from gate / stress (%) | Pixels near 0.55, gate / stress (%) | Pixels that flip (%): gate / stress only | Out of range: any feature / `axis_t` (%) |
|---|---|---|---|---|---|---|---|---|
| 20240514 | compat | 11.9 / 17.8 | 7.7 / 13.2 | — | 1.8 / 13.6 | 4.0 / 26.8 | 22.1: 8.8 / 13.3 | 1.0 / 0.1 |
| 20240514 | fixed | 8.8 / 7.3 | 6.3 / 4.5 | 7.3 / 7.2 | 2.1 / 4.5 | 4.1 / 9.5 | 13.9: 6.9 / 7.0 | 0.6 / 0.1 |
| 20250518 | compat | 36.6 / 29.3 | 33.1 / 22.2 | — | 16.2 / 8.7 | 32.6 / 29.7 | 63.3: 61.5 / 1.8 | 81.6 / 77.5 |
| 20250518 | fixed | 43.8 / 23.4 | 43.9 / 23.8 | 0.1 / 15.1 | 17.6 / 0.1 | 33.3 / 10.1 | 64.1: 64.1 / 0.0 | 58.2 / 51.4 |
| 20230509 (next worst) | fixed | 3.3 / 1.8 | 3.1 / 1.8 | 0.5 / 0.3 | 3.7 / 0.0 | 7.3 / 0.1 | 13.0: 12.9 / 0.1 | 2.1 / 0.0 |
| June dates | fixed | ≤ 1.9 / ≤ 1.2 | ≤ 1.3 / ≤ 0.9 | ≤ 1.9 / ≤ 1.1 | ≤ 0.4 / ≤ 5.7 | ≤ 0.7 / ≤ 9.1 | ≤ 6.6 | ≤ 2.3 / 0.0 |
| The other 7 dates | both | ≤ 0.1 | ≤ 0.1 | ≤ 0.1 | 0 | 0 | ≤ 0.2 | ≤ 1.8 / ≤ 0.2 |

Notes on the columns:

- **Pixels near 0.55, gate:** share of all in-mask pixels whose gate confidence is within ±0.05 of 0.55.
- **Pixels near 0.55, stress:** share of in-mask pixels that the gate routes to the stress model and whose stress confidence is within ±0.05 of 0.55.
- **Pixels that flip:** in-mask pixels whose final label is not the same in all five seeds. They are split into flips where the gate's own decision (HPY / stressed / UNK) changes, and flips where only the stress model's decision changes.

Per-seed class fractions at 1,000 trees, as HPY / PLE / UNK %. Paper values from Table 5: 20240514 PLE 22.8 and UNK 18.2; 20250518 HPY / PLE / UNK 44.4 / 23.4 / 31.9.

| Date | Mode | Seed 0 | Seed 1 | Seed 2 | Seed 3 | Seed 4 |
|---|---|---|---|---|---|---|
| 20240514 | compat | 60.5 / 29.0 / 10.4 | 60.2 / 24.8 / 15.0 | 57.9 / 29.9 / 12.2 | 54.2 / 22.2 / 23.6 | 60.2 / 22.8 / 17.0 |
| 20240514 | fixed | 59.8 / 34.6 / 5.5 | 58.8 / 34.3 / 6.9 | 57.8 / 38.8 / 3.4 | 57.3 / 38.2 / 4.6 | 59.6 / 32.5 / 7.9 |
| 20250518 | compat | 67.2 / 15.7 / 17.0 | 38.9 / 23.0 / 37.8 | 41.9 / 39.6 / 18.1 | 67.5 / 16.1 / 16.1 | 54.8 / 6.6 / 38.3 |
| 20250518 | fixed | 80.5 / 1.9 / 17.5 | 39.7 / 40.2 / 19.8 | 44.2 / 45.8 / 9.6 | 48.8 / 36.1 / 14.9 | 59.9 / 6.4 / 33.5 |
| 20250518 | fixed, rows held at seed 0 (bootstrap only, 250 trees) | 81.0 / 1.9 / 17.0 | 71.9 / 1.9 / 26.1 | 87.1 / 1.8 / 11.0 | 81.1 / 1.9 / 17.0 | 77.1 / 1.9 / 20.8 |

## Answers

### Q1. Where UNK comes from

**20240514: mostly the stress model.**
- **Source:** in compat, 13.6% of in-mask pixels are UNK because the stress model is unsure, against 1.8% from the gate.
- **Near the threshold:** 27% of in-mask pixels are routed to the stress model and sit within ±0.05 of 0.55. That is about three times any other date (June dates ≤ 11%, stable dates 0%).
- **Fixed units:** the stress UNK falls to 4.5% and the near-threshold share to 9.5%. This is the Py6 finding that the units fix moves 20240514 pixels from UNK to PLE.

**20250518: almost entirely the gate.**
- **Source:** 16–18% of in-mask pixels are UNK because the gate is unsure.
- **Near the threshold:** a third of all in-mask pixels (33%) have gate confidence within ±0.05 of 0.55. The next highest date is 20230509 at 7%; the stable dates are at 0%.

### Q2. Which model disagrees between seeds

**20240514: both, about equally.**
- **Flips:** 14–22% of pixels flip. Half are stress-only flips (PLE ↔ UNK at the 0.55 stress cut) and half are gate flips (HPY ↔ UNK).
- **Class pairs:** PLE–UNK 64–66%, HPY–UNK 28–32%, HPY–PLE 1–7%.

**20250518: the gate.**
- **Flips:** 64% of all in-mask pixels flip between seeds, and essentially every flip is a gate flip (64.1 of 64.1% in fixed, 61.5 of 63.3% in compat).
- **Class pairs (fixed):** PLE–UNK 37%, HPY–PLE 33%, HPY–UNK 30%. The gate moves a large block of the reef between HPY and stressed, and the stress model then calls almost all of that block PLE.
- **Not the stress model:** its scores barely move between seeds (cross-seed SD 0.016–0.025, against 0.105–0.118 for the gate).

### Q3. Out of domain

**20240514 is not out of domain.**
- **Features:** only 0.6–1.0% of pixels have any of the 21 features outside the training range, and 0.1% on `axis_t`.
- **Flips:** none of the flipping pixels is out of range.

**20250518 is out of domain on `axis_t`.**
- **`axis_t`:** 77.5% of pixels are beyond the training range in compat, which matches the 77.6% reported in the paper, and 51.4% in fixed.
- **Other features:** `norm_G2` (15%) and `norm_Y` (11%) come next. Everything else is below 1%.

**But being out of range does not pick out the flipping pixels:**
- **Flip rates:** P(flip | out of range) = 63–65%, against 63–64% for pixels inside the range. Flips are 59–81% out of range, about the same as the date as a whole.
- **What this means:** the whole scene sits on the gate's boundary, not just its out-of-range part.
- **Why `axis_t` still matters:** a tree cannot extrapolate. Every pixel beyond the training range lands in the same edge leaves, so the label for that block is set by the few training rows at the edge of the `axis_t` range. Those rows change with the row draw (see Q4).

### Q4. Would more trees fix it?

**20240514: partly. It is bootstrap noise in the stress model.**
- **More trees:** going from 250 to 1,000 trees halves the cross-seed spread of the stress score (SD 0.020 → 0.009). The ranges shrink: compat PLE 11.9 → 7.7 pp and UNK 17.8 → 13.2 pp; fixed PLE 8.8 → 6.3 pp and UNK 7.3 → 4.5 pp.
- **Bootstrap-only runs:** they give almost the full spread (fixed PLE / UNK 7.3 / 7.2 pp, against 8.8 / 7.3 pp with the rows varying). So the variance comes from the bootstrap, and more trees help.

**20250518: no. It is the gate's training-row draw.**
- **More trees:** at 1,000 trees the range is unchanged (fixed PLE 43.8 → 43.9 pp; compat 36.6 → 33.1 pp), and the cross-seed SD of the gate score stays at 0.105.
- **Bootstrap-only runs:** with the 3,000 + 3,000 gate rows held fixed, PLE is 1.8–1.9% in every seed (range 0.1 pp). Only HPY ↔ UNK moves (UNK range 15 pp).
- **What this means:** which labelled rows the gate draws decides whether this scene is HPY or stressed.

**Per-tree vote spread:**
- **Gate:** trees disagree more on flipping pixels than on stable ones (per-tree SD of the gate score 0.39–0.40 against 0.29–0.32 on 20250518, and 0.29–0.32 against 0.21–0.22 on 20240514; 0.5 would be an even split).
- **Stress model:** its trees disagree almost equally everywhere (SD about 0.36–0.41). That only matters for pixels near 0.55.

### Q5. Scene conditions

Scene statistics are medians over the baseline coral mask, in fixed units, from `run_py8_scene.py`. Sun angle, cloud and haze come from the PlanetScope metadata summary.

- **Both dates:** 0% cloud and 0% haze, and the two highest sun elevations in the series (20240514 70.1°, 20250518 72.4°). The wave diagnostic calls both CALM.
- **20240514: strong glint.**
  - **Glint:** the glint correction removes 15.1% of raw blue over the reef, similar to the June dates (9.5–18.1%) and 20230629 (14.1%, the hazy, TURBULENT date), against ≤ 7.6% on all other dates.
  - **Wave flags:** the wave diagnostic flags 6 directional (AI) sectors but no roughness, which fits glint rather than waves.
  - **Brightness:** 0.43, near the middle of the series.
  - **Compat clipping:** compat mode clips none of its NIR.
  - **Summary:** it looks like a sun-glint scene whose spectra fall between HPY and PLE after the correction. The stress model is unsure about them, not out of range.
- **20250518: almost no glint, but the most green-heavy spectrum of the series.**
  - **Glint:** the correction removes 1.4% of blue.
  - **Spectrum:** Green1 0.076 against Blue 0.044. 2·Blue/(Red + Green1) is 0.80, the lowest of all 14 dates; the closest is 20230509 at 0.82, the most PLE-heavy baseline date. That shift puts `axis_t` beyond the training range.
  - **Compat clipping:** compat mode floors 6.9% of its NIR pixels at 1e-6. This is modest; 20240307 is at 84%. The pattern is the same in both modes, so NIR clipping is not the cause.

### Q6. Bottom line

**Why these two dates:**
- **20240514** is a glint-heavy scene. Many reef pixels land close to the stress model's 0.55 cut. Which side they fall on depends on bootstrap noise in the 250-tree stress forest.
- **20250518** has an unusual, green-shifted spectrum that sits outside the training range on `axis_t`. A third of the reef lies right on the gate's HPY / stressed boundary, and where that boundary falls for this scene depends on which 6,000 training rows the gate happens to draw. The date has no stable answer under this training design: the labelled data do not pin it down.

**Which remedy works:**

| Remedy | 20240514 | 20250518 |
|---|---|---|
| More trees (1,000) | Helps. PLE range 8.8 → 6.3 pp, UNK 7.3 → 4.5 pp (fixed); compat UNK 17.8 → 13.2 pp | No effect (PLE 43.8 → 43.9 pp) |
| Averaging scores over seeds | Helps, about the same as more trees (2-seed averages: PLE range 5.5 pp in fixed) | Gives one reproducible number (5-seed average, fixed: HPY / PLE / UNK 49.1 / 32.8 / 17.9). It cannot make the answer less row-dependent: 2-seed averages still span 42 pp of PLE |
| Changing the 0.55 threshold | At 0.5, UNK disappears and the PLE range drops to 1.7 pp (fixed) or 5.5 pp (compat). This changes the method | No help. At 0.5 the HPY / PLE range grows to 48 pp; at 0.6 the PLE range falls only because PLE moves into UNK (UNK up to 43–51%) |
| Reporting a range across seeds | Appropriate | **The only honest option** |

**Recommendation** (no method change):
- **20240514:** report the 5-seed mean with its range, and note that it is threshold-sensitive because of glint.
- **20250518:** report the range across seeds (fixed: PLE 2–46%, HPY 40–81%) together with the out-of-range `axis_t` share, and say the date is outside what the training data constrain. Treat it as indeterminate, not as a measured composition.
- **The other 12 dates:** these move ≤ 3.3 pp between seeds and need no change.
