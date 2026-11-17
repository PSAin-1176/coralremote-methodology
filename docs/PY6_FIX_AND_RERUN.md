# Py6: fix the units (P3, P4), audit the labels (P5), re-run (P6) and compare (P7)

Status 16 Nov 2026. Branch `python-port`, local only. The manuscript and the paper-of-record files are untouched.

## P7 decision: **material under the plan's Table 5 criterion; not material for OA, F1 or the headline narrative**
| Criterion (plan §4, P7) | Result | Material? |
|---|---|---|
| Stage 2 OA moves > 1 pp | −0.39 pp (0.9709 vs 0.9748); −0.19 pp against the same retraining in paper units | No |
| Any class F1 moves > 0.02 | max 0.015 (PLE 0.835 vs 0.850) | No |
| June BLD / September HPY narrative | Dominant class unchanged on every date in every seed: 20240613 and 20240616 BLD, 20240620 PLE (as in the paper), 20240926 HPY. 20250518 is HPY in 4/5 seeds; one seed flips in each mode, so this is seed noise | No |
| Any Table 5 cell moves > 2 pp | **Yes.** Beyond seed noise (same sign in all 5 paired seeds): 20240620 BLD **+4.7 pp** (37.2 → 41.4), PLE −3.8; 20240514 PLE **+9.1 pp** and UNK −8.8 pp; 20240616 BLD −2.6 pp on average (4 of 5 seeds) | **Yes** |

"Material" here is the plan's own trigger: any Table 5 cell moving more than 2 pp. What follows is the author's call under plan §5: with R1 under review and a material change, tell the editor in round 2. The changes are second-order:
- OA, kappa and F1 are unchanged within seed noise.
- No date's dominant class changes.

The fix moves pixels between PLE and BLD on 20240620, and from UNK to PLE on 20240514.

## P3: fixed stack
- `units.FIXED` = Source glint stack × 10,000, i.e. true reflectance.
- **Check, all 14 dates:** glint NIR equals raw NIR to within one float32 step (Hedley leaves NIR unchanged), with relative error ≤ 9.8e-8.
- **Check, 8 dates:** the surviving Phase 3 `GLINT_<date>.tif` files agree with the fixed stack to the same tolerance.
- The fixed reflectance passes `assert_reflectance`.
- The paper-of-record stack isn't modified. The conversion happens at load time.

## P4: floor
- `FIXED.floor` = **1e-4 reflectance**, applied by clipping (`max`). This keeps the paper's rule and just changes the threshold.
- **Old floor (paper):** the 1e-6 floor on the double-scaled stack equals 0.01 reflectance. Within the baseline mask it clipped, for example:
  - NIR on 84% of pixels on 20240307, 66% on 20240926, 23% on 20230404, 20240317 and 20241111;
  - red-edge on 60% on 20240307.
- **New floor:** clips at most 0.18% of any band on any date, so no band is flagged under the plan's 1% rule.
- **Negatives:** glint-corrected values are ≤ 0.01% negative, so clipping and excluding them are equivalent; clipping is kept.
- The per-date, per-band table is printed by `scripts/run_py6_units.py`.

## P5: labels: no relabelling needed
The strict and unified label polygons come from this chain:
- `Run_Autolabel_HealthyPerDate` and `Run_Autolabel_BleachedPerDate` (superpixels + seed statistics);
- the HPU and BPU polygons;
- `Build_StrickLabel_Dataset`, filtered by `highconf_*_mask_strict`.

These scripts read the double-scaled glint stack and the `SourceStack_AllTestedFeatures` stack. That stack is double-scaled too: brightness is about 3e-6, and DII is shifted by the constant −0.18·ln(1e4) = −1.66.

All of their decision thresholds are relative, so they don't depend on scale:
- z-scores against seed medians and MADs;
- percentile stretches (p1–p99);
- overlap fractions and keep-fractions.

The only absolute constants are small guards. Their measured effect on the paper data:

| Guard | Value on the double-scaled data | Binds? |
|---|---|---|
| `max(hi−lo, 1e-6)` in `robustScale` (brightness p1–p99) | 1.18e-6 to 4.7e-6 on every labelled date (interior mask) | No. Lowest headroom is 18%, on 20240408; it binds only on 20241111, which wasn't labelled |
| `max(MAD, 1e-9)` | brightness MAD ≥ 1.5e-7 | No |
| `max(2·blue, 1e-8)`, `max(blue, 1e-8)` | blue ≥ 2e-7, except a few blue = 0 pixels on 20240620 | Only on those pixels, which are outliers in any units |
| DII log floor (1e-6 on stack values) | B or R < 1e-6 on ≤ 0.72% of mask pixels (20240307) | Marginally |

So the frozen label polygons would come out the same in correct units. They stay as the training and test labels. The autolabel scripts weren't ported: superpixel and k-means steps aren't reproducible, and the audit shows no reason to regenerate.

## P6: re-run
`scripts/run_py6_rerun.py` runs the same procedure in both units:
- MATLAB-rule CART, 250 trees, MinLeaf 5, 0.55 UNK rule;
- curated PLE and the block split, recomputed in each mode;
- seeds 0–4 with the same sampling stream.

Comparing the two modes seed by seed separates the fix from retraining noise; Py4 showed that 20240514 and 20250518 vary 17–30 pp between seeds. Maps are in the cache folder, outside the repo.

**Deterministic steps that change with the fix:**

| | Paper units | Fixed |
|---|---|---|
| Healthy / bleached axis thresholds | 0.7568 / 0.8945 | 0.7647 / 0.8927 |
| Curated PLE polygons kept | 88 / 923 | 81 / 923 |
| PLE train / test pixels | 1,290 / 831 | 1,163 / 801 |
| HPY, BLD pixels | same | same |

**Test metrics** (mean [range] over 5 seeds):

| | Paper (frozen) | Retrain, paper units | Retrain, fixed | Fixed, on the paper's test pixels |
|---|---|---|---|---|
| OA | 0.9748 | 0.9728 [0.9702, 0.9759] | 0.9709 [0.9670, 0.9730] | 0.9706 |
| κ | 0.9344 | 0.9295 | 0.9246 | 0.9242 |
| F1 HPY / PLE / BLD | 0.985 / 0.850 / 0.976 | 0.983 / 0.845 / 0.976 | 0.981 / 0.835 / 0.976 | 0.981 / 0.836 / 0.974 |

**Table 5 cells that move.** Values are %; "fix effect" = fixed − paper units, paired by seed: mean [min, max].

| Date | Class | Paper | Paper units, retrain mean | Fixed mean | Fix effect |
|---|---|---|---|---|---|
| 20230509 | HPY | 41.97 | 40.32 | 38.87 | −1.45 [−2.39, −0.58] |
| 20240514 | PLE | 22.79 | 25.88 | 34.99 | **+9.11 [+1.65, +22.31]** |
| 20240514 | UNK | 18.22 | 15.35 | 6.61 | **−8.75 [−24.67, −0.29]** |
| 20240616 | BLD | 93.74 | 92.93 | 90.37 | −2.56 [−4.17, +0.11] |
| 20240620 | PLE | 55.62 | 56.34 | 52.59 | **−3.76 [−4.33, −3.45]** |
| 20240620 | BLD | 37.23 | 36.70 | 41.38 | **+4.68 [+4.39, +4.90]** |
| 20250518 | HPY / PLE / UNK | 44.4 / 23.4 / 31.9 | 54.1 / 20.7 / 24.9 | 55.2 / 26.9 / 17.7 | +1.1 / +6.2 / −7.3; each spans about ±20 pp across seeds (noise) |

Other dates:
- every other cell moves < 2 pp;
- the calm baseline and recovery dates are unchanged (≤ 0.2 pp);
- 20240613 BLD 93.97 vs 94.46.

Per-pixel agreement, fixed vs paper units for the same seed, averaged over all 14 dates, is 95.7%. The lowest dates are 20250518 (73.9%), 20240514 (88.1%) and 20240620 (91.2%).

**Reading.** The old floor clipped NIR on 23–84% of mask pixels on five dates (most on 20240307 and 20240926) and red-edge on up to 60% (20240307). Restoring them changes the curated PLE set and sharpens the stress model's PLE/BLD boundary:
- 20240620 reads more bleached (BLD +4.7 pp);
- 20240514 gives fewer low-confidence (UNK) pixels, which go to PLE instead.

Py8 shows that 20240514 is not out of domain: its UNK pixels sit near the stress model's 0.55 cut ([PY8](PY8_SEED_SENSITIVITY.md)). 20250518 stays seed-unstable in both units, so its Table 5 row should be reported as a mean over seeds in any re-run.

## Reproduce
```bash
cd coralremote-python
PYTHONPATH=. python scripts/run_py6_units.py        # P3/P4 checks, about 1 min
PYTHONPATH=. python scripts/run_py6_rerun.py        # P6, about 9 min (5 seeds x 2 modes)
PYTHONPATH=. python scripts/run_py6_compare.py      # P7 tables and criteria
```
