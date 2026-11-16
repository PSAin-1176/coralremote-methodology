# Py5: baseline coral mask and wave diagnostic

Status 16 Nov 2026. Both reproduce the paper of record exactly.

## Baseline coral mask (`coralremote/mask.py`)
The provenance was traced to two MATLAB stages:

1. **Consensus votes** (`Baseline_Coral_LogRatio.m`, A3 run `Baseline_Coral/run_20260401_103612`):
   - inputs are the Phase 9 Stage 1 class maps for 20230404, 20230509 and 20240307 (`Phase9_Classify/run_20260401_102734`);
   - each pixel gets one vote per date classified CORAL, zeroed outside the reef;
   - tier 2 = ≥ 2/3 votes, tier 3 = 3/3; `mask0` = tier ≥ 2 = **21,460 px** (the number in the submitted manuscript).
2. **Morphology on 6 Apr** (scripts in `done/baseline_mask/Autolabel`):
   - `Clean_CoralMask_Conservative`: drop objects < 12 px, fill holes ≤ 12 px, one majority pass limited to the boundary zone, then repeat the object and hole steps → 21,961 px;
   - `Peel_CoralMask_OnePixel`: erode with a 3×3 cross → 20,273 px;
   - `Finalize_CoralMask_MinorFixes`: drop objects < 4 px, fill holes ≤ 40 px → **20,330 px** = `Asset/baseline_coral_mask.mat`.

| Array saved in `Asset/baseline_coral_mask.mat` | Python | Differing px |
|---|---|---|
| `vote_map`, `tier_map` | from the 3 class maps | 0 |
| `mask0` | 21,460 | 0 |
| `maskClean` | 21,961 | 0 |
| `maskPeel` | 20,273 | 0 |
| `maskFinal` = `baseline_coral_mask` | 20,330 | 0 |

The MATLAB conventions reproduced are:
- `bwareaopen` uses 8-connectivity;
- `imfill 'holes'` uses a 4-connected background;
- `bwperim` uses 4-connectivity;
- `bwmorph majority` sets a pixel when ≥ 5 of its 3×3 neighbourhood are set;
- `imerode` pads with foreground.

**Provenance note.** The Phase 9 maps came from an earlier Stage 1 model, `Phase_Stage1/run_20260401_102559` (OOB 0.1001, UNK threshold 0.55). That isn't the paper's reported Stage 1 model (`run_20260401_170150`, OOB 0.1004), and the 102559 model itself didn't survive. Both models were trained on true-reflectance A3 features (1 Apr), so the mask isn't affected by the Stage 2 double-scaling bug. **§2 of the plan changes from "probably no" to confirmed.** A rebuild of the mask with the 170150 model hasn't been run. It needs the UDM2 masks on the F: drive.

## Wave FFT diagnostic (`coralremote/wave.py`)
This is a port of `Diag_WaveFFT_ThreeClass.m`:
- NDWI ocean tiles of 128 × 128 px with ≥ 80% ocean;
- per band, a Hann-windowed 2-D FFT and a 6–96 m annulus (radial bins 4–64);
- WBP (roughness) and a 12-sector AI (directionality) per band;
- cross-date median + 1.5 MAD flags;
- verdict WAVE if ≥ 4 roughness and ≥ 6 AI flags, TURBULENT if ≥ 4 roughness flags, otherwise CALM.

The original input was the Phase 1 co-registered tifs on `C:\SAT`, which are lost. The Source raw stack holds the same co-registered reflectance on the analysis grid.

| Date | Tiles | Roughness flags | AI flags | Python | Paper |
|---|---|---|---|---|---|
| 20230404 | 76 | 0 | 0 | CALM | CALM |
| 20230509 | 76 | 0 | 0 | CALM | CALM |
| 20230519 | 75 | 8 | 7 | WAVE | WAVE |
| 20230629 | **38 (LOW)** | 7 | 3 | TURBULENT | TURBULENT (LOW, 38 tiles) |
| 20240307 | 76 | 0 | 0 | CALM | CALM |
| 20240317 | 75 | 0 | 0 | CALM | CALM |
| 20240408 | 60 | 2 | 0 | CALM | CALM |
| 20240514 | 76 | 0 | 6 | CALM | CALM (directional AI but no roughness) |
| 20240613 | 76 | 6 | 8 | WAVE | WAVE |
| 20240616 | 76 | 8 | 8 | WAVE | WAVE |
| 20240620 | 76 | 8 | 8 | WAVE | WAVE |
| 20240926 | 76 | 0 | 0 | CALM | CALM |
| 20241111 | 76 | 0 | 0 | CALM | CALM |
| 20250518 | 76 | 0 | 0 | CALM | CALM |

All 14 verdicts and confidence tiers are identical, with 13 dates at HIGH. The paper's numeric AI and WBP values weren't saved, so verdicts, tile counts and confidence are the comparison.

**Units.** The diagnostic can't be affected by a scale error. NDWI, the ocean-pixel test and AI are scale-free, and WBP is compared only with cross-date median + MAD thresholds. Scaling every date by the same factor leaves all flags unchanged. **§2 changes from "probably no" to confirmed.**

**Table 1 wording (for the author).** The manuscript's Table 1 sea-state labels don't follow these diagnostic classes one to one:
- 20240514 is "Low-moderate turbulent" in the table, CALM here;
- 20230629 is "Moderate wave" in the table, TURBULENT (LOW) here.

The table was presumably written from visual review, so this isn't a port issue. It may be worth one sentence in the methods.

## Reproduce
```bash
cd E:/CoralRemote/coralremote-python
PYTHONPATH=. python scripts/run_py5.py      # ~10 s; exits 0 on full parity
```
