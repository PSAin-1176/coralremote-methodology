# Units, floors and the scale check

## Convention
- **Every array is surface reflectance on a 0–1 scale.** PlanetScope SR GeoTIFFs store reflectance × 10,000 (`units.SR_SCALE`). Divide once, on reading the GeoTIFF, and never again.
- `units.assert_reflectance(A, name)` runs on load. It checks:
  - the finite fraction;
  - the median visible reflectance is 0.005–0.3;
  - the median NIR is 0.001–0.1;
  - the maximum is below 1.5.
- A failure raises `UnitError`. A stack that is 10,000× too small or too large fails at once.
- **Floor:** log and normalised-spectrum features clip reflectance at **1e-4** (`units.FIXED.floor`), i.e. 0.01% reflectance. Hedley glint correction can produce small negatives. Clipping them is equivalent to excluding them here, because they are ≤ 0.01% of coral-mask pixels on every date. On the 14 Koh Man Nai dates the floor raises at most 0.18% of any band within the coral mask (`scripts/run_py6_units.py`).

## Modes
| Mode | Stack multiplier | Floor | Use |
|---|---|---|---|
| `fixed` (default for new runs) | × 10,000 on the paper-of-record Source stack, i.e. true reflectance | 1e-4 | All new work |
| `compat` | × 1 (stack left double-scaled) | 1e-6 (= 0.01 reflectance in true units) | Only to reproduce the AJSTR paper-of-record numbers |

The Source glint stack and the 18-feature Source stack were built on 5 Apr 2026, and both divided by `SR_SCALE` a second time. The raw stack wasn't double-scaled. `fixed` mode corrects this when the stack is loaded; the stored files are left as they are.

## Check that a correction is right
Hedley glint correction leaves NIR unchanged, so glint-corrected NIR must equal raw NIR. In `fixed` mode it does, to within one float32 rounding step, on all 14 dates (`tests/test_fixed_units.py`).

## Methods text (for future papers)
> Surface reflectance was obtained by dividing the PlanetScope SR product by 10,000. Before the logarithmic and normalised-spectrum features were computed, reflectance was clipped at a floor of 10⁻⁴ to keep logarithms finite. The floor affected no more than 0.18% of values in any band on any date within the coral mask.

For a correction note on the AJSTR paper, if the author decides one is needed:

> In the analysis reported, the reflectance stack used for Stage 2 was inadvertently scaled by an additional factor of 10⁻⁴. The numerical floor of 10⁻⁶ applied before the logarithmic features therefore acted as a floor of 0.01 reflectance, mainly on near-infrared values. Re-running the identical procedure on correctly scaled reflectance changed overall accuracy by −0.4 percentage points and left the dominant class of every date unchanged. Individual class fractions changed by up to about 9 percentage points (seed-averaged): on 20240514 the uncertain fraction fell from 18.2% to 6.6% in favour of pale; 20240616 bleached changed from 93.7% to 90.4%; and 20240620 bleached from 37.2% to 41.4%.
