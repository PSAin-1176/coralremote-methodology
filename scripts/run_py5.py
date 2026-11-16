"""Py5: rebuild the baseline coral mask from the Phase 9 class maps and run the
wave FFT diagnostic; compare both with the paper of record."""
import numpy as np
import scipy.io as sio

from coralremote import mask, paths, wave
from coralremote.preprocess import raw_stack

m = sio.loadmat(str(paths.BASELINE_MASK), squeeze_me=True)
r = mask.build_baseline_mask(mask.load_phase9_class_maps(), m["reef_mask"] > 0)
print("Baseline mask chain (Python px / pixels differing from the saved MATLAB array):")
for k in ("vote_map", "tier_map"):
    print(f"  {k:10s} differing px {int((r[k] != m[k]).sum())}")
for k in ("mask0", "maskClean", "maskPeel", "maskFinal"):
    print(f"  {k:10s} {int(r[k].sum()):6d} px, differing {int((r[k] != (m[k] > 0)).sum())}")

s = raw_stack()
out, thr = wave.run_wave_diagnostic(s.get, s.dates)
print("\nWave FFT diagnostic:  date  tiles  Nrough  Ndir  verdict  confidence  (paper)")
ok = 0
for d, x in out.items():
    ok += x["verdict"] == wave.PAPER_VERDICTS[d]
    print(f"  {d}  {x['n_tiles_used']:3d}  {x['n_rough_bands']}/8  {x['n_bands_flag']}/8  "
          f"{x['verdict']:9s}  {x['confidence']:6s}  ({wave.PAPER_VERDICTS[d]})")
print(f"verdicts identical: {ok}/{len(out)}")
raise SystemExit(0 if ok == len(out) else 1)
