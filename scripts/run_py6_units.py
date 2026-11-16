"""Py6 P3/P4 checks: (1) the fixed stack (Source glint stack x SR_SCALE) gives
true reflectance: Hedley leaves NIR unchanged, so glint NIR must equal raw NIR;
the surviving Phase 3 GLINT tifs must agree. (2) per-band, per-date fraction
of in-mask pixels that the floor raises, compat (1e-6 on the double-scaled
stack = 0.01 reflectance) vs fixed (1e-4 reflectance), plus negatives."""
import json
import sys

import numpy as np
import rasterio

from coralremote import io, paths
from coralremote.features import prepare_spectra
from coralremote.preprocess import raw_stack
from coralremote.units import COMPAT, FIXED, assert_reflectance

g = io.load_glint_stack()
r = raw_stack()
m = io.load_baseline_mask()[0]
out = {}
print("date      NIR glint*1e4 vs raw: maxrel  exact%   GLINT.tif maxrel | floored % in mask (B1..B8): compat | fixed | negative")
for d in g.dates:
    G = prepare_spectra(g.get(d), FIXED)
    R = r.get(d)
    ok = np.isfinite(G[..., 7]) & np.isfinite(R[..., 7]) & (R[..., 7] > 0)
    rel = float((np.abs(G[..., 7].astype(np.float64) - R[..., 7])[ok] / R[..., 7][ok]).max())
    exact = float((G[..., 7][ok] == R[..., 7][ok]).mean() * 100)
    p = paths.PHASE3_GLINT_DIR / f"GLINT_{d}.tif"
    tif = None
    if p.exists():
        with rasterio.open(p) as ds:
            T = np.moveaxis(ds.read(), 0, -1).astype(np.float64)
        o = np.isfinite(T) & np.isfinite(G) & (T != 0)
        tif = float((np.abs(T - G)[o] / np.abs(T[o])).max())
    X = g.get(d)[m]
    X = X[np.isfinite(X).all(1)]
    assert_reflectance(prepare_spectra(X, FIXED), f"glint {d} (fixed)")
    fc = (prepare_spectra(X, COMPAT) < COMPAT.floor).mean(0) * 100
    ff = (prepare_spectra(X, FIXED) < FIXED.floor).mean(0) * 100
    neg = (X < 0).mean(0) * 100
    out[d] = {"nir_maxrel": rel, "nir_exact_pct": exact, "tif_maxrel": tif,
              "floored_compat_pct": fc.tolist(), "floored_fixed_pct": ff.tolist(), "negative_pct": neg.tolist()}
    print(f"{d}  {rel:9.1e} {exact:6.1f}   {'-' if tif is None else f'{tif:9.1e}'} | "
          + " ".join(f"{v:4.1f}" for v in fc) + " | " + " ".join(f"{v:4.2f}" for v in ff)
          + " | " + " ".join(f"{v:4.2f}" for v in neg))
if len(sys.argv) > 1:
    json.dump(out, open(sys.argv[1], "w"), indent=1)
