"""Regression tile (plan P8): a small fixed set of baseline-mask pixels on all
14 dates, used for fast unit, range and scale-invariance tests.

The tile is cut from the PlanetScope-derived glint stack, so it is written
to the local cache folder only and never committed (Planet data-sharing
terms). The repository keeps only aggregate golden numbers
(``tests/golden/regression_tile.json``).
"""
from __future__ import annotations

import numpy as np

from . import io, paths
from .features import build_hybrid_features, compute_feature_base, prepare_spectra
from .matlab_compat import find_colmajor
from .units import COMPAT, FIXED

N_PIXELS = 400


def tile_indices(n=N_PIXELS):
    """Every k-th baseline-mask pixel in MATLAB (column-major) order."""
    support = io.load_baseline_mask()[0]
    r, c = find_colmajor(support)
    step = len(r) // n
    sel = np.arange(n) * step
    return r[sel], c[sel]


def load_tile(n=N_PIXELS):
    """(dates, (n_dates, n, 8) stored stack values); cached as .npz."""
    f = paths.cache_dir() / f"regression_tile_{n}.npz"
    if f.exists():
        z = np.load(f)
        return list(z["dates"]), z["X"]
    g = io.load_glint_stack()
    r, c = tile_indices(n)
    X = np.stack([g.get(d)[r, c, :] for d in g.dates]).astype(np.float32)
    np.savez_compressed(f, dates=np.array(g.dates), X=X)
    return list(g.dates), X


def tile_summary(dates, X, gate=None, stress=None, mu=None):
    """Aggregate numbers stored as golden values (no pixel values)."""
    out = {}
    for i, d in enumerate(dates):
        ok = np.isfinite(X[i]).all(1)
        Xf = prepare_spectra(X[i][ok], FIXED)
        Fb = compute_feature_base(Xf, FIXED)
        rec = {"n_valid": int(ok.sum()),
               "mean_reflectance": np.round(Xf.mean(0).astype(np.float64), 7).tolist(),
               "floored_fixed": np.round((Xf < FIXED.floor).mean(0), 6).tolist(),
               "floored_compat": np.round((prepare_spectra(X[i][ok], COMPAT) < COMPAT.floor).mean(0), 6).tolist(),
               "mean_feature_base_fixed": np.round(Fb.mean(0).astype(np.float64), 6).tolist()}
        if gate is not None:
            from .stage2 import Stage2Options, hierarchical_predict
            mu_h, mu_b = mu
            F = build_hybrid_features(compute_feature_base(X[i][ok], COMPAT), mu_h, mu_b)
            p, _, _ = hierarchical_predict(F, gate, stress, Stage2Options(mode=COMPAT))
            rec["paper_forest_counts"] = np.bincount(p, minlength=5)[1:].tolist()
        out[d] = rec
    return out
