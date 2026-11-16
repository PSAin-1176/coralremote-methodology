"""Py6 (P3/P4): fixed units and floor."""
import numpy as np
import pytest

from coralremote import io, paths
from coralremote.features import compute_feature_base, prepare_spectra
from coralremote.units import COMPAT, FIXED


def test_norm_spectrum_scale_invariant_when_floor_does_not_bind():
    rng = np.random.default_rng(0)
    X = rng.uniform(0.002, 0.2, (200, 8)).astype(np.float32)
    a = compute_feature_base(X, FIXED)
    b = compute_feature_base(X * np.float32(10), FIXED)
    assert np.allclose(a[:, :8], b[:, :8], rtol=1e-6)               # normalised spectrum
    assert np.allclose(b[:, 8:16] - a[:, 8:16], np.log(10), atol=1e-5)  # log spectrum shifts by log 10


def test_compat_floor_equals_one_percent_reflectance():
    true = np.full((1, 8), 0.005, np.float32)                       # 0.5 % reflectance
    stored = true / np.float32(10000)                              # double-scaled stack
    fc = compute_feature_base(prepare_spectra(stored, COMPAT), COMPAT)
    ff = compute_feature_base(prepare_spectra(stored, FIXED), FIXED)
    assert np.allclose(np.exp(fc[0, 8:16]) * 1e4, 0.01, rtol=1e-5)  # clipped to 0.01
    assert np.allclose(np.exp(ff[0, 8:16]), 0.005, rtol=1e-5)       # kept


@pytest.mark.data
@pytest.mark.skipif(not (paths.GLINT_STACK.exists() and paths.RAW_STACK.exists()), reason="stacks not available")
@pytest.mark.parametrize("date", ["20240307", "20240620", "20241111"])
def test_fixed_stack_nir_equals_raw_and_floor_rarely_binds(date):
    from coralremote.preprocess import raw_stack
    G = prepare_spectra(io.load_glint_stack().get(date), FIXED)
    R = raw_stack().get(date)
    ok = np.isfinite(G[..., 7]) & np.isfinite(R[..., 7]) & (R[..., 7] > 0)
    assert (np.abs(G[..., 7][ok] - R[..., 7][ok]) <= np.spacing(R[..., 7][ok])).all()   # Hedley leaves NIR unchanged; 1 ulp
    m = io.load_baseline_mask()[0]
    X = G[m]
    X = X[np.isfinite(X).all(1)]
    assert ((X < FIXED.floor).mean(0) < 0.01).all()                 # P4: no band floored on > 1 % of pixels
