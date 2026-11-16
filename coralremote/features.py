"""Stage 2 spectral features (21 per pixel), single precision like MATLAB.

Port of ``computeFeatureBase``, ``projectToAxis``, ``spectralAngle`` and
``buildHybridFeatures`` from ``Run_Stage2_Hierarchical_CuratedHoldout_QuickTest.m``.

Feature order (columns):
  0-7   normalised spectrum  X / sum(X)
  8-15  log spectrum         log(X)
  16    brightness           sum(X)
  17    axis position t      projection on the healthy->bleached axis
  18    axis residual        distance from that axis
  19    spectral angle to the healthy mean
  20    spectral angle to the bleached mean

All arithmetic is float32 and row sums are accumulated left to right, as
MATLAB does for ``sum(X, 2)`` on a single matrix. ``log`` and ``acos`` are
evaluated in double and rounded to single: against the 9,318 MATLAB training
rows saved in the forests this leaves 1-ulp differences (<=1.1e-7 relative)
in 451 rows (numpy's float32 versions: 1,525 rows); all other features are
bit-identical.
"""
from __future__ import annotations

import numpy as np

from .units import COMPAT, FIXED, get_mode

F32 = np.float32


def _log32(X):
    return np.log(np.asarray(X, np.float64)).astype(F32)


def rowsum(X):
    """Sequential float32 row sum (MATLAB ``sum(X,2)`` on single)."""
    X = np.asarray(X, F32)
    acc = X[:, 0].copy()
    for j in range(1, X.shape[1]):
        acc += X[:, j]
    return acc


def prepare_spectra(X, mode=FIXED):
    """Bring stored stack values to the feature scale of ``mode``.

    compat: values untouched (double-scaled stack, as in the paper).
    fixed:  multiplied by SR_SCALE, i.e. true reflectance.
    """
    mode = get_mode(mode)
    X = np.asarray(X, F32)
    if mode.stack_multiplier != 1.0:
        X = X * F32(mode.stack_multiplier)
    return X


def compute_feature_base(X, mode=FIXED):
    """[normSpec(8), logSpec(8), brightness(1)] from (n, 8) spectra that are
    already on the ``mode`` scale (see :func:`prepare_spectra`)."""
    mode = get_mode(mode)
    floor = F32(mode.floor)
    X = np.fmax(np.asarray(X, F32), floor)          # MATLAB max() ignores NaN
    s = rowsum(X)
    norm = X / np.fmax(s, floor)[:, None]
    logs = _log32(X)
    return np.concatenate([norm, logs, s[:, None]], axis=1).astype(F32, copy=False)


def floored_fraction(X, mode=FIXED):
    """Per-band fraction of values that the floor of ``mode`` raises."""
    mode = get_mode(mode)
    X = np.asarray(X, F32)
    return (X < F32(mode.floor)).mean(axis=0)


def axis_vector(mu_h, mu_b):
    mu_h = np.asarray(mu_h, F32).ravel()
    mu_b = np.asarray(mu_b, F32).ravel()
    v = (mu_b - mu_h).astype(F32)
    vn = rowsum((v * v)[None, :])[0]
    vnorm2 = F32(max(vn, F32(1e-8)))
    return v, vnorm2


def project_to_axis(X, mu_h, v, vnorm2):
    X = np.asarray(X, F32)
    mu_h = np.asarray(mu_h, F32)
    centered = X - mu_h
    t = rowsum(centered * v) / F32(vnorm2)
    proj = mu_h + t[:, None] * v
    d = X - proj
    residual = np.sqrt(rowsum(d * d))
    return t.astype(F32), residual.astype(F32)


def spectral_angle(X, ref):
    X = np.asarray(X, F32)
    ref = np.asarray(ref, F32).ravel()
    num = rowsum(X * ref)
    den = np.sqrt(rowsum(X * X)) * np.sqrt(rowsum((ref * ref)[None, :])[0])
    cosang = np.fmax(np.fmin(num / np.fmax(den, F32(1e-12)), F32(1)), F32(-1))
    return np.arccos(cosang.astype(np.float64)).astype(F32)


def build_hybrid_features(Fbase, mu_h, mu_b, v=None, vnorm2=None):
    """Append axis position, residual and the two spectral angles (21 cols)."""
    if v is None:
        v, vnorm2 = axis_vector(mu_h, mu_b)
    spec = Fbase[:, :8]
    t, r = project_to_axis(spec, mu_h, v, vnorm2)
    a_h = spectral_angle(spec, mu_h)
    a_b = spectral_angle(spec, mu_b)
    return np.concatenate([Fbase, t[:, None], r[:, None], a_h[:, None], a_b[:, None]],
                          axis=1).astype(F32, copy=False)


FEATURE_NAMES = (
    [f"norm_{b}" for b in ("C", "B", "G1", "G2", "Y", "R", "RE", "NIR")]
    + [f"log_{b}" for b in ("C", "B", "G1", "G2", "Y", "R", "RE", "NIR")]
    + ["brightness", "axis_t", "axis_residual", "angle_healthy", "angle_bleached"]
)
