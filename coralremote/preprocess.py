"""Deterministic preprocessing: rebuild of the A3 Phase 3-5 products that
Stage 1 used (the original ``FEATURES_<date>.tif`` files were lost with
``C:/SAT``).

Chain, from ``Phase3_LogRatioNormalize.m``, ``Phase4_WaveHandling.m``,
``Phase5_FeatureStack_LogRatio.m`` and ``buildFeatureStack_LogRatio.m``,
with the A3 config (``pipeline_logratio_config.mat``):

* Phase 3: Hedley glint-corrected reflectance (true 0-1 scale; the
  surviving ``GLINT_<date>.tif`` files, else the Source glint stack x SR_SCALE);
  ``APPLY_NORMALIZATION = false``, so NORM = glint reflectance; NaN where any
  raw band is <= 0 or non-finite, or where the prepared valid mask is false.
  The prepared valid mask is the UDM2 ``clear`` band (verified: on 20230629
  it removes exactly the 120 light-haze test pixels the paper excluded),
  read from the raw scenes in ``paths.RAW_PLANET_DIR`` and placed on the
  grid by nearest-neighbour resampling.
* Phase 4: wave pixels NaN where raw NIR > 0.07 inside the reef mask.
* Phase 5: B/G1, G1/G2, B/R, DII_BR = log B - 0.82 log R,
  DII_GR = log G1 - 0.80 log R (EPS = 1e-6 floors); NaN outside the reef mask.

Verified against the 17,372 training rows saved in the paper Stage 1
forest: pixel counts identical; ratio features bit-identical where the
GLINT tif survives (<= 2.5e-7 relative otherwise); DII within one float32
ulp of the log terms (<= 7.2e-7 absolute).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from . import io, paths
from .units import SR_SCALE

F32 = np.float32
A3_EPS = 1e-6
K_BR, K_GR = 0.82, 0.80
NIR_WAVE_THRESHOLD = 0.07
A3_FEATURES = ("B_over_G1", "G1_over_G2", "B_over_R", "DII_BR", "DII_GR")


def _log32(x):
    return np.log(np.asarray(x, np.float64)).astype(F32)


def a3_features(A):
    """Port of ``buildFeatureStack_LogRatio`` (mode A3) on (..., 8) reflectance."""
    A = np.asarray(A, F32)
    eps = F32(A3_EPS)
    B, G1, G2, R = A[..., 1], A[..., 2], A[..., 3], A[..., 5]
    return np.stack([
        B / np.fmax(G1, eps),
        G1 / np.fmax(G2, eps),
        B / np.fmax(R, eps),
        _log32(np.fmax(B, eps)) - F32(K_BR) * _log32(np.fmax(R, eps)),
        _log32(np.fmax(G1, eps)) - F32(K_GR) * _log32(np.fmax(R, eps)),
    ], axis=-1).astype(F32)


@lru_cache(maxsize=None)
def reef_mask():
    import rasterio
    with rasterio.open(paths.REEF_MASK_TIF) as ds:
        return ds.read(1) > 0


@lru_cache(maxsize=None)
def raw_stack():
    return io.StackFile(paths.RAW_STACK, "stack_raw", on_disk_scale=1.0)


def glint_reflectance(date, prefer_tif=True):
    """Phase 3 glint-corrected reflectance (H, W, 8), true 0-1 scale.

    Uses the original ``GLINT_<date>.tif`` when it survives, else undoes the
    stack builder's extra division (x SR_SCALE; may differ in the last bit)."""
    p = paths.PHASE3_GLINT_DIR / f"GLINT_{date}.tif"
    if prefer_tif and p.exists():
        import rasterio
        with rasterio.open(p) as ds:
            return np.moveaxis(ds.read(), 0, -1).astype(F32)
    return (io.load_glint_stack().get(date) * F32(SR_SCALE)).astype(F32)


@lru_cache(maxsize=None)
def udm2_clear(date):
    """UDM2 'clear' band on the analysis grid, or None if the raw scene is absent."""
    import glob
    import warnings
    import rasterio
    from rasterio.warp import Resampling, reproject
    hits = sorted(glob.glob(str(paths.RAW_PLANET_DIR / f"{date}_*_udm2_clip.tif")))
    if not hits:
        warnings.warn(f"UDM2 for {date} not found under {paths.RAW_PLANET_DIR}; validity mask uses raw > 0 only")
        return None
    ref = io.load_glint_stack().ref
    out = np.zeros((ref.height, ref.width), np.uint8)
    with rasterio.open(hits[0]) as ds:
        reproject(ds.read(1), out, src_transform=ds.transform, src_crs=ds.crs,
                  dst_transform=rasterio.Affine(ref.cell_x, 0, ref.x_min, 0, -ref.cell_y, ref.y_max),
                  dst_crs=f"EPSG:{ref.epsg}", resampling=Resampling.nearest)
    return out == 1


def valid_and_wave(date, use_udm=True):
    """Phase 3 validity and Phase 4 wave masks from the raw stack (+ UDM2 clear)."""
    raw = raw_stack().get(date)
    valid = np.all(np.isfinite(raw) & (raw > 0), axis=2)
    if use_udm:
        clear = udm2_clear(date)
        if clear is not None:
            valid &= clear
    wave = (raw[:, :, 7] > F32(NIR_WAVE_THRESHOLD)) & reef_mask() & valid
    return valid, wave


@lru_cache(maxsize=16)
def phase5_features(date, prefer_tif=True):
    """(H, W, 5) A3 features exactly as the paper's ``FEATURES_<date>.tif``
    (NaN outside reef, invalid or wave pixels)."""
    A = glint_reflectance(date, prefer_tif).copy()
    valid, wave = valid_and_wave(date)
    A[~valid | wave] = np.nan
    F = a3_features(A)
    F[~reef_mask()] = np.nan
    return F
