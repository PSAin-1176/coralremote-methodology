"""Image-only wave diagnostic (Py5): port of ``Diag_WaveFFT_ThreeClass.m``
(NDWI ocean tiles, 2-D FFT, roughness + anisotropy, 3-class verdict).

Input is Phase 1 co-registered raw surface reflectance on a true 0-1 scale
(``single(tif) / SR_SCALE`` in MATLAB). The Phase 1 tifs were lost with
``C:/SAT``; the Source raw stack (``SourceStack_CoregRawReflectance.mat``)
holds the same co-registered reflectance on the analysis grid.

Per date: 128 x 128 tiles with >= 80% NDWI > 0.05 ocean pixels; per tile and
band, non-ocean pixels are filled with the ocean mean, the tile is
mean-centred, Hann-windowed and FFT'd. On the 6-96 m wavelength annulus:
WEF = annulus / total (DC excluded) power, WBP = mean annulus power,
AI = 1 - min/max power over 12 angular sectors. Medians over tiles give one
value per band; across dates, a band is flagged when AI (WBP) exceeds the
cross-date median + 1.5 MAD. Verdict: WAVE if >= 4 roughness and >= 6 AI
flags, TURBULENT if >= 4 roughness flags, else CALM. Confidence by the
number of ocean tiles (HIGH >= 60, MEDIUM >= 40, else LOW).
"""
from __future__ import annotations

import math
import warnings

import numpy as np

GSD_M = 3.0
NDWI_THR = 0.05
TILE_PX = 128
TILE_OCEAN_FRAC = 0.80
WAVE_LAMBDA_MIN_M, WAVE_LAMBDA_MAX_M = 6, 96
N_SECTORS = 12
AI_K = 1.5
WBP_K = 1.5
ROUGH_MIN_BANDS = 4
WAVE_MIN_FLAGS = 6
CONF_HIGH_MIN_TILES, CONF_MED_MIN_TILES = 60, 40
GREEN1, NIR = 2, 7                    # 0-based band indices (BANDS.GREEN1 = 3, NIR = 8)
BAND_NAMES = ("Coastal", "Blue", "Green1", "Green2", "Yellow", "Red", "RedEdge", "NIR")

# Paper verdicts (Code/Wave_Diagnostic_Methodology_and_Results_Note.txt and
# Date_Coral_Wave_Status_Summary_20260406.txt: 9 CALM, 1 TURBULENT, 4 WAVE;
# 13 dates HIGH support, 20230629 LOW with 38 tiles).
PAPER_VERDICTS = {
    "20230404": "CALM", "20230509": "CALM", "20230519": "WAVE", "20230629": "TURBULENT",
    "20240307": "CALM", "20240317": "CALM", "20240408": "CALM", "20240514": "CALM",
    "20240613": "WAVE", "20240616": "WAVE", "20240620": "WAVE", "20240926": "CALM",
    "20241111": "CALM", "20250518": "CALM",
}


def _hann(n):
    """MATLAB ``hann(n)`` (symmetric)."""
    k = np.arange(n)
    return 0.5 * (1 - np.cos(2 * np.pi * k / (n - 1)))


def tile_geometry(tile_px=TILE_PX):
    rb_min = max(1, math.floor(tile_px * GSD_M / WAVE_LAMBDA_MAX_M))
    rb_max = min(tile_px // 2, math.ceil(tile_px * GSD_M / WAVE_LAMBDA_MIN_M))
    c = tile_px // 2                                  # 0-based of MATLAB floor(T/2)+1
    yy, xx = np.mgrid[0:tile_px, 0:tile_px]
    rr = np.sqrt((xx - c) ** 2.0 + (yy - c) ** 2.0)
    theta = np.mod(np.arctan2(yy - c, xx - c), 2 * np.pi)
    sector = np.minimum(N_SECTORS, np.floor(theta / (2 * np.pi) * N_SECTORS) + 1).astype(int)
    ann = (rr >= rb_min) & (rr <= rb_max)
    return {"rb_min": rb_min, "rb_max": rb_max, "hann": np.outer(_hann(tile_px), _hann(tile_px)),
            "ann": ann, "dc_excl": rr > 0, "sector_masks": [ann & (sector == s) for s in range(1, N_SECTORS + 1)]}


def ocean_mask(A):
    g = A[:, :, GREEN1].astype(np.float64)
    n = A[:, :, NIR].astype(np.float64)
    den = g + n
    ndwi = np.zeros(g.shape)
    v = den > 0                                       # NaN compares False, as in MATLAB
    ndwi[v] = (g[v] - n[v]) / den[v]
    return ndwi > NDWI_THR


def date_metrics(A, geom=None, tile_px=TILE_PX):
    """Per-band tile-median WEF, WBP and AI for one (H, W, 8) image."""
    geom = geom or tile_geometry(tile_px)
    H, W, nb = A.shape
    ocean = ocean_mask(A)
    rows = range(0, H - tile_px + 1, tile_px)
    cols = range(0, W - tile_px + 1, tile_px)
    tiles = [(r, c) for r in rows for c in cols
             if ocean[r:r + tile_px, c:c + tile_px].mean() >= TILE_OCEAN_FRAC]
    n = len(tiles)
    out = {"n_tiles_used": n}
    if n == 0:
        return out
    wef = np.full((n, nb), np.nan)
    wbp = np.full((n, nb), np.nan)
    ai = np.full((n, nb), np.nan)
    ann, dc, hann = geom["ann"], geom["dc_excl"], geom["hann"]
    n_ann = max(int(ann.sum()), 1)
    for t, (r, c) in enumerate(tiles):
        ot = ocean[r:r + tile_px, c:c + tile_px]
        for b in range(nb):
            tile = A[r:r + tile_px, c:c + tile_px, b].astype(np.float64)
            vals = tile[ot]
            vals = vals[np.isfinite(vals) & (vals > 0)]
            if vals.size < 10:
                continue
            tf = tile.copy()
            tf[~ot] = vals.mean()
            tf = tf - tf.mean()
            P = np.abs(np.fft.fftshift(np.fft.fft2(tf * hann))) ** 2
            total, wave_pow = P[dc].sum(), P[ann].sum()
            if total > 0:
                wef[t, b] = wave_pow / total
            wbp[t, b] = wave_pow / n_ann
            sp = np.array([P[m].sum() for m in geom["sector_masks"]])
            smax, smin = sp.max(), sp.min()
            if np.isfinite(smax) and smax > 0 and np.isfinite(smin):
                ai[t, b] = 1 - smin / smax
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out.update(band_wef=np.nanmedian(wef, 0), band_wbp=np.nanmedian(wbp, 0), band_ai=np.nanmedian(ai, 0))
    return out


def _med_mad(col):
    v = col[np.isfinite(col)]
    if v.size < 2:
        return np.nan, np.nan
    m = np.median(v)
    return m, np.median(np.abs(v - m))


def classify(results):
    """Cross-date thresholds and verdicts; ``results`` maps date -> date_metrics.
    Returns (per-date dict, thresholds dict)."""
    dates = [d for d in results if results[d]["n_tiles_used"] > 0]
    ai_mat = np.array([results[d]["band_ai"] for d in dates])
    wbp_mat = np.array([results[d]["band_wbp"] for d in dates])
    nb = ai_mat.shape[1]
    ai_med, ai_mad, wbp_med, wbp_mad = (np.full(nb, np.nan) for _ in range(4))
    for b in range(nb):
        ai_med[b], ai_mad[b] = _med_mad(ai_mat[:, b])
        wbp_med[b], wbp_mad[b] = _med_mad(wbp_mat[:, b])
    ai_thr = ai_med + AI_K * ai_mad
    wbp_thr = wbp_med + WBP_K * wbp_mad
    out = {}
    for i, d in enumerate(dates):
        ai_flag = np.isfinite(ai_mat[i]) & (ai_mat[i] > ai_thr)
        wbp_flag = np.isfinite(wbp_mat[i]) & np.isfinite(wbp_thr) & (wbp_mat[i] > wbp_thr)
        nflag, nrough = int(ai_flag.sum()), int(wbp_flag.sum())
        if nrough >= ROUGH_MIN_BANDS and nflag >= WAVE_MIN_FLAGS:
            verdict = "WAVE"
        elif nrough >= ROUGH_MIN_BANDS:
            verdict = "TURBULENT"
        else:
            verdict = "CALM"
        nt = results[d]["n_tiles_used"]
        conf = "HIGH" if nt >= CONF_HIGH_MIN_TILES else "MEDIUM" if nt >= CONF_MED_MIN_TILES else "LOW"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out[d] = {"n_tiles_used": nt, "n_rough_bands": nrough, "n_bands_flag": nflag,
                      "verdict": verdict, "confidence": conf,
                      "ai_median": float(np.nanmedian(ai_mat[i])),
                      "wbp_median": float(np.nanmedian(wbp_mat[i])),
                      "wef_median": float(np.nanmedian(results[d]["band_wef"])),
                      "band_ai_flag": ai_flag, "band_wbp_flag": wbp_flag}
    thr = {"ai_med": ai_med, "ai_mad": ai_mad, "ai_thr": ai_thr,
           "wbp_med": wbp_med, "wbp_mad": wbp_mad, "wbp_thr": wbp_thr}
    return out, thr


def run_wave_diagnostic(get_image, dates):
    """``get_image(date)`` -> (H, W, 8) true-scale raw reflectance."""
    geom = tile_geometry()
    res = {d: date_metrics(np.asarray(get_image(d), np.float32), geom) for d in dates}
    return classify(res)
