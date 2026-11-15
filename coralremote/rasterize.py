"""GeoJSON label polygons -> pixel masks, reproducing the MATLAB Stage 2
rasteriser exactly (``rasterizeSingleFeature`` + ``extractRings`` +
``latlon2utm47`` + ``poly2mask``).

Two MATLAB quirks are reproduced on purpose (compat behaviour):

1. Only the **outer ring** of each polygon is rasterised; holes are filled.
2. ``extractRings`` relies on how ``jsondecode`` nests arrays. A polygon that
   ``jsondecode`` returns as a numeric R x N x 2 array *inside a cell*
   (MultiPolygon whose parts differ in shape, or a bare ``Polygon``) is
   indexed as ``ring(:,2)``/``ring(:,1)`` on a 3-D array, which yields
   lon/lon of single vertices, so that polygon contributes no pixels.
   None of the paper's label files contain such geometries (all are
   single-part MultiPolygons), but the rule is kept for faithfulness.

Vertices outside the grid are dropped *before* ``poly2mask``, as in MATLAB.
"""
from __future__ import annotations

import warnings

import numpy as np

from .matlab_compat import mround, poly2mask_indices


def latlon2utm47(lat, lon):
    """Port of ``Code/utils/latlon2utm47.m`` (WGS84 -> UTM 47N), same
    operation order so that results are bit-identical in float64."""
    lat = np.asarray(lat, np.float64)
    lon = np.asarray(lon, np.float64)
    a = 6378137.0
    f = 1 / 298.257223563
    k0 = 0.9996
    E0 = 500000
    lon0 = 99
    e2 = 2 * f - f ** 2
    e4 = e2 ** 2
    e6 = e2 ** 3
    d2r = np.pi / 180
    lr = d2r * lat
    wr = d2r * lon
    dw = wr - d2r * lon0
    N_ = a / np.sqrt(1 - e2 * np.sin(lr) ** 2)
    T = np.tan(lr) ** 2
    C = e2 / (1 - e2) * np.cos(lr) ** 2
    A = np.cos(lr) * dw
    M = a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lr
             - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * np.sin(2 * lr)
             + (15 * e4 / 256 + 45 * e6 / 1024) * np.sin(4 * lr)
             - (35 * e6 / 3072) * np.sin(6 * lr))
    x = k0 * N_ * (A + (1 - T + C) * A ** 3 / 6
                   + (5 - 18 * T + T ** 2 + 72 * C - 58 * e2 / (1 - e2)) * A ** 5 / 120) + E0
    y = k0 * (M + N_ * np.tan(lr) * (A ** 2 / 2
                                      + (5 - T + 9 * C + 4 * C ** 2) * A ** 4 / 24
                                      + (61 - 58 * T + T ** 2 + 600 * C - 330 * e2 / (1 - e2)) * A ** 6 / 720))
    return x, y


def _uniform(arr):
    """True if ``jsondecode`` would turn this nested list into a numeric array."""
    try:
        a = np.asarray(arr, dtype=np.float64)
    except (ValueError, TypeError):
        return False
    return a.dtype != object


def extract_rings(geom):
    """Port of ``extractRings.m`` under ``jsondecode`` nesting rules.

    Returns a list of (N, 2) [lon, lat] arrays that MATLAB would rasterise.
    """
    rings = []
    gtype = geom.get("type")
    c = geom.get("coordinates")
    if gtype == "Polygon":
        if _uniform(c):
            # numeric R x N x 2 -> ring(:,1)/ring(:,2) pick single vertices: no area
            warnings.warn("bare Polygon decoded as numeric array: contributes no pixels (MATLAB quirk)")
            return []
        rings.append(np.asarray(c[0], np.float64))
    elif gtype == "MultiPolygon":
        if _uniform(c):                       # numeric P x R x N x 2
            a = np.asarray(c, np.float64)
            for p in range(a.shape[0]):
                rings.append(a[p, 0])
            return rings
        for poly in c:                        # cell array of polygons
            if _uniform(poly):
                warnings.warn("MultiPolygon part decoded as numeric array: contributes no pixels (MATLAB quirk)")
                continue
            r = np.asarray(poly[0], np.float64)
            if r.ndim == 2 and r.shape[1] >= 2:
                rings.append(r)
    return rings


def feature_pixels(feature, ref, rule="stage2"):
    """0-based (rows, cols) MATLAB would set for one label feature.

    rule='stage2': ``round((x - xmin)/cw) + 1``; vertices outside the grid
        are dropped (``rasterizeSingleFeature`` in the Stage 2 script).
    rule='stage1': ``round((x - xmin)/cw + 0.5)`` clamped to the grid; rings
        with < 3 vertices skipped (local ``rasterizeGeoJSONLabels`` in
        ``Train_Stage1_LogRatio.m``). This is half a pixel different from
        the Stage 2 rule; both are kept as-is for parity.
    """
    H, W = ref.height, ref.width
    rr_all, cc_all = [], []
    for ring in extract_rings(feature["geometry"]):
        if rule == "stage1" and ring.shape[0] < 3:
            continue
        xu, yu = latlon2utm47(ring[:, 1], ring[:, 0])
        if rule == "stage1":
            cw = (ref.x_min + W * ref.cell_x - ref.x_min) / W
            ch = (ref.y_max - (ref.y_max - H * ref.cell_y)) / H
            col = np.clip(mround((xu - ref.x_min) / cw + 0.5), 1, W)
            row = np.clip(mround((ref.y_max - yu) / ch + 0.5), 1, H)
            rr, cc = poly2mask_indices(col, row, H, W)
            rr_all.append(rr)
            cc_all.append(cc)
            continue
        col = mround((xu - ref.x_min) / ref.cell_x) + 1
        row = mround((ref.y_max - yu) / ref.cell_y) + 1
        valid = (col >= 1) & (col <= W) & (row >= 1) & (row <= H)
        if not valid.any():
            continue
        rr, cc = poly2mask_indices(col[valid], row[valid], H, W)
        rr_all.append(rr)
        cc_all.append(cc)
    if not rr_all:
        return np.zeros(0, int), np.zeros(0, int)
    return np.concatenate(rr_all), np.concatenate(cc_all)


def rasterize_date_labels(features, ref, class_to_id, date):
    """Port of ``rasterizeDateLabels``: union mask and first-writer-wins
    label image for the features of one date."""
    H, W = ref.height, ref.width
    mask = np.zeros((H, W), bool)
    label = np.zeros((H, W), np.uint8)
    for f in features:
        p = f["properties"]
        if "date" not in p or "Class" not in p:
            continue
        if str(p["date"]) != date:
            continue
        cls = str(p["Class"])
        if cls not in class_to_id:
            continue
        rr, cc = feature_pixels(f, ref)
        if rr.size == 0:
            continue
        free = label[rr, cc] == 0
        label[rr[free], cc[free]] = class_to_id[cls]
        mask[rr, cc] = True
    return mask, label


def rasterize_labels_stage1(features, ref, class_map):
    """Port of the Stage 1 ``rasterizeGeoJSONLabels``: all features (no
    date filter), first writer wins, ``class_map`` e.g. {'HC':1,'SND':2,'RBL':3}."""
    H, W = ref.height, ref.width
    mask = np.zeros((H, W), bool)
    label = np.zeros((H, W), np.uint8)
    for f in features:
        p = f["properties"]
        cls = p.get("Class", p.get("class"))
        if cls not in class_map:
            continue
        rr, cc = feature_pixels(f, ref, rule="stage1")
        if rr.size == 0:
            continue
        free = label[rr, cc] == 0
        label[rr[free], cc[free]] = class_map[cls]
        mask[rr, cc] = True
    return mask, label
