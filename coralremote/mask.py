"""Baseline coral mask (Py5): consensus votes from the Stage 1 class maps,
then the 6 Apr morphology chain that produced ``Asset/baseline_coral_mask.mat``.

Chain (MATLAB sources):

1. ``Baseline_Coral_LogRatio.m`` (A3 run 20260401_103612): CORAL votes over
   the Phase 9 ``CLASS_<date>.tif`` maps for 20230404, 20230509 and 20240307,
   zeroed outside the reef; tier 1 = >= 1 vote, 2 = >= ceil(0.6 N) votes,
   3 = N votes; ``mask0 = tier >= 2`` (21,460 px).
2. ``Clean_CoralMask_Conservative.m`` (MinObjectPx 12, MaxHolePx 12, one
   boundary-only majority pass) -> ``maskClean`` (21,961 px).
3. ``Peel_CoralMask_OnePixel.m``: erosion with ``strel('disk', 1, 0)``
   (3x3 cross) -> ``maskPeel`` (20,273 px).
4. ``Finalize_CoralMask_MinorFixes.m`` (MinObjectPx 4, MaxHolePx 40)
   -> ``maskFinal`` = the paper's baseline mask (20,330 px).

MATLAB defaults reproduced: ``bwareaopen`` 8-connected; ``imfill 'holes'``
4-connected background; ``bwperim`` 4-connected; ``bwmorph majority`` = 1 if
>= 5 of the 3x3 neighbourhood are 1 (outside = 0); ``imerode`` pads with 1.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi

from . import paths

BASELINE_DATES = ("20230404", "20230509", "20240307")
CORAL_ID = 1
MOD_FRAC = 0.6
PHASE9_DIR = paths.A3_ROOT / "Phase9_Classify" / "run_20260401_102734"

_CONN8 = np.ones((3, 3), bool)
_CONN4 = ndi.generate_binary_structure(2, 1)


# ---- MATLAB morphology primitives -------------------------------------------
def bwareaopen(bw, p, conn=8):
    """Remove connected components with fewer than ``p`` pixels."""
    lab, n = ndi.label(bw, _CONN8 if conn == 8 else _CONN4)
    if n == 0:
        return np.zeros_like(bw, bool)
    size = np.bincount(lab.ravel())
    keep = size >= p
    keep[0] = False
    return keep[lab]


def imfill_holes(bw):
    """``imfill(bw, 'holes')`` (2-D, 4-connected background)."""
    return ndi.binary_fill_holes(bw, _CONN4)


def bwperim(bw):
    """``bwperim`` with 4-connectivity: foreground pixels with a 4-neighbour
    background pixel (pixels beyond the image count as background)."""
    bw = np.asarray(bw, bool)
    return bw & ~ndi.binary_erosion(bw, _CONN4, border_value=0)


def majority(bw):
    """One pass of ``bwmorph(bw, 'majority')``."""
    s = ndi.correlate(np.asarray(bw, np.uint8), np.ones((3, 3), np.uint8), mode="constant", cval=0)
    return s >= 5


def imerode(bw, se):
    """Binary ``imerode`` (outside the image counts as foreground)."""
    return ndi.binary_erosion(bw, se, border_value=1)


def disk_strel(r):
    """``strel('disk', r, 0)``: pixels with Euclidean distance <= r."""
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def fill_small_holes(bw, max_hole_px):
    filled = imfill_holes(bw)
    holes = filled & ~bw
    return bw | (holes & ~bwareaopen(holes, max_hole_px + 1))


# ---- chain ------------------------------------------------------------------
def consensus(class_maps, reef, coral_id=CORAL_ID, mod_frac=MOD_FRAC):
    """Vote and tier maps (``Baseline_Coral_LogRatio.m``); returns
    ``vote, tier, mask0``."""
    n = len(class_maps)
    vote = np.zeros(reef.shape, np.uint8)
    for c in class_maps:
        vote += (np.asarray(c) == coral_id).astype(np.uint8)
    vote[~reef] = 0
    mod = max(1, math.ceil(n * mod_frac))
    tier = np.zeros(reef.shape, np.uint8)
    tier[reef & (vote >= 1)] = 1
    tier[reef & (vote >= mod)] = 2
    tier[reef & (vote >= n)] = 3
    return vote, tier, tier >= 2


def clean_conservative(mask0, min_object_px=12, max_hole_px=12, boundary_dilate_px=1, majority_passes=1):
    m = fill_small_holes(bwareaopen(mask0, min_object_px), max_hole_px)
    for _ in range(majority_passes):
        zone = ndi.binary_dilation(bwperim(m), np.ones((2 * boundary_dilate_px + 1,) * 2, bool))
        maj = majority(m)
        m = m.copy()
        m[zone] = maj[zone]
    return fill_small_holes(bwareaopen(m, min_object_px), max_hole_px)


def peel(mask, radius_px=1):
    return imerode(mask, disk_strel(radius_px))


def finalize(mask, min_object_px=4, max_hole_px=40):
    return fill_small_holes(bwareaopen(mask, min_object_px), max_hole_px)


def build_baseline_mask(class_maps, reef):
    """Full chain; returns a dict with every intermediate the MATLAB scripts saved."""
    vote, tier, mask0 = consensus(class_maps, reef)
    clean = clean_conservative(mask0)
    pl = peel(clean)
    final = finalize(pl)
    return {"vote_map": vote, "tier_map": tier, "mask0": mask0, "maskClean": clean,
            "maskPeel": pl, "maskFinal": final}


def load_phase9_class_maps(dates=BASELINE_DATES, folder=PHASE9_DIR):
    import rasterio
    out = []
    for d in dates:
        with rasterio.open(folder / f"CLASS_{d}.tif") as ds:
            out.append(ds.read(1))
    return out
