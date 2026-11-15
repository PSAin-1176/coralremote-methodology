"""Bit-faithful ports of the MATLAB built-ins whose exact behaviour the
paper-of-record results depend on.

Everything here is deterministic and documented against the MATLAB source
semantics, so that the Python port reproduces MATLAB outputs exactly:

* ``mround``      - MATLAB ``round`` (half away from zero; numpy rounds half to even).
* ``mquantile``   - MATLAB ``quantile`` default method (midpoint plotting positions,
                    linear interpolation, clamped at the ends).
* ``poly2mask``   - MATLAB ``poly2mask`` (5x5 sub-pixel edge list + parity scan).
* ``fortran_nonzero`` - ``find`` order (column-major linear indices).
"""
from __future__ import annotations

import numpy as np


def mround(x):
    """MATLAB ``round``: halves are rounded away from zero."""
    x = np.asarray(x, dtype=np.float64)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def mquantile(x, p):
    """MATLAB ``quantile(x, p)`` for a vector ``x`` (default 'exact' method).

    Sorted values sit at cumulative probabilities (i - 0.5)/n; values between
    are linearly interpolated, values outside are clamped to min / max.
    NaNs are ignored, as in MATLAB.
    """
    x = np.asarray(x).ravel()
    x = np.sort(x[~np.isnan(x)])
    n = x.size
    if n == 0:
        return np.nan
    out_dtype = np.float32 if x.dtype == np.float32 else np.float64
    xd = x.astype(np.float64)
    pos = n * float(p) + 0.5          # 1-based fractional index
    if pos < 1:
        return out_dtype(xd[0])
    if pos >= n:
        return out_dtype(xd[-1])
    lo = int(np.floor(pos))
    frac = pos - lo
    return out_dtype(xd[lo - 1] + frac * (xd[lo] - xd[lo - 1]))


def _intline(x1, x2, y1, y2):
    """Port of images/private/intline.m (integer line between two points)."""
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    if dx == 0 and dy == 0:
        return np.array([x1], dtype=np.float64), np.array([y1], dtype=np.float64)
    flip = False
    if dx >= dy:
        if x1 > x2:
            x1, x2, y1, y2 = x2, x1, y2, y1
            flip = True
        m = (y2 - y1) / (x2 - x1)
        x = np.arange(x1, x2 + 1, dtype=np.float64)
        y = mround(y1 + m * (x - x1))
    else:
        if y1 > y2:
            x1, x2, y1, y2 = x2, x1, y2, y1
            flip = True
        m = (x2 - x1) / (y2 - y1)
        y = np.arange(y1, y2 + 1, dtype=np.float64)
        x = mround(x1 + m * (y - y1))
    if flip:
        x = x[::-1]
        y = y[::-1]
    return x, y


def _poly2edgelist(x, y, scale=5):
    x = mround(scale * (np.asarray(x, np.float64) - 0.5) + 1)
    y = mround(scale * (np.asarray(y, np.float64) - 0.5) + 1)
    xs, ys = [], []
    for k in range(len(x) - 1):
        a, b = _intline(x[k], x[k + 1], y[k], y[k + 1])
        xs.append(a)
        ys.append(b)
    if not xs:
        return np.zeros(0), np.zeros(0)
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    d = np.diff(x)
    ei = np.nonzero(d)[0]
    xe = x[ei].copy()
    xe[d[ei] < 0] -= 1
    ye = np.minimum(y[ei], y[ei + 1])
    return xe, ye


def poly2mask_indices(x, y, M, N, scale=5):
    """Return (rows, cols) 0-based of the pixels MATLAB ``poly2mask(x,y,M,N)``
    sets true. ``x`` are 1-based column coords, ``y`` 1-based row coords."""
    x = np.asarray(x, np.float64).ravel()
    y = np.asarray(y, np.float64).ravel()
    if x.size == 0:
        return np.zeros(0, int), np.zeros(0, int)
    if x[-1] != x[0] or y[-1] != y[0]:
        x = np.append(x, x[0])
        y = np.append(y, y[0])
    xe, ye = _poly2edgelist(x, y, scale)
    shift = (scale - 1) / 2
    xe = (xe + shift) / scale
    keep = xe == np.floor(xe)
    xe = xe[keep]
    ye = ye[keep]
    ye = np.ceil((ye + shift) / scale)
    bad = (xe < 1) | (xe > N) | (ye > M)
    xe = xe[~bad].astype(np.int64)
    ye = np.maximum(ye[~bad], 1).astype(np.int64)
    if xe.size == 0:
        return np.zeros(0, int), np.zeros(0, int)
    c0, c1 = xe.min(), xe.max()
    r1 = ye.max()
    # Parity scan only over the bounding columns, rows 1..r1 (rows above the
    # first edge in a column have even parity, so starting at row 1 is exact).
    S = np.zeros((r1, c1 - c0 + 1), dtype=np.int64)
    np.add.at(S, (ye - 1, xe - c0), 1)
    BW = (np.cumsum(S, axis=0) % 2) == 1
    BW = BW[:M]                        # drop row M+1
    rr, cc = np.nonzero(BW)
    return rr, cc + (c0 - 1)


def poly2mask(x, y, M, N):
    """Full boolean mask equivalent of MATLAB ``poly2mask``."""
    bw = np.zeros((M, N), dtype=bool)
    rr, cc = poly2mask_indices(x, y, M, N)
    bw[rr, cc] = True
    return bw


def find_colmajor(mask):
    """0-based (row, col) of true pixels in MATLAB ``find`` (column-major) order."""
    cols, rows = np.nonzero(np.asarray(mask).T)
    return rows, cols
