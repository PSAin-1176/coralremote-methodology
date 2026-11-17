"""Readers for the frozen MATLAB-era products.

* v7.3 (HDF5) stacks via ``h5py``. HDF5 stores MATLAB arrays transposed, so
  ``stack_glint`` (MATLAB H x W x 8 x N) appears as (N, 8, W, H) and is
  returned here as per-date (H, W, 8) float32.
* v5 ``.mat`` files via ``scipy.io``.
* GeoJSON label files via ``json``.

All loaders are read-only.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np
import scipy.io as sio

from . import paths
from .units import SR_SCALE, assert_reflectance


@dataclass(frozen=True)
class RasterRef:
    """Subset of MATLAB ``MapCellsReference`` the pipeline uses."""
    x_min: float        # XWorldLimits(1), outer edge of the first column
    y_max: float        # YWorldLimits(2), outer edge of the first row
    cell_x: float       # CellExtentInWorldX
    cell_y: float       # CellExtentInWorldY
    height: int
    width: int
    epsg: int = 32647

    @property
    def transform(self):
        """Affine (rasterio-style) transform tuple."""
        return (self.cell_x, 0.0, self.x_min, 0.0, -self.cell_y, self.y_max)


def _h5_str(f, ref):
    return "".join(map(chr, f[ref][()].ravel()))


def _read_mapcells_ref(f) -> RasterRef:
    """Decode the MapCellsReference property struct stored under #refs#."""
    for key in f["#refs#"]:
        g = f["#refs#"][key]
        if isinstance(g, h5py.Group) and "TiePointWorld" in g:
            tw = g["TiePointWorld"][()].ravel()
            ti = g["TiePointIntrinsic"][()].ravel()
            num = g["DeltaNumerator"][()].ravel()
            den = g["DeltaDenominator"][()].ravel()
            size = g["RasterSize"][()].ravel()
            dx, dy = num[0] / den[0], num[1] / den[1]   # dy < 0 (rows go south)
            x_min = tw[0] - (ti[0] - 0.5) * dx
            y_max = tw[1] - (ti[1] - 0.5) * dy
            return RasterRef(float(x_min), float(y_max), float(abs(dx)), float(abs(dy)),
                             int(size[0]), int(size[1]))
    raise KeyError("no MapCellsReference property struct found")


# --------------------------------------------------------------------------
# derived-array cache: keyed on the source file's identity and validated
def cache_path_for(src, var):
    """Cache file for variable ``var`` of ``src``. The name includes a hash of
    the resolved source path, so same-named files in different folders never
    share an entry."""
    key = hashlib.sha1(str(Path(src).resolve()).lower().encode()).hexdigest()[:12]
    return paths.cache_dir() / f"{Path(src).stem}__{var}__{key}.npy"


def cache_meta(src, var, shape, dates):
    st = Path(src).stat()
    return {"source": str(Path(src).resolve()), "var": var, "size": int(st.st_size),
            "mtime_ns": int(st.st_mtime_ns), "shape": [int(x) for x in shape], "dates": list(map(str, dates))}


def _meta_path(cache_file):
    return Path(cache_file).with_suffix(".meta.json")


def write_cache_meta(cache_file, meta):
    _meta_path(cache_file).write_text(json.dumps(meta, indent=1))


def cache_is_valid(cache_file, meta):
    """True only if the cache file exists, its sidecar matches ``meta``
    (source path, size, modification time, shape, dates) and the stored array
    has that shape."""
    cache_file, mp = Path(cache_file), _meta_path(cache_file)
    if not (cache_file.exists() and mp.exists()):
        return False
    try:
        if json.loads(mp.read_text()) != meta:
            return False
        return list(np.load(cache_file, mmap_mode="r").shape) == meta["shape"]
    except (OSError, ValueError):
        return False


class StackFile:
    """Lazy reader for an (H, W, 8, N) v7.3 stack.

    ``on_disk_scale`` is the factor that turns stored values into
    reflectance (10,000 for the double-scaled paper stack, 1 otherwise);
    it is recorded so callers can check units explicitly.
    """

    def __init__(self, path, var, on_disk_scale=1.0, cache=True):
        self.path = Path(path)
        self.var = var
        self.on_disk_scale = on_disk_scale
        with h5py.File(self.path, "r") as f:
            self.dates = [_h5_str(f, r) for r in f["dates_col"][()].ravel()]
            self.band_names = ([_h5_str(f, r) for r in f["rawBandNames"][()].ravel()]
                               if "rawBandNames" in f else None)
            self.ref = _read_mapcells_ref(f)
            self.shape_h5 = f[var].shape
        self._use_cache = cache
        self._arr = None

    def _expected_shape(self):
        n, b, w, h = self.shape_h5
        return [int(n), int(h), int(w), int(b)]

    def _cache_path(self):
        return cache_path_for(self.path, self.var)

    def _load_all(self):
        if self._arr is not None:
            return self._arr
        cp = self._cache_path() if self._use_cache else None
        meta = cache_meta(self.path, self.var, self._expected_shape(), self.dates)
        if cp is not None and cache_is_valid(cp, meta):
            self._arr = np.load(cp, mmap_mode="r")
            return self._arr
        with h5py.File(self.path, "r") as f:
            a = f[self.var][()]                   # (N, 8, W, H)
        a = np.ascontiguousarray(a.transpose(0, 3, 2, 1))   # (N, H, W, 8)
        if cp is not None:
            np.save(cp, a)
            write_cache_meta(cp, meta)
            a = np.load(cp, mmap_mode="r")
        self._arr = a
        return a

    def date_index(self, date):
        return self.dates.index(str(date))

    def get(self, date_or_index) -> np.ndarray:
        """(H, W, 8) float32 array for one date, exactly as stored."""
        i = date_or_index if isinstance(date_or_index, int) else self.date_index(date_or_index)
        return np.asarray(self._load_all()[i], dtype=np.float32)

    def check_units(self, i, mask=None):
        return assert_reflectance(self.get(i), f"{self.path.name}[{self.dates[i]}]",
                                  mask=mask, scale=self.on_disk_scale)


@lru_cache(maxsize=None)
def load_glint_stack(path=None) -> StackFile:
    """Paper-of-record Stage 2 input (double-scaled: reflectance / SR_SCALE)."""
    return StackFile(path or paths.GLINT_STACK, "stack_glint", on_disk_scale=SR_SCALE)


def load_baseline_mask(path=None):
    """Baseline coral mask (H x W bool) plus (H, W)."""
    m = sio.loadmat(str(path or paths.BASELINE_MASK),
                    variable_names=["baseline_coral_mask", "H", "W"])
    mask = m["baseline_coral_mask"].astype(bool)
    H, W = int(m["H"].ravel()[0]), int(m["W"].ravel()[0])
    assert mask.shape == (H, W)
    return mask, H, W


@lru_cache(maxsize=None)
def load_geojson(path):
    """Feature list of a GeoJSON FeatureCollection (order preserved)."""
    with open(path, "r", encoding="utf-8") as fh:
        return tuple(json.load(fh)["features"])


def load_labels(name):
    return load_geojson(str(paths.LABELS[name]))


def _plain(v):
    """Unwrap scipy.io MATLAB structs / cells / scalars into Python values."""
    if isinstance(v, np.ndarray):
        if v.dtype.names:
            if v.size == 1:
                return {n: _plain(v[n].ravel()[0]) for n in v.dtype.names}
            return [{n: _plain(e[n]) for n in v.dtype.names} for e in v.ravel()]
        if v.dtype == object:
            return [_plain(e) for e in v.ravel()]
        if v.dtype.kind in "US":
            return str(v.ravel()[0]) if v.size else ""
        if v.size == 1:
            return v.ravel()[0]
        return v
    return v


def load_stage2_output(path=None):
    """Plain (non-object) variables of the paper Stage 2 output .mat."""
    names = ["muHealthy", "muBleached", "gatePredTest", "stressPredTest", "metricsTest",
             "dateSummary", "pleStats", "pleSplitStats", "healthyAnchorRecall",
             "bleachedAnchorRecall"]
    m = sio.loadmat(str(path or paths.STAGE2_OUTPUT), variable_names=names)
    out = {k: _plain(m[k]) for k in names if k in m}
    out["muHealthy"] = np.asarray(m["muHealthy"], np.float32).ravel()
    out["muBleached"] = np.asarray(m["muBleached"], np.float32).ravel()
    out["gatePredTest"] = m["gatePredTest"].ravel()
    out["stressPredTest"] = m["stressPredTest"].ravel()
    return out


def read_date_summary_csv(path=None):
    """Paper Table 5 source: list of (date, HPY, PLE, BLD, UNK) as float64."""
    rows = []
    with open(path or paths.STAGE2_DATE_SUMMARY, "r", encoding="utf-8") as fh:
        header = fh.readline().strip().split(",")
        for line in fh:
            if line.strip():
                parts = line.strip().split(",")
                rows.append({header[0]: parts[0],
                             **{h: float(x) for h, x in zip(header[1:], parts[1:])}})
    return rows
