"""Unit contract (patch plan P2).

Convention: every array handed between pipeline stages is surface
reflectance on a 0-1 scale. Loaders state the scale of what is on disk and
convert; ``assert_reflectance`` checks the result at every boundary.

The paper-of-record Stage 2 stack (``SourceStack_HedleyGlintReflectance.mat``)
was divided by SR_SCALE twice, so it holds reflectance / 10,000. In
``compat`` mode the pipeline keeps those values (and the 1e-6 floor) to
reproduce the paper; in ``fixed`` mode it multiplies by SR_SCALE first.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SR_SCALE = 10_000.0

#: Band order of every 8-band stack (PlanetScope SuperDove).
BANDS = ("Coastal", "Blue", "Green1", "Green2", "Yellow", "Red", "RedEdge", "NIR")
VISIBLE = slice(0, 6)
NIR = 7


class UnitError(ValueError):
    """Raised when an array violates the reflectance contract."""


@dataclass(frozen=True)
class Mode:
    """Unit handling for Stage 2 features.

    compat: keep the double-scaled stack and the bare 1e-6 floor (paper of record).
    fixed:  true reflectance and an explicit, documented floor (plan P3/P4).
    """
    name: str
    stack_multiplier: float
    floor: float

    @property
    def is_compat(self) -> bool:
        return self.name == "compat"


COMPAT = Mode("compat", 1.0, 1e-6)
FIXED = Mode("fixed", SR_SCALE, 1e-4)


def get_mode(mode) -> Mode:
    if isinstance(mode, Mode):
        return mode
    return {"compat": COMPAT, "fixed": FIXED}[str(mode)]


def reflectance_stats(A, mask=None):
    """Median visible and NIR reflectance, max and finite fraction of an
    (..., 8) array, optionally restricted to ``mask`` (H x W)."""
    A = np.asarray(A)
    X = A[mask] if mask is not None else A.reshape(-1, A.shape[-1])
    fin = np.isfinite(X).all(axis=1)
    Xf = X[fin]
    if Xf.size == 0:
        return dict(finite_frac=0.0, median_vis=np.nan, median_nir=np.nan, max=np.nan)
    return dict(
        finite_frac=float(fin.mean()),
        median_vis=float(np.median(Xf[:, VISIBLE])),
        median_nir=float(np.median(Xf[:, NIR])),
        max=float(Xf.max()),
    )


def assert_reflectance(A, name="array", mask=None, *, scale=1.0,
                       vis_range=(0.005, 0.3), nir_range=(0.001, 0.1),
                       max_value=1.5, min_finite=0.5):
    """Check that ``A * scale`` looks like 0-1 surface reflectance.

    ``scale`` lets a caller check a stack whose on-disk scale is declared
    (e.g. ``scale=SR_SCALE`` for the double-scaled paper stack). Raises
    :class:`UnitError` with the offending statistic; returns the stats.
    """
    s = reflectance_stats(A, mask)
    s = {k: (v * scale if k != "finite_frac" else v) for k, v in s.items()}
    problems = []
    if s["finite_frac"] < min_finite:
        problems.append(f"finite fraction {s['finite_frac']:.3f} < {min_finite}")
    if not (vis_range[0] <= s["median_vis"] <= vis_range[1]):
        problems.append(f"median visible {s['median_vis']:.3g} outside {vis_range}")
    if not (nir_range[0] <= s["median_nir"] <= nir_range[1]):
        problems.append(f"median NIR {s['median_nir']:.3g} outside {nir_range}")
    if not (s["max"] < max_value):
        problems.append(f"max {s['max']:.3g} >= {max_value}")
    if problems:
        raise UnitError(f"{name}: not reflectance (scale={scale}): " + "; ".join(problems))
    return s
