"""Locations of the frozen paper-of-record inputs and outputs.

The data root defaults to ``E:/CoralRemote`` and can be overridden with the
``CORALREMOTE_ROOT`` environment variable. Nothing under ``Asset``,
``Source``, ``Results`` or ``docs/Mar6`` is ever written by this package.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("CORALREMOTE_ROOT", "E:/CoralRemote"))
ASSET = ROOT / "Asset"
SOURCE = ROOT / "Source"
RESULTS = ROOT / "Results"

READ_ONLY_DIRS = (ASSET, SOURCE, RESULTS, ROOT / "docs" / "Mar6")

GLINT_STACK = SOURCE / "SourceStack_HedleyGlintReflectance.mat"
RAW_STACK = SOURCE / "SourceStack_CoregRawReflectance.mat"
BASELINE_MASK = ASSET / "baseline_coral_mask.mat"

STAGE2_RUN = RESULTS / "stage2_hierarchical_curated_holdout_quicktest" / "run_20260408_192522"
STAGE2_OUTPUT = STAGE2_RUN / "hierarchical_curated_holdout_output.mat"
STAGE2_DATE_SUMMARY = STAGE2_RUN / "hierarchical_curated_holdout_date_summary.csv"
STAGE2_REPORT = STAGE2_RUN / "hierarchical_curated_holdout_report.txt"

LABELS = {
    "strict_train": ASSET / "stricklabel_Training_SpatialDefault.geojson",
    "strict_test": ASSET / "stricklabel_Testing_SpatialDefault.geojson",
    "strict_hold": ASSET / "stricklabel_Holdout_20240514.geojson",
    "unified_train": ASSET / "Unified_Training_SpatialDefault.geojson",
    "unified_test": ASSET / "Unified_Testing_SpatialDefault.geojson",
    "unified_hold": ASSET / "Unified_Holdout_20240514.geojson",
}


def assert_writable(path) -> Path:
    """Refuse to write inside the protected paper-of-record folders."""
    p = Path(path).resolve()
    for d in READ_ONLY_DIRS:
        try:
            p.relative_to(d.resolve())
        except ValueError:
            continue
        raise PermissionError(f"{p} is inside read-only folder {d}")
    return p


def cache_dir() -> Path:
    d = Path(os.environ.get("CORALREMOTE_CACHE", ROOT / "coralremote-python-cache"))
    assert_writable(d)
    d.mkdir(parents=True, exist_ok=True)
    return d
