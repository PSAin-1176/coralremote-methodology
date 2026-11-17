"""Locations of the frozen paper-of-record inputs and outputs.

Machine-specific locations come from, in order: environment variables
(``CORALREMOTE_ROOT``, ``CORALREMOTE_A3_ROOT``, ``CORALREMOTE_RAW_PLANET``,
``CORALREMOTE_CACHE``, ``CORALREMOTE_ACA_BATHY``), then a git-ignored ``paths.local.json`` at the repo
root (keys ``root``, ``a3_root``, ``raw_planet``, ``cache``, ``aca_bathymetry``; see
``paths.example.json``), then defaults relative to the data root.
Nothing under ``Asset``, ``Source``, ``Results``, ``docs/Mar6``, the A3
outputs or the raw scenes is ever written by this package.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

_LOCAL_FILE = Path(__file__).resolve().parent.parent / "paths.local.json"
_LOCAL = json.loads(_LOCAL_FILE.read_text(encoding="utf-8")) if _LOCAL_FILE.exists() else {}


def _setting(env, key, default):
    return Path(os.environ.get(env) or _LOCAL.get(key) or default)


ROOT = _setting("CORALREMOTE_ROOT", "root", "data")
ASSET = ROOT / "Asset"
SOURCE = ROOT / "Source"
RESULTS = ROOT / "Results"

READ_ONLY_DIRS = (ASSET, SOURCE, RESULTS, ROOT / "docs" / "Mar6")

GLINT_STACK = SOURCE / "SourceStack_HedleyGlintReflectance.mat"
RAW_STACK = SOURCE / "SourceStack_CoregRawReflectance.mat"
BASELINE_MASK = ASSET / "baseline_coral_mask.mat"
# Allen Coral Atlas bathymetry GeoTIFF (only for the R1 mask-depth summary)
ACA_BATHY = _setting("CORALREMOTE_ACA_BATHY", "aca_bathymetry", ROOT / "Allen" / "bathymetry_0.tif")

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
    d = _setting("CORALREMOTE_CACHE", "cache", ROOT / "coralremote-python-cache")
    assert_writable(d)
    d.mkdir(parents=True, exist_ok=True)
    return d

# Paper Stage 1 run and the A3 pipeline outputs that survive (a copy of the
# original Output_LogRatio/A3 folder, made 5 Apr 2026).
A3_ROOT = _setting("CORALREMOTE_A3_ROOT", "a3_root", ROOT / "A3")
STAGE1_RUN = A3_ROOT / "Phase_Stage1" / "run_20260401_170150"
STAGE1_MODEL = STAGE1_RUN / "rf_stage1_model.mat"
STAGE1_REPORT = STAGE1_RUN / "stage1_report.txt"
PHASE3_GLINT_DIR = A3_ROOT / "Phase3_LogRatio" / "run_20260401_102508"
REEF_MASK_TIF = A3_ROOT / "Phase2_Masking" / "run_20260401_102458" / "reef_mask.tif"
STAGE1_LABELS = {"train": ASSET / "TrainingStage1.geojson", "test": ASSET / "TestingStage1.geojson"}
READ_ONLY_DIRS = READ_ONLY_DIRS + (A3_ROOT,)
# Raw PlanetScope scenes + UDM2 (removable backup drive; may be absent).
RAW_PLANET_DIR = _setting("CORALREMOTE_RAW_PLANET", "raw_planet", ROOT / "RawPlanet")
READ_ONLY_DIRS = READ_ONLY_DIRS + (RAW_PLANET_DIR,)
