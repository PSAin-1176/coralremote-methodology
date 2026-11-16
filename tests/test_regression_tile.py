"""P8 regression tests on a 400-pixel tile, all 14 dates (< 10 s).

Checks: the unit contract, expected feature ranges, scale invariance, and
golden aggregates (reflectance means, floored fractions, feature means, and
the paper forests' class counts). The tile itself stays in the local cache
(Planet terms); only aggregates are in ``tests/golden``.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from coralremote import paths
from coralremote.features import compute_feature_base, prepare_spectra
from coralremote.units import COMPAT, FIXED, assert_reflectance

pytestmark = [pytest.mark.data,
              pytest.mark.skipif(not (paths.GLINT_STACK.exists() and paths.BASELINE_MASK.exists()),
                                 reason="Source stack / baseline mask not available")]
GOLDEN = json.loads((Path(__file__).parent / "golden" / "regression_tile.json").read_text())


@pytest.fixture(scope="module")
def tile():
    from coralremote import regression
    return regression.load_tile(GOLDEN["n_pixels"])


def test_tile_dates_and_unit_contract(tile):
    dates, X = tile
    assert dates == list(GOLDEN["dates"])
    for i, d in enumerate(dates):
        Xf = prepare_spectra(X[i][np.isfinite(X[i]).all(1)], FIXED)
        assert_reflectance(Xf, f"tile {d}")
        assert np.nanmedian(X[i]) < 1e-3                       # stored stack really is double-scaled


def test_feature_ranges(tile):
    dates, X = tile
    for i in range(len(dates)):
        F = compute_feature_base(prepare_spectra(X[i], FIXED), FIXED)
        F = F[np.isfinite(F).all(1)]
        norm, logs, bright = F[:, :8], F[:, 8:16], F[:, 16]
        assert np.allclose(norm.sum(1), 1, atol=1e-5) and (norm >= 0).all() and (norm <= 1).all()
        assert (logs >= np.log(np.float32(FIXED.floor)) - 1e-6).all() and (logs < np.log(1.5)).all()
        assert (bright >= 8 * FIXED.floor * 0.999).all() and (bright < 8 * 1.5).all()


def test_scale_invariance_where_floor_does_not_bind(tile):
    dates, X = tile
    for i in range(len(dates)):
        Xf = prepare_spectra(X[i], FIXED)
        keep = np.isfinite(Xf).all(1) & (Xf > 10 * FIXED.floor).all(1)
        a = compute_feature_base(Xf[keep], FIXED)
        b = compute_feature_base(Xf[keep] * np.float32(10), FIXED)
        assert np.allclose(a[:, :8], b[:, :8], rtol=1e-5, atol=1e-7)
        # compat and fixed agree on the normalised spectrum wherever neither floor binds
        keepc = keep & (Xf > 0.01).all(1)
        c = compute_feature_base(prepare_spectra(X[i][keepc], COMPAT), COMPAT)
        f = compute_feature_base(prepare_spectra(X[i][keepc], FIXED), FIXED)
        assert np.allclose(c[:, :8], f[:, :8], rtol=1e-5, atol=1e-7)


def test_golden_aggregates(tile):
    from coralremote import regression
    dates, X = tile
    got = regression.tile_summary(dates, X)
    for d in dates:
        g, w = got[d], GOLDEN["dates"][d]
        assert g["n_valid"] == w["n_valid"]
        for k in ("floored_fixed", "floored_compat"):
            assert g[k] == w[k], (d, k)
        np.testing.assert_allclose(g["mean_reflectance"], w["mean_reflectance"], rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(g["mean_feature_base_fixed"], w["mean_feature_base_fixed"], rtol=1e-5, atol=1e-6)


@pytest.mark.skipif(not paths.STAGE2_OUTPUT.exists(), reason="paper Stage 2 run not available")
def test_paper_forest_counts_on_tile(tile):
    from coralremote import io, regression
    from coralremote.forest import load_treebagger
    dates, X = tile
    ref = io.load_stage2_output()
    got = regression.tile_summary(dates, X, load_treebagger(paths.STAGE2_OUTPUT, "rfGate"),
                                  load_treebagger(paths.STAGE2_OUTPUT, "rfStress"),
                                  (ref["muHealthy"], ref["muBleached"]))
    assert {d: got[d]["paper_forest_counts"] for d in dates} == \
        {d: GOLDEN["dates"][d]["paper_forest_counts"] for d in dates}
