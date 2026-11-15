"""Stage 1 parity against the paper run (needs the OneDrive A3 outputs and
the raw scenes on the backup drive for UDM2; skipped otherwise)."""
import numpy as np
import pytest

from coralremote import paths

pytestmark = pytest.mark.data
needs = pytest.mark.skipif(not (paths.STAGE1_MODEL.exists() and paths.RAW_PLANET_DIR.exists()
                                and paths.RAW_STACK.exists()), reason="Stage 1 inputs not available")


@pytest.fixture(scope="module")
def s1():
    from matio import load_from_mat
    from coralremote import stage1
    from coralremote.forest import load_treebagger
    mat = load_from_mat(str(paths.STAGE1_MODEL))
    forest, props = load_treebagger(None, "rf_model", return_object=True, mat=mat)
    return mat, forest, props, stage1.prepare_inputs()


@needs
def test_training_rows_reproduced(s1):
    mat, forest, props, inp = s1
    PX = np.asarray(props["PrivateX"], np.float32)
    Y = np.array([int(str(np.asarray(y).ravel()[0])) for y in np.asarray(props["Y"]).ravel()])
    assert inp["Xtrain"].shape == PX.shape and np.array_equal(inp["Ytrain"], Y)
    assert np.abs(inp["Xtrain"] - PX).max() < 1e-6      # ratios exact or 1 ulp; DII 1 ulp of log


@needs
def test_frozen_oob_secondary_and_transfer(s1):
    from coralremote import stage1
    mat, forest, props, inp = s1
    assert round(stage1.oob_error(forest, inp["Xtrain"], inp["Ytrain"], props["OOBIndices"]), 4) == 0.1004
    r = stage1.run_frozen(forest, inp)
    assert np.array_equal(r["secondary"]["C"], np.asarray(mat["secondary_eval"]["C"].ravel()[0], float))
    for row, q in zip(r["transfer"], mat["transfer_results"].ravel()):
        want = [float(np.asarray(q[k]).ravel()[0]) for k in ("n", "OA", "kappa", "coral_prec", "coral_rec", "coral_f1")]
        got = [row[k] for k in ("n", "OA", "kappa", "coral_prec", "coral_rec", "coral_f1")]
        assert np.allclose(got, want, atol=1e-12, equal_nan=True), row["date"]
