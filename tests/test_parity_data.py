"""Parity against the frozen paper-of-record run (needs the data under
CORALREMOTE_ROOT; skipped otherwise)."""
import numpy as np
import pytest

from coralremote import io, paths, stage2

pytestmark = pytest.mark.data
needs_data = pytest.mark.skipif(not paths.available(paths.GLINT_STACK), reason="frozen data not available")


@pytest.fixture(scope="module")
def inputs():
    return stage2.prepare_frozen_inputs(stage2.Stage2Options(mode="compat"))


@needs_data
def test_stack_and_mask_load():
    s = io.load_glint_stack()
    mask, H, W = io.load_baseline_mask()
    assert (H, W) == (1366, 1078) and mask.sum() == 20330
    assert (s.ref.x_min, s.ref.y_max, s.ref.cell_x) == (790485.0, 1399065.0, 3.0)
    assert len(s.dates) == 14 and s.get(0).shape == (1366, 1078, 8)


@needs_data
def test_glint_stack_is_double_scaled():
    """Documents the bug: the paper stack only passes the unit check once
    its extra /SR_SCALE is declared."""
    from coralremote.units import UnitError, assert_reflectance
    s = io.load_glint_stack()
    mask, _, _ = io.load_baseline_mask()
    i = s.date_index("20230404")
    with pytest.raises(UnitError):
        assert_reflectance(s.get(i), "glint", mask=mask)
    s.check_units(i, mask=mask)


@needs_data
def test_axis_means_bit_identical(inputs):
    ref = io.load_stage2_output()
    assert np.array_equal(inputs["mu_h"], ref["muHealthy"])
    assert np.array_equal(inputs["mu_b"], ref["muBleached"])


@needs_data
def test_thresholds_match_report(inputs):
    t, ps = inputs["thresholds"], inputs["ple_stats"]
    got = [t["hThr"], t["bThr"], t["rAnchor"], ps["rPleCut"], ps["hKeep"], ps["bKeep"]]
    want = [0.756828, 0.894454, 0.116670, 0.095561, 0.736184, 0.915098]
    assert [f"{g:.6f}" for g in got] == [f"{w:.6f}" for w in want]


@needs_data
def test_ple_curation_and_split(inputs):
    ref = io.load_stage2_output()
    ps, ss = inputs["ple_stats"], inputs["split_stats"]
    assert (ps["nTotalPolygons"], ps["nKeptPolygons"]) == (923, 88)
    assert np.array_equal(ps["keepMask"], np.asarray(ref["pleStats"]["keepMask"]).astype(bool).ravel())
    assert (ss["nTrainPolygons"], ss["nTestPolygons"]) == (51, 37)


@needs_data
def test_sample_counts(inputs):
    assert inputs["n_train"] == {"HPY": 9883, "BLD": 2028, "PLE": 1290}
    assert [int((inputs["Ytest"] == k).sum()) for k in (1, 2, 3)] == [7068, 831, 1279]
    assert len(inputs["anchors"]["healthy"][1]) == 1262
    assert len(inputs["anchors"]["bleached"][1]) == 434


@pytest.fixture(scope="module")
def frozen(inputs):
    from coralremote.forest import load_treebagger
    gate = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
    stress = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
    return stage2.run_frozen(gate, stress, inputs)


@needs_data
def test_frozen_test_predictions_identical(frozen):
    ref = io.load_stage2_output()
    assert np.array_equal(frozen["gate_pred"], ref["gatePredTest"])
    assert np.array_equal(frozen["stress_pred"], ref["stressPredTest"])


@needs_data
def test_frozen_confusion_matrix_and_anchors(frozen):
    from coralremote.evaluate import parse_report_confusion
    C_ref = parse_report_confusion(paths.STAGE2_REPORT.read_text())
    assert np.array_equal(frozen["metrics"]["C"], C_ref)
    assert f"{frozen['metrics']['OA']:.6f}" == "0.974841"
    assert f"{frozen['metrics']['kappa']:.6f}" == "0.934418"
    assert (round(frozen["healthy_anchor_recall"], 4), frozen["n_healthy_anchor"]) == (0.8613, 1262)
    assert (frozen["bleached_anchor_recall"], frozen["n_bleached_anchor"]) == (1.0, 434)


@needs_data
def test_frozen_table5_all_56_cells(frozen):
    from coralremote.evaluate import compare_date_summary
    n, bad, maxd, _ = compare_date_summary(frozen["date_summary"], io.read_date_summary_csv(), tol=1e-9)
    assert (n, bad) == (56, 0)


@needs_data
@pytest.mark.slow
def test_py4_training_parity_one_seed(inputs, frozen):
    """Py4: MATLAB-rule CART trained on the paper's rows (seed 1, ~75 s)."""
    from coralremote.train import train_stage2
    g, s, ts = train_stage2(inputs, 1)
    assert len(ts["ygate"]) == 6000 and len(ts["ystress"]) == 3318
    r = stage2.run_frozen(g, s, inputs, return_maps=True)
    assert abs(r["metrics"]["OA"] - frozen["metrics"]["OA"]) < 0.01
    ref_maps = stage2.run_frozen(*_frozen_forests(), inputs, return_maps=True)["maps"]
    sup = inputs["support"]
    agree = np.mean([(r["maps"][d][sup] == ref_maps[d][sup]).mean() for d in ref_maps])
    assert agree >= 0.97


def _frozen_forests():
    from coralremote.forest import load_treebagger
    return (load_treebagger(paths.STAGE2_OUTPUT, "rfGate"), load_treebagger(paths.STAGE2_OUTPUT, "rfStress"))
