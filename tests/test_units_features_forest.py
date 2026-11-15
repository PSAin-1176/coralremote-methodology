import numpy as np
import pytest

from coralremote import features as fe
from coralremote.forest import Forest, tree_from_matlab
from coralremote.units import COMPAT, FIXED, SR_SCALE, UnitError, assert_reflectance


def _fake_reflectance(n=500, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.uniform(0.02, 0.15, (n, 8)).astype(np.float32)
    X[:, 7] = rng.uniform(0.002, 0.05, n)
    return X


def test_assert_reflectance_accepts_true_reflectance():
    assert_reflectance(_fake_reflectance()[None], "ok")


def test_assert_reflectance_rejects_double_scaled():
    X = _fake_reflectance() / SR_SCALE
    with pytest.raises(UnitError):
        assert_reflectance(X[None], "double-scaled")
    assert_reflectance(X[None], "declared", scale=SR_SCALE)   # passes once scale is declared


def test_compat_floor_binds_on_double_scaled_values():
    R = _fake_reflectance()
    frac = fe.floored_fraction(R / SR_SCALE, COMPAT)
    # the 1e-6 floor on double-scaled data clips everything below 0.01 reflectance
    assert np.isclose(frac[7], (R[:, 7] < 0.01).mean(), atol=0.01) and frac[7] > 0.1
    assert fe.floored_fraction(fe.prepare_spectra(X, FIXED), FIXED).max() == 0


def test_normalised_spectrum_scale_invariant_without_floor():
    X = _fake_reflectance()
    a = fe.compute_feature_base(X, FIXED)[:, :8]
    b = fe.compute_feature_base(X * 10, FIXED)[:, :8]
    assert np.allclose(a, b, rtol=1e-6)


def test_hybrid_features_shape_and_axis():
    X = _fake_reflectance()
    Fb = fe.compute_feature_base(X, FIXED)
    mu_h, mu_b = Fb[:250, :8].mean(0), Fb[250:, :8].mean(0) + 0.01
    F = fe.build_hybrid_features(Fb, mu_h, mu_b)
    assert F.shape == (len(X), 21) and F.dtype == np.float32
    v, vn = fe.axis_vector(mu_h, mu_b)
    t, r = fe.project_to_axis(np.stack([mu_h, mu_b]), mu_h, v, vn)
    assert np.allclose(t, [0, 1], atol=1e-5) and np.allclose(r, 0, atol=1e-5)


def test_forest_evaluator_matlab_conventions():
    # node1: x0 < 0.5 -> node2 (leaf, class 1) else node3 (leaf, class 2)
    t1 = tree_from_matlab([1, 0, 0], [0.5, 0, 0], [[2, 3], [0, 0], [0, 0]],
                          [[0.5, 0.5], [1.0, 0.0], [0.2, 0.8]])
    t2 = tree_from_matlab([0], [0], [[0, 0]], [[0.6, 0.4]])
    f = Forest([t1, t2], ("1", "2"))
    lab, sc = f.predict(np.array([[0.4], [0.5], [0.9]]))
    # 0.5 is not < 0.5, so it goes right, as in MATLAB
    assert list(lab) == [1, 2, 2]
    assert np.allclose(sc, [[0.8, 0.2], [0.4, 0.6], [0.4, 0.6]])
    # ties go to the first class (MATLAB max)
    tie = Forest([tree_from_matlab([0], [0], [[0, 0]], [[0.5, 0.5]])], ("1", "2"))
    assert tie.predict(np.zeros((1, 1)))[0][0] == 1
