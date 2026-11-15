import numpy as np

from coralremote.cart import grow_forest, grow_tree


def test_tie_goes_to_first_predictor_and_midpoint_cut():
    # columns 0 and 1 separate the classes equally well -> MATLAB picks column 0
    X = np.array([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [10, 10], [11, 11], [12, 12], [13, 13], [14, 14]],
                 np.float32)
    y = np.array([1] * 5 + [2] * 5)
    t = grow_tree(X, y)
    assert t.cut_var[0] == 0 and t.cut_point[0] == 7.0
    assert (t.cut_var[1:] == -1).all()                    # both children pure
    assert np.allclose(t.class_prob[1:], [[1, 0], [0, 1]])


def test_min_leaf_and_min_parent():
    X = np.arange(12, dtype=np.float32)[:, None]
    y = np.array([1, 2] * 6)
    t = grow_tree(X, y, min_leaf=5, min_parent=10)
    leaves = t.cut_var == -1
    # leaf sizes can't be checked from the arrays directly; re-route the data
    sizes = np.bincount(t.leaf_index(X), minlength=len(t.cut_var))[leaves]
    assert sizes.min() >= 5


def test_forest_predicts_separable_data():
    rng = np.random.default_rng(0)
    X = np.r_[rng.normal(0, 1, (100, 3)), rng.normal(4, 1, (100, 3))].astype(np.float32)
    y = np.r_[np.ones(100, int), np.full(100, 2)]
    f = grow_forest(X, y, n_trees=20, seed=1)
    lab, sc = f.predict(X)
    assert (lab == y).mean() > 0.97 and np.allclose(sc.sum(1), 1)
