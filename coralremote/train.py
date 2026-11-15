"""Py4: scikit-learn port of the Stage 2 training (gate + stress split).

TreeBagger settings read from the saved paper forests (``ModelParams``):
Gini ('gdi'), MinParent 10, MinLeaf 5, all 21 predictors at every split,
no pruning, no surrogates, bootstrap of size n (FBoot 1), empirical class
priors with equal observation weights, symmetric cost. The scikit-learn
equivalent is ``RandomForestClassifier(n_estimators=250, criterion='gini',
min_samples_split=10, min_samples_leaf=5, max_features=None, bootstrap=True)``,
except for split tie-breaking; the default engine is therefore the
MATLAB-rule CART in :mod:`coralremote.cart`.

TreeBagger and scikit-learn draw bootstrap samples and break split ties
differently, so trained forests agree statistically, not bit for bit
(patch plan §7.6). The frozen evaluator (Py3) is the exact check.

Training rows
-------------
MATLAB downsamples with ``rng(seed); randperm(n, k)`` (gate: 3000 healthy,
seed 11; 3000 stressed, seed 17; stress model keeps all 1290 PLE + 2028 BLD).
MATLAB's ``randperm`` stream is not reproducible in numpy, so ``compat``
mode recovers the exact rows (and their order) from the training matrix
saved inside the paper forests (``PrivateX``). Other modes draw with numpy.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier

from . import paths
from .features import build_hybrid_features, compute_feature_base

# Columns that are bit-identical between Python and MATLAB features
# (normalised spectrum, brightness, axis t, residual); log/acos may differ by 1 ulp.
EXACT_COLS = list(range(8)) + [16, 17, 18]

TREEBAGGER_PARAMS = dict(n_estimators=250, criterion="gini", min_samples_split=10,
                         min_samples_leaf=5, max_features=None, bootstrap=True,
                         max_samples=None, class_weight=None)


class SkForest:
    """scikit-learn forest with the ``predict -> (labels, scores)`` interface
    of :class:`coralremote.forest.Forest`. Labels are 1 / 2 like MATLAB."""

    def __init__(self, model):
        self.model = model
        self.class_names = tuple(str(c) for c in model.classes_)

    def predict(self, X):
        score = self.model.predict_proba(np.asarray(X, np.float32))
        labels = self.model.classes_[np.argmax(score, axis=1)]
        return labels.astype(int), score

    @property
    def oob_error(self):
        return 1.0 - self.model.oob_score_


def fit_forest(F, y, seed, n_jobs=-1, **overrides):
    params = {**TREEBAGGER_PARAMS, **overrides}
    rf = RandomForestClassifier(random_state=seed, n_jobs=n_jobs, oob_score=True, **params)
    rf.fit(np.asarray(F, np.float32), np.asarray(y).astype(int))
    return SkForest(rf)


def recover_rows(pool_F, saved_F):
    """Indices into ``pool_F`` of each row of ``saved_F`` (exact match on
    :data:`EXACT_COLS`). Raises if any row is missing or ambiguous."""
    keys = {}
    for i, r in enumerate(np.ascontiguousarray(pool_F[:, EXACT_COLS])):
        keys.setdefault(r.tobytes(), []).append(i)
    idx = []
    for r in np.ascontiguousarray(np.asarray(saved_F, np.float32)[:, EXACT_COLS]):
        hits = keys.get(r.tobytes())
        if not hits:
            raise ValueError("saved training row not found in Python pool")
        idx.append(hits[0])      # duplicates are identical pixels, so any one will do
    return np.asarray(idx)


def build_training_sets(inputs, mode_rows="compat", seed=0, n_gate=3000, n_stress=2500):
    """Gate and stress training matrices (21 features) and labels.

    ``mode_rows='compat'``: exact paper rows recovered from the frozen forests.
    ``mode_rows='numpy'``: same sizes, drawn with ``numpy`` (seeded).
    """
    mu_h, mu_b, v, vn = inputs["mu_h"], inputs["mu_b"], inputs["v"], inputs["vnorm2"]
    mode = inputs["view"].mode
    Xtr = inputs["Xtrain"]

    def F(X):
        return build_hybrid_features(compute_feature_base(X, mode), mu_h, mu_b, v, vn)

    FH = F(Xtr["HPY"])
    FS = F(np.concatenate([Xtr["PLE"], Xtr["BLD"]]))      # stressed = PLE then BLD
    FP, FB = F(Xtr["PLE"]), F(Xtr["BLD"])

    if mode_rows == "compat":
        from .forest import load_treebagger
        _, gp = load_treebagger(paths.STAGE2_OUTPUT, "rfGate", return_object=True)
        _, sp = load_treebagger(paths.STAGE2_OUTPUT, "rfStress", return_object=True)
        PXg = np.asarray(gp["PrivateX"], np.float32)
        nh = int(min(n_gate, len(FH)))
        ih = recover_rows(FH, PXg[:nh])
        is_ = recover_rows(FS, PXg[nh:])
        Fg = np.concatenate([FH[ih], FS[is_]])
        yg = np.r_[np.ones(len(ih), int), np.full(len(is_), 2)]
        PXs = np.asarray(sp["PrivateX"], np.float32)
        Fst = np.concatenate([FP, FB])
        ist = recover_rows(Fst, PXs)
        ys_all = np.r_[np.ones(len(FP), int), np.full(len(FB), 2)]
        Fs, ys = Fst[ist], ys_all[ist]
        rows = dict(gate_healthy=ih, gate_stressed=is_, stress=ist)
    else:
        rng = np.random.default_rng(seed)

        def down(Fm, k):
            return Fm if len(Fm) <= k else Fm[rng.choice(len(Fm), k, replace=False)]

        Fh_, Fs_ = down(FH, n_gate), down(FS, n_gate)
        Fg = np.concatenate([Fh_, Fs_])
        yg = np.r_[np.ones(len(Fh_), int), np.full(len(Fs_), 2)]
        Fp_, Fb_ = down(FP, n_stress), down(FB, n_stress)
        Fs = np.concatenate([Fp_, Fb_])
        ys = np.r_[np.ones(len(Fp_), int), np.full(len(Fb_), 2)]
        rows = None
    return dict(Fgate=Fg, ygate=yg, Fstress=Fs, ystress=ys, rows=rows)


ENGINES = ("matlab_cart", "sklearn")


def train_stage2(inputs, seed, mode_rows="compat", engine="matlab_cart", n_jobs=-1, **overrides):
    """Train gate and stress forests (250 trees, MinLeaf 5, all predictors).

    engine='matlab_cart' (default): :mod:`coralremote.cart`, MATLAB split
        rules (ties to the first predictor). Tree sizes match the paper
        forests (gate 73.3 vs 73.2-73.6 nodes, stress 127.4 vs 126.6-128.3).
    engine='sklearn': scikit-learn RandomForestClassifier. Matches on the
        test set but breaks split ties at random, which shifts the
        out-of-distribution date 20250518 by ~30 pp (HPY -> PLE); see
        docs/PY4_TRAINING_PARITY.md.
    """
    ts = build_training_sets(inputs, mode_rows, seed)
    if engine == "sklearn":
        gate = fit_forest(ts["Fgate"], ts["ygate"], seed, n_jobs, **overrides)
        stress = fit_forest(ts["Fstress"], ts["ystress"], seed + 1, n_jobs, **overrides)
    elif engine == "matlab_cart":
        from .cart import grow_forest
        kw = dict(n_trees=overrides.get("n_estimators", 250),
                  min_leaf=overrides.get("min_samples_leaf", 5),
                  min_parent=overrides.get("min_samples_split", 10))
        gate = grow_forest(ts["Fgate"], ts["ygate"], seed=seed, **kw)
        stress = grow_forest(ts["Fstress"], ts["ystress"], seed=seed + 1000, **kw)
    else:
        raise ValueError(f"engine must be one of {ENGINES}")
    return gate, stress, ts
