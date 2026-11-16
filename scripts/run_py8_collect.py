"""Py8 step 1: train the Stage 2 forests with 1,000 trees per seed (both
units, seeds 0-4) and store, for every in-mask pixel of every date, the gate
and stress class-1 scores and their per-tree spread after the first 250 trees
and after all 1,000.

The bootstrap stream is drawn tree by tree, so the first 250 trees of a
1,000-tree forest are exactly the 250-tree forest of the same seed (checked
against the cached Py6 maps). Output goes to the cache folder only:
``py8/{mode}_seed{seed}.npz`` and ``py8/features_{mode}.npz``.
"""
import argparse
import time

import numpy as np

from coralremote import paths, stage2
from coralremote.cart import grow_forest
from coralremote.features import build_hybrid_features, compute_feature_base
from coralremote.matlab_compat import find_colmajor
from coralremote.train import build_training_sets

ap = argparse.ArgumentParser()
ap.add_argument("--modes", nargs="+", default=["compat", "fixed"])
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--trees", type=int, default=1000)
ap.add_argument("--split", type=int, default=250)
ap.add_argument("--row-seed", type=int, default=None,
                help="hold the training-row draw fixed at this seed (only the bootstrap varies)")
ap.add_argument("--tag", default="")
a = ap.parse_args()

out_dir = paths.cache_dir() / "py8"
out_dir.mkdir(exist_ok=True)
t0 = time.time()


def tree_stats(forest, F, split):
    """Running class-1 mean and std of the leaf probabilities over trees."""
    s = np.zeros(len(F))
    ss = np.zeros(len(F))
    res = {}
    for k, t in enumerate(forest.trees, 1):
        p = t.predict_proba(F)[:, 0]
        s += p
        ss += p * p
        if k in (split, len(forest.trees)):
            m = s / k
            res[k] = (m.astype(np.float32), np.sqrt(np.maximum(ss / k - m * m, 0)).astype(np.float32))
    return res


for mode in a.modes:
    inp = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=mode))
    view, support = inp["view"], inp["support"]
    feats, coords = {}, {}
    for d in view.dates:
        A = view.get(d)
        valid = support & np.isfinite(A).all(axis=2)
        r, c = find_colmajor(valid)
        feats[d] = build_hybrid_features(compute_feature_base(A[r, c, :], view.mode),
                                         inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])
        coords[d] = np.stack([r, c])
    if not a.tag:
        np.savez_compressed(out_dir / f"features_{mode}.npz",
                            **{f"F_{d}": feats[d] for d in feats}, **{f"rc_{d}": coords[d] for d in coords})
    print(f"{mode}: features ready [{time.time() - t0:.0f}s]", flush=True)
    for seed in a.seeds:
        ts = build_training_sets(inp, "numpy", seed if a.row_seed is None else a.row_seed)
        gate = grow_forest(ts["Fgate"], ts["ygate"], n_trees=a.trees, seed=seed)
        stress = grow_forest(ts["Fstress"], ts["ystress"], n_trees=a.trees, seed=seed + 1000)
        rec = {"train_min": np.nanmin(np.r_[ts["Fgate"], ts["Fstress"]], 0),
               "train_max": np.nanmax(np.r_[ts["Fgate"], ts["Fstress"]], 0),
               "gate_min": np.nanmin(ts["Fgate"], 0), "gate_max": np.nanmax(ts["Fgate"], 0),
               "stress_min": np.nanmin(ts["Fstress"], 0), "stress_max": np.nanmax(ts["Fstress"], 0)}
        for d, F in feats.items():
            for name, f in (("g", gate), ("s", stress)):
                st = tree_stats(f, F, a.split)
                for k, (m, sd) in st.items():
                    rec[f"{name}{k}_{d}"] = m
                    rec[f"{name}{k}sd_{d}"] = sd
        np.savez_compressed(out_dir / f"{mode}{a.tag}_seed{seed}.npz", **rec)
        print(f"{mode} seed {seed} done [{time.time() - t0:.0f}s]", flush=True)
