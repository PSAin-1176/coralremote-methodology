"""A small CART that follows MATLAB ``fitctree`` / ``TreeBagger`` rules as
used by the paper (classification, 'gdi', 'allsplits', all predictors,
no surrogates, no pruning):

* every predictor, every cut between consecutive distinct sorted values;
* both children need >= MinLeaf observations, the node >= MinParent;
* the split with the largest Gini decrease wins; ties go to the lowest
  predictor index, then the lowest cut (MATLAB ``max`` returns the first);
* cut point = midpoint of the two neighbouring values in single precision
  (the saved ``CutPoint`` arrays are single); ``x < cut`` goes left;
* leaf/node class probabilities are the (equal-weight) class proportions,
  as in the saved forests (root ClassProb = bootstrap class shares).

TreeBagger: each tree is grown on a bootstrap sample of size n drawn with
replacement (numpy RNG here; MATLAB's stream can't be reproduced).
Used to separate implementation effects from genuine model instability
when comparing with scikit-learn (Py4).
"""
from __future__ import annotations

import numpy as np

from .forest import Forest, Tree


def _best_split(Xn, yn, min_leaf):
    """Best (gain, feature, cut) for one node; None if no valid split."""
    m, p = Xn.shape
    order = np.argsort(Xn, axis=0, kind="stable")
    xs = np.take_along_axis(Xn, order, axis=0)
    c1 = np.cumsum(yn[order] == 1, axis=0)[:-1].astype(np.float64)   # class-1 count on the left, left size i+1
    nl = np.arange(1, m, dtype=np.float64)[:, None]
    nr = m - nl
    c1_tot = float((yn == 1).sum())
    c1r = c1_tot - c1
    gl = 1 - (c1 / nl) ** 2 - ((nl - c1) / nl) ** 2
    gr = 1 - (c1r / nr) ** 2 - ((nr - c1r) / nr) ** 2
    p1 = c1_tot / m
    g0 = 1 - p1 ** 2 - (1 - p1) ** 2
    gain = g0 - (nl * gl + nr * gr) / m
    valid = (xs[:-1] < xs[1:]) & (nl >= min_leaf) & (nr >= min_leaf)
    gain = np.where(valid, gain, -np.inf)
    pos = np.argmax(gain, axis=0)                 # first best position per feature
    best_per_f = gain[pos, np.arange(p)]
    j = int(np.argmax(best_per_f))                # first best feature
    g = best_per_f[j]
    if not np.isfinite(g) or g <= 1e-12:
        return None
    i = pos[j]
    a, b = xs[i, j], xs[i + 1, j]
    cut = np.float32((np.float32(a) + np.float32(b)) / np.float32(2))
    return g, j, cut


def grow_tree(X, y, min_leaf=5, min_parent=10):
    """Grow one tree on (X, y); y in {1, 2}. Returns a :class:`Tree`."""
    X = np.asarray(X, np.float32)
    y = np.asarray(y)
    cut_var, cut_pt, children, prob = [], [], [], []

    def new_node(idx):
        k = len(cut_var)
        cut_var.append(-1)
        cut_pt.append(0.0)
        children.append([-1, -1])
        n1 = float((y[idx] == 1).sum())
        prob.append([n1 / len(idx), 1 - n1 / len(idx)])
        return k

    stack = [(new_node(np.arange(len(y))), np.arange(len(y)))]
    while stack:
        k, idx = stack.pop(0)                    # breadth-first like MATLAB node numbering
        yn = y[idx]
        if len(idx) < min_parent or (yn == yn[0]).all():
            continue
        s = _best_split(X[idx], yn, min_leaf)
        if s is None:
            continue
        _, j, cut = s
        left = idx[X[idx, j] < cut]
        right = idx[X[idx, j] >= cut]
        cut_var[k], cut_pt[k] = j, float(cut)
        kl = new_node(left)
        kr = new_node(right)
        children[k] = [kl, kr]
        stack.append((kl, left))
        stack.append((kr, right))
    return Tree(np.asarray(cut_var, np.int64), np.asarray(cut_pt, np.float64),
                np.asarray(children, np.int64), np.asarray(prob, np.float64))


def grow_forest(X, y, n_trees=250, min_leaf=5, min_parent=10, seed=0):
    rng = np.random.default_rng(seed)
    n = len(y)
    trees = []
    for _ in range(n_trees):
        b = rng.integers(0, n, n)
        trees.append(grow_tree(X[b], np.asarray(y)[b], min_leaf, min_parent))
    return Forest(trees, ("1", "2"))
