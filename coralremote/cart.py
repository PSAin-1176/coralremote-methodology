"""A small CART that follows MATLAB ``fitctree`` / ``TreeBagger`` rules as
used by the paper (classification, 'gdi', 'allsplits', all predictors,
no surrogates, no pruning):

* every predictor, every cut between consecutive distinct sorted values;
* both children need >= MinLeaf observations, the node >= MinParent;
* the split with the largest (weighted) Gini decrease wins; ties go to the
  lowest predictor index, then the lowest cut (MATLAB ``max`` returns the first);
* cut point = midpoint of the two neighbouring values in single precision
  (the saved ``CutPoint`` arrays are single); ``x < cut`` goes left;
* node class probabilities are weighted class shares. Weights are 1 for a
  symmetric cost (Stage 2). With a non-0/1 cost (Stage 1) MATLAB folds the
  cost into the prior, so class k is weighted by the row sum of Cost
  (saved Stage 1 trees: root ClassProb ~ [0.36 0.28 0.36] vs empirical
  [0.22 0.34 0.44]).

TreeBagger: each tree is grown on a bootstrap sample of size n drawn with
replacement (numpy RNG here; MATLAB's stream can't be reproduced).
"""
from __future__ import annotations

import numpy as np

from .forest import Forest, Tree


def _best_split(Xn, Yoh, wn, min_leaf):
    """Best (gain, feature, cut) for one node; None if no valid split.

    ``Yoh``: (m, K) one-hot classes; ``wn``: (m,) observation weights.
    Gain = P(t)g(t) - P(L)g(L) - P(R)g(R) with weighted Gini (MATLAB 'gdi');
    MinLeaf counts observations, not weight.
    """
    m, p = Xn.shape
    order = np.argsort(Xn, axis=0, kind="stable")
    xs = np.take_along_axis(Xn, order, axis=0)
    Wm = Yoh * wn[:, None]                          # (m, K) weighted class mass
    cw = np.cumsum(Wm[order], axis=0)[:-1]          # (m-1, p, K) mass left of each cut
    tot = Wm.sum(axis=0)                            # (K,)
    T = tot.sum()
    wl = cw.sum(axis=2)
    wr = T - wl
    cr = tot - cw
    gl = 1 - np.sum((cw / wl[..., None]) ** 2, axis=2)
    gr = 1 - np.sum((cr / np.maximum(wr, 1e-300)[..., None]) ** 2, axis=2)
    g0 = 1 - np.sum((tot / T) ** 2)
    gain = g0 * T - wl * gl - wr * gr
    nl = np.arange(1, m)[:, None]
    valid = (xs[:-1] < xs[1:]) & (nl >= min_leaf) & (m - nl >= min_leaf)
    gain = np.where(valid, gain, -np.inf)
    pos = np.argmax(gain, axis=0)                   # first best position per feature
    best_per_f = gain[pos, np.arange(p)]
    j = int(np.argmax(best_per_f))                  # first best feature
    g = best_per_f[j]
    if not np.isfinite(g) or g <= 1e-12 * T:
        return None
    i = pos[j]
    a, b = xs[i, j], xs[i + 1, j]
    cut = np.float32((np.float32(a) + np.float32(b)) / np.float32(2))
    return g, j, cut


def grow_tree(X, y, min_leaf=5, min_parent=10, classes=(1, 2), class_weight=None):
    """Grow one tree on (X, y); ``class_weight[k]`` weights every observation
    of ``classes[k]``. Returns a :class:`Tree`."""
    X = np.asarray(X, np.float32)
    y = np.asarray(y)
    Yoh = (y[:, None] == np.asarray(classes)[None, :]).astype(np.float64)
    cwt = np.ones(len(classes)) if class_weight is None else np.asarray(class_weight, np.float64)
    w = Yoh @ cwt
    cut_var, cut_pt, children, prob = [], [], [], []

    def new_node(idx):
        k = len(cut_var)
        cut_var.append(-1)
        cut_pt.append(0.0)
        children.append([-1, -1])
        mass = (Yoh[idx] * w[idx, None]).sum(axis=0)
        prob.append(mass / mass.sum())
        return k

    root = np.arange(len(y))
    queue = [(new_node(root), root)]
    while queue:
        k, idx = queue.pop(0)                       # breadth-first like MATLAB node numbering
        yn = y[idx]
        if len(idx) < min_parent or (yn == yn[0]).all():
            continue
        s = _best_split(X[idx], Yoh[idx], w[idx], min_leaf)
        if s is None:
            continue
        _, j, cut = s
        left = idx[X[idx, j] < cut]
        right = idx[X[idx, j] >= cut]
        cut_var[k], cut_pt[k] = j, float(cut)
        kl = new_node(left)
        kr = new_node(right)
        children[k] = [kl, kr]
        queue.append((kl, left))
        queue.append((kr, right))
    return Tree(np.asarray(cut_var, np.int64), np.asarray(cut_pt, np.float64),
                np.asarray(children, np.int64), np.asarray(prob, np.float64))


def grow_forest(X, y, n_trees=250, min_leaf=5, min_parent=10, seed=0, classes=(1, 2),
                cost=None, return_oob=False):
    """TreeBagger-style forest (bootstrap of size n per tree).

    With a ``cost`` matrix, class k is weighted by ``cost[k].sum()`` during
    growth and prediction uses the cost-prior rule in :class:`Forest`.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    n = len(y)
    cw = None if cost is None else np.asarray(cost, np.float64).sum(axis=1)
    trees = []
    oob = np.zeros((n, n_trees), bool)
    for t in range(n_trees):
        b = rng.integers(0, n, n)
        oob[:, t] = np.bincount(b, minlength=n) == 0
        trees.append(grow_tree(X[b], y[b], min_leaf, min_parent, classes, cw))
    prior = np.array([(y == c).mean() for c in classes])
    f = Forest(trees, tuple(str(c) for c in classes),
               None if cost is None else np.asarray(cost, np.float64),
               prior if cost is not None else None)
    return (f, oob) if return_oob else f
