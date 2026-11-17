"""Numpy evaluator for frozen MATLAB ``TreeBagger`` classification forests.

A forest is a list of :class:`Tree` objects holding the arrays MATLAB's
``CompactClassificationTree`` stores (1-based in MATLAB, 0-based here):

* ``cut_var``    (n_nodes,)  predictor index of the split, -1 for a leaf
* ``cut_point``  (n_nodes,)  split value; ``x < cut_point`` goes left
* ``children``   (n_nodes, 2) left/right child, -1 for a leaf
* ``class_prob`` (n_nodes, n_classes) class posterior at each node

``predict`` returns MATLAB's ``[label, score]``: the score is the mean of
the leaf class probabilities over trees (equal tree weights), and the label
is the class with the largest score (first one on ties, like ``max``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Tree:
    cut_var: np.ndarray
    cut_point: np.ndarray
    children: np.ndarray
    class_prob: np.ndarray

    def leaf_index(self, X):
        X = np.asarray(X)
        node = np.zeros(X.shape[0], np.int64)
        rows = np.arange(X.shape[0])
        active = self.cut_var[node] >= 0
        while active.any():
            idx = rows[active]
            nd = node[idx]
            xv = X[idx, self.cut_var[nd]].astype(np.float64)
            cp = self.cut_point[nd]
            # MATLAB: x < cut -> left; NaN predictors stop at the node
            go_left = xv < cp
            nan = np.isnan(xv)
            nxt = np.where(go_left, self.children[nd, 0], self.children[nd, 1])
            nxt = np.where(nan, nd, nxt)
            node[idx] = nxt
            active = np.zeros_like(active)
            active[idx] = (~nan) & (self.cut_var[nxt] >= 0)
        return node

    def predict_proba(self, X):
        return self.class_prob[self.leaf_index(X)]


@dataclass
class Forest:
    trees: list
    class_names: tuple          # e.g. ('1', '2'); labels returned as ints when numeric
    cost: np.ndarray = None     # (K, K) misclassification cost, rows = true class
    prior: np.ndarray = None    # (K,) class prior given to TreeBagger (empirical)

    def predict(self, X):
        """MATLAB ``predict``: mean leaf probabilities over trees.

        With a non-0/1 ``Cost``, MATLAB grows the trees under cost-adjusted
        priors (prior_j * sum_k Cost(j,k), normalised), so the stored leaf
        probabilities already include the cost. The label is the expected-cost
        minimiser after mapping the scores back to the original prior; for the
        paper's cost matrices this equals the argmax of the stored scores.
        Verified on Stage 1: OOB 0.1004 and the secondary-test confusion
        matrix are reproduced (applying the cost a second time gives 0.1041).
        """
        X = np.asarray(X)
        score = np.zeros((X.shape[0], len(self.class_names)), np.float64)
        for t in self.trees:
            score += t.predict_proba(X)
        score /= len(self.trees)
        k = self._decide(score)
        labels = np.array([int(c) for c in self.class_names])[k]
        return labels, score

    def _decide(self, score):
        if self.cost is None or self.prior is None:
            return np.argmax(score, axis=1)
        C = np.asarray(self.cost, np.float64)
        P = np.asarray(self.prior, np.float64).ravel()
        padj = P * C.sum(axis=1)
        padj = padj / padj.sum()
        s = score * (P / padj)
        return np.argmin(s @ C, axis=1)

    def save_npz(self, path):
        """Compact, MATLAB-free copy of the forest (for tests / Py4 comparison)."""
        d = {"class_names": np.array(self.class_names), "n_trees": len(self.trees)}
        for i, t in enumerate(self.trees):
            d[f"t{i}_cut_var"] = t.cut_var
            d[f"t{i}_cut_point"] = t.cut_point
            d[f"t{i}_children"] = t.children
            d[f"t{i}_class_prob"] = t.class_prob
        from . import paths
        np.savez_compressed(paths.assert_writable(path), **d)

    @classmethod
    def load_npz(cls, path):
        z = np.load(path, allow_pickle=False)
        trees = [Tree(z[f"t{i}_cut_var"], z[f"t{i}_cut_point"], z[f"t{i}_children"],
                      z[f"t{i}_class_prob"]) for i in range(int(z["n_trees"]))]
        return cls(trees, tuple(str(c) for c in z["class_names"]))


def _str(o):
    a = np.asarray(o).ravel()
    return str(a[0]) if a.size else ""


def _classlabel(o):
    """Class names from a ``classreg.learning.internal.ClassLabel`` object."""
    L = o.properties["L"] if hasattr(o, "properties") else o
    return tuple(_str(c) for c in np.asarray(list(L), dtype=object).ravel())


def load_treebagger(path, var, return_object=False, mat=None):
    """Read a saved MATLAB ``TreeBagger`` (classification) with ``mat-io``.

    Uses ``Compact.Trees{k}.Impl`` (TreeImpl: CutVar, CutPoint, Children,
    ClassProb) and ``Compact.ClassNames``. Surrogate splits and categorical
    cuts are not supported (the paper forests use neither; checked here).
    """
    from matio import load_from_mat

    obj = (mat or load_from_mat(str(path), variable_names=[var]))[var]
    props = obj.properties
    compact = props["Compact"].properties
    class_names = tuple(_str(c) for c in np.asarray(compact["ClassNames"]).ravel())
    trees = []
    for t in np.asarray(compact["Trees"]).ravel():
        impl = t.properties["Impl"].properties
        if np.asarray(impl["SurrCutVar"]).size:
            raise NotImplementedError("surrogate splits")
        cs = t.properties["ClassSummary"]
        tree_classes = _classlabel(cs["ClassNames"].ravel()[0])
        nonzero = _classlabel(cs["NonzeroProbClasses"].ravel()[0])
        if tree_classes != class_names:
            raise ValueError(f"tree class order {tree_classes} != forest {class_names}")
        prob = np.asarray(impl["ClassProb"], np.float64)
        if nonzero != class_names:      # ClassProb columns follow NonzeroProbClasses
            full = np.zeros((prob.shape[0], len(class_names)))
            for j, c in enumerate(nonzero):
                full[:, class_names.index(c)] = prob[:, j]
            prob = full
        trees.append(tree_from_matlab(impl["CutVar"], impl["CutPoint"], impl["Children"], prob))
    cost = np.asarray(props["ClassSummary"]["Cost"].ravel()[0], np.float64)
    prior = np.asarray(props["ClassSummary"]["Prior"].ravel()[0], np.float64).ravel()
    forest = Forest(trees, class_names, cost if cost.size else None, prior if prior.size else None)
    return (forest, props) if return_object else forest


def tree_from_matlab(cut_predictor_idx, cut_point, children, class_prob):
    """Build a :class:`Tree` from MATLAB 1-based arrays.

    ``cut_predictor_idx``: 1-based predictor index per node, 0 for leaves.
    ``children``: (n, 2) 1-based child node ids, 0 for leaves.
    """
    cv = np.asarray(cut_predictor_idx, np.int64).ravel() - 1
    ch = np.asarray(children, np.int64).reshape(-1, 2) - 1
    leaf = ch[:, 0] < 0
    cv[leaf] = -1
    return Tree(cv, np.asarray(cut_point, np.float64).ravel(), ch,
                np.asarray(class_prob, np.float64).reshape(len(cv), -1))
