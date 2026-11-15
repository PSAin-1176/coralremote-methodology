"""Stage 1 substrate classifier (CORAL / SND / RBL), port of
``Code/Train_Stage1_LogRatio.m`` (paper run ``run_20260401_170150``).

Inputs are the A3 Phase 5 features rebuilt by :mod:`coralremote.preprocess`.
Stage 1 is independent of the Stage 2 double-scaling bug: it read the
Phase 5 stacks built from true reflectance (1 Apr), not the 5 Apr stack.
"""
from __future__ import annotations

import numpy as np

from . import io, paths
from .matlab_compat import find_colmajor
from .preprocess import phase5_features
from .rasterize import rasterize_labels_stage1

CLASS_NAMES = ("CORAL", "SND", "RBL")
POLY_CLASS_MAP = {"HC": 1, "SND": 2, "RBL": 3}
CORAL_DATES = ("20230404", "20230509", "20240307", "20241111")
SND_RBL_DATES = ("20230404", "20230509", "20230519", "20230629", "20240307", "20240317")
CORAL_TEST_DATE = "20240408"
SNDRBL_TEST_DATE = "20230519"
TRANSFER_TEST_DATES = ("20230404", "20230509", "20230519", "20230629",
                       "20240307", "20240317", "20240408", "20241111")
RF_COST = np.array([[0, 2, 2], [1, 0, 1], [1, 1, 0]], float)


def _sample(A, mask, label=None):
    v = mask & np.isfinite(A).all(axis=2)
    if label is not None:
        v &= label > 0
    r, c = find_colmajor(v)
    return A[r, c], (label[r, c] if label is not None else None)


def prepare_inputs(prefer_tif=True):
    ref = io.load_glint_stack().ref
    tr_mask, tr_lab = rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS["train"])), ref, POLY_CLASS_MAP)
    te_mask, te_lab = rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS["test"])), ref, POLY_CLASS_MAP)

    def F(d):
        return phase5_features(d, prefer_tif)

    Xs, Ys = [], []
    for d in CORAL_DATES:
        X, _ = _sample(F(d), tr_mask & (tr_lab == 1))
        Xs.append(X); Ys.append(np.full(len(X), 1, np.uint8))
    for d in SND_RBL_DATES:
        for k in (2, 3):
            X, _ = _sample(F(d), tr_mask & (tr_lab == k))
            Xs.append(X); Ys.append(np.full(len(X), k, np.uint8))
    Xtrain, Ytrain = np.concatenate(Xs), np.concatenate(Ys)

    Xc, Yc = _sample(F(CORAL_TEST_DATE), te_mask & (te_lab == 1), te_lab)
    Xs_, Ys_ = _sample(F(SNDRBL_TEST_DATE), te_mask & ((te_lab == 2) | (te_lab == 3)), te_lab)
    transfer = {d: _sample(F(d), te_mask, te_lab) for d in TRANSFER_TEST_DATES}
    return dict(ref=ref, train_mask=tr_mask, train_labels=tr_lab, test_mask=te_mask, test_labels=te_lab,
                Xtrain=Xtrain, Ytrain=Ytrain,
                Xtest=np.concatenate([Xc, Xs_]), Ytest=np.concatenate([Yc, Ys_]), transfer=transfer)


def confusion_metrics(ytrue, ypred):
    """Port of the Stage 1 ``computeConfusionMetrics`` (classes present only)."""
    ytrue = np.asarray(ytrue).ravel(); ypred = np.asarray(ypred).ravel()
    ids = np.unique(np.r_[ytrue, ypred]); ids = ids[ids > 0]
    C = np.array([[np.count_nonzero((ytrue == r) & (ypred == c)) for c in ids] for r in ids], float)
    tot = C.sum()
    oa = np.trace(C) / tot
    pe = np.sum((C.sum(0) / tot) * (C.sum(1) / tot))
    kappa = (oa - pe) / max(1 - pe, 1e-9)
    tp = np.diag(C); total = C.sum(1)
    prec = tp / np.maximum(C.sum(0), 1); rec = tp / np.maximum(total, 1)
    f1 = 2 * prec * rec / np.maximum(prec + rec, 1e-9)
    return dict(ids=ids, C=C, OA=oa, kappa=kappa, prec=prec, rec=rec, f1=f1, total=total, tp=tp)


def oob_error(forest, X, y, oob_indices):
    """Ensemble OOB misclassification rate from per-tree OOB masks
    (MATLAB ``oobError(B,'Mode','ensemble')`` with equal weights)."""
    oob = np.asarray(oob_indices, bool)
    score = np.zeros((len(X), len(forest.class_names)))
    cnt = np.zeros(len(X))
    for t, tree in enumerate(forest.trees):
        m = oob[:, t]
        if m.any():
            score[m] += tree.predict_proba(X[m])
            cnt[m] += 1
    has = cnt > 0
    score[has] /= cnt[has, None]
    labels = np.array([int(c) for c in forest.class_names])
    pred = labels[forest._decide(score)]
    return float(np.mean(pred[has] != np.asarray(y).ravel()[has]))


def run_frozen(forest, inputs):
    lab, _ = forest.predict(inputs["Xtest"])
    sec = confusion_metrics(inputs["Ytest"], lab)
    rows = []
    for d, (X, Y) in inputs["transfer"].items():
        p, _ = forest.predict(X)
        m = confusion_metrics(Y, p)
        ci = np.nonzero(m["ids"] == 1)[0]
        cp, cr, cf = ((m["prec"][ci[0]], m["rec"][ci[0]], m["f1"][ci[0]]) if ci.size else (np.nan,) * 3)
        rows.append(dict(date=d, n=len(Y), OA=m["OA"], kappa=m["kappa"], coral_prec=cp, coral_rec=cr, coral_f1=cf))
    return dict(secondary=sec, transfer=rows)
