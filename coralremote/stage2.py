"""Stage 2 hierarchical classifier (gate HPY vs stressed, then PLE vs BLD).

Port of ``Code/Run_Stage2_Hierarchical_CuratedHoldout_QuickTest.m``.
This module holds the deterministic parts: sample extraction, PLE curation,
the spatial-block split, hierarchical prediction with the 0.55 UNK rule and
whole-date classification. Training (Py4) is not here yet; the paper's own
forests are evaluated through :mod:`coralremote.forest`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import io
from .features import (axis_vector, build_hybrid_features, compute_feature_base,
                       prepare_spectra, project_to_axis, rowsum)
from .matlab_compat import find_colmajor, mquantile, mround
from .rasterize import rasterize_date_labels
from .units import COMPAT, FIXED, get_mode

HPY, PLE, BLD, UNK = 1, 2, 3, 4
CLASS_NAMES = ("HPY", "PLE", "BLD")


@dataclass
class Stage2Options:
    gate_threshold: float = 0.55
    stress_threshold: float = 0.55
    healthy_quantile: float = 0.80
    bleached_quantile: float = 0.10
    anchor_residual_quantile: float = 0.95
    ple_polygon_residual_quantile: float = 0.95
    ple_band_expand_frac: float = 0.15
    ple_train_blocks: tuple = (2,)
    ple_test_blocks: tuple = (7,)
    healthy_anchor_date: str = "20230509"
    bleached_anchor_date: str = "20240613"
    interpretation_date: str = "20240514"
    mode: object = field(default=FIXED)


# --------------------------------------------------------------------------
# label filtering
def _train_eligible(p):
    te = False
    if "train_eligible" in p:
        te = bool(p["train_eligible"])
    if "TrainEligible" in p:
        te = bool(p["TrainEligible"])
    return te


def filter_features(features, keep_classes, require_train_eligible):
    out = []
    for f in features:
        p = f["properties"]
        if "Class" not in p or str(p["Class"]) not in keep_classes:
            continue
        if require_train_eligible and not _train_eligible(p):
            continue
        out.append(f)
    return out


def filter_date_features(features, date):
    return [f for f in features if "date" in f["properties"] and str(f["properties"]["date"]) == date]


def _spatial_block(p):
    if "spatial_block" in p:
        return float(p["spatial_block"])
    if "SpatialBlock" in p:
        return float(p["SpatialBlock"])
    return np.nan


# --------------------------------------------------------------------------
# sampling
class StackView:
    """Per-date (H, W, 8) spectra on the ``mode`` scale, cached in memory."""

    def __init__(self, stack, mode=FIXED):
        self.stack = stack
        self.mode = get_mode(mode)
        self.dates = list(stack.dates)
        self.ref = stack.ref
        self._cache = {}

    def get(self, date):
        if date not in self._cache:
            self._cache[date] = prepare_spectra(self.stack.get(date), self.mode)
        return self._cache[date]


def build_raw_spectral_split(view, support, features, class_to_id, return_dates=False):
    """Port of ``buildRawSpectralSplit``: dates in order of first appearance,
    pixels in column-major order within each date. With ``return_dates`` a
    third array gives each row's acquisition date."""
    dates = []
    for f in features:
        d = str(f["properties"]["date"])
        if d and d not in dates:
            dates.append(d)
    Xs, Ys, Ds = [], [], []
    for d in dates:
        mask, label = rasterize_date_labels(features, view.ref, class_to_id, d)
        mask &= support
        A = view.get(d)
        valid = mask & np.isfinite(A).all(axis=2) & (label > 0)
        r, c = find_colmajor(valid)
        Xs.append(A[r, c, :])
        Ys.append(label[r, c])
        Ds.append(np.full(r.size, d))
    if not Xs:
        empty = (np.zeros((0, 8), np.float32), np.zeros(0, np.uint8))
        return empty + (np.zeros(0, "<U8"),) if return_dates else empty
    out = (np.concatenate(Xs), np.concatenate(Ys))
    return out + (np.concatenate(Ds),) if return_dates else out


def mean_omitnan_f32(X):
    """Column mean of a single matrix (MATLAB ``mean(X,1,'omitnan')``).

    MATLAB sums a single column with four interleaved float32 accumulators
    (row i goes to lane i mod 4), then adds the lanes in order. Reproducing
    that gives the saved ``muHealthy`` / ``muBleached`` bit for bit; plain
    sequential or pairwise sums differ by up to 4e-7.
    """
    X = np.asarray(X, np.float32)
    ok = ~np.isnan(X)
    s = np.where(ok, X, 0).astype(np.float32)
    n, p = s.shape
    lanes = np.zeros((4, p), np.float32)
    full = n - n % 4
    for i in range(0, full, 4):
        lanes += s[i:i + 4]
    for k, i in enumerate(range(full, n)):
        lanes[k] += s[i]
    tot = lanes[0].copy()
    for k in range(1, 4):
        tot += lanes[k]
    return (tot / ok.sum(axis=0).astype(np.float32)).astype(np.float32)


# --------------------------------------------------------------------------
# PLE curation and spatial split
def score_one_feature(f, view, support, mu_h, v, vnorm2):
    p = f["properties"]
    H, W = support.shape
    if not all(k in p for k in ("date", "Class", "centroid_row", "centroid_col")):
        return None
    d = str(p["date"])
    if d not in view.dates:
        return None
    row = int(max(1, min(H, mround(float(p["centroid_row"])))))
    col = int(max(1, min(W, mround(float(p["centroid_col"])))))
    if not support[row - 1, col - 1]:
        return None
    pix = view.get(d)[row - 1, col - 1, :][None, :]
    if not np.isfinite(pix).all():
        return None
    F = compute_feature_base(pix, view.mode)
    t, r = project_to_axis(F[:, :8], mu_h, v, vnorm2)
    return float(t[0]), float(r[0])


def select_middle_pale(features, view, support, mu_h, v, vnorm2, h_thr, b_thr, r_anchor, opts):
    n = len(features)
    meanT = np.full(n, np.nan, np.float32)
    meanR = np.full(n, np.nan, np.float32)
    npix = np.zeros(n, int)
    for i, f in enumerate(features):
        s = score_one_feature(f, view, support, mu_h, v, vnorm2)
        if s is not None:
            meanT[i], meanR[i] = s
            npix[i] = 1
    valid = (npix > 0) & np.isfinite(meanT) & np.isfinite(meanR)
    f32 = np.float32           # MATLAB: single op double -> single
    r_ple = f32(min(f32(r_anchor), f32(mquantile(meanR[valid], opts.ple_polygon_residual_quantile))))
    band = f32(max(f32(b_thr) - f32(h_thr), f32(1e-6)))
    h_keep = f32(f32(h_thr) - f32(opts.ple_band_expand_frac) * band)
    b_keep = f32(f32(b_thr) + f32(opts.ple_band_expand_frac) * band)
    keep = valid & (meanT > h_keep) & (meanT < b_keep) & (meanR <= r_ple)
    stats = dict(rPleCut=r_ple, nTotalPolygons=n, nKeptPolygons=int(keep.sum()),
                 nRejectedPolygons=int(n - keep.sum()), nKeptPixels=int(npix[keep].sum()),
                 nRejectedPixels=int(npix[~keep].sum()), hKeep=h_keep, bKeep=b_keep,
                 keepMask=keep, meanT=meanT, meanR=meanR)
    return [f for f, k in zip(features, keep) if k], stats


def split_by_block(features, train_blocks, test_blocks):
    blocks = np.array([_spatial_block(f["properties"]) for f in features])
    tr = np.isin(blocks, train_blocks)
    te = np.isin(blocks, test_blocks)
    stats = dict(trainMask=tr, testMask=te, nTrainPolygons=int(tr.sum()), nTestPolygons=int(te.sum()))
    return ([f for f, k in zip(features, tr) if k], [f for f, k in zip(features, te) if k], stats)


# --------------------------------------------------------------------------
# prediction
def hierarchical_predict(F, gate, stress, opts):
    """Gate (1 HPY / 2 stressed) then stress split (1 PLE / 2 BLD); UNK (4)
    when the gate's max score, or the stress max score for gate-stressed
    pixels, is below its threshold."""
    gate_pred, gate_score = gate.predict(F)
    stress_pred, stress_score = stress.predict(F)
    gate_conf = gate_score.max(axis=1)
    stress_conf = stress_score.max(axis=1)
    pred = np.zeros(F.shape[0], np.uint8)
    pred[gate_pred == 1] = HPY
    pred[(gate_pred == 2) & (stress_pred == 1)] = PLE
    pred[(gate_pred == 2) & (stress_pred == 2)] = BLD
    pred[gate_conf < opts.gate_threshold] = UNK
    pred[(gate_pred == 2) & (stress_conf < opts.stress_threshold)] = UNK
    return pred, gate_pred.astype(np.uint8), stress_pred.astype(np.uint8)


def classify_all_dates(view, support, mu_h, mu_b, gate, stress, opts, return_maps=False):
    v, vnorm2 = axis_vector(mu_h, mu_b)
    rows, maps = [], {}
    for d in view.dates:
        A = view.get(d)
        valid = support & np.isfinite(A).all(axis=2)
        r, c = find_colmajor(valid)
        F = build_hybrid_features(compute_feature_base(A[r, c, :], view.mode), mu_h, mu_b, v, vnorm2)
        pred, _, _ = hierarchical_predict(F, gate, stress, opts)
        n = max(pred.size, 1)
        rows.append(dict(date=d,
                         HPY=100 * np.count_nonzero(pred == HPY) / n,
                         PLE=100 * np.count_nonzero(pred == PLE) / n,
                         BLD=100 * np.count_nonzero(pred == BLD) / n,
                         UNK=100 * np.count_nonzero(pred == UNK) / n))
        if return_maps:
            m = np.zeros(support.shape, np.uint8)
            m[r, c] = pred
            maps[d] = m
    return (rows, maps) if return_maps else rows


# --------------------------------------------------------------------------
# frozen-model run (Py3)
def prepare_frozen_inputs(opts=None, stack=None, drop_train_dates=()):
    """Everything up to the trained forests that is deterministic: axis
    statistics, PLE curation, split, test pixels (in MATLAB order).

    ``drop_train_dates``: training polygons (HPY, BLD, PLE) on these dates
    are left out before anything is computed, so the axis means, thresholds
    and PLE curation are re-derived without them (R1 leave-date-out check).
    Test polygons and anchors are unchanged. ``Dtrain`` gives each training
    row's date."""
    opts = opts or Stage2Options()
    stack = stack or io.load_glint_stack()
    support, H, W = io.load_baseline_mask()
    assert (H, W) == (stack.ref.height, stack.ref.width)
    view = StackView(stack, opts.mode)

    u_train, u_test = io.load_labels("unified_train"), io.load_labels("unified_test")
    s_train, s_test = io.load_labels("strict_train"), io.load_labels("strict_test")

    train_h = filter_features(u_train, {"HPY"}, True)
    train_b = filter_features(s_train, {"BLD"}, True)
    train_p_all = filter_features(u_train, {"PLE"}, True)
    if drop_train_dates:
        drop = {str(d) for d in drop_train_dates}

        def keep(fs):
            return [f for f in fs if str(f["properties"].get("date")) not in drop]
        train_h, train_b, train_p_all = keep(train_h), keep(train_b), keep(train_p_all)

    Xh, _, Dh = build_raw_spectral_split(view, support, train_h, {"HPY": HPY}, return_dates=True)
    Xb, _, Db = build_raw_spectral_split(view, support, train_b, {"BLD": BLD}, return_dates=True)
    Fh = compute_feature_base(Xh, view.mode)
    Fb = compute_feature_base(Xb, view.mode)
    mu_h = mean_omitnan_f32(Fh[:, :8])
    mu_b = mean_omitnan_f32(Fb[:, :8])
    v, vnorm2 = axis_vector(mu_h, mu_b)
    tH, rH = project_to_axis(Fh[:, :8], mu_h, v, vnorm2)
    tB, rB = project_to_axis(Fb[:, :8], mu_h, v, vnorm2)
    h_thr = mquantile(tH, opts.healthy_quantile)
    b_thr = mquantile(tB, opts.bleached_quantile)
    r_anchor = mquantile(np.concatenate([rH, rB]), opts.anchor_residual_quantile)

    curated, ple_stats = select_middle_pale(train_p_all, view, support, mu_h, v, vnorm2,
                                            h_thr, b_thr, r_anchor, opts)
    train_p, test_p, split_stats = split_by_block(curated, opts.ple_train_blocks, opts.ple_test_blocks)
    Xp, _, Dp = build_raw_spectral_split(view, support, train_p, {"PLE": PLE}, return_dates=True)

    test_feats = (filter_features(u_test, {"HPY"}, False) + test_p
                  + filter_features(s_test, {"BLD"}, False))
    Xtest, Ytest, Dtest = build_raw_spectral_split(view, support, test_feats,
                                                   {"HPY": HPY, "PLE": PLE, "BLD": BLD}, return_dates=True)
    anchor_h = filter_date_features(filter_features(u_test, {"HPY"}, False), opts.healthy_anchor_date)
    anchor_b = filter_date_features(filter_features(s_test, {"BLD"}, False), opts.bleached_anchor_date)
    XaH, YaH = build_raw_spectral_split(view, support, anchor_h, {"HPY": HPY})
    XaB, YaB = build_raw_spectral_split(view, support, anchor_b, {"BLD": BLD})

    return dict(opts=opts, view=view, support=support,
                mu_h=mu_h, mu_b=mu_b, v=v, vnorm2=vnorm2,
                thresholds=dict(hThr=h_thr, bThr=b_thr, rAnchor=r_anchor),
                ple_stats=ple_stats, split_stats=split_stats,
                n_train=dict(HPY=len(Xh), BLD=len(Xb), PLE=len(Xp)),
                Xtrain=dict(HPY=Xh, BLD=Xb, PLE=Xp),
                Dtrain=dict(HPY=Dh, BLD=Db, PLE=Dp),
                Xtest=Xtest, Ytest=Ytest, Dtest=Dtest,
                anchors=dict(healthy=(XaH, YaH), bleached=(XaB, YaB)))


def run_frozen(gate, stress, inputs=None, mu_h=None, mu_b=None, opts=None, return_maps=False):
    """Evaluate two forests (the paper's frozen ones, or Py4-trained ones)
    on Python-built features.

    ``mu_h``/``mu_b`` default to the Python-recomputed axis means; pass the
    saved MATLAB values to isolate the evaluator from the mean computation.
    """
    from .evaluate import compute_metrics3

    inputs = inputs or prepare_frozen_inputs(opts)
    opts = inputs["opts"]
    view, support = inputs["view"], inputs["support"]
    mu_h = inputs["mu_h"] if mu_h is None else np.asarray(mu_h, np.float32)
    mu_b = inputs["mu_b"] if mu_b is None else np.asarray(mu_b, np.float32)
    v, vnorm2 = axis_vector(mu_h, mu_b)

    def feats(X):
        return build_hybrid_features(compute_feature_base(X, view.mode), mu_h, mu_b, v, vnorm2)

    pred, gate_pred, stress_pred = hierarchical_predict(feats(inputs["Xtest"]), gate, stress, opts)
    metrics = compute_metrics3(inputs["Ytest"], pred)
    XaH, YaH = inputs["anchors"]["healthy"]
    XaB, YaB = inputs["anchors"]["bleached"]
    pH, _, _ = hierarchical_predict(feats(XaH), gate, stress, opts)
    pB, _, _ = hierarchical_predict(feats(XaB), gate, stress, opts)
    maps = None
    if return_maps:
        rows, maps = classify_all_dates(view, support, mu_h, mu_b, gate, stress, opts, return_maps=True)
    else:
        rows = classify_all_dates(view, support, mu_h, mu_b, gate, stress, opts)
    return dict(maps=maps, metrics=metrics, pred=pred, gate_pred=gate_pred, stress_pred=stress_pred,
                healthy_anchor_recall=float(np.mean(pH == YaH)), n_healthy_anchor=int(YaH.size),
                bleached_anchor_recall=float(np.mean(pB == YaB)), n_bleached_anchor=int(YaB.size),
                date_summary=rows)
