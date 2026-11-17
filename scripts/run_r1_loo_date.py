"""R1 (Reviewer D): leave-one-June-date-out test of the Stage 2 BLD fractions.

Every BLD training pixel comes from the three June 2024 dates, which are
also the wave-affected dates. This script retrains Stage 2 with one June
date's training pixels removed (all classes, both models), classifies all
mask pixels on every date, and compares with a baseline trained on the same
seeds with nothing removed.

Settings (each for every mode and seed):
  baseline        nothing removed (the Py6 retraining recipe)
  rows_D          training rows from date D removed; axis means, thresholds
                  and PLE curation stay as in the paper (the literal request)
  strict_D        date D's training polygons removed before anything is
                  computed, so the axis means, thresholds and PLE curation
                  are re-derived without D as well
  alljune         all three June dates removed. No BLD training pixels are
                  left, so only the gate (HPY vs stressed, stressed = PLE
                  from 20230509) can be trained; results are gate-level only

Recipe: MATLAB-rule CART, 250 trees, MinLeaf 5, MinParent 10, all 21
predictors, gate 3,000 HPY + up to 3,000 stressed (numpy draw), stress model
all PLE + BLD (cap 2,500 each), UNK if max score < 0.55 in either model.

Aggregates only go to the JSON (cache folder); no imagery or spectra.
"""
import argparse
import json
import time

import numpy as np

from coralremote import paths, stage2
from coralremote.cart import grow_forest
from coralremote.evaluate import compute_metrics3
from coralremote.features import build_hybrid_features, compute_feature_base
from coralremote.forest import load_treebagger
from coralremote.matlab_compat import find_colmajor
from coralremote.train import build_training_sets, train_stage2
from coralremote.units import COMPAT, FIXED

JUNE = ("20240613", "20240616", "20240620")
EXTRA = ("20240514",)
NAMES = ("HPY", "PLE", "BLD", "UNK")
THR = 0.55

ap = argparse.ArgumentParser()
ap.add_argument("--modes", nargs="+", default=["compat", "fixed"])
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--settings", nargs="+",
                default=["baseline"] + [f"rows_{d}" for d in JUNE] + [f"strict_{d}" for d in JUNE] + ["alljune"])
ap.add_argument("--json", default=str(paths.cache_dir() / "r1_loo_date.json"))
a = ap.parse_args()
t0 = time.time()


def counts_by_date(D):
    u, c = np.unique(D, return_counts=True)
    return {str(k): int(v) for k, v in zip(u, c)}


def drop_rows(inp, dates):
    """Copy of the inputs with the training rows on ``dates`` removed."""
    new = dict(inp)
    keep = {k: ~np.isin(D, list(dates)) for k, D in inp["Dtrain"].items()}
    new["Xtrain"] = {k: X[keep[k]] for k, X in inp["Xtrain"].items()}
    new["Dtrain"] = {k: D[keep[k]] for k, D in inp["Dtrain"].items()}
    new["n_train"] = {k: int(keep[k].sum()) for k in keep}
    return new


def date_features(inp):
    """Per-date (features, n) for all valid mask pixels, in this input's axis."""
    view, support = inp["view"], inp["support"]
    out = {}
    for d in view.dates:
        A = view.get(d)
        r, c = find_colmajor(support & np.isfinite(A).all(axis=2))
        out[d] = build_hybrid_features(compute_feature_base(A[r, c, :], view.mode),
                                       inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])
    return out


def feats(inp, X):
    return build_hybrid_features(compute_feature_base(X, inp["view"].mode),
                                 inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])


def gate_level(gate, F):
    """Gate decision alone: 1 HPY, 2 stressed, 4 UNK (gate max score < 0.55)."""
    lab, sc = gate.predict(F)
    return np.where(sc.max(1) < THR, 4, lab)


def frac(p, codes=(1, 2, 3, 4)):
    n = max(p.size, 1)
    return [round(100 * float(np.count_nonzero(p == k)) / n, 3) for k in codes]


def gate_rows_by_date(inp, seed, ts):
    """Dates of the stressed rows the gate drew (replays the numpy draw and
    checks it against the training matrix)."""
    rng = np.random.default_rng(seed)
    nH = len(inp["Xtrain"]["HPY"])
    ih = rng.choice(nH, 3000, replace=False) if nH > 3000 else np.arange(nH)
    Ds = np.r_[inp["Dtrain"]["PLE"], inp["Dtrain"]["BLD"]]
    Cs = np.r_[np.full(len(inp["Dtrain"]["PLE"]), "PLE"), np.full(len(inp["Dtrain"]["BLD"]), "BLD")]
    nS = len(Ds)
    is_ = rng.choice(nS, 3000, replace=False) if nS > 3000 else np.arange(nS)
    Xs = np.concatenate([inp["Xtrain"]["PLE"], inp["Xtrain"]["BLD"]])
    chk = feats(inp, Xs[is_])
    assert np.array_equal(chk, ts["Fgate"][len(ih):], equal_nan=True), "gate draw replay mismatch"
    out = {}
    for cls in ("PLE", "BLD"):
        sel = Cs[is_] == cls
        out[cls] = counts_by_date(Ds[is_][sel])
    return out


def evaluate(inp, gate, stress, F_dates):
    res = {"dates": {}, "test": {}}
    for d, F in F_dates.items():
        g = gate_level(gate, F)
        row = {"n": int(len(F)), "gate_HPY_STR_UNK": frac(g, (1, 2, 4))}
        if stress is not None:
            p, _, _ = stage2.hierarchical_predict(F, gate, stress, inp["opts"])
            row["HPY_PLE_BLD_UNK"] = frac(p)
        res["dates"][d] = row
    Ft = feats(inp, inp["Xtest"])
    Y, Dt = inp["Ytest"], inp["Dtest"]
    if stress is not None:
        p, _, _ = stage2.hierarchical_predict(Ft, gate, stress, inp["opts"])
        m = compute_metrics3(Y, p)
        res["test"]["OA"] = float(m["OA"])
        res["test"]["f1"] = [float(x) for x in m["f1"][:3]]
    else:
        p = gate_level(gate, Ft)
    # per June date: what the test BLD polygons on that date are called
    for d in JUNE + EXTRA:
        sel = (Dt == d) & (Y == 3)
        if sel.any():
            res["test"][f"BLD_test_{d}"] = {"n": int(sel.sum()),
                                            "pred": frac(p[sel]) if stress is not None else frac(p[sel], (1, 2, 4))}
    return res


out = {"meta": {"seeds": a.seeds, "settings": a.settings, "threshold": THR}, "frozen_paper": None, "modes": {}}
if "compat" in a.modes:
    base_c = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=COMPAT))
    fg = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
    fs = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
    fr = evaluate(base_c, fg, fs, date_features(base_c))
    fr["test_metrics_check"] = stage2.run_frozen(fg, fs, base_c)["metrics"]["OA"]
    out["frozen_paper"] = fr
    print(f"frozen paper forests evaluated [{time.time() - t0:.0f}s]", flush=True)

for mode in a.modes:
    m = COMPAT if mode == "compat" else FIXED
    base = base_c if mode == "compat" else stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=m))
    F_base = date_features(base)
    mo = out["modes"][mode] = {"pool_by_date": {k: counts_by_date(D) for k, D in base["Dtrain"].items()},
                               "settings": {}}
    for setting in a.settings:
        if setting == "baseline":
            inp, F_dates = base, F_base
        elif setting.startswith("rows_"):
            inp, F_dates = drop_rows(base, [setting[5:]]), F_base
        elif setting.startswith("strict_"):
            inp = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=m), drop_train_dates=[setting[7:]])
            F_dates = date_features(inp)
        elif setting == "alljune":
            inp, F_dates = drop_rows(base, JUNE), F_base
        else:
            raise ValueError(setting)
        so = mo["settings"][setting] = {
            "pool_by_date": {k: counts_by_date(D) for k, D in inp["Dtrain"].items()},
            "n_train": {k: int(len(X)) for k, X in inp["Xtrain"].items()},
            "runs": []}
        if setting.startswith("strict_"):
            so["ple_kept_polygons"] = int(inp["ple_stats"]["nKeptPolygons"])
            so["axis_shift"] = {"mu_h": float(np.abs(inp["mu_h"] - base["mu_h"]).max()),
                                "mu_b": float(np.abs(inp["mu_b"] - base["mu_b"]).max())}
        for seed in a.seeds:
            if len(inp["Xtrain"]["BLD"]) == 0:
                # gate only: HPY vs stressed (PLE) with the same numpy draw as the full recipe
                ts = build_training_sets(inp, "numpy", seed)
                gate = grow_forest(ts["Fgate"], ts["ygate"], n_trees=250, seed=seed)
                stress = None
            else:
                gate, stress, ts = train_stage2(inp, seed, "numpy", "matlab_cart")
            r = evaluate(inp, gate, stress, F_dates)
            r["seed"] = seed
            r["gate_rows"] = {"HPY": int((ts["ygate"] == 1).sum()), "stressed": int((ts["ygate"] == 2).sum()),
                              "stressed_by_date": gate_rows_by_date(inp, seed, ts)}
            r["stress_rows"] = {"PLE": int((ts["ystress"] == 1).sum()), "BLD": int((ts["ystress"] == 2).sum())}
            so["runs"].append(r)
            show = {d: (r["dates"][d].get("HPY_PLE_BLD_UNK") or r["dates"][d]["gate_HPY_STR_UNK"]) for d in JUNE + EXTRA}
            print(f"{mode} {setting} seed {seed}: " + " | ".join(f"{d[4:]} {v}" for d, v in show.items())
                  + f" [{time.time() - t0:.0f}s]", flush=True)
            with open(a.json, "w") as fh:
                json.dump(out, fh)
print("done", flush=True)
