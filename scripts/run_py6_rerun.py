"""Py6 / P6: re-run Stage 2 on the fixed units (true reflectance, 1e-4 floor)
and, as the like-for-like reference, on the paper's units (compat), with the
same procedure, seeds and sampling stream (MATLAB-rule CART, 250 trees,
MinLeaf 5, 0.55 UNK rule, curated PLE). Every deterministic step (axis
means, thresholds, PLE curation, block split) is recomputed in each mode.

Each fixed-mode forest is also scored on the paper's own test pixels
(compat test set, features in fixed units), so test metrics can be compared
on identical pixels. Results go to JSON; per-date class maps go to the
cache folder (outside the repo).
"""
import argparse
import json
import time

import numpy as np

from coralremote import paths, stage2
from coralremote.evaluate import compute_metrics3
from coralremote.forest import load_treebagger
from coralremote.stage2 import hierarchical_predict
from coralremote.features import build_hybrid_features, compute_feature_base
from coralremote.train import train_stage2
from coralremote.units import COMPAT, FIXED, SR_SCALE

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--modes", nargs="+", default=["compat", "fixed"])
ap.add_argument("--json", default=str(paths.cache_dir() / "py6_rerun.json"))
a = ap.parse_args()

t0 = time.time()
inputs = {"compat": stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=COMPAT))}
if "fixed" in a.modes:
    inputs["fixed"] = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=FIXED))

fg = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
fs = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
paper = stage2.run_frozen(fg, fs, inputs["compat"])


def summarize(r):
    m = r["metrics"]
    return dict(OA=m["OA"], kappa=m["kappa"], f1=m["f1"][:3].tolist(), C=np.asarray(m["C"]).tolist(),
                anchors=[r["healthy_anchor_recall"], r["bleached_anchor_recall"]],
                table5=[{k: row[k] for k in ("date", "HPY", "PLE", "BLD", "UNK")} for row in r["date_summary"]])


def det(inp):
    return dict(n_train=inp["n_train"], n_test=np.bincount(inp["Ytest"], minlength=4)[1:].tolist(),
                thresholds={k: float(v) for k, v in inp["thresholds"].items()},
                ple_kept=int(inp["ple_stats"]["nKeptPolygons"]), ple_total=int(inp["ple_stats"]["nTotalPolygons"]),
                ple_split=[int(inp["split_stats"]["nTrainPolygons"]), int(inp["split_stats"]["nTestPolygons"])])


out = {"paper": summarize(paper), "deterministic": {m: det(inputs[m]) for m in inputs}, "runs": {}}
print(f"paper OA {paper['metrics']['OA']:.4f} [{time.time() - t0:.0f}s]")
cdir = paths.cache_dir() / "py6_maps"
cdir.mkdir(exist_ok=True)
Xpaper_fixed = inputs["compat"]["Xtest"] * np.float32(SR_SCALE)
for mode in a.modes:
    inp = inputs[mode]
    out["runs"][mode] = []
    for seed in a.seeds:
        g, s, _ = train_stage2(inp, seed, "numpy", "matlab_cart")
        r = stage2.run_frozen(g, s, inp, return_maps=True)
        res = summarize(r)
        res["seed"] = seed
        if mode == "fixed":
            F = build_hybrid_features(compute_feature_base(Xpaper_fixed, FIXED),
                                      inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])
            p, _, _ = hierarchical_predict(F, g, s, inp["opts"])
            mp = compute_metrics3(inputs["compat"]["Ytest"], p)
            res["on_paper_test"] = dict(OA=mp["OA"], kappa=mp["kappa"], f1=mp["f1"][:3].tolist())
        np.savez_compressed(cdir / f"{mode}_seed{seed}.npz", **{d: m.astype(np.uint8) for d, m in r["maps"].items()})
        out["runs"][mode].append(res)
        print(f"{mode} seed {seed}: OA {res['OA']:.4f} kappa {res['kappa']:.4f} F1 {np.round(res['f1'], 3).tolist()}"
              + (f" | on paper test OA {res['on_paper_test']['OA']:.4f}" if mode == "fixed" else "")
              + f" [{time.time() - t0:.0f}s]", flush=True)
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
print("saved", a.json)
