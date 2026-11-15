"""Py4 for Stage 1: retrain the substrate forest (MATLAB-rule CART with
cost-adjusted weights) on the paper's training pixels and compare with the
frozen paper forest (report metrics + per-pixel reef maps on 8 dates)."""
import argparse
import time

import numpy as np
from matio import load_from_mat

from coralremote import paths, stage1
from coralremote.cart import grow_forest
from coralremote.forest import load_treebagger
from coralremote.matlab_compat import find_colmajor
from coralremote.preprocess import phase5_features, reef_mask

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
a = ap.parse_args()

mat = load_from_mat(str(paths.STAGE1_MODEL))
frozen = load_treebagger(None, "rf_model", mat=mat)
inp = stage1.prepare_inputs()
reef = reef_mask()


def maps(f):
    out = {}
    for d in stage1.TRANSFER_TEST_DATES:
        F = phase5_features(d)
        r, c = find_colmajor(reef & np.isfinite(F).all(2))
        out[d] = f.predict(F[r, c])[0]
    return out


ref_maps = maps(frozen)
ref = stage1.run_frozen(frozen, inp)
print(f"paper: OOB 0.1004 secondary OA {ref['secondary']['OA']:.4f} kappa {ref['secondary']['kappa']:.4f} "
      f"F1 {np.round(ref['secondary']['f1'], 3).tolist()}")
for seed in a.seeds:
    t = time.time()
    f, oob = grow_forest(inp["Xtrain"], inp["Ytrain"], n_trees=200, seed=seed, classes=(1, 2, 3),
                         cost=stage1.RF_COST, return_oob=True)
    nodes = np.mean([len(x.cut_var) for x in f.trees])
    r = stage1.run_frozen(f, inp)
    m = maps(f)
    agree = {d: float((m[d] == ref_maps[d]).mean()) for d in m}
    s = r["secondary"]
    print(f"seed {seed}: grow {time.time() - t:.0f}s nodes {nodes:.0f} (MATLAB 1540) "
          f"OOB {stage1.oob_error(f, inp['Xtrain'], inp['Ytrain'], oob):.4f} "
          f"secondary OA {s['OA']:.4f} kappa {s['kappa']:.4f} F1 {np.round(s['f1'], 3).tolist()} | "
          f"reef-map agreement {np.mean(list(agree.values())) * 100:.2f}% "
          f"(min {min(agree.values()) * 100:.2f}% {min(agree, key=agree.get)})")
    print("   transfer OA py(paper):", " ".join(f"{x['date'][2:]} {x['OA']:.3f}({y['OA']:.3f})"
                                         for x, y in zip(r["transfer"], ref["transfer"])))
    print("   transfer CORAL F1 py(paper):", " ".join(f"{x['coral_f1']:.3f}({y['coral_f1']:.3f})"
                                               for x, y in zip(r["transfer"], ref["transfer"])))
