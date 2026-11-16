"""Py6 / P7: compare the fixed-units re-run with the paper of record.

Material change (patch plan section 4, P7) if any Table 5 cell moves > 2 pp,
Stage 2 OA > 1 pp, any class F1 > 0.02, or the June BLD / September HPY
narrative changes. Two comparisons:

A. fixed (mean of seeds) vs the paper as published (frozen forests);
B. fixed vs compat (same retraining, same seeds) = effect of the fix net of
   retraining noise; per-seed paired differences.
"""
import argparse
import json

import numpy as np

from coralremote import paths

ap = argparse.ArgumentParser()
ap.add_argument("--json", default=str(paths.cache_dir() / "py6_rerun.json"))
a = ap.parse_args()
R = json.load(open(a.json))
paper, runs = R["paper"], R["runs"]
CL = ("HPY", "PLE", "BLD", "UNK")
JUNE = ("20240613", "20240616", "20240620")
SEPT = ("20240926",)


def t5(run):
    return {row["date"]: np.array([row[k] for k in CL]) for row in run["table5"]}


P = t5(paper)
dates = list(P)
M = {m: {d: np.array([t5(r)[d] for r in runs[m]]) for d in dates} for m in runs}   # (seeds, 4)

print("Deterministic steps:", json.dumps(R["deterministic"], indent=None))
print("\nTable 5 (%): paper | compat retrain mean [range] | fixed mean [range] | fixed-paper | fixed-compat (paired mean, range)")
flags_a, flags_b = [], []
for d in dates:
    for i, k in enumerate(CL):
        c, f = M["compat"][d][:, i], M["fixed"][d][:, i]
        da, db = f.mean() - P[d][i], f - c
        flags_a += [(d, k, da)] if abs(da) > 2 else []
        flags_b += [(d, k, db.mean())] if abs(db.mean()) > 2 else []
        print(f"{d} {k}: {P[d][i]:6.2f} | {c.mean():6.2f} [{c.min():6.2f},{c.max():6.2f}] | "
              f"{f.mean():6.2f} [{f.min():6.2f},{f.max():6.2f}] | {da:+6.2f} | {db.mean():+6.2f} [{db.min():+6.2f},{db.max():+6.2f}]")


def mstat(m, key, idx=None):
    v = [r[key] if idx is None else r[key][idx] for r in runs[m]]
    return np.mean(v), np.min(v), np.max(v)


print("\nMetrics: paper | compat mean [range] | fixed mean [range] | fixed on the paper's test pixels")
oa_p = paper["OA"]
for name, key, idx in [("OA", "OA", None), ("kappa", "kappa", None), ("F1 HPY", "f1", 0), ("F1 PLE", "f1", 1), ("F1 BLD", "f1", 2)]:
    pv = paper[key] if idx is None else paper[key][idx]
    c, f = mstat("compat", key, idx), mstat("fixed", key, idx)
    op = [r["on_paper_test"][key] if idx is None else r["on_paper_test"][key][idx] for r in runs["fixed"]]
    print(f"{name:7s}: {pv:.4f} | {c[0]:.4f} [{c[1]:.4f},{c[2]:.4f}] | {f[0]:.4f} [{f[1]:.4f},{f[2]:.4f}] | {np.mean(op):.4f}")
print("anchors (healthy, bleached recall): paper", np.round(paper["anchors"], 4),
      "| fixed mean", np.round(np.mean([r["anchors"] for r in runs["fixed"]], 0), 4))


def dominant(tab, d):
    return CL[int(np.argmax(tab[d]))]


print("\nNarrative: dominant class, paper | every fixed seed (June, September and any date that differs)")
narr_diff = []
for d in dates:
    fs = [CL[int(np.argmax(x))] for x in M["fixed"][d]]
    differs = any(f != dominant(P, d) for f in fs)
    if differs:
        narr_diff.append(d)
    if differs or d in JUNE + SEPT:
        print(d, dominant(P, d), "|", fs)
print("(compat retrain, 20250518:", [CL[int(np.argmax(x))] for x in M["compat"]["20250518"]], ")")

print("\nSystematic fix effects (paired fixed - compat, same sign in every seed, |mean| > 1 pp):")
for d in dates:
    for i, k in enumerate(CL):
        db = M["fixed"][d][:, i] - M["compat"][d][:, i]
        if abs(db.mean()) > 1 and (np.all(db > 0) or np.all(db < 0)):
            print(f"  {d} {k}: {db.mean():+.2f} pp [{db.min():+.2f}, {db.max():+.2f}]")

oa_f = mstat("fixed", "OA")[0]
f1_d = max(abs(mstat("fixed", "f1", i)[0] - paper["f1"][i]) for i in range(3))
narr = not narr_diff
print(f"\nP7 criteria (A, fixed mean vs paper): cells > 2 pp: {len(flags_a)} {[(d, k, round(x, 2)) for d, k, x in flags_a]}; "
      f"|dOA| {abs(oa_f - oa_p) * 100:.2f} pp; max |dF1| {f1_d:.4f}; dominant class unchanged on every date and seed: {narr} {narr_diff}")
print(f"P7 criteria (B, fixed vs compat retrain, paired mean): cells > 2 pp: {len(flags_b)} "
      f"{[(d, k, round(x, 2)) for d, k, x in flags_b]}")

# per-pixel map agreement fixed vs compat for the same seed
cdir = paths.cache_dir() / "py6_maps"
from coralremote import io
sup = io.load_baseline_mask()[0]
agree = []
for r in runs["fixed"]:
    s = r["seed"]
    A, B = np.load(cdir / f"compat_seed{s}.npz"), np.load(cdir / f"fixed_seed{s}.npz")
    agree.append({d: float((A[d][sup] == B[d][sup]).mean()) for d in A.files})
mean_by_date = {d: np.mean([x[d] for x in agree]) for d in agree[0]}
print("\nPixel agreement fixed vs compat (same seed), mean over seeds:",
      " ".join(f"{d[2:]} {v * 100:.1f}" for d, v in mean_by_date.items()),
      f"| overall {np.mean(list(mean_by_date.values())) * 100:.2f}%")
