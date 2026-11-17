"""R1 add-on g: classify the PLE candidate polygons that the curated split
did not use, with the frozen paper Stage 2 forests (compat units).

Candidates = train-eligible PLE polygons in the unified training labels
(the set ``select_middle_pale`` scores). Used = curated and in the train
block (2) or test block (7). Each polygon is classified on its own date;
results are pixel-weighted per group, plus polygon counts by majority class.
Groups follow the curation outcome and the polygon's centroid position on
the healthy (t = 0) to bleached (t = 1) axis. Aggregates only.
"""
import json

import numpy as np

from coralremote import io, paths, stage2
from coralremote.features import build_hybrid_features, compute_feature_base
from coralremote.forest import load_treebagger
from coralremote.rasterize import feature_pixels
from coralremote.units import COMPAT

CLS = ("HPY", "PLE", "BLD", "UNK")
opts = stage2.Stage2Options(mode=COMPAT)
inp = stage2.prepare_frozen_inputs(opts)
view, support = inp["view"], inp["support"]
fg = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
fs = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")

cands = stage2.filter_features(io.load_labels("unified_train"), {"PLE"}, True)
ps = inp["ple_stats"]
keep, meanT, meanR = ps["keepMask"], ps["meanT"], ps["meanR"]
blocks = np.array([stage2._spatial_block(f["properties"]) for f in cands])
used = keep & np.isin(blocks, list(opts.ple_train_blocks) + list(opts.ple_test_blocks))


def reason(i):
    if used[i]:
        return "used (curated, block 2 or 7)"
    if keep[i]:
        return "curated, other block (unused)"
    if not np.isfinite(meanT[i]):
        return "rejected: centroid outside mask / no data"
    if meanT[i] <= ps["hKeep"]:
        return "rejected: t below band (healthy side)"
    if meanT[i] >= ps["bKeep"]:
        return "rejected: t above band (bleached side)"
    return "rejected: residual above cut"


def tbin(t):
    if not np.isfinite(t):
        return "no centroid value"
    for lo, hi in ((-np.inf, 0), (0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0), (1.0, np.inf)):
        if lo <= t < hi:
            return f"t {lo:g} to {hi:g}"


# Stage 2 training pixels per date (HPY, BLD, used PLE), to flag overlap
train_sets = (stage2.filter_features(io.load_labels("unified_train"), {"HPY"}, True)
              + stage2.filter_features(io.load_labels("strict_train"), {"BLD"}, True)
              + [f for f, u in zip(cands, used) if u])
train_mask = {}
for d in view.dates:
    m, _ = stage2.rasterize_date_labels(train_sets, view.ref, {"HPY": 1, "BLD": 3, "PLE": 2}, d)
    train_mask[d] = m & support
rows = []
for i, f in enumerate(cands):
    d = str(f["properties"]["date"])
    if d not in view.dates:
        continue
    rr, cc = feature_pixels(f, view.ref)
    if rr.size == 0:
        continue
    A = view.get(d)
    ok = support[rr, cc] & np.isfinite(A[rr, cc]).all(1)
    if not ok.any():
        continue
    F = build_hybrid_features(compute_feature_base(A[rr[ok], cc[ok]], view.mode),
                              inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])
    p, _, _ = stage2.hierarchical_predict(F, fg, fs, opts)
    cnt = np.bincount(p, minlength=5)[1:]
    rows.append(dict(i=i, date=d, reason=reason(i), tbin=tbin(meanT[i]), unused=not used[i],
                     counts=cnt, major=CLS[int(np.argmax(cnt))],
                     n_train_overlap=int(train_mask[d][rr[ok], cc[ok]].sum())))

out = {"n_candidates": len(cands), "n_with_mask_pixels": len(rows),
       "n_used": int(used.sum()), "n_curated": int(keep.sum()),
       "hKeep": float(ps["hKeep"]), "bKeep": float(ps["bKeep"]), "rPleCut": float(ps["rPleCut"])}


def summarise(sel):
    c = np.sum([r["counts"] for r in sel], 0) if sel else np.zeros(4)
    n = int(c.sum())
    maj = {k: sum(r["major"] == k for r in sel) for k in CLS}
    ov = sum(r["n_train_overlap"] for r in sel)
    return {"polygons": len(sel), "pixels": n, "pct_pixels_in_training": round(100 * ov / max(n, 1), 1),
            "pct": [round(100 * float(x) / max(n, 1), 1) for x in c], "majority_polygons": maj}


groups = {}
groups["ALL unused candidates"] = summarise([r for r in rows if r["unused"]])
groups["used (for reference)"] = summarise([r for r in rows if not r["unused"]])
for g in sorted({r["reason"] for r in rows}):
    groups[g] = summarise([r for r in rows if r["reason"] == g])
for g in sorted({r["tbin"] for r in rows if r["unused"]}, key=lambda s: (s.startswith("no"), float(s.split()[1]) if not s.startswith("no") else 0)):
    groups[f"unused, {g}"] = summarise([r for r in rows if r["unused"] and r["tbin"] == g])
for d in sorted({r["date"] for r in rows if r["unused"]}):
    groups[f"unused, date {d}"] = summarise([r for r in rows if r["unused"] and r["date"] == d])
out["groups"] = groups
json.dump(out, open(paths.cache_dir() / "r1_uncurated_ple.json", "w"), indent=1)

print({k: v for k, v in out.items() if k != "groups"})
print("| Group | Polygons | Pixels | % pixels also Stage 2 training pixels | HPY | PLE | BLD | UNK | Polygons by majority class (HPY/PLE/BLD/UNK) |")
print("|---|---|---|---|---|---|---|---|---|")
for k, v in groups.items():
    m = v["majority_polygons"]
    print(f"| {k} | {v['polygons']} | {v['pixels']} | {v['pct_pixels_in_training']:.1f} | " + " | ".join(f"{x:.1f}" for x in v["pct"])
          + f" | {m['HPY']}/{m['PLE']}/{m['BLD']}/{m['UNK']} |")
