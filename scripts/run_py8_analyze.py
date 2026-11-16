"""Py8 step 2: why do the 20240514 and 20250518 PLE/UNK fractions swing
between retraining seeds? Reads the cache written by run_py8_collect.py.

Q1  where UNK comes from (gate vs stress) and how many pixels sit near 0.55;
Q2  per-pixel cross-seed flips: which model disagrees, between which classes;
Q3  out-of-range features and whether flips coincide with them;
Q4  250 vs 1,000 trees; per-tree spread; bootstrap-only (rows fixed) runs.
Writes a JSON summary next to the cache files.
"""
import argparse
import json
from itertools import combinations

import numpy as np

from coralremote import paths
from coralremote.features import FEATURE_NAMES

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
a = ap.parse_args()
D = paths.cache_dir() / "py8"
UNSTABLE = ("20240514", "20250518")
NAMES = {1: "HPY", 2: "PLE", 3: "BLD", 4: "UNK"}
THR = 0.55


def decide(g, s):
    """Hierarchical rule from class-1 scores: returns pred and the UNK source."""
    gate1 = g >= 0.5                                  # ties go to the first class, like max()
    gconf = np.maximum(g, 1 - g)
    sconf = np.maximum(s, 1 - s)
    pred = np.where(gate1, 1, np.where(s >= 0.5, 2, 3))
    unk_gate = gconf < THR
    unk_stress = ~unk_gate & ~gate1 & (sconf < THR)
    pred = np.where(unk_gate | unk_stress, 4, pred)
    return pred, unk_gate, unk_stress, gate1


def gate_state(g):
    """Gate decision on its own: 0 HPY, 1 stressed, 2 gate-UNK."""
    return np.where(np.maximum(g, 1 - g) < THR, 2, np.where(g >= 0.5, 0, 1))


summary = {}
for mode in ("compat", "fixed"):
    feats = np.load(D / f"features_{mode}.npz")
    runs = [np.load(D / f"{mode}_seed{s}.npz") for s in a.seeds]
    dates = sorted(k[2:] for k in feats.files if k.startswith("F_"))
    tmin = np.min([r["train_min"] for r in runs], 0)
    tmax = np.max([r["train_max"] for r in runs], 0)
    out = {}
    # Py6 check: the first 250 trees reproduce the cached 250-tree maps
    check = []
    for s, r in zip(a.seeds, runs):
        m = np.load(paths.cache_dir() / "py6_maps" / f"{mode}_seed{s}.npz")
        same = []
        for d in dates:
            rr, cc = feats[f"rc_{d}"]
            same.append(np.array_equal(decide(r[f"g250_{d}"], r[f"s250_{d}"])[0], m[d][rr, cc]))
        check.append(all(same))
    out["first250_equals_py6_maps"] = bool(all(check))

    for d in dates:
        F = feats[f"F_{d}"]
        n = len(F)
        rec = {"n": n}
        for T in (250, 1000):
            P, UG, US, G1, GS, SS = [], [], [], [], [], []
            for r in runs:
                g, s = r[f"g{T}_{d}"], r[f"s{T}_{d}"]
                p, ug, us, g1 = decide(g, s)
                P.append(p); UG.append(ug); US.append(us); G1.append(g1); GS.append(g); SS.append(s)
            P = np.array(P); GS = np.array(GS); SS = np.array(SS)
            frac = {NAMES[k]: (100 * (P == k).mean(1)).tolist() for k in NAMES}
            rec[f"T{T}"] = {
                "fractions": frac,
                "range_pp": {k: float(np.ptp(v)) for k, v in frac.items()},
                "sd_pp": {k: float(np.std(v, ddof=1)) for k, v in frac.items()},
                # Q1 sources of UNK (% of in-mask pixels, mean over seeds)
                "unk_from_gate": float(100 * np.mean([u.mean() for u in UG])),
                "unk_from_stress": float(100 * np.mean([u.mean() for u in US])),
                "near_gate_thr": float(100 * np.mean(np.abs(np.maximum(GS, 1 - GS) - THR) <= 0.05)),
                "near_stress_thr_routed": float(100 * np.mean([(~g1 & (np.abs(np.maximum(s, 1 - s) - THR) <= 0.05)).mean()
                                                            for g1, s in zip(G1, SS)])),
                "cross_seed_sd_gate_score": float(np.std(GS, 0, ddof=1).mean()),
                "cross_seed_sd_stress_score": float(np.std(SS, 0, ddof=1).mean()),
                "tree_sd_gate": float(np.mean([r[f"g{T}sd_{d}"].mean() for r in runs])),
                "tree_sd_stress": float(np.mean([r[f"s{T}sd_{d}"].mean() for r in runs])),
            }
            if T == 250:
                # Q2 flips: pixels whose final label is not the same in every seed
                flip = (P != P[0]).any(0)
                gst = np.array([gate_state(g) for g in GS])
                gate_flip = (gst != gst[0]).any(0)
                pairs = {}
                for i, j in combinations(range(len(P)), 2):
                    dif = P[i] != P[j]
                    for x, y in zip(P[i][dif], P[j][dif]):
                        k = "-".join(sorted((NAMES[x], NAMES[y])))
                        pairs[k] = pairs.get(k, 0) + 1
                tot = sum(pairs.values()) or 1
                rec["flip_pct"] = float(100 * flip.mean())
                rec["flip_gate_pct"] = float(100 * (flip & gate_flip).mean())
                rec["flip_stress_only_pct"] = float(100 * (flip & ~gate_flip).mean())
                rec["flip_pairs_pct"] = {k: round(100 * v / tot, 1) for k, v in sorted(pairs.items(), key=lambda kv: -kv[1])}
                # Q3 out of range (training range over all seeds' training matrices)
                oor = (F < tmin) | (F > tmax)
                anyo = oor.any(1)
                rec["oor_any_pct"] = float(100 * anyo.mean())
                rec["oor_axis_t_pct"] = float(100 * oor[:, FEATURE_NAMES.index("axis_t")].mean())
                top = np.argsort(-oor.mean(0))[:5]
                rec["oor_top_features"] = {FEATURE_NAMES[i]: round(100 * float(oor[:, i].mean()), 1) for i in top}
                rec["flip_given_oor"] = float(100 * flip[anyo].mean()) if anyo.any() else None
                rec["flip_given_in_range"] = float(100 * flip[~anyo].mean()) if (~anyo).any() else None
                rec["oor_among_flips"] = float(100 * anyo[flip].mean()) if flip.any() else None
                rec["tree_sd_gate_flip_vs_stable"] = [float(np.mean([r[f"g{T}sd_{d}"][flip].mean() for r in runs])) if flip.any() else None,
                                                      float(np.mean([r[f"g{T}sd_{d}"][~flip].mean() for r in runs]))]
                rec["tree_sd_stress_flip_vs_stable"] = [float(np.mean([r[f"s{T}sd_{d}"][flip].mean() for r in runs])) if flip.any() else None,
                                                        float(np.mean([r[f"s{T}sd_{d}"][~flip].mean() for r in runs]))]
        # Q4b bootstrap only (rows fixed at seed 0; fixed units, 250 trees)
        if mode == "fixed":
            bo = [runs[0]] + [np.load(D / f"fixed_rows0_seed{s}.npz") for s in a.seeds[1:]
                              if (D / f"fixed_rows0_seed{s}.npz").exists()]
            if len(bo) > 1:
                P = np.array([decide(r[f"g250_{d}"], r[f"s250_{d}"])[0] for r in bo])
                fr = {NAMES[k]: (100 * (P == k).mean(1)).tolist() for k in NAMES}
                rec["bootstrap_only_T250"] = {"n_runs": len(bo), "fractions": fr,
                                              "range_pp": {k: float(np.ptp(v)) for k, v in fr.items()}}
        out[d] = rec
    summary[mode] = out

json.dump(summary, open(D / "py8_summary.json", "w"), indent=1)
for mode, out in summary.items():
    print(f"\n=== {mode}  (first 250 trees = Py6 maps: {out['first250_equals_py6_maps']})")
    print("date      PLE/UNK range 250T | 1000T | boot-only 250T | UNK gate/stress % | near thr gate/stress % | flip% (gate / stress-only) | OOR any / axis_t %")
    for d, r in out.items():
        if d == "first250_equals_py6_maps":
            continue
        a250, a1k = r["T250"]["range_pp"], r["T1000"]["range_pp"]
        bo = r.get("bootstrap_only_T250", {}).get("range_pp")
        t = r["T250"]
        print(f"{d}  {a250['PLE']:5.1f}/{a250['UNK']:5.1f} | {a1k['PLE']:5.1f}/{a1k['UNK']:5.1f} | "
              + (f"{bo['PLE']:5.1f}/{bo['UNK']:5.1f}" if bo else "    -     ")
              + f" | {t['unk_from_gate']:5.1f}/{t['unk_from_stress']:5.1f} | {t['near_gate_thr']:5.1f}/{t['near_stress_thr_routed']:5.1f}"
              + f" | {r['flip_pct']:5.1f} ({r['flip_gate_pct']:4.1f}/{r['flip_stress_only_pct']:4.1f}) | {r['oor_any_pct']:5.1f}/{r['oor_axis_t_pct']:5.1f}")
    for d in UNSTABLE:
        r = out[d]
        print(f"\n{mode} {d}: flip pairs {r['flip_pairs_pct']}")
        print(f"   OOR top features {r['oor_top_features']}; P(flip|OOR) {r['flip_given_oor']}, "
              f"P(flip|in range) {r['flip_given_in_range']}, OOR among flips {r['oor_among_flips']}")
        print(f"   per-tree SD gate flip/stable {np.round(r['tree_sd_gate_flip_vs_stable'], 3).tolist()}, "
              f"stress {np.round(r['tree_sd_stress_flip_vs_stable'], 3).tolist()}; "
              f"cross-seed SD of scores gate/stress 250T {r['T250']['cross_seed_sd_gate_score']:.3f}/{r['T250']['cross_seed_sd_stress_score']:.3f}, "
              f"1000T {r['T1000']['cross_seed_sd_gate_score']:.3f}/{r['T1000']['cross_seed_sd_stress_score']:.3f}")
