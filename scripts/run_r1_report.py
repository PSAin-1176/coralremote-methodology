"""Markdown tables for docs/R1_LOO_DATE.md from the leave-date-out JSON
(run_r1_loo_date.py). Prints to stdout."""
import argparse
import json

import numpy as np

from coralremote import paths

ap = argparse.ArgumentParser()
ap.add_argument("--json", default=str(paths.cache_dir() / "r1_loo_date.json"))
a = ap.parse_args()
R = json.load(open(a.json))
JUNE = ("20240613", "20240616", "20240620")
CLS = ("HPY", "PLE", "BLD", "UNK")


def stats(vals):
    v = np.asarray(vals, float)
    return v.mean(0), v.min(0), v.max(0)


def fmt(m, lo, hi):
    return f"{m:.1f} [{lo:.1f}–{hi:.1f}]"


def runs_frac(setting, mode, date, key="HPY_PLE_BLD_UNK"):
    return [r["dates"][date][key] for r in R["modes"][mode]["settings"][setting]["runs"]]


fz = R["frozen_paper"]
for mode in R["modes"]:
    S = R["modes"][mode]["settings"]
    print(f"\n### {mode}: held-out date classified (mean [min–max] over {len(S['baseline']['runs'])} seeds, % of mask pixels)\n")
    print("| Held-out date | Class | Frozen paper model | Baseline (nothing removed) | Rows of D removed | Δ | D removed incl. axis + PLE curation | Δ |")
    print("|---|---|---|---|---|---|---|---|")
    for d in JUNE:
        b = stats(runs_frac("baseline", mode, d))
        rw = stats(runs_frac(f"rows_{d}", mode, d)) if f"rows_{d}" in S else None
        st = stats(runs_frac(f"strict_{d}", mode, d)) if f"strict_{d}" in S else None
        for k, c in enumerate(CLS):
            row = [d if k == 0 else "", c, f"{fz['dates'][d]['HPY_PLE_BLD_UNK'][k]:.1f}" if mode == "compat" else "–",
                   fmt(b[0][k], b[1][k], b[2][k])]
            for x in (rw, st):
                row += [fmt(x[0][k], x[1][k], x[2][k]), f"{x[0][k] - b[0][k]:+.1f}"] if x else ["–", "–"]
            print("| " + " | ".join(row) + " |")

    print(f"\n### {mode}: BLD on every June date when one date is removed (mean over seeds)\n")
    print("| Setting | " + " | ".join(f"BLD {d}" for d in JUNE) + " | PLE 20240514 | UNK 20240514 |")
    print("|---|" + "---|" * (len(JUNE) + 2))
    for s in [x for x in S if x != "alljune"]:
        m = {d: np.mean([r["dates"][d]["HPY_PLE_BLD_UNK"] for r in S[s]["runs"]], 0) for d in JUNE + ("20240514",)}
        print(f"| {s} | " + " | ".join(f"{m[d][2]:.1f}" for d in JUNE)
              + f" | {m['20240514'][1]:.1f} | {m['20240514'][3]:.1f} |")

    print(f"\n### {mode}: test BLD polygons on the held-out date (prediction shares, mean over seeds)\n")
    print("| Date | n test BLD px | Frozen paper BLD recall | Baseline BLD / PLE / UNK | Rows of D removed BLD / PLE / UNK | Strict BLD / PLE / UNK |")
    print("|---|---|---|---|---|---|")
    for d in JUNE:
        key = f"BLD_test_{d}"

        def t(s):
            if s not in S:
                return "–"
            v = np.mean([r["test"][key]["pred"] for r in S[s]["runs"]], 0)
            return f"{v[2]:.1f} / {v[1]:.1f} / {v[3]:.1f}"
        n = S["baseline"]["runs"][0]["test"][key]["n"]
        fzr = f"{fz['test'][key]['pred'][2]:.1f}" if mode == "compat" else "–"
        print(f"| {d} | {n} | {fzr} | {t('baseline')} | {t(f'rows_{d}')} | {t(f'strict_{d}')} |")

    print(f"\n### {mode}: test-set OA (mean [min–max])\n")
    print("| Setting | OA | F1 HPY / PLE / BLD |")
    print("|---|---|---|")
    for s in [x for x in S if x != "alljune"]:
        oa = [r["test"]["OA"] for r in S[s]["runs"]]
        f1 = np.mean([r["test"]["f1"] for r in S[s]["runs"]], 0)
        print(f"| {s} | {np.mean(oa):.4f} [{min(oa):.4f}–{max(oa):.4f}] | {f1[0]:.3f} / {f1[1]:.3f} / {f1[2]:.3f} |")

    if "alljune" in S:
        print(f"\n### {mode}: all three June dates removed (gate only; % HPY / stressed / UNK, mean [min–max])\n")
        print("| Date | Baseline gate | All June removed | Δ stressed |")
        print("|---|---|---|---|")
        for d in JUNE + ("20240514",):
            b = stats(runs_frac("baseline", mode, d, "gate_HPY_STR_UNK"))
            x = stats(runs_frac("alljune", mode, d, "gate_HPY_STR_UNK"))
            print(f"| {d} | " + " / ".join(fmt(b[0][k], b[1][k], b[2][k]) for k in range(3)) + " | "
                  + " / ".join(fmt(x[0][k], x[1][k], x[2][k]) for k in range(3)) + f" | {x[0][1] - b[0][1]:+.1f} |")
        print("\nTest BLD polygons under the all-June gate (HPY / stressed / UNK %):",
              {d: np.round(np.mean([r["test"][f"BLD_test_{d}"]["pred"] for r in S["alljune"]["runs"]], 0), 1).tolist()
               for d in JUNE})

    print(f"\n### {mode}: training pixels by date\n")
    print("| Setting | BLD pool by date | PLE pool by date | Gate stressed rows drawn (seed 0): BLD by date; PLE by date | Stress model rows PLE / BLD |")
    print("|---|---|---|---|---|")
    for s in S:
        pb = S[s]["pool_by_date"]
        r0 = S[s]["runs"][0]
        gr = r0["gate_rows"]["stressed_by_date"]
        print(f"| {s} | {pb['BLD'] or 'none'} (n={sum(pb['BLD'].values())}) | {pb['PLE']} (n={sum(pb['PLE'].values())}) | "
              f"{gr['BLD'] or 'none'}; {gr['PLE']} (HPY {r0['gate_rows']['HPY']}, stressed {r0['gate_rows']['stressed']}) | "
              f"{r0['stress_rows']['PLE']} / {r0['stress_rows']['BLD']} |")
    for s in S:
        if "axis_shift" in S[s]:
            print(f"\n{mode} {s}: axis means max shift mu_h {S[s]['axis_shift']['mu_h']:.2e}, "
                  f"mu_b {S[s]['axis_shift']['mu_b']:.2e}; PLE polygons kept {S[s]['ple_kept_polygons']}")

print("\n### Frozen paper model (Table 5 reference), % HPY / PLE / BLD / UNK\n")
for d in JUNE + ("20240514",):
    print(f"- {d}: " + " / ".join(f"{x:.2f}" for x in fz["dates"][d]["HPY_PLE_BLD_UNK"]))
print("frozen test OA check:", fz.get("test_metrics_check"))
