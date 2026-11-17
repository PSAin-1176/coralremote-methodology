"""Py4: train the Stage 2 forests with scikit-learn on the paper's training
pixels and compare with the paper of record (frozen forests = paper maps).

Done-when (patch plan §7.4): pixel agreement >= 97 % with the paper maps,
every Table 5 cell within 1 pp, OA within 1 pp."""
import argparse, json, time
import numpy as np
from coralremote import io, paths, stage2
from coralremote.forest import load_treebagger
from coralremote.train import train_stage2

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--rows", default="compat", choices=["compat", "numpy"])
ap.add_argument("--json", default=None)
ap.add_argument("--engine", default="matlab_cart", choices=["matlab_cart", "sklearn"])
a = ap.parse_args()
if a.json:
    paths.assert_writable(a.json)

t0 = time.time()
inp = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode="compat"))
fg = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
fs = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
ref = stage2.run_frozen(fg, fs, inp, return_maps=True)
support = inp["support"]
print(f"paper (frozen) OA {ref['metrics']['OA']:.4f}  [{time.time()-t0:.0f}s]")

out = []
for seed in a.seeds:
    g, s, ts = train_stage2(inp, seed, a.rows, a.engine)
    r = stage2.run_frozen(g, s, inp, return_maps=True)
    agree = {d: float((r["maps"][d][support] == ref["maps"][d][support]).mean()) for d in ref["maps"]}
    all_agree = float(np.mean([agree[d] for d in agree]))
    cells = [(row["date"], k, row[k] - q[k]) for row, q in zip(r["date_summary"], ref["date_summary"])
             for k in ("HPY", "PLE", "BLD", "UNK")]
    maxcell = max(cells, key=lambda c: abs(c[2]))
    m, m0 = r["metrics"], ref["metrics"]
    res = dict(seed=seed, OA=m["OA"], dOA=m["OA"] - m0["OA"], kappa=m["kappa"],
               f1=dict(zip(("HPY", "PLE", "BLD"), m["f1"][:3].round(4).tolist())),
               oob_gate=getattr(g, 'oob_error', None), oob_stress=getattr(s, 'oob_error', None),
               pixel_agreement=all_agree, min_date_agreement=min(agree.values()),
               worst_date=min(agree, key=agree.get), max_cell=maxcell,
               n_cells_over_1pp=sum(abs(c[2]) > 1 for c in cells),
               anchors=(round(r["healthy_anchor_recall"], 4), round(r["bleached_anchor_recall"], 4)),
               table5=r["date_summary"], C=m["C"].tolist(), per_date_agreement=agree)
    out.append(res)
    print(f"seed {seed}: OA {m['OA']:.4f} (d {res['dOA']*100:+.2f} pp) kappa {m['kappa']:.4f} F1 {res['f1']} "
          f"| agree {all_agree*100:.2f}% "
          f"(min {res['min_date_agreement']*100:.2f}% {res['worst_date']}) | max cell {maxcell[0]} {maxcell[1]} "
          f"{maxcell[2]:+.2f} pp, cells>1pp {res['n_cells_over_1pp']} | anchors {res['anchors']} [{time.time()-t0:.0f}s]")
if a.json:
    with open(a.json, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
