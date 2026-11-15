"""Py3: evaluate the paper's frozen Stage 2 forests on Python-built features
and compare with the paper-of-record outputs (Table 5, confusion matrix)."""
import sys, time
import numpy as np
from coralremote import io, paths, stage2, features as fe
from coralremote.evaluate import compare_date_summary, format_metric_block, parse_report_confusion
from coralremote.forest import load_treebagger

t0 = time.time()
gate, gprops = load_treebagger(paths.STAGE2_OUTPUT, "rfGate", return_object=True)
stress, sprops = load_treebagger(paths.STAGE2_OUTPUT, "rfStress", return_object=True)
print(f"loaded forests: gate {len(gate.trees)} trees {gate.class_names}, stress {len(stress.trees)} trees ({time.time()-t0:.1f}s)")

inp = stage2.prepare_frozen_inputs()
ref = io.load_stage2_output()

# 1. training features: every PrivateX row must be an exact row of the Python features
def feats(X):
    return fe.build_hybrid_features(fe.compute_feature_base(X), inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])
Xtr = inp["Xtrain"]
pool_gate = np.concatenate([feats(Xtr["HPY"]), feats(np.concatenate([Xtr["PLE"], Xtr["BLD"]]))])
pool_stress = np.concatenate([feats(Xtr["PLE"]), feats(Xtr["BLD"])])
for name, props, pool in (("gate", gprops, pool_gate), ("stress", sprops, pool_stress)):
    PX = np.asarray(props["PrivateX"], np.float32)
    s = {r.tobytes() for r in pool}
    hit = sum(r.tobytes() in s for r in PX)
    print(f"{name} PrivateX rows found bit-exact in Python features: {hit}/{len(PX)}")

res = stage2.run_frozen(gate, stress, inp)
gp, sp = ref["gatePredTest"], ref["stressPredTest"]
print("gatePredTest identical:", np.array_equal(res["gate_pred"], gp), int((res["gate_pred"] != gp).sum()), "diffs")
print("stressPredTest identical:", np.array_equal(res["stress_pred"], sp), int((res["stress_pred"] != sp).sum()), "diffs")

C_ref = parse_report_confusion(paths.STAGE2_REPORT.read_text())
print("\nSpatial test metrics (Python):\n" + format_metric_block(res["metrics"]))
print("confusion matrix identical:", np.array_equal(res["metrics"]["C"], C_ref))
print(f"anchors: HPY {res['healthy_anchor_recall']:.4f} (n={res['n_healthy_anchor']}), "
      f"BLD {res['bleached_anchor_recall']:.4f} (n={res['n_bleached_anchor']})")

rows_ref = io.read_date_summary_csv()
n, bad, maxd, det = compare_date_summary(res["date_summary"], rows_ref, tol=1e-9)
print(f"\nTable 5: {n} cells, {bad} differ, max |diff| = {maxd:.3g} pp")
for d in det: print("  ", d)
print("date      HPY      PLE      BLD      UNK")
for r in res["date_summary"]:
    print(f"{r['date']} {r['HPY']:8.4f} {r['PLE']:8.4f} {r['BLD']:8.4f} {r['UNK']:8.4f}")
print(f"total {time.time()-t0:.1f}s")
sys.exit(0 if (bad == 0 and np.array_equal(res["metrics"]["C"], C_ref)) else 1)
