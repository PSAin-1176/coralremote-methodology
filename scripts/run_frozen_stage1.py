"""Py3 for Stage 1: evaluate the paper's frozen Stage 1 forest on rebuilt
Phase 5 features and compare with the saved run (rf_stage1_model.mat)."""
import numpy as np
from matio import load_from_mat
from coralremote import paths, stage1
from coralremote.forest import load_treebagger

mat = load_from_mat(str(paths.STAGE1_MODEL))
forest, props = load_treebagger(None, "rf_model", return_object=True, mat=mat)
print("trees", len(forest.trees), "classes", forest.class_names, "cost\n", forest.cost)
inp = stage1.prepare_inputs()
PX = np.asarray(props["PrivateX"], np.float32)
print("train rows", inp["Xtrain"].shape, "saved", PX.shape)
eq = inp["Xtrain"] == PX
print("train cols exact frac", eq.mean(0).round(4).tolist(), "max abs diff", np.abs(inp["Xtrain"] - PX).max(0).tolist())
Ysaved = np.array([int(str(np.asarray(y).ravel()[0])) for y in np.asarray(props["Y"]).ravel()])
print("labels identical:", np.array_equal(inp["Ytrain"], Ysaved))

oob_py = stage1.oob_error(forest, inp["Xtrain"], inp["Ytrain"], props["OOBIndices"])
oob_saved = stage1.oob_error(forest, PX, Ysaved, props["OOBIndices"])
print(f"OOB error: {oob_py:.4f} on Python features, {oob_saved:.4f} on saved rows; report {float(np.asarray(mat['oob_err']).ravel()[0]):.4f}")

r = stage1.run_frozen(forest, inp)
se = mat["secondary_eval"]
C_ref = np.asarray(se["C"].ravel()[0], float)
s = r["secondary"]
print("secondary test n", int(s["C"].sum()), "C identical:", np.array_equal(s["C"], C_ref))
print(s["C"].astype(int)); print(f"OA {s['OA']:.4f} kappa {s['kappa']:.4f} F1 {np.round(s['f1'], 3).tolist()}")
tr_ref = mat["transfer_results"].ravel()
print("date        n   OA     kappa  C_P    C_R    C_F1   | saved")
allsame = True
for row, q in zip(r["transfer"], tr_ref):
    qv = [float(np.asarray(q[k]).ravel()[0]) for k in ("n", "OA", "kappa", "coral_prec", "coral_rec", "coral_f1")]
    pv = [row[k] for k in ("n", "OA", "kappa", "coral_prec", "coral_rec", "coral_f1")]
    same = np.allclose(pv, qv, atol=1e-12, equal_nan=True); allsame &= same
    print(row["date"], f"{pv[0]:4d} " + " ".join(f"{x:.3f}" for x in pv[1:]), "|", "identical" if same else qv)
print("transfer table identical:", allsame)
