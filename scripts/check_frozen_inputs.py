"""Run the deterministic Stage 2 steps and compare with the paper report."""
import time, numpy as np
from coralremote import io, stage2
t=time.time()
inp = stage2.prepare_frozen_inputs()
print("time", round(time.time()-t,1))
ref = io.load_stage2_output()
print("thresholds", {k: float(v) for k,v in inp["thresholds"].items()})
ps = inp["ple_stats"]; print("rPleCut", float(ps["rPleCut"]), "hKeep", float(ps["hKeep"]), "bKeep", float(ps["bKeep"]))
print("PLE total/kept", ps["nTotalPolygons"], ps["nKeptPolygons"])
print("keepMask equal:", np.array_equal(ps["keepMask"], np.asarray(ref["pleStats"]["keepMask"]).astype(bool).ravel()))
ss = inp["split_stats"]; print("split", ss["nTrainPolygons"], ss["nTestPolygons"])
print("mu_h diff", np.abs(inp["mu_h"]-ref["muHealthy"]).max(), "mu_b diff", np.abs(inp["mu_b"]-ref["muBleached"]).max())
print("mu_h ulp-equal", np.array_equal(inp["mu_h"], ref["muHealthy"]), np.array_equal(inp["mu_b"], ref["muBleached"]))
print("n_train", inp["n_train"])
print("test support", [int((inp["Ytest"]==k).sum()) for k in (1,2,3)], "len", len(inp["Ytest"]))
print("anchors", len(inp["anchors"]["healthy"][1]), len(inp["anchors"]["bleached"][1]))
