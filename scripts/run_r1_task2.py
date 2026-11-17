"""R1 Task 2: existing numbers for the response letter (no new modelling).

a) wave / surface-state diagnostic for all 14 dates (Py5 port, 14/14 match)
c) Stage 2 training and test counts per class (paper rows recovered from the
   frozen forests)
d) the 20,330-pixel mask vs the reef polygon: Allen Coral Atlas depth and
   distance from the shoreline (median and IQR)
e) 20240514 gate-stressed pixels in the Stage 2 feature space (axis position
   and residual) vs June BLD and 20230509 PLE training pixels

The shoreline is not on disk, so land is taken from the imagery: median
NDWI (Green1 vs NIR) over the CALM dates < 0, keeping connected regions of
at least 1,000 pixels. Aggregates only; nothing raster-shaped is written.
"""
import json

import numpy as np
import scipy.io as sio
from scipy import ndimage

from coralremote import io, paths, stage2, wave
from coralremote.features import FEATURE_NAMES, build_hybrid_features, compute_feature_base
from coralremote.forest import load_treebagger
from coralremote.preprocess import raw_stack
from coralremote.train import build_training_sets
from coralremote.units import COMPAT

ACA_BATHY = str(paths.ACA_BATHY)
OUT = paths.cache_dir() / "r1_task2.json"
res = {}


def q(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None
    p25, p50, p75 = np.percentile(x, [25, 50, 75])
    return {"n": int(x.size), "median": round(float(p50), 4), "q25": round(float(p25), 4), "q75": round(float(p75), 4)}


# ---------------------------------------------------------------- a) wave
s = raw_stack()
wv, thr = wave.run_wave_diagnostic(s.get, s.dates)
res["wave"] = {"dates": {d: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in x.items()}
                         for d, x in wv.items()},
               "thresholds": {k: np.asarray(v).tolist() for k, v in thr.items()},
               "band_names": list(wave.BAND_NAMES),
               "matches_paper": sum(wv[d]["verdict"] == wave.PAPER_VERDICTS[d] for d in wv)}
print("a) wave: date tiles AI(median over bands) WBP(median) roughness-flags AI-flags verdict confidence")
for d, x in wv.items():
    print(f"   {d} {x['n_tiles_used']:3d} {x['ai_median']:.4f} {x['wbp_median']:.3e} "
          f"{x['n_rough_bands']}/8 {x['n_bands_flag']}/8 {x['verdict']} {x['confidence']}")
print("   thresholds per band (AI med/MAD/thr | WBP med/MAD/thr):")
for b, name in enumerate(wave.BAND_NAMES):
    print(f"   {name:8s} {thr['ai_med'][b]:.4f}/{thr['ai_mad'][b]:.4f}/{thr['ai_thr'][b]:.4f} | "
          f"{thr['wbp_med'][b]:.3e}/{thr['wbp_mad'][b]:.3e}/{thr['wbp_thr'][b]:.3e}")

# ---------------------------------------------------------------- c) counts
inp = stage2.prepare_frozen_inputs(stage2.Stage2Options(mode=COMPAT))
ts = build_training_sets(inp, "compat")
nP = len(inp["Xtrain"]["PLE"])
gs, st = ts["rows"]["gate_stressed"], ts["rows"]["stress"]
yt = inp["Ytest"]
res["counts"] = {
    "gate": {"HPY": int(len(ts["rows"]["gate_healthy"])), "PLE": int((gs < nP).sum()), "BLD": int((gs >= nP).sum())},
    "stress": {"PLE": int((st < nP).sum()), "BLD": int((st >= nP).sum())},
    "pools": inp["n_train"],
    "test": {"HPY": int((yt == 1).sum()), "PLE": int((yt == 2).sum()), "BLD": int((yt == 3).sum())},
    "anchors": {"healthy_20230509": int(inp["anchors"]["healthy"][1].size),
                "bleached_20240613": int(inp["anchors"]["bleached"][1].size)},
}
print("c)", res["counts"])

# ---------------------------------------------------------------- d) mask vs reef polygon
m = sio.loadmat(str(paths.BASELINE_MASK), squeeze_me=True)
final = m["maskFinal"] > 0
reef = m["reef_mask"] > 0
ref = s.ref
import rasterio
from rasterio.warp import Resampling, reproject

depth = np.full((ref.height, ref.width), np.nan, np.float32)
with rasterio.open(ACA_BATHY) as src:
    a = src.read(1).astype(np.float32)
    a[a == src.nodata] = np.nan
    reproject(a, depth, src_transform=src.transform, src_crs=src.crs, src_nodata=np.nan,
              dst_transform=rasterio.Affine(*ref.transform), dst_crs=f"EPSG:{ref.epsg}",
              dst_nodata=np.nan, resampling=Resampling.nearest)
depth_m = depth / 100.0                       # ACA depth is in cm, positive down

calm = [d for d in s.dates if wave.PAPER_VERDICTS[d] == "CALM"]
nd = []
for d in calm:
    A = s.get(d)
    g, n = A[..., wave.GREEN1].astype(np.float32), A[..., wave.NIR].astype(np.float32)
    nd.append((g - n) / (g + n))
with np.errstate(invalid="ignore"):
    ndwi = np.nanmedian(np.stack(nd), 0)
land = ndwi < 0
lab, nl = ndimage.label(land)
sizes = ndimage.sum(land, lab, range(1, nl + 1))
land = np.isin(lab, 1 + np.flatnonzero(sizes >= 1000))
dist_m = ndimage.distance_transform_edt(~land) * ref.cell_x

res["mask"] = {"n_mask": int(final.sum()), "n_reef_polygon": int(reef.sum()),
               "mask_inside_polygon": int((final & reef).sum()),
               "land_px": int(land.sum()), "mask_on_land": int((final & land).sum()),
               "calm_dates_for_land": calm}
for name, sel in (("mask", final), ("reef_polygon", reef), ("polygon_outside_mask", reef & ~final)):
    res["mask"][name] = {"depth_m": q(depth_m[sel]), "depth_coverage_pct": round(100 * float(np.isfinite(depth_m[sel]).mean()), 1),
                         "dist_shore_m": q(dist_m[sel])}
print("d)", json.dumps(res["mask"], indent=1))

# ---------------------------------------------------------------- e) 20240514 in feature space
fg = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
fsm = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
it, ir = FEATURE_NAMES.index("axis_t"), FEATURE_NAMES.index("axis_residual")


def F(X):
    return build_hybrid_features(compute_feature_base(X, inp["view"].mode), inp["mu_h"], inp["mu_b"], inp["v"], inp["vnorm2"])


view, support = inp["view"], inp["support"]
A = view.get("20240514")
from coralremote.matlab_compat import find_colmajor
r, c = find_colmajor(support & np.isfinite(A).all(axis=2))
F14 = F(A[r, c, :])
pred, gp, _ = stage2.hierarchical_predict(F14, fg, fsm, inp["opts"])
groups = {
    "20240514 gate-stressed (all)": F14[gp == 2],
    "20240514 classified PLE": F14[pred == 2],
    "20240514 classified UNK": F14[pred == 4],
    "20240514 classified HPY": F14[pred == 1],
}
Db, Dp, Dh = inp["Dtrain"]["BLD"], inp["Dtrain"]["PLE"], inp["Dtrain"]["HPY"]
FB, FP, FH = F(inp["Xtrain"]["BLD"]), F(inp["Xtrain"]["PLE"]), F(inp["Xtrain"]["HPY"])
groups["June BLD training (all 3 dates)"] = FB
for d in ("20240613", "20240616", "20240620"):
    groups[f"BLD training {d}"] = FB[Db == d]
groups["PLE training 20230509"] = FP[Dp == "20230509"]
groups["PLE training 20240620"] = FP[Dp == "20240620"]
groups["HPY training (all dates)"] = FH
res["feature_space"] = {k: {"axis_t": q(v[:, it]), "axis_residual": q(v[:, ir])} for k, v in groups.items()}
print("e) group: n | axis_t median [IQR] | residual median [IQR]")
for k, v in res["feature_space"].items():
    t, rr = v["axis_t"], v["axis_residual"]
    print(f"   {k:34s} {t['n']:6d} | {t['median']:.3f} [{t['q25']:.3f}, {t['q75']:.3f}] | "
          f"{rr['median']:.4f} [{rr['q25']:.4f}, {rr['q75']:.4f}]")

json.dump(res, open(OUT, "w"), indent=1, default=float)
print("written", OUT)
