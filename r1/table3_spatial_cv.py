"""Table 3 spatial cross-validation column and the earlier exploratory test column (Supplement S2).

Rebuilds the controlled Stage 1 feature-family ablation that supplies the spatial-CV macro-F1 of Table 3
(0.741, 0.867, 0.860, 0.871 for R, R+B, R+G, R+BG) and the earlier exploratory test OA that Supplement S2
reports as replaced (70.0, 96.9, 95.0, 96.7%). Ported from the exploratory ablation of the revision
workspace; the design and computations are unchanged:
  - training rows: Stage 1 training polygons, coral on five dates and sand/rubble on seven (20,596 rows);
    earlier exploratory test: coral 20241111, sand/rubble 20230519 (905 rows);
  - three StratifiedGroupKFold folds on 32-pixel tiles (shuffle, random_state 42), all dates of a tile kept
    together;
  - per fold, the pooled StableSand OLS pair is re-estimated with that fold's validation tiles excluded
    (no calibration sand falls in any validation tile, so it equals the pooled pair);
  - 200-tree MATLAB-rule forests, equal costs, seeds 0 to 2; macro-F1 is the mean of the three class F1.
The output layout is the one r1/table3_retest.py reads (set its REVISION_OLD_RUN to this TABLE3_CV_OUT).

Environment (required unless marked):
  CORALREMOTE_ROOT, CORALREMOTE_A3_ROOT, REVISION_UDM, CORALREMOTE_STABLESAND   as for revision_runner.py
  TABLE3_CV_OUT   (optional) output folder, default r1/output/table3_cv
  CV_REFERENCE    (optional) folder of a previous ablation run; every fit's input hashes and metrics are
                  compared with the file of the same name there
  CV_WORKERS      (optional) parallel processes, default 3

Run from the repository root (about 10 to 15 minutes with three workers):  python r1/table3_spatial_cv.py
"""
import concurrent.futures as cf
import csv
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.argv = sys.argv[:1]
sys.path.insert(0, str(HERE))
import revision_runner as r                       # noqa: E402  (repository root on sys.path; UDM2 and roots)
import numpy as np                                # noqa: E402
from coralremote import io, paths, stage1, preprocess as pp   # noqa: E402
from coralremote.cart import grow_forest          # noqa: E402
from coralremote.evaluate import compute_metrics3  # noqa: E402
from coralremote.matlab_compat import find_colmajor  # noqa: E402
from coralremote.rasterize import rasterize_labels_stage1  # noqa: E402

OUT = r.guard_out(Path(os.environ.get('TABLE3_CV_OUT') or HERE / 'output' / 'table3_cv'))
ABL = OUT / 'stage1_ablation'
SEEDS = (0, 1, 2)
FEATURE_SETS = {'R': [0, 1, 2], 'R+B': [0, 1, 2, 3], 'R+G': [0, 1, 2, 4], 'R+BG': [0, 1, 2, 3, 4]}
ORIGINAL = (0.82, 0.80)                            # coefficients built into the stored Phase 5 contrasts
# Dates of the exploratory ablation design (as recorded with its saved settings).
CORAL_DATES = ('20230404', '20230509', '20240307', '20240317', '20240408')
SNDRBL_DATES = ('20230404', '20230509', '20230519', '20230629', '20240307', '20240317', '20240408')
TEST_PLAN = (('20241111', (1,)), ('20230519', (2, 3)))


def calibration(exclude_groups=()):
    """Pooled StableSand OLS pair on the three calm dates, with the given 32-pixel tiles excluded."""
    from rasterio.features import rasterize
    from rasterio.warp import transform_geom
    from affine import Affine
    ref = io.load_glint_stack().ref
    rr_, cc_ = np.indices((ref.height, ref.width)); groups = (rr_ // 32) * 10000 + cc_ // 32
    geo = json.loads(r.stablesand_path().read_text())
    m = rasterize([(transform_geom('EPSG:4326', 'EPSG:32647', f['geometry']), 1) for f in geo['features']],
                  out_shape=(ref.height, ref.width), transform=Affine(*ref.transform), dtype='uint8') > 0
    m &= ~np.isin(groups, exclude_groups)
    xx, ys = [], [[], []]
    for d in r.BASE:
        A = pp.glint_reflectance(d); nir = pp.raw_stack().get(d)[:, :, 7]
        ok = m & np.isfinite(A[:, :, [1, 2, 5]]).all(2) & (nir > .005) & (nir < .12)
        v = [np.log(np.maximum(A[:, :, i][ok].astype(float), 1e-4)) for i in (5, 1, 2)]
        xx.append(v[0] - v[0].mean())
        for j, y0 in enumerate(v[1:]):
            ys[j].append(y0 - y0.mean())
    x = np.concatenate(xx)
    return [float(x @ np.concatenate(yl) / (x @ x)) for yl in ys]


def inputs():
    dest = ABL / 'inputs.npz'
    if dest.exists():
        return
    from sklearn.model_selection import StratifiedGroupKFold
    ref = io.load_glint_stack().ref
    tr, tl = rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['train'])), ref, stage1.POLY_CLASS_MAP)
    te, el = rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['test'])), ref, stage1.POLY_CLASS_MAP)
    assert not np.any(tr & te), 'Stage 1 training/test polygons overlap'

    def collect(training):
        X, Y, L, K, G = [], [], [], [], []
        plan = ([(d, (1,)) for d in CORAL_DATES] + [(d, (2, 3)) for d in SNDRBL_DATES]) if training else TEST_PLAN
        for d, classes in plan:
            A = pp.phase5_features(d); lab = tl if training else el; ma = tr if training else te
            rows, cols = find_colmajor(ma & np.isin(lab, classes) & np.isfinite(A).all(2))
            X.append(A[rows, cols]); Y.append(lab[rows, cols])
            L.append(pp._log32(np.maximum(pp.glint_reflectance(d)[rows, cols, 5], np.float32(1e-6))))
            G.append((rows // 32) * 10000 + cols // 32); K.extend(f'{d}:{a}:{b}:{int(c)}' for a, b, c in zip(rows, cols, lab[rows, cols]))
        return np.concatenate(X), np.concatenate(Y), np.concatenate(L), np.array(K), np.concatenate(G)
    X, Y, L, K, G = collect(True); Xt, Yt, Lt, Kt, Gt = collect(False)
    folds = np.full(len(Y), -1, int)
    for i, (_, va) in enumerate(StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=42).split(X, Y, G)):
        folds[va] = i
    ks, foldinfo = {}, {}
    for fold in range(3):
        va = folds == fold
        assert set(Y[va]) == {1, 2, 3} and set(Y[~va]) == {1, 2, 3} and not set(G[va]) & set(G[~va])
        ks[str(fold)] = {'original': list(ORIGINAL), 'stable_OLS': calibration(np.unique(G[va]))}
        foldinfo[str(fold)] = dict(train_n=int((~va).sum()), validation_n=int(va.sum()), coefficients=ks[str(fold)])
    ks['historical'] = {'original': list(ORIGINAL), 'stable_OLS': calibration()}
    ABL.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, X=X, Y=Y, L=L, K=K, G=G, Xtest=Xt, Ytest=Yt, Ltest=Lt, Ktest=Kt, Gtest=Gt, fold=folds)
    r.save(ABL / 'design.json', dict(training_n=len(Y), historical_test_n=len(Yt), training_classes=np.bincount(Y, minlength=4)[1:],
                                     folds=foldinfo, coefficients=ks, feature_sets=FEATURE_SETS, n_trees=200, min_leaf=5, min_parent=10,
                                     cost=None, training_keys_sha=r.ah(K), historical_test_keys_sha=r.ah(Kt)))
    r.log(f'Ablation design fixed: {len(Y)} training and {len(Yt)} exploratory test rows; 3 spatial folds')


def changed(X, L, k):
    """Replace the stored 0.82/0.80 contrasts by those of coefficient pair k."""
    A = X.copy(); A[:, 3] += np.float32(ORIGINAL[0] - k[0]) * L; A[:, 4] += np.float32(ORIGINAL[1] - k[1]) * L
    return A


def job(task):
    scenario, feature, fold, seed = task
    tag = f'{scenario}_{feature.replace("+", "_")}_{fold}_seed{seed}'; dest = ABL / f'{tag}.json'
    if dest.exists():
        return tag + ' cached'
    t = time.time(); z = np.load(ABL / 'inputs.npz'); design = r.load(ABL / 'design.json')
    k = design['coefficients'][str(fold)][scenario]; cols = FEATURE_SETS[feature]
    if fold == 'historical':
        tr = np.ones(len(z['Y']), bool); Xv, Yv, Lv = z['Xtest'], z['Ytest'], z['Ltest']
    else:
        tr = z['fold'] != int(fold); Xv, Yv, Lv = z['X'][~tr], z['Y'][~tr], z['L'][~tr]
    X = changed(z['X'][tr], z['L'][tr], k)[:, cols]; V = changed(Xv, Lv, k)[:, cols]
    f = grow_forest(X, z['Y'][tr], n_trees=200, seed=seed, classes=(1, 2, 3)); p, _ = f.predict(V)
    metrics = compute_metrics3(Yv, p)
    r.save(dest, dict(scenario=scenario, feature_set=feature, fold=fold, seed=seed, coefficients=k, seconds=time.time() - t,
                      n_train=int(tr.sum()), n_test=len(p), metrics=metrics, macro_f1=float(np.mean(metrics['f1'][:3])),
                      train_input_sha=r.ah(X), test_input_sha=r.ah(V)))
    return f'{tag}: {time.time() - t:.1f}s OA={metrics["OA"]:.4f}'


def summary():
    rows = []
    for feature in FEATURE_SETS:
        scen = 'original' if feature == 'R' else 'stable_OLS'      # R does not use the coefficients
        cv = [r.load(ABL / f'{scen}_{feature.replace("+", "_")}_{fold}_seed{s}.json') for fold in '012' for s in SEEDS]
        hist = [r.load(ABL / f'{scen}_{feature.replace("+", "_")}_historical_seed{s}.json') for s in SEEDS]
        rows.append(dict(scenario='stable_OLS', feature_set=feature,
                         fold_seed_mean_macro_f1=float(np.mean([x['macro_f1'] for x in cv])),
                         fold_seed_min_macro_f1=float(min(x['macro_f1'] for x in cv)),
                         fold_seed_max_macro_f1=float(max(x['macro_f1'] for x in cv)),
                         mean_historical_OA=float(np.mean([x['metrics']['OA'] for x in hist])),
                         mean_historical_macro_f1=float(np.mean([x['macro_f1'] for x in hist]))))
    (OUT / 'report').mkdir(parents=True, exist_ok=True)
    r.csvsave(OUT / 'report' / 'ablation_summary.csv', rows)
    for x in rows:
        print(f"{x['feature_set']:5s} spatial-CV macro-F1 {x['fold_seed_mean_macro_f1']:.3f}   exploratory test OA {100 * x['mean_historical_OA']:.1f}%")


def compare(ref):
    ref = Path(ref); bad = []
    for f in sorted(ABL.glob('*_seed*.json')):
        g = ref / 'stage1_ablation' / f.name
        if not g.exists():
            bad.append(f'{f.name}: no reference'); continue
        a, b = r.load(f), r.load(g)
        same = a['train_input_sha'] == b['train_input_sha'] and a['test_input_sha'] == b['test_input_sha'] and a['metrics']['C'] == b['metrics']['C']
        if not same:
            bad.append(f'{f.name}: differs')
    old = {x['feature_set']: x for x in csv.DictReader((ref / 'report' / 'ablation_summary.csv').open(encoding='utf-8')) if x['scenario'] == 'stable_OLS'}
    for x in csv.DictReader((OUT / 'report' / 'ablation_summary.csv').open(encoding='utf-8')):
        for k in ('fold_seed_mean_macro_f1', 'mean_historical_OA'):
            if abs(float(x[k]) - float(old[x['feature_set']][k])) > 1e-12:
                bad.append(f"summary {x['feature_set']} {k} differs")
    print('reference comparison:', 'all fits and summaries identical' if not bad else bad)
    return not bad


if __name__ == '__main__':
    inputs()
    tasks = [(('original' if f == 'R' else 'stable_OLS'), f, fold, s) for fold in ('historical', '0', '1', '2') for s in SEEDS for f in FEATURE_SETS]
    with cf.ProcessPoolExecutor(max_workers=int(os.environ.get('CV_WORKERS', '3'))) as pool:
        for res in pool.map(job, tasks):
            r.log(res)
    summary()
    if os.environ.get('CV_REFERENCE'):
        sys.exit(0 if compare(os.environ['CV_REFERENCE']) else 1)
