"""Summary tables and context statistics of the R1 manuscript, built from a completed revision run.

Reads the outputs of r1/revision_runner.py (all policies, seeds, omissions and LOO; folder REVISION_RUN)
and writes CSV/JSON summaries to TABLES_OUT (default r1/output/tables). No model is refitted; forests
are only reloaded. Ported from the post-processing scripts of the original revision workspace; the
computations are unchanged.

Environment (required unless marked):
  CORALREMOTE_ROOT, CORALREMOTE_A3_ROOT, REVISION_UDM   as for revision_runner.py
  REVISION_RUN       output folder of the completed revision run
  CORALREMOTE_ACA_BATHY  Allen Coral Atlas bathymetry GeoTIFF (centimetres), for the mask depth summary
  TABLES_OUT         (optional) output folder
  REFERENCE_REPORT   (optional) folder of previously written tables; if set, every output is compared
                     with the file of the same name there (numbers to 1e-9) and the result is printed

Run from the repository root:  python r1/make_tables.py

Outputs and the manuscript items they support:
  candidate_Table2.csv, OOB_and_drawn_training_counts.csv     Table 2 (training draw, OOB errors)
  reference_polygons.json                                      Table 2 caption (polygon counts)
  candidate_Table4A.csv, stage1_accuracy.csv                   Table 4A
  Stage1_transfer_metrics.csv                                  Section 3.2.1 (0.90 on 20240307, 0.035 on 20240317)
  candidate_Table4B.csv, stage2_evaluation_and_training_counts.csv   Table 4B; S4 common (8,854) and
                                                               location-disjoint (8,320) subsets
  candidate_Table5.csv, paired_stage2_dates.csv                Table 5; S3 seed and policy ranges
  paired_masks.csv                                             S3 mask counts, IoU and area change
  calm_date_common_support.csv                                 S1 common-pixel table
  June_omission_results.csv                                    S4 date omission
  LOO_polygons.csv                                             Section 3.2.1 leave-one-polygon-out recall
  DHW_correlations.csv                                         Section 3.5.4 and S5
  context_statistics.json                                      mask size and depth (Section 3.5.6), haze on
                                                               20240616, pale polygons, anchor checks,
                                                               May 2025 axis domain (S5), uncurated pale
                                                               candidates (S4), NIR-threshold shares (S4)
"""
import csv
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _need(name):
    v = os.environ.get(name)
    if not v:
        raise RuntimeError(f'{name} is not set (see the docstring of r1/make_tables.py)')
    return v


RUN = Path(_need('REVISION_RUN'))
os.environ['REVISION_OUT'] = str(RUN)            # the runner resolves its output folder from this
sys.argv = sys.argv[:1]                           # the runner reads --config from argv; none is used here
sys.path.insert(0, str(HERE))
import revision_runner as r                       # noqa: E402  (puts the repository root on sys.path)
import numpy as np                                # noqa: E402
from scipy.stats import pearsonr                  # noqa: E402
from coralremote.train import build_training_sets  # noqa: E402

P = r.OUT
OUT = r.guard_out(Path(os.environ.get('TABLES_OUT') or HERE / 'output' / 'tables'))
OUT.mkdir(parents=True, exist_ok=True)
DHW = [0, .91, 2.46, 2.46, 0, 0, .33, 7.57, 9.71, 9.71, 9.71, 0, 0, 1.09]   # Table 1, in stack date order


def Cmat(y, p):
    C = np.zeros((4, 4), np.int64)
    np.add.at(C, (y.astype(int) - 1, p.astype(int) - 1), 1)
    return C


def row(j, d):
    return next(x for x in j['table'] if x['date'] == d)


def stage2_view():
    return r.stage2.StackView(r.io.load_glint_stack(), 'fixed')


# ------------------------------------------------------------------ summaries (no refitting)
def summarize():
    cal = r.load(P / 'calibration/calibration.json')
    primary = r.load(P / 'manifest/selection.json')['selected_candidate']
    prep = r.load(P / 'stage2/prepared.json')
    common, clean = set(prep['common_keys']), set(prep['location_disjoint_common_keys'])
    models, date_rows, metricrows = {}, [], []
    for policy in ('pooled', 'per_date'):
        for seed in r.SEEDS:
            tag = f'{policy}_seed{seed}'
            m = r.load(P / f'stage2/{tag}.json'); models[tag] = m
            z = np.load(P / f'stage2/{tag}_predictions.npz'); inp = r.undump(P / f'stage2/{tag}_inputs.pkl')
            er = dict(policy=policy, seed=seed, n_scenario=len(z['test_y']), OA_core=m['metrics']['OA'],
                      accuracy_with_UNK=m['accuracy_with_UNK'], PLE_F1=m['metrics']['f1'][1])
            for label, keys in [('common', common), ('location_disjoint', clean)]:
                ix = np.array([k in keys for k in z['test_keys']]); assert int(ix.sum()) == len(keys)
                C = Cmat(z['test_y'][ix], z['test_pred'][ix]); met = r.compute_metrics3(z['test_y'][ix], z['test_pred'][ix])
                er[label + '_n'] = len(keys); er[label + '_accuracy'] = float(np.trace(C) / C.sum()); er[label + '_PLE_F1'] = met['f1'][1]
            ts = build_training_sets({**inp, 'view': stage2_view()}, mode_rows='numpy', seed=seed)
            er.update(gate_n=len(ts['ygate']), stress_n=len(ts['ystress']), gate_HPY_n=int(np.sum(ts['ygate'] == 1)),
                      gate_stressed_n=int(np.sum(ts['ygate'] == 2)), stress_PLE_n=int(np.sum(ts['ystress'] == 1)),
                      stress_BLD_n=int(np.sum(ts['ystress'] == 2)))
            for k, v in ts.items():
                if isinstance(v, np.ndarray):
                    assert r.ah(v) == m['feature_hashes'][k], f'training feature hash differs: {tag} {k}'
            metricrows.append(er)
            for d in m['table']:
                valid = z[d['date']] > 0; counts = np.bincount(z[d['date']][valid], minlength=5)[1:]
                base = row(models[f'pooled_seed{seed}'], d['date'])
                date_rows.append(dict(policy=policy, seed=seed, date=d['date'], n=int(valid.sum()), **{k: d[k] for k in r.NAMES},
                                      **{k + '_count': int(counts[i]) for i, k in enumerate(r.NAMES)},
                                      **{k + '_delta_pooled_pp': d[k] - base[k] for k in r.NAMES}))
    r.csvsave(OUT / 'paired_stage2_dates.csv', date_rows); r.csvsave(OUT / 'stage2_evaluation_and_training_counts.csv', metricrows)
    stage1, maskrows = [], []
    for policy in ('pooled', 'per_date'):
        for seed in r.SEEDS:
            z = r.load(P / f'stage1/{policy}_accuracy_seed{seed}.json')
            stage1.append(dict(policy=policy, seed=seed, OA=z['metrics']['OA'], macro_F1=z['macro_f1'], n=z['n_test'],
                               **{f'{cl}_{metric}': z['metrics'][metric][i] for i, cl in enumerate(('CORAL', 'SND', 'RBL')) for metric in ('prec', 'rec', 'f1')}))
            a = np.load(P / f'stage1/pooled_mask_seed{seed}_masks.npz')['maskFinal']; b = np.load(P / f'stage1/{policy}_mask_seed{seed}_masks.npz')['maskFinal']
            inter = int((a & b).sum()); union = int((a | b).sum())
            maskrows.append(dict(policy=policy, seed=seed, n=int(b.sum()), pooled_n=int(a.sum()), area_m2=int(b.sum()) * 9, iou=inter / union,
                                 lost=int((a & ~b).sum()), gained=int((b & ~a).sum()), area_delta_pct=100 * (b.sum() - a.sum()) / a.sum()))
    r.csvsave(OUT / 'stage1_accuracy.csv', stage1); r.csvsave(OUT / 'paired_masks.csv', maskrows)
    comparison = []
    for d in r.BASE:
        a = cal['per_date'][d]; b = cal['common_support'][d]
        comparison.append(dict(date=d, own_support_n=a['n'], common_support_n=b['n'], own_BR=a['pair'][0], own_GR=a['pair'][1],
                               common_BR=b['pair'][0], common_GR=b['pair'][1], own_R2_BR=a['R2'][0], own_R2_GR=a['R2'][1],
                               common_R2_BR=b['R2'][0], common_R2_GR=b['R2'][1]))
    r.csvsave(OUT / 'calm_date_common_support.csv', comparison)
    omissions = []
    for held in r.JUNE:
        for seed in r.SEEDS:
            tag = f'{primary}_omit{held}_seed{seed}'; m = r.load(P / f'omissions/{tag}.json'); base = models[f'{primary}_seed{seed}']
            d = row(m, held); q = row(base, held); z = np.load(P / f'omissions/{tag}_predictions.npz')
            ix = np.array([k.startswith(held + ':') for k in z['test_keys']]) & (z['test_y'] == 3)
            omissions.append(dict(held_date=held, seed=seed, held_BLD=d['BLD'], held_PLE=d['PLE'], held_HPY=d['HPY'], held_UNK=d['UNK'],
                                  held_stressed=d['PLE'] + d['BLD'], full_BLD=q['BLD'], delta_BLD_pp=d['BLD'] - q['BLD'],
                                  held_BLD_label_n=int(ix.sum()), held_BLD_label_recall=float(np.mean(z['test_pred'][ix] == 3)),
                                  n_test=m['n_test'], OA_core=m['metrics']['OA'], accuracy_with_UNK=m['accuracy_with_UNK']))
    r.csvsave(OUT / 'June_omission_results.csv', omissions)
    loofits = [r.load(f) for f in sorted((P / 'loo').glob('polygon*_seed*.json')) if not f.name.endswith('.settings.json')]
    r.csvsave(OUT / 'LOO_polygons.csv', loofits)
    dates = list(r.io.load_glint_stack().dates); dhw = dict(zip(dates, DHW))
    ten = [d for d in dates if d not in ('20240317', '20240408', '20241111', '20250518')]; seven = [d for d in ten if d not in r.JUNE]
    allmodels = {**models, **{f'{primary}_omit{d}_seed{s}': r.load(P / f'omissions/{primary}_omit{d}_seed{s}.json') for d in r.JUNE for s in r.SEEDS}}
    correlations = []
    for tag, m in allmodels.items():
        for label, dd in [('interpreted10', ten), ('all14', dates), ('interpreted_without_June7', seven)]:
            correlations.append(dict(run=tag, date_set=label, n=len(dd), pearson_r=float(pearsonr([dhw[d] for d in dd], [row(m, d)['BLD'] for d in dd]).statistic)))
    r.csvsave(OUT / 'DHW_correlations.csv', correlations)
    producer = models[f'{primary}_seed0']; table5 = []
    for d in dates:
        q = row(producer, d); rs = [row(models[f'{primary}_seed{s}'], d) for s in r.SEEDS]
        zz = next(x for x in date_rows if x['policy'] == primary and x['seed'] == 0 and x['date'] == d)
        table5.append(dict(date=d, producer=f'{primary}_seed0', n=zz['n'], **{k: round(q[k], 2) for k in r.NAMES},
                           **{k + '_count': zz[k + '_count'] for k in r.NAMES},
                           **{k + '_replicate_mean': float(np.mean([z[k] for z in rs])) for k in r.NAMES},
                           **{k + '_replicate_min': float(min(z[k] for z in rs)) for k in r.NAMES},
                           **{k + '_replicate_max': float(max(z[k] for z in rs)) for k in r.NAMES}))
    r.csvsave(OUT / 'candidate_Table5.csv', table5)
    return primary


# ------------------------------------------------------------------ OOB, Table 2 draw, Stage 1 transfer
def oob(f, X, y, seed):
    """Out-of-bag error from the exact NumPy bootstrap stream of a stored forest (no refitting)."""
    rng = np.random.default_rng(seed); bags = np.zeros((len(y), len(f.trees)), bool)
    for i in range(len(f.trees)):
        idx = rng.integers(0, len(y), len(y)); bags[:, i] = np.bincount(idx, minlength=len(y)) == 0
    assert bags.any(axis=1).all()
    return r.stage1.oob_error(f, X, y, bags)


def oob_and_transfer(policy):
    z = np.load(P / 'stage1/inputs.npz'); tr = r.zrows(z, 'accuracy'); X = r.apply_rows(tr, policy); rows = []
    for seed in r.SEEDS:
        f = r.read_forest(P / f'stage1/{policy}_accuracy_seed{seed}_forest.npz'); err = oob(f, X, tr['Y'], seed)
        inp = r.undump(P / f'stage2/{policy}_seed{seed}_inputs.pkl'); inp['view'] = stage2_view()
        ts = build_training_sets(inp, mode_rows='numpy', seed=seed)
        g = r.read_forest(P / f'stage2/{policy}_seed{seed}_gate.npz'); s = r.read_forest(P / f'stage2/{policy}_seed{seed}_stress.npz')
        ge = oob(g, ts['Fgate'], ts['ygate'], seed); se = oob(s, ts['Fstress'], ts['ystress'], seed + 1000)
        rng = np.random.default_rng(seed); nh = len(inp['Xtrain']['HPY']); ns = len(inp['Xtrain']['PLE']) + len(inp['Xtrain']['BLD'])
        if nh > 3000:
            rng.choice(nh, 3000, replace=False)          # advances the stream exactly as training did
        ix = rng.choice(ns, 3000, replace=False) if ns > 3000 else np.arange(ns)
        npale = len(inp['Xtrain']['PLE']); classpool = np.r_[np.full(npale, 2), np.full(len(inp['Xtrain']['BLD']), 3)]
        raw = np.concatenate([inp['Xtrain']['PLE'], inp['Xtrain']['BLD']])
        F = r.fe.build_hybrid_features(r.fe.compute_feature_base(raw, 'fixed'), inp['mu_h'], inp['mu_b'], inp['v'], inp['vnorm2'])
        assert np.array_equal(F[ix], ts['Fgate'][min(3000, nh):])
        rows.append(dict(seed=seed, Stage1_OOB_error=err, Stage2_gate_OOB_error=ge, Stage2_stress_OOB_error=se, gate_HPY_n=min(3000, nh),
                         gate_PLE_n=int(np.sum(classpool[ix] == 2)), gate_BLD_n=int(np.sum(classpool[ix] == 3)),
                         stress_PLE_n=npale, stress_BLD_n=len(inp['Xtrain']['BLD'])))
    r.csvsave(OUT / 'OOB_and_drawn_training_counts.csv', rows)
    draw = rows[0]
    r.csvsave(OUT / 'candidate_Table2.csv', [
        dict(model='gate', HPY=int(draw['gate_HPY_n']), PLE=int(draw['gate_PLE_n']), BLD=int(draw['gate_BLD_n']), n=6000),
        dict(model='stress split', HPY=0, PLE=int(draw['stress_PLE_n']), BLD=int(draw['stress_BLD_n']),
             n=int(draw['stress_PLE_n']) + int(draw['stress_BLD_n']))])
    ref = r.io.load_glint_stack().ref
    te, lab = r.rasterize_labels_stage1(r.io.load_geojson(str(r.paths.STAGE1_LABELS['test'])), ref, r.stage1.POLY_CLASS_MAP)
    transfer = []
    for d in r.stage1.TRANSFER_TEST_DATES:
        A = r.date_features(d, policy); rr_, cc_ = r.find_colmajor(te & np.isfinite(A).all(2) & (lab > 0)); Y = lab[rr_, cc_]
        for seed in r.SEEDS:
            f = r.read_forest(P / f'stage1/{policy}_accuracy_seed{seed}_forest.npz'); pred = f.predict(A[rr_, cc_])[0]; met = r.compute_metrics3(Y, pred)
            transfer.append(dict(date=d, seed=seed, n=len(Y), OA=met['OA'], CORAL_precision=met['prec'][0], CORAL_recall=met['rec'][0],
                                 CORAL_F1=met['f1'][0], date_in_CORAL_training=d in r.stage1.CORAL_DATES,
                                 date_in_SNDRBL_training=d in r.stage1.SND_RBL_DATES))
    r.csvsave(OUT / 'Stage1_transfer_metrics.csv', transfer)


def tables4(policy):
    producer = r.load(P / f'stage2/{policy}_seed0.json'); s1 = r.load(P / f'stage1/{policy}_accuracy_seed0.json')
    r.csvsave(OUT / 'candidate_Table4A.csv', [dict(class_name=cl, precision=s1['metrics']['prec'][i], recall=s1['metrics']['rec'][i],
                                                   F1=s1['metrics']['f1'][i], n=int(s1['metrics']['total'][i])) for i, cl in enumerate(('CORAL', 'SND', 'RBL'))])
    r.csvsave(OUT / 'candidate_Table4B.csv', [dict(class_name=cl, precision=producer['metrics']['prec'][i], recall=producer['metrics']['rec'][i],
                                                   F1=producer['metrics']['f1'][i], n=int(producer['metrics']['total'][i])) for i, cl in enumerate(('HPY', 'PLE', 'BLD'))])


# ------------------------------------------------------------------ context statistics
def context(policy):
    import rasterio
    from rasterio.warp import reproject, Resampling
    x = r.undump(P / f'stage2/{policy}_seed0_inputs.pkl'); support = x['support']; reef = r.pp.reef_mask(); ref = r.io.load_glint_stack().ref
    with rasterio.open(_need('CORALREMOTE_ACA_BATHY')) as src:
        dest = np.full(support.shape, np.nan, float)
        reproject(src.read(1).astype(float), dest, src_transform=src.transform, src_crs=src.crs, src_nodata=src.nodata,
                  dst_transform=rasterio.Affine(*ref.transform), dst_crs=f'EPSG:{ref.epsg}', dst_nodata=np.nan, resampling=Resampling.bilinear)
    depth = dest / 100

    def q(a):
        return dict(zip(('q25', 'median', 'q75'), np.nanpercentile(a, [25, 50, 75]).tolist()))
    udm = np.load(r.UDM / 'udm2_20240616.npz')['udm2']; haze = float(np.mean(udm[3][support] == 1) * 100)
    audit = r.load(P / 'stage2/sample_audit.json')[f'{policy}_seed0']; ts = audit['split_stats']
    stats = dict(mask_n=int(support.sum()), mask_area_m2=int(support.sum()) * 9, reef_fraction_percent=float(support.sum() / reef.sum() * 100),
                 depth_mask=q(depth[support]), depth_other_reef=q(depth[reef & ~support]), depth_reef=q(depth[reef]),
                 haze_20240616_percent=haze, PLE_curated_polygons=audit['ple_stats']['nKeptPolygons'],
                 PLE_train_polygons=ts['nTrainPolygons'], PLE_test_polygons=ts['nTestPolygons'])
    g = r.read_forest(P / f'stage2/{policy}_seed0_gate.npz'); s = r.read_forest(P / f'stage2/{policy}_seed0_stress.npz')
    for key, label in [('healthy', 'anchor_healthy_recall'), ('bleached', 'anchor_bleached_recall')]:
        X, Y = x['anchors'][key]
        F = r.fe.build_hybrid_features(r.fe.compute_feature_base(X, 'fixed'), x['mu_h'], x['mu_b'], x['v'], x['vnorm2'])
        pred = r.stage2.hierarchical_predict(F, g, s, x['opts'])[0]; stats[label] = float(np.mean(pred == Y)); stats[label + '_n'] = len(Y)
    res = []
    for seed in r.SEEDS:
        z = r.undump(P / f'stage2/{policy}_seed{seed}_inputs.pkl'); tr = np.concatenate(list(z['Xtrain'].values()))
        F = r.fe.compute_feature_base(tr, 'fixed'); t, _ = r.fe.project_to_axis(F[:, :8], z['mu_h'], z['v'], z['vnorm2'])
        A = stage2_view().get('20250518'); rows, cols = np.where(z['support'] & np.isfinite(A).all(2))
        G = r.fe.compute_feature_base(A[rows, cols], 'fixed'); u, _ = r.fe.project_to_axis(G[:, :8], z['mu_h'], z['v'], z['vnorm2'])
        res.append(dict(seed=seed, axis_training_min=float(t.min()), axis_training_max=float(t.max()),
                        outside_axis_training_range_percent=float(100 * np.mean((u < t.min()) | (u > t.max())))))
    stats['May2025_axis_domain'] = res
    allp = r.stage2.filter_features(r.io.load_labels('unified_train'), {'PLE'}, True)
    unused = [f for f, k in zip(allp, audit['ple_stats']['keepMask']) if not k]
    maps = np.load(P / f'stage2/{policy}_seed0_predictions.npz'); rows = []
    for d in r.io.load_glint_stack().dates:
        mm, _ = r.stage2.rasterize_date_labels(unused, ref, {'PLE': 2}, d); pr = maps[d][mm & support & (maps[d] > 0)]
        if len(pr):
            rows.append(dict(date=d, n=len(pr), **{name: float(100 * np.mean(pr == i)) for i, name in enumerate(r.NAMES, 1)}))
    stats['uncurated_PLE_candidate_predictions'] = rows
    wave = []
    for d in r.io.load_glint_stack().dates:
        _, w = r.pp.valid_and_wave(d); wave.append(dict(date=d, removed_percent=float(100 * np.sum(w & reef) / reef.sum())))
    stats['Stage1_wave_removed'] = wave
    r.save(OUT / 'context_statistics.json', stats)
    u_train, u_test = r.io.load_labels('unified_train'), r.io.load_labels('unified_test')
    s_train, s_test = r.io.load_labels('strict_train'), r.io.load_labels('strict_test')
    ff = r.stage2.filter_features
    r.save(OUT / 'reference_polygons.json', dict(
        HPY_train=len(ff(u_train, {'HPY'}, True)), HPY_test=len(ff(u_test, {'HPY'}, False)),
        PLE_train=ts['nTrainPolygons'], PLE_test=ts['nTestPolygons'],
        BLD_train=len(ff(s_train, {'BLD'}, True)), BLD_test=len(ff(s_test, {'BLD'}, False))))


# ------------------------------------------------------------------ comparison with stored tables
def _cells_equal(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b)))
    except (TypeError, ValueError):
        return False


def _same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return all(k in b and _same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return _cells_equal(a, b)


def compare(ref_dir):
    ref_dir = Path(ref_dir); report = []
    for f in sorted(OUT.iterdir()):
        g = ref_dir / f.name
        if f.name == 'context_statistics.json':
            g = ref_dir.parent / 'document_candidates' / 'candidate_context_statistics.json'
        if not g.exists():
            report.append((f.name, 'no reference')); continue
        if f.suffix == '.csv':
            a = list(csv.DictReader(f.open(encoding='utf-8'))); b = list(csv.DictReader(g.open(encoding='utf-8')))
            key = lambda x: json.dumps(x, sort_keys=True)
            if 'polygon_id' in (a[0] if a else {}):
                a.sort(key=key); b.sort(key=key)
            ok = len(a) == len(b) and all(set(x) <= set(y) and all(_cells_equal(x[k], y[k]) for k in x) for x, y in zip(a, b))
        else:
            ok = _same(json.loads(f.read_text(encoding='utf-8')), json.loads(g.read_text(encoding='utf-8')))
        report.append((f.name, 'match' if ok else 'DIFFERS'))
    for name, status in report:
        print(f'{status:12s} {name}')
    return all(s in ('match', 'no reference') for _, s in report)


if __name__ == '__main__':
    policy = summarize()
    oob_and_transfer(policy)
    tables4(policy)
    context(policy)
    r.log(f'Tables written to {OUT}')
    if os.environ.get('REFERENCE_REPORT'):
        sys.exit(0 if compare(os.environ['REFERENCE_REPORT']) else 1)
