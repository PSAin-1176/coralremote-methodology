"""Table 3 feature-family ablation re-run on the Table 4A held-out design.

Reuses, without modification, the Stage 1 accuracy-model design of the
revision run (the output folder of r1/revision_runner.py):
  * training rows and 905-pixel test rows from stage1/inputs.npz
    (accuracy_* and test_*), built by revision_runner.collect();
  * pooled DII coefficients from calibration/calibration.json, applied by
    revision_runner.apply_rows(..., 'pooled');
  * forest = coralremote.cart.grow_forest(n_trees=200, min_leaf=5,
    min_parent=10, classes=(1,2,3), cost=stage1.RF_COST, seed in 0,1,2);
  * metrics = coralremote.evaluate.compute_metrics3.
Only the predictor columns change between feature families.

Sanity gate: the five-feature family must reproduce the saved Table 4A
predictions (stage1/pooled_accuracy_seed{0,1,2}_evaluation.npz) exactly and
the metrics in report/stage1_accuracy.csv. If not, no results are written.

Writes only metrics, counts, configs and hashes (no imagery, no reflectance
or feature values) to TABLE3_OUT (default r1/output/table3_retest/).

Locations come from environment variables (all required unless marked):
  CORALREMOTE_ROOT  data root (also needed by the runner and the package)
  REVISION_RUN      output folder of the completed revision run (stage1/, calibration/, report/, manifest/)
  REVISION_OLD_RUN  output folder of the earlier exploratory ablation run (report/ablation_summary.csv,
                    stage1_ablation/), read for the unchanged spatial-CV and exploratory OA columns
  TABLE3_OUT        (optional) output folder
Run from the repository root:  python r1/table3_retest.py
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
def _need(name):
    v = os.environ.get(name)
    if not v:
        raise RuntimeError(f'{name} is not set (see the docstring of r1/table3_retest.py)')
    return v
import sys
sys.dont_write_bytecode = True            # never write __pycache__ into read-only trees
import csv, json, hashlib, time, subprocess
import concurrent.futures as cf
from datetime import datetime
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
RUNNER_DIR = SCRIPT.parent                        # r1/
PORT = RUNNER_DIR.parent                          # repository root (holds coralremote/)
_need('CORALREMOTE_ROOT')
RUN = Path(_need('REVISION_RUN'))
OLD = Path(_need('REVISION_OLD_RUN'))
OUT = Path(os.environ.get('TABLE3_OUT') or RUNNER_DIR / 'output' / 'table3_retest')
os.environ['REVISION_OUT'] = str(RUN)             # the runner reads calibration/ from its output folder

sys.path.insert(0, str(RUNNER_DIR))
import revision_runner as rr                      # noqa: E402  (inserts the repository root on sys.path)
import coralremote                                # noqa: E402
from coralremote import io, paths, stage1, preprocess as pp   # noqa: E402
from coralremote.cart import grow_forest          # noqa: E402
from coralremote.evaluate import compute_metrics3 # noqa: E402
from coralremote.rasterize import rasterize_labels_stage1     # noqa: E402

SEEDS = (0, 1, 2)
N_TREES = 200
MIN_LEAF = 5
MIN_PARENT = 10
CLASSES = (1, 2, 3)
CLASS_NAMES = ('CORAL', 'SND', 'RBL')
FEATURE_NAMES = tuple(pp.A3_FEATURES)
# Same four families as candidate Table 3 (old ablation design.json feature_sets).
FAMILIES = {'R': [0, 1, 2], 'R+B': [0, 1, 2, 3], 'R+G': [0, 1, 2, 4], 'R+BG': [0, 1, 2, 3, 4]}
FAMILY_LABEL = {'R': 'Three ratios', 'R+B': 'Ratios + BR', 'R+G': 'Ratios + GR', 'R+BG': 'Ratios + both'}
PAPER_FAMILY = 'R+BG'


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def plain(x):
    if isinstance(x, np.ndarray):
        return plain(x.tolist())
    if isinstance(x, np.generic):
        return plain(x.item())
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, Path):
        return str(x)
    return x


def save_json(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(plain(x), indent=2), encoding='utf-8')


def load_design_arrays():
    z = np.load(RUN / 'stage1' / 'inputs.npz')
    train = rr.zrows(z, 'accuracy')
    test = rr.zrows(z, 'test')
    X = rr.apply_rows(train, 'pooled')
    V = rr.apply_rows(test, 'pooled')
    return z, train, test, X, V


def fit_job(task):
    family, seed = task
    t0 = time.time()
    _, train, test, X, V = load_design_arrays()
    cols = FAMILIES[family]
    Xf = np.ascontiguousarray(X[:, cols])
    Vf = np.ascontiguousarray(V[:, cols])
    f = grow_forest(Xf, train['Y'], n_trees=N_TREES, min_leaf=MIN_LEAF, min_parent=MIN_PARENT,
                    seed=seed, classes=CLASSES, cost=stage1.RF_COST)
    pred = f.predict(Vf)[0]
    met = compute_metrics3(test['Y'], pred)
    return dict(family=family, seed=seed, seconds=time.time() - t0, pred=np.asarray(pred),
                metrics=met, train_sha=rr.ah(Xf), test_sha=rr.ah(Vf),
                pred_sha=hashlib.sha256(np.ascontiguousarray(np.asarray(pred, np.int64)).tobytes()).hexdigest())


def overlap_checks(z, train, test):
    ref = io.load_glint_stack().ref
    te_mask, te_lab = rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['test'])), ref, stage1.POLY_CLASS_MAP)
    tr_mask = z['train_mask']
    R, C, D, Y = train['R'], train['C'], train['D'].astype(str), train['Y']
    tR, tC, tD, tY = test['R'], test['C'], test['D'].astype(str), test['Y']
    W = ref.width
    trkey = R.astype(np.int64) * W + C
    tekey = tR.astype(np.int64) * W + tC
    tile = lambda r, c: (r // 32) * 10000 + (c // 32)
    res = dict(
        n_train=int(len(Y)), n_test=int(len(tY)),
        train_class_counts={n: int((Y == k).sum()) for k, n in zip(CLASSES, CLASS_NAMES)},
        test_class_counts={n: int((tY == k).sum()) for k, n in zip(CLASSES, CLASS_NAMES)},
        train_rows_by_date_class={f'{d}_{n}': int(((D == d) & (Y == k)).sum())
                                  for d in sorted(set(D)) for k, n in zip(CLASSES, CLASS_NAMES)
                                  if ((D == d) & (Y == k)).any()},
        test_rows_by_date_class={f'{d}_{n}': int(((tD == d) & (tY == k)).sum())
                                 for d in sorted(set(tD)) for k, n in zip(CLASSES, CLASS_NAMES)
                                 if ((tD == d) & (tY == k)).any()},
        train_rows_inside_test_polygons=int(te_mask[R, C].sum()),
        train_rows_coral_20240408=int(((D == '20240408') & (Y == 1)).sum()),
        train_rows_any_class_20240408=int((D == '20240408').sum()),
        train_pixel_locations_equal_to_test_pixel_locations=int(np.isin(np.unique(trkey), np.unique(tekey)).sum()),
        train_and_test_polygon_masks_intersect_pixels=int((tr_mask & te_mask).sum()),
        test_rows_inside_test_polygons=int(te_mask[tR, tC].sum()),
        test_rows_label_matches_test_raster=int((te_lab[tR, tC] == tY).sum()),
        test_rows_inside_training_polygons=int(tr_mask[tR, tC].sum()),
        informational_train_rows_sharing_a_32px_tile_with_a_test_pixel=int(np.isin(tile(R, C), np.unique(tile(tR, tC))).sum()),
        informational_train_rows_SND_RBL_on_test_date_20230519_other_locations=int(((D == '20230519') & (Y > 1)).sum()),
    )
    res['pass'] = (res['train_rows_inside_test_polygons'] == 0 and res['train_rows_coral_20240408'] == 0
                   and res['train_pixel_locations_equal_to_test_pixel_locations'] == 0
                   and res['train_and_test_polygon_masks_intersect_pixels'] == 0
                   and res['test_rows_inside_test_polygons'] == len(tY)
                   and res['test_rows_label_matches_test_raster'] == len(tY)
                   and res['n_train'] == 17372 and res['n_test'] == 905
                   and res['test_class_counts'] == dict(CORAL=200, SND=322, RBL=383))
    return res


def read_stage1_accuracy_csv():
    rows = list(csv.DictReader((RUN / 'report' / 'stage1_accuracy.csv').open(encoding='utf-8')))
    return {int(r['seed']): r for r in rows if r['policy'] == 'pooled'}


def sanity_gate(results, test):
    ref_csv = read_stage1_accuracy_csv()
    out = dict(family=PAPER_FAMILY, seeds={}, criteria='predictions identical to saved Table 4A evaluation arrays; '
               'OA, macro-F1 and per-class P/R/F1 equal report/stage1_accuracy.csv within 1e-12')
    ok = True
    for s in SEEDS:
        r = results[(PAPER_FAMILY, s)]
        ev = np.load(RUN / 'stage1' / f'pooled_accuracy_seed{s}_evaluation.npz')
        saved = json.loads((RUN / 'stage1' / f'pooled_accuracy_seed{s}.json').read_text(encoding='utf-8'))
        m = r['metrics']
        same_rows = (np.array_equal(ev['y'], test['Y']) and np.array_equal(ev['row'], test['R'])
                     and np.array_equal(ev['col'], test['C']) and np.array_equal(ev['date'].astype(str), test['D'].astype(str)))
        pred_equal = bool(np.array_equal(ev['pred'], r['pred']))
        c = ref_csv[s]
        diffs = dict(OA=abs(m['OA'] - float(c['OA'])), macro_F1=abs(float(np.mean(m['f1'][:3])) - float(c['macro_F1'])))
        for i, n in enumerate(CLASS_NAMES):
            for key, col in (('prec', 'prec'), ('rec', 'rec'), ('f1', 'f1')):
                diffs[f'{n}_{key}'] = abs(m[key][i] - float(c[f'{n}_{col}']))
        feat_hash = dict(train=r['train_sha'] == saved['train_features_sha'], test=r['test_sha'] == saved['test_features_sha'])
        seed_ok = same_rows and pred_equal and max(diffs.values()) <= 1e-12 and all(feat_hash.values())
        ok &= seed_ok
        out['seeds'][s] = dict(pass_=seed_ok, test_rows_identical=same_rows, predictions_identical=pred_equal,
                               n_prediction_differences=int(np.sum(ev['pred'] != r['pred'])),
                               feature_hashes_match_saved_run=feat_hash, max_abs_metric_difference=max(diffs.values()),
                               OA=m['OA'], OA_reference=float(c['OA']), CORAL_recall=m['rec'][0],
                               CORAL_recall_reference=float(c['CORAL_rec']), kappa=m['kappa'],
                               kappa_reference_saved_json=saved['metrics']['kappa'])
    out['pass'] = bool(ok)
    return out


def old_table3():
    """Existing candidate Table 3 values, read from the old ablation outputs (unchanged)."""
    summ = {r['feature_set']: r for r in csv.DictReader((OLD / 'report' / 'ablation_summary.csv').open(encoding='utf-8'))
            if r['scenario'] == 'stable_OLS'}
    out = {}
    for fam in FAMILIES:
        per = []
        for s in SEEDS:
            # R does not depend on the coefficients; the old run stored it under 'original' only.
            name = f'original_R_historical_seed{s}.json' if fam == 'R' else f'stable_OLS_{fam.replace("+", "_")}_historical_seed{s}.json'
            per.append(json.loads((OLD / 'stage1_ablation' / name).read_text(encoding='utf-8'))['metrics']['OA'])
        out[fam] = dict(spatial_cv_macro_f1=float(summ[fam]['fold_seed_mean_macro_f1']),
                        old_exploratory_test_OA_mean=float(summ[fam]['mean_historical_OA']),
                        old_exploratory_test_OA_seeds=per)
    return out


def module_hashes():
    mods = {}
    for name, m in sorted(sys.modules.items()):
        f = getattr(m, '__file__', None)
        if f and (name == 'coralremote' or name.startswith('coralremote.') or name == 'revision_runner'):
            mods[name] = dict(path=str(Path(f).resolve()), sha256=sha(f))
    return mods


def git_hash(repo):
    try:
        h = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'coralremote'], capture_output=True, text=True).stdout.strip()
        return dict(commit=h, coralremote_package_dirty=bool(dirty))
    except Exception as e:  # pragma: no cover
        return dict(error=str(e))


def main():
    t0 = time.time()
    rr.guard_out(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    assert Path(coralremote.__file__).resolve().parent == (PORT / 'coralremote').resolve(), coralremote.__file__
    assert FEATURE_NAMES == ('B/G1', 'G1/G2', 'B/R', 'DII_BR', 'DII_GR') or len(FEATURE_NAMES) == 5, FEATURE_NAMES
    cal = json.loads((RUN / 'calibration' / 'calibration.json').read_text(encoding='utf-8'))
    pooled_pair = cal['pooled']['pair']
    assert abs(pooled_pair[0] - 0.380987039988) < 1e-11 and abs(pooled_pair[1] - 0.374286191215) < 1e-11, pooled_pair

    z, train, test, X, V = load_design_arrays()
    overlap = overlap_checks(z, train, test)
    save_json(OUT / 'overlap_checks.json', overlap)
    if not overlap['pass']:
        raise SystemExit(f'Design/overlap check failed: {overlap}')

    tasks = [(f, s) for f in FAMILIES for s in SEEDS]
    results = {}
    with cf.ProcessPoolExecutor(max_workers=int(os.environ.get('T3_WORKERS', '6'))) as pool:
        for r in pool.map(fit_job, tasks):
            results[(r['family'], r['seed'])] = r
            print(f"{datetime.now().isoformat(timespec='seconds')} {r['family']} seed{r['seed']}: "
                  f"OA={r['metrics']['OA']:.4f} ({r['seconds']:.0f}s)", flush=True)

    gate = sanity_gate(results, test)
    save_json(OUT / 'sanity_gate.json', gate)
    if not gate['pass']:
        raise SystemExit('SANITY GATE FAILED: five-feature family does not reproduce Table 4A; no results written.')

    old = old_table3()
    per_seed_rows, summary_rows = [], []
    for fam, cols in FAMILIES.items():
        for s in SEEDS:
            m = results[(fam, s)]['metrics']
            row = dict(family=fam, features=';'.join(FEATURE_NAMES[c] for c in cols), n_features=len(cols), seed=s,
                       OA_pct=100 * m['OA'], kappa=m['kappa'], macro_F1=float(np.mean(m['f1'][:3])),
                       n_test=int(m['total'][:3].sum()), n_correct=int(m['tp'][:3].sum()), n_train=len(train['Y']))
            for i, n in enumerate(CLASS_NAMES):
                row.update({f'{n}_P': m['prec'][i], f'{n}_R': m['rec'][i], f'{n}_F1': m['f1'][i], f'{n}_n': int(m['total'][i])})
            row['confusion_3x3_rows_true'] = json.dumps(m['C'][:3, :3].tolist())
            row['pred_sha256'] = results[(fam, s)]['pred_sha']
            per_seed_rows.append(row)
        seeds = [r for r in per_seed_rows if r['family'] == fam]
        s0 = seeds[0]
        srow = dict(family=fam, label=FAMILY_LABEL[fam], features=s0['features'], n_features=len(cols))
        for k in ['OA_pct', 'kappa', 'macro_F1'] + [f'{n}_{q}' for n in CLASS_NAMES for q in ('P', 'R', 'F1')]:
            vals = [r[k] for r in seeds]
            srow[f'{k}_seed0'] = vals[0]
            srow[f'{k}_min'] = min(vals)
            srow[f'{k}_max'] = max(vals)
            srow[f'{k}_mean'] = float(np.mean(vals))
        for n in CLASS_NAMES:
            srow[f'{n}_n'] = s0[f'{n}_n']
        srow['n_test'] = s0['n_test']
        srow['n_train'] = s0['n_train']
        srow['spatial_cv_macro_F1_unchanged'] = old[fam]['spatial_cv_macro_f1']
        srow['old_exploratory_test_OA_pct_mean'] = 100 * old[fam]['old_exploratory_test_OA_mean']
        srow['old_exploratory_test_OA_pct_seed0'] = 100 * old[fam]['old_exploratory_test_OA_seeds'][0]
        summary_rows.append(srow)

    def rank(key):
        order = sorted(summary_rows, key=lambda r: -r[key])
        return [r['family'] for r in order]
    ranking = dict(new_OA_seed0=rank('OA_pct_seed0'), new_OA_mean=rank('OA_pct_mean'),
                   old_exploratory_OA_mean=rank('old_exploratory_test_OA_pct_mean'),
                   spatial_cv_macro_F1=rank('spatial_cv_macro_F1_unchanged'))

    def write_csv(p, rows):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with p.open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(plain(rows))
    write_csv(OUT / 'results.csv', summary_rows)
    write_csv(OUT / 'results_per_seed.csv', per_seed_rows)

    inputs = [RUN / 'stage1' / 'inputs.npz', RUN / 'calibration' / 'calibration.json', RUN / 'report' / 'stage1_accuracy.csv',
              RUN / 'manifest' / 'design.json', RUNNER_DIR / 'revision_runner.py',
              *[RUN / 'stage1' / f'pooled_accuracy_seed{s}{suf}' for s in SEEDS for suf in ('.json', '_evaluation.npz', '_forest.settings.json')],
              paths.STAGE1_LABELS['train'], paths.STAGE1_LABELS['test'], paths.GLINT_STACK,
              OLD / 'report' / 'ablation_summary.csv', OLD / 'stage1_ablation' / 'design.json',
              *[OLD / 'stage1_ablation' / (f'original_R_historical_seed{s}.json' if f == 'R' else f'stable_OLS_{f.replace("+", "_")}_historical_seed{s}.json')
                for f in FAMILIES for s in SEEDS]]
    mods = module_hashes()
    frozen = {Path(h['path']).name: h['sha256'] for h in json.loads((RUN / 'manifest' / 'design.json').read_text(encoding='utf-8'))['input_hashes']
              if '\\coralremote-python\\coralremote\\' in h['path'] or '/coralremote-python/coralremote/' in h['path']}
    mod_vs_frozen = {k: (frozen.get(Path(v['path']).name) == v['sha256']) for k, v in mods.items() if k.startswith('coralremote')}
    design = dict(
        created=datetime.now().astimezone().isoformat(timespec='seconds'),
        purpose='Table 3 feature-family ablation recomputed on the Table 4A Stage 1 held-out design',
        script=dict(path=str(SCRIPT), sha256=sha(SCRIPT)),
        python=sys.version, numpy=np.__version__,
        port_git=git_hash(PORT), imported_modules=mods,
        imported_port_modules_match_revision_run_freeze_hashes=mod_vs_frozen,
        inputs=[dict(path=str(p), sha256=sha(p), bytes=Path(p).stat().st_size) for p in inputs],
        design=dict(training='stage1/inputs.npz accuracy_* rows (17,372; CORAL from 20230404, 20230509, 20240307, 20241111; '
                             'SND/RBL from 20230404, 20230509, 20230519, 20230629, 20240307, 20240317)',
                    test='stage1/inputs.npz test_* rows (905; CORAL 200 from 20240408, SND 322 and RBL 383 from 20230519; Stage 1 test polygons)',
                    features=FEATURE_NAMES, families={k: [FEATURE_NAMES[c] for c in v] for k, v in FAMILIES.items()},
                    coefficients=dict(policy='pooled', pair=pooled_pair, applied_by='revision_runner.apply_rows -> adjust'),
                    forest=dict(function='coralremote.cart.grow_forest', n_trees=N_TREES, min_leaf=MIN_LEAF, min_parent=MIN_PARENT,
                                predictors_per_split='all', classes=CLASSES, cost=stage1.RF_COST,
                                prior='empirical class shares of the training rows (set inside grow_forest)'),
                    metrics='coralremote.evaluate.compute_metrics3', seeds=SEEDS,
                    varied='predictor columns only; bootstrap draws per seed are identical across families'),
        sanity_gate_pass=gate['pass'], overlap_checks_pass=overlap['pass'], ranking=ranking,
        runtime_seconds=time.time() - t0)
    save_json(OUT / 'design.json', design)
    print(json.dumps(plain(dict(gate=gate['pass'], ranking=ranking)), indent=1))
    for r in summary_rows:
        print(f"{r['family']:5s} OA seed0 {r['OA_pct_seed0']:.1f}% range {r['OA_pct_min']:.1f}-{r['OA_pct_max']:.1f} "
              f"kappa0 {r['kappa_seed0']:.3f} CORAL R0 {r['CORAL_R_seed0']:.3f} | old {r['old_exploratory_test_OA_pct_mean']:.1f} "
              f"CV {r['spatial_cv_macro_F1_unchanged']:.3f}")


if __name__ == '__main__':
    main()
