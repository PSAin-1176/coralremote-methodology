"""T13--T19 isolated revised-analysis candidate, local release entry point.

Prepared GLINT/raw stacks and cached coregistered UDM2 are licensed inputs.
The historical manuscript, source code, labels and outputs are read-only.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
import sys,json,hashlib,time,pickle,csv,inspect,shutil,traceback,argparse,re
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import concurrent.futures as cf
HERE=Path(__file__).resolve().parent
REPO=HERE.parent
def _config_root(argv):
    # Fallback for CORALREMOTE_ROOT: a literal (non-placeholder) data_root in the --config file.
    for i,x in enumerate(argv):
        if x=='--config' and i+1<len(argv):
            try:v=json.loads(Path(argv[i+1]).read_text(encoding='utf-8')).get('data_root')
            except (OSError,ValueError):return None
            return v if v and '${' not in v else None
    return None
_root=os.environ.get('CORALREMOTE_ROOT') or _config_root(sys.argv)
if not _root:raise RuntimeError('CORALREMOTE_ROOT is not set. Set it to the folder holding Source/ and Asset/ (or give a literal data_root in the configuration).')
ROOT=Path(_root);os.environ.setdefault('CORALREMOTE_ROOT',str(ROOT))
OUT=Path(os.environ.get('REVISION_OUT',str(HERE/'output')))
UDM=Path(os.environ.get('REVISION_UDM',str(ROOT/'coralremote-python-cache/udm2')))
sys.path.insert(0,str(REPO))
import numpy as np
from coralremote import io,paths,preprocess as pp,stage1,stage2,mask,features as fe
from coralremote.cart import grow_forest
from coralremote.forest import Forest
from coralremote.matlab_compat import find_colmajor
from coralremote.evaluate import compute_metrics3
from coralremote.train import train_stage2
from coralremote.rasterize import rasterize_labels_stage1,feature_pixels
BASE=('20230404','20230509','20240307')
NEEDED=('20230404','20230509','20230519','20230629','20240307','20240317','20240408','20241111')
JUNE=('20240613','20240616','20240620')
SEEDS=(0,1,2);NAMES=('HPY','PLE','BLD','UNK')
WORKERS=int(os.environ.get('REVISION_WORKERS','8'))
def stablesand_path():
    v=os.environ.get('CORALREMOTE_STABLESAND')
    if not v:raise RuntimeError('CORALREMOTE_STABLESAND is not set. Set it to the StableSand.geojson file identified in the configuration.')
    return Path(v)
def expand(s):
    # Resolve ${NAME} placeholders in configuration paths.
    known={'CORALREMOTE_ROOT':str(ROOT),'CORALREMOTE_A3_ROOT':str(paths.A3_ROOT),'REVISION_UDM':str(UDM)}
    def one(m):
        n=m.group(1);v=known.get(n) or os.environ.get(n)
        if not v:raise RuntimeError(f'Configuration placeholder ${{{n}}} is not set')
        return v
    return re.sub(r'\$\{([A-Za-z0-9_]+)\}',one,s)
pp.udm2_clear=lambda d:np.load(UDM/f'udm2_{d}.npz')['udm2'][0]==1

def plain(x):
    if isinstance(x,np.ndarray):return plain(x.tolist())
    if isinstance(x,np.generic):return plain(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    if isinstance(x,Path):return str(x)
    return x
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(plain(x),indent=2,allow_nan=False),encoding='utf-8')
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def ah(x):return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def dump(p,x):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with Path(p).open('wb') as f:pickle.dump(x,f,protocol=5)
def undump(p):
    with Path(p).open('rb') as f:return pickle.load(f)
def csvsave(p,rows):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(plain(rows))
def log(s):print(datetime.now(ZoneInfo('Asia/Bangkok')).isoformat(timespec='seconds'),s,flush=True)
def parallel(fn,tasks,workers=WORKERS):
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        fs={pool.submit(fn,t):t for t in tasks}
        for f in cf.as_completed(fs):
            try:log(f.result())
            except Exception:
                save(OUT/f'failures/{time.time_ns()}.json',dict(task=str(fs[f]),traceback=traceback.format_exc()));raise

def freeze():
    p=OUT/'manifest/design.json'
    if p.exists():return
    for d in ('manifest','calibration','stage1','stage2','loo','omissions','release','report','manuscript_audit'): (OUT/d).mkdir(parents=True,exist_ok=True)
    source=stablesand_path()
    plan=os.environ.get('REVISION_PLAN');plan=Path(plan) if plan and Path(plan).exists() else None
    required=[paths.RAW_STACK,paths.GLINT_STACK,paths.REEF_MASK_TIF,source,*paths.STAGE1_LABELS.values(),*paths.LABELS.values()]
    required+=list(paths.PHASE3_GLINT_DIR.glob('GLINT_*.tif'))+[UDM/f'udm2_{d}.npz' for d in NEEDED]
    required+=list((REPO/'coralremote').glob('*.py'))
    save(p,dict(created_Bangkok=datetime.now(ZoneInfo('Asia/Bangkok')).isoformat(),plan=str(plan) if plan else None,plan_sha256=sha(plan) if plan else None,root=str(ROOT),a3_root=str(paths.A3_ROOT),udm_root=str(UDM),input_hashes=[dict(path=str(x),sha256=sha(x),bytes=x.stat().st_size) for x in required],sand_reference=str(source),estimator='OLS with intercept, within-date centring for pooled fit',calibration_dates=BASE,additional_calm_dates=['20240317','20240408'],Stage1_required_dates=NEEDED,raw_NIR_strict_interval=[.005,.12],calibration_floor_reflectance=1e-4,Stage1_feature_floor_reflectance=1e-6,Stage2_corrected_floor_reflectance=1e-4,unit_multiplier_stored_to_reflectance=10000,seeds=SEEDS,workers=WORKERS,threads_per_worker=1,features=list(pp.A3_FEATURES),Stage1_mask_model=dict(trees=200,min_leaf=5,min_parent=10,cost=None,coral_dates=BASE,threshold=.55),Stage1_accuracy_model=dict(trees=200,min_leaf=5,min_parent=10,cost=stage1.RF_COST,coral_dates=stage1.CORAL_DATES,test_coral_date=stage1.CORAL_TEST_DATE),Stage2=dict(trees_per_model=250,min_leaf=5,min_parent=10,threshold=.55,mode='fixed',engine='matlab_cart',sampling='numpy'),eligibility=dict(n_min=50,occupied_sand_polygons_min=4,log_red_std_min=.05,R2_min=.5,CI_lower_positive=True,relative_95CI_width_max=1.,bootstrap_units='first-writer sand polygon IDs, all dates jointly resampled; spatial independence not established',bootstrap_replicates=1000,bootstrap_seed=20261118,fallback='eligible three-baseline pooled pair; label every fallback'),selection_rule=dict(per_date_requires_all_required_dates_eligible_without_fallback=True,common_support_baseline_contrasts_require_both_band_CIs_exclude_zero=True,spatial_fold_macroF1_guard_drop_max=.01,guard_rationale='Prospective pragmatic complexity guard against a >=1pp loss, not a scientific equivalence margin.',otherwise='pooled candidate; evidence for per-date adoption inconclusive',independent_of='June20 ordering, final historical test scores, desired temporal story'),evaluation=dict(cv='3 fixed StratifiedGroupKFold folds, 32pixel tiles, shuffle seed42; all dates of each tile kept together',target_scene_calibration='Sand from target acquisition permitted; habitat labels and validation tiles excluded from calibration.',final_test='Reused historical spatial/date test, not untouched confirmation',no_universal_equivalence_margin=True),LOO=dict(polygons='36 training HC geometries',trees=50,cost=[[0,3,3],[1,0,1],[1,1,0]],seeds=SEEDS,baseline='New explicit nanmedian of the three baseline-date feature arrays; original Phase6 construction not fully recovered; do not label exact historical replay',exclude='Held polygon geometry from CORAL training on every date; sand and rubble from revised baseline composite'),revised_candidate='Local candidate, not adoption into original manuscript; no upload/push'))
    log('Updated design frozen before new paired outcomes')

def sand_geometry():
    from rasterio.features import rasterize
    from rasterio.warp import transform_geom
    from affine import Affine
    ref=io.load_glint_stack().ref;geo=load(load(OUT/'manifest/design.json')['sand_reference'])
    shapes=[(transform_geom('EPSG:4326','EPSG:32647',f['geometry']),i+1) for i,f in reversed(list(enumerate(geo['features'])))]
    own=rasterize(shapes,out_shape=(ref.height,ref.width),transform=Affine(*ref.transform),dtype='int16')
    _,labels=rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['test'])),ref,stage1.POLY_CLASS_MAP)
    rr,cc=np.indices(own.shape);tiles=rr//32*10000+cc//32;exclude=np.unique(tiles[labels>0]);own[np.isin(tiles,exclude)]=0
    return own,tiles,exclude

def raw_calibration(d,own):
    A=pp.glint_reflectance(d);nir=pp.raw_stack().get(d)[:,:,7]
    valid=(own>0)&np.isfinite(A[:,:,[1,2,5]]).all(2)&(nir>.005)&(nir<.12)
    x=np.log(np.maximum(A[:,:,5][valid].astype(float),1e-4));y=np.stack([np.log(np.maximum(A[:,:,i][valid].astype(float),1e-4)) for i in (1,2)],axis=1)
    return valid,x,y,own[valid]

def moments(x,y,ids,npoly):
    return np.stack([np.bincount(ids,weights=z,minlength=npoly+1)[1:] for z in (np.ones(len(x)),x,y[:,0],y[:,1],x*x,x*y[:,0],x*y[:,1],y[:,0]**2,y[:,1]**2)],axis=1)
def estimate(M,weights=None):
    # Dates x polygons x moments -> within-date-centred pooled moments.
    if weights is None:weights=np.ones(M.shape[1])
    z=(M*weights[None,:,None]).sum(1);n=z[:,0];good=n>1;z=z[good];n=n[good]
    vx=np.sum(z[:,4]-z[:,1]**2/n);cov=np.sum(z[:,5:7]-z[:,1,None]*z[:,2:4]/n[:,None],axis=0);vy=np.sum(z[:,7:9]-z[:,2:4]**2/n[:,None],axis=0)
    return cov/vx, cov**2/(vx*vy),float(np.sqrt(vx/n.sum())),int(n.sum())
def bootstrap(M):
    rng=np.random.default_rng(20261118);n=M.shape[1];res=[]
    for _ in range(1000):
        w=np.bincount(rng.integers(0,n,n),minlength=n)
        try:k,_,_,_=estimate(M,w);res.append(k)
        except (FloatingPointError,ZeroDivisionError):res.append([np.nan,np.nan])
    return np.array(res)
def fit_summary(M):
    k,r2,std,n=estimate(M);bs=bootstrap(M);ci=np.nanquantile(bs,[.025,.975],axis=0);polys=int((M[:,:,0].sum(0)>0).sum())
    reasons=[]
    if n<50:reasons.append('n<50')
    if polys<4:reasons.append('fewer than4 occupied polygons')
    if std<.05:reasons.append('log-red std<0.05')
    if np.any(r2<.5):reasons.append('R2<0.50')
    if np.any(ci[0]<=0):reasons.append('nonpositive CI bound')
    if np.any((ci[1]-ci[0])/k>1):reasons.append('relative CI width>1')
    return dict(pair=k,R2=r2,log_red_std=std,n=n,occupied_polygons=polys,CI95=ci,eligible=not reasons,failure_reasons=reasons),bs

def calibrate():
    dest=OUT/'calibration/calibration.json'
    if dest.exists():return
    own,tiles,excluded=sand_geometry();npoly=int(own.max());raw={d:raw_calibration(d,own) for d in NEEDED};M={d:moments(v[1],v[2],v[3],npoly) for d,v in raw.items()}
    pool,bs=fit_summary(np.stack([M[d] for d in BASE]));assert pool['eligible'],'Shared fallback fails prespecified eligibility'
    fits={};boots={};mapping={}
    for d in NEEDED:
        fits[d],boots[d]=fit_summary(M[d][None]);mapping[d]=fits[d]['pair'] if fits[d]['eligible'] else pool['pair']
    common=np.logical_and.reduce([raw[d][0] for d in BASE]);CM={}
    for d in BASE:
        A=pp.glint_reflectance(d);x=np.log(np.maximum(A[:,:,5][common].astype(float),1e-4));y=np.stack([np.log(np.maximum(A[:,:,i][common].astype(float),1e-4)) for i in (1,2)],axis=1);CM[d]=moments(x,y,own[common],npoly)
    commons={};cboots={}
    for d in BASE:commons[d],cboots[d]=fit_summary(CM[d][None])
    contrasts=[]
    for i,d in enumerate(BASE):
        for e in BASE[i+1:]:
            delta=cboots[e]-cboots[d];ci=np.nanquantile(delta,[.025,.975],axis=0);contrasts.append(dict(date_a=d,date_b=e,difference=np.asarray(commons[e]['pair'])-commons[d]['pair'],CI95=ci,both_bands_exclude_zero=bool(np.all((ci[0]>0)|(ci[1]<0)))))
    drop={d:fit_summary(np.stack([M[e] for e in BASE if e!=d]))[0] for d in BASE}
    fallback=[d for d in NEEDED if not fits[d]['eligible']]
    # Fixed-pixel feature changes and verified date-to-row mapping in next phase.
    save(dest,dict(pooled=pool,per_date=fits,per_date_applied=mapping,fallback_dates=fallback,common_support=commons,common_support_n=int(common.sum()),common_support_contrasts=contrasts,leave_one_calibration_date_out=drop,spatial_test_tiles_excluded=excluded,sand_polygons=npoly,units='true-reflectance log bands; floor1e-4',note='Bootstrap intervals use polygon resampling; polygons may be spatially correlated. Eligibility gates are operational criteria, not physical validation.'))
    np.savez_compressed(OUT/'calibration/calibration_pixels.npz',owner=own,excluded_test_tiles=excluded,common_support=common,**{f'valid_{d}':v[0] for d,v in raw.items()})
    rows=[]
    for d,f in fits.items():
        for j,b in enumerate(('BR','GR')):rows.append(dict(date=d,band=b,slope=f['pair'][j],applied=mapping[d][j],n=f['n'],occupied_polygons=f['occupied_polygons'],log_red_std=f['log_red_std'],R2=f['R2'][j],lower95=f['CI95'][0][j],upper95=f['CI95'][1][j],eligible=f['eligible'],fallback_reason=';'.join(f['failure_reasons'])))
    csvsave(OUT/'calibration/date_coefficients.csv',rows);log(f'Calibration complete: pooled {pool["pair"]}, per-date fallbacks {fallback}')

def policy_pair(policy,d):
    cal=load(OUT/'calibration/calibration.json');return cal['pooled']['pair'] if policy=='pooled' else cal['per_date_applied'][d]
def adjust(F,L,k):
    F=F.copy();F[:,3]+=np.float32(.82-k[0])*L;F[:,4]+=np.float32(.80-k[1])*L;return F
def date_features(d,policy):
    F=pp.phase5_features(d);A=pp.glint_reflectance(d);L=pp._log32(np.maximum(A[:,:,5],np.float32(1e-6)));return adjust(F.reshape(-1,5),L.ravel(),policy_pair(policy,d)).reshape(F.shape)

def collect(coral_dates):
    ref=io.load_glint_stack().ref;tr,tl=rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['train'])),ref,stage1.POLY_CLASS_MAP);te,el=rasterize_labels_stage1(io.load_geojson(str(paths.STAGE1_LABELS['test'])),ref,stage1.POLY_CLASS_MAP)
    assert not np.any(tr&te)
    def rows(plans,ma,lab):
        arrays={k:[] for k in ('F','L','D','R','C','Y')}
        for d,k in plans:
            F=pp.phase5_features(d);r,c=find_colmajor(ma&(lab==k)&np.isfinite(F).all(2));L=pp._log32(np.maximum(pp.glint_reflectance(d)[r,c,5],np.float32(1e-6)))
            for n,v in dict(F=F[r,c],L=L,D=np.full(len(r),d),R=r,C=c,Y=lab[r,c]).items():arrays[n].append(v)
        return {k:np.concatenate(v) for k,v in arrays.items()}
    train=rows([(d,1) for d in coral_dates]+[(d,k) for d in stage1.SND_RBL_DATES for k in (2,3)],tr,tl);test=rows([(stage1.CORAL_TEST_DATE,1),(stage1.SNDRBL_TEST_DATE,2),(stage1.SNDRBL_TEST_DATE,3)],te,el)
    return train,test,tr,tl

def apply_rows(z,policy,mapping=None):
    X=z['F'].copy()
    for d in np.unique(z['D']):
        ix=z['D']==d;k=(mapping[d] if mapping else policy_pair(policy,d));X[ix]=adjust(z['F'][ix],z['L'][ix],k)
    assert np.array_equal(X[:,:3],z['F'][:,:3]);return X

def prepare_stage1():
    if (OUT/'stage1/prepared.json').exists():return
    from sklearn.model_selection import StratifiedGroupKFold
    a,t,tr,tl=collect(stage1.CORAL_DATES);b,_,_,_=collect(BASE)
    assert len(a['Y'])==17372 and len(b['Y'])==16434 and len(t['Y'])==905
    groups=a['R']//32*10000+a['C']//32;folds=np.full(len(groups),-1,int)
    for i,(_,ix) in enumerate(StratifiedGroupKFold(3,shuffle=True,random_state=42).split(a['F'],a['Y'],groups)):folds[ix]=i
    np.savez_compressed(OUT/'stage1/inputs.npz',**{f'{name}_{k}':v for name,z in [('accuracy',a),('test',t),('mask',b)] for k,v in z.items()},folds=folds,groups=groups,train_mask=tr,train_labels=tl)
    feature_changes=[]
    for d in NEEDED:
        ix=a['D']==d
        if not ix.any():continue
        aa=apply_rows({k:v[ix] for k,v in a.items()},'pooled');bb=apply_rows({k:v[ix] for k,v in a.items()},'per_date');dif=bb-aa
        feature_changes.append(dict(date=d,n=int(ix.sum()),ratio_columns_bit_exact=np.array_equal(aa[:,:3],bb[:,:3]),DII_BR_change_std=float(dif[:,3].std()),DII_GR_change_std=float(dif[:,4].std()),max_abs_difference=float(np.max(np.abs(dif)))))
    csvsave(OUT/'stage1/feature_changes.csv',feature_changes)
    save(OUT/'stage1/prepared.json',dict(accuracy_n=len(a['Y']),mask_n=len(b['Y']),test_n=len(t['Y']),folds=[dict(fold=i,train_n=int((folds!=i).sum()),validation_n=int((folds==i).sum()),validation_groups=np.unique(groups[folds==i])) for i in range(3)],training_dates={str(d):int(np.sum(a['D']==d)) for d in np.unique(a['D'])},date_mapping_verified=True,ratio_identity=True))
    log('Stage1 rows/date mapping prepared; all ratio columns identical')

def zrows(z,prefix):return {k:z[f'{prefix}_{k}'] for k in ('F','L','D','R','C','Y')}
def save_forest(f,p):
    f.save_npz(p);save(Path(p).with_suffix('.settings.json'),dict(cost=f.cost,prior=f.prior))
def read_forest(p):
    f=Forest.load_npz(p);extra=load(Path(p).with_suffix('.settings.json'))
    f.cost=None if extra['cost'] is None else np.asarray(extra['cost']);f.prior=None if extra['prior'] is None else np.asarray(extra['prior']);return f

def stage1_job(task):
    policy,kind,seed=task;tag=f'{policy}_{kind}_seed{seed}';dest=OUT/f'stage1/{tag}.json'
    if dest.exists():return tag+' reused'
    started=time.time();z=np.load(OUT/'stage1/inputs.npz');train=zrows(z,'mask' if kind=='mask' else 'accuracy');test=zrows(z,'test');mapping=None
    if kind.startswith('cv'):
        fold=int(kind[-1]);va=z['folds']==fold;test={k:v[va] for k,v in train.items()};train={k:v[~va] for k,v in train.items()}
        # Calibration sand must exclude the held spatial tiles as well as final-test tiles.
        own,tiles,_=sand_geometry();exclude=np.unique(z['groups'][va]);own[np.isin(tiles,exclude)]=0;npoly=int(own.max());MM={d:moments(*raw_calibration(d,own)[1:],npoly) for d in NEEDED};pool,_=fit_summary(np.stack([MM[d] for d in BASE]));assert pool['eligible'];mapping={}
        for d in NEEDED:
            fit,_=fit_summary(MM[d][None]);mapping[d]=fit['pair'] if policy=='per_date' and fit['eligible'] else pool['pair']
    X=apply_rows(train,policy,mapping);V=apply_rows(test,policy,mapping)
    cost=None if kind=='mask' else stage1.RF_COST
    f=grow_forest(X,train['Y'],n_trees=200,seed=seed,classes=(1,2,3),cost=cost)
    p=f.predict(V)[0];met=compute_metrics3(test['Y'],p)
    result=dict(policy=policy,kind=kind,seed=seed,seconds=time.time()-started,n_train=len(X),n_test=len(V),metrics=met,macro_f1=float(np.mean(met['f1'][:3])),train_features_sha=ah(X),test_features_sha=ah(V),train_dates=train['D'],date_policy_mapping=mapping,ratio_columns_sha=ah(X[:,:3]))
    np.savez_compressed(OUT/f'stage1/{tag}_evaluation.npz',pred=p,y=test['Y'],date=test['D'],row=test['R'],col=test['C'])
    save_forest(f,OUT/f'stage1/{tag}_forest.npz')
    if kind=='mask':
        reef=pp.reef_mask();maps=[]
        for d in BASE:
            A=date_features(d,policy);r,c=np.where(reef&np.isfinite(A).all(2));pr,sc=f.predict(A[r,c]);pr[sc.max(1)<.55]=4;cm=np.zeros(reef.shape,np.uint8);cm[r,c]=pr;maps.append(cm)
        chain=mask.build_baseline_mask(maps,reef);np.savez_compressed(OUT/f'stage1/{tag}_masks.npz',**chain,**{f'class_{d}':a for d,a in zip(BASE,maps)});result['mask_n']=int(chain['maskFinal'].sum())
    save(dest,result);return f'{tag}: {time.time()-started:.1f}s macroF1={result["macro_f1"]:.4f}'

def select():
    c=load(OUT/'calibration/calibration.json');cv={}
    for p in ('pooled','per_date'):
        rr=[load(OUT/f'stage1/{p}_cv{i}_seed{s}.json') for i in range(3) for s in SEEDS if (OUT/f'stage1/{p}_cv{i}_seed{s}.json').exists()]
        cv[p]=dict(n_completed=len(rr),mean_macro_f1=None if not rr else float(np.mean([r['macro_f1'] for r in rr])))
    reliable=not c['fallback_dates'];common=any(r['both_bands_exclude_zero'] for r in c['common_support_contrasts']);complete=all(z['n_completed']==9 for z in cv.values());guard=complete and cv['per_date']['mean_macro_f1']>=cv['pooled']['mean_macro_f1']-.01
    # The failure of the all-date eligibility requirement suffices to reject
    # adoption without waiting for or selecting on downstream maps.
    choice='per_date' if reliable and common and guard else 'pooled'
    save(OUT/'manifest/selection.json',dict(selected_candidate=choice,documented_Bangkok=datetime.now(ZoneInfo('Asia/Bangkok')).isoformat(),all_dates_eligible=reliable,common_support_evidence=common,spatial_cv_guard_pass=guard if complete else None,cv=cv,comparison_complete=complete,fallback_dates=c['fallback_dates'],decision='Per-date adoption inconclusive under frozen reliability criteria; pooled is the documented candidate.' if choice=='pooled' else 'Per-date criteria met, without considering June20 ordering or historical final-test scores.',scope='Local candidate selection for tests; original manuscript not replaced.'))
    log(f'Candidate rule {choice}; completed CV {complete}; per-date fallbacks {c["fallback_dates"]}')

def sample_keys(view,support,fs,cm):
    out=[]
    for d in dict.fromkeys(str(f['properties']['date']) for f in fs):
        mm,ll=stage2.rasterize_date_labels(fs,view.ref,cm,d);r,c=find_colmajor(mm&support&np.isfinite(view.get(d)).all(2)&(ll>0))
        out.extend(f'{d}:{ri}:{ci}:{int(y)}' for ri,ci,y in zip(r,c,ll[r,c]))
    return np.array(out,dtype='<U40')

def prepare_support(support,drop=()):
    oldmask=io.load_baseline_mask;oldsample=stage2.build_raw_spectral_split;captured=[]
    def wrapper(view,m,fs,cm,return_dates=False):
        r=oldsample(view,m,fs,cm,return_dates);captured.append((r[0],sample_keys(view,m,fs,cm),fs));return r
    io.load_baseline_mask=lambda *a,**kw:(support,*support.shape);stage2.build_raw_spectral_split=wrapper
    try:z=stage2.prepare_frozen_inputs(stage2.Stage2Options(mode='fixed'),drop_train_dates=drop)
    finally:io.load_baseline_mask=oldmask;stage2.build_raw_spectral_split=oldsample
    def keys(X):return next(k for a,k,_ in captured if a is X)
    z['test_keys']=keys(z['Xtest']);z['train_keys']={c:keys(X) for c,X in z['Xtrain'].items()}
    # exact date/pixel overlap is forbidden; cross-date location reuse reported.
    tr={k.rsplit(':',1)[0] for kk in z['train_keys'].values() for k in kk};te={k.rsplit(':',1)[0] for k in z['test_keys']};assert not tr&te
    assert not set(drop)&set(np.concatenate(list(z['Dtrain'].values())))
    z['drop_train_dates']=drop;z['training_polygons']={c:next(fs for a,k,fs in captured if a is X) for c,X in z['Xtrain'].items()}
    z['test_polygons']=next(fs for a,k,fs in captured if a is z['Xtest'])
    return z

def unit_checks():
    if (OUT/'manifest/unit_checks.json').exists():return
    from coralremote.units import UnitError
    from types import SimpleNamespace
    values=np.array([[.02]*8,[1e-5]*8,[1e-4]*8,[np.nan]+[.02]*7],np.float32).reshape(1,4,8)
    class Probe:
        on_disk_scale=10000;dates=['probe'];ref=SimpleNamespace(height=1,width=4)
        def get(self,d):return values/np.float32(10000)
    view=stage2.StackView(Probe(),'fixed');A=view.get('probe');assert np.allclose(A,values,atol=1e-9,rtol=1e-6,equal_nan=True)
    F=fe.compute_feature_base(A[0,:3],'fixed');expected=np.maximum(values[0,:3],np.float32(1e-4))
    assert np.allclose(F[:,8:16],np.log(expected.astype(float)),atol=2e-6,rtol=0)
    assert np.allclose(F[:,16],expected.sum(1),atol=1e-8,rtol=1e-6)
    probe=Probe();probe.on_disk_scale=1;rejected=False
    try:stage2.StackView(probe,'fixed')
    except UnitError:rejected=True
    assert rejected,'Already-reflectance input was accepted for double conversion'
    class Certain:
        def predict(self,X):return np.ones(len(X),int),np.tile([1.,0.],(len(X),1))
    rows,maps=stage2.classify_all_dates(view,np.ones((1,4),bool),np.ones(8,np.float32)/8,np.arange(1,9,dtype=np.float32)/36,Certain(),Certain(),stage2.Stage2Options(mode='fixed'),return_maps=True)
    assert maps['probe'][0,3]==0 and np.count_nonzero(maps['probe'])==3
    product=[]
    for d in io.load_glint_stack().dates:
        actual=stage2.StackView(io.load_glint_stack(),'fixed').get(d);ref=pp.glint_reflectance(d);reef=pp.reef_mask();fin=reef&np.isfinite(actual).all(2)&np.isfinite(ref).all(2)
        product.append(dict(date=d,max_abs_vs_independent_GLINT=float(np.max(np.abs(actual[fin].astype(float)-ref[fin]))),n=int(fin.sum()),floor_reflectance=1e-4))
    save(OUT/'manifest/unit_checks.json',dict(known_positive=True,below_floor=True,at_floor=True,invalid_row_excluded=True,already_scaled_input_rejected=True,conversion_count=1,log_tolerance_absolute=2e-6,product_comparisons=product))
    log('T14 corrected-scale boundary and invalid-row checks passed')

def prepare_stage2_inputs():
    if (OUT/'stage2/prepared.json').exists():return
    unit_checks();common=None;audit={}
    for policy in ('pooled','per_date'):
        for seed in SEEDS:
            tag=f'{policy}_seed{seed}';support=np.load(OUT/f'stage1/{policy}_mask_seed{seed}_masks.npz')['maskFinal'];z=prepare_support(support)
            assert set(z['Ytest'])=={1,2,3} and all(z['n_train'][c]>0 for c in ('HPY','PLE','BLD'))
            kk=set(z['test_keys']);common=kk if common is None else common&kk
            audit[tag]=dict(n_train=z['n_train'],n_test=len(z['Ytest']),mu_h=z['mu_h'],mu_b=z['mu_b'],thresholds=z['thresholds'],ple_stats=z['ple_stats'],split_stats=z['split_stats'],train_keys=z['train_keys'],test_keys=z['test_keys'],training_polygons={c:[f['properties'] for f in fs] for c,fs in z['training_polygons'].items()},test_polygons=[f['properties'] for f in z['test_polygons']])
            z.pop('view');z.pop('training_polygons');z.pop('test_polygons');dump(OUT/f'stage2/{tag}_inputs.pkl',z)
            log(f'Prepared Stage2 {tag}: train {z["n_train"]}, test {len(z["Ytest"])}')
    common=sorted(common);coords=lambda k:tuple(k.split(':')[1:3]);trcoords={coords(k) for a in audit.values() for kk in a['train_keys'].values() for k in kk};clean=[k for k in common if coords(k) not in trcoords]
    assert {k.rsplit(':',1)[1] for k in common}=={'1','2','3'} and {k.rsplit(':',1)[1] for k in clean}=={'1','2','3'}
    save(OUT/'stage2/sample_audit.json',audit);save(OUT/'stage2/prepared.json',dict(common_keys=common,location_disjoint_common_keys=clean,n_common=len(common),n_location_disjoint=len(clean),fixed_before_fitting=True))

def fit2(tag,folder,seed):
    dest=folder/f'{tag}.json'
    if dest.exists():return tag+' reused'
    t=time.time();z=undump(folder/f'{tag}_inputs.pkl');z['view']=stage2.StackView(io.load_glint_stack(),'fixed',check_mask=z['support'])
    if not all(z['n_train'][c]>0 for c in ('HPY','PLE','BLD')):
        save(dest,dict(tag=tag,seed=seed,evaluable=False,n_train=z['n_train'],reason='Required class absent; no retuning'));return tag+' UNEVALUABLE'
    g,s,ts=train_stage2(z,seed,mode_rows='numpy',engine='matlab_cart');r=stage2.run_frozen(g,s,z,return_maps=True)
    save_forest(g,folder/f'{tag}_gate.npz');save_forest(s,folder/f'{tag}_stress.npz')
    np.savez_compressed(folder/f'{tag}_predictions.npz',**r['maps'],test_pred=r['pred'],test_y=z['Ytest'],test_keys=z['test_keys'],gate_pred=r['gate_pred'],stress_pred=r['stress_pred'])
    def F(X):return fe.build_hybrid_features(fe.compute_feature_base(X,'fixed'),z['mu_h'],z['mu_b'],z['v'],z['vnorm2'])
    gp=g.predict(F(z['Xtest']))[0];sp=s.predict(F(z['Xtest']))[0];st=z['Ytest']!=1
    result=dict(tag=tag,seed=seed,evaluable=True,seconds=time.time()-t,metrics=r['metrics'],accuracy_with_UNK=float(np.mean(r['pred']==z['Ytest'])),n_train=z['n_train'],n_test=len(z['Ytest']),table=r['date_summary'],gate_metrics=compute_metrics3(np.where(z['Ytest']==1,1,2),gp),stress_metrics=compute_metrics3(z['Ytest'][st]-1,sp[st]),feature_hashes={k:ah(v) for k,v in ts.items() if isinstance(v,np.ndarray)},drop_train_dates=z['drop_train_dates'])
    save(dest,result);return f'{tag}: {time.time()-t:.1f}s accuracy={result["accuracy_with_UNK"]:.4f}'
def stage2_job(task):
    policy,seed=task;return fit2(f'{policy}_seed{seed}',OUT/'stage2',seed)
def run_stage2():
    prepare_stage2_inputs();parallel(stage2_job,[(p,s) for p in ('pooled','per_date') for s in SEEDS],WORKERS)

def run_omissions():
    policy=load(OUT/'manifest/selection.json')['selected_candidate'];audit={}
    for seed in SEEDS:
        support=np.load(OUT/f'stage1/{policy}_mask_seed{seed}_masks.npz')['maskFinal']
        for d in JUNE:
            tag=f'{policy}_omit{d}_seed{seed}';p=OUT/f'omissions/{tag}_inputs.pkl'
            if p.exists():continue
            z=prepare_support(support,(d,));audit[tag]=dict(heldout_date=d,n_train=z['n_train'],training_dates={c:dict(zip(*np.unique(dd,return_counts=True))) for c,dd in z['Dtrain'].items()},test_keys=z['test_keys'],train_keys=z['train_keys'],mu_h=z['mu_h'],mu_b=z['mu_b'],thresholds=z['thresholds'],ple_stats=z['ple_stats'],split_stats=z['split_stats'],held_date_absent_from_training=True,cached_references_reused=False)
            z.pop('view');z.pop('training_polygons');z.pop('test_polygons');dump(p,z)
    if audit:save(OUT/'omissions/input_audit.json',audit)
    parallel(omission_job,[(policy,d,s) for d in JUNE for s in SEEDS],WORKERS)
def omission_job(task):
    p,d,s=task;return fit2(f'{p}_omit{d}_seed{s}',OUT/'omissions',s)

def prepare_loo():
    if (OUT/'loo/inputs.npz').exists():return
    policy=load(OUT/'manifest/selection.json')['selected_candidate'];z=np.load(OUT/'stage1/inputs.npz');train=zrows(z,'accuracy');X=apply_rows(train,policy)
    # Explicit revised composite, not an assertion about unrecovered Phase6 code.
    with np.errstate(all='ignore'):
        import warnings
        with warnings.catch_warnings():warnings.simplefilter('ignore',RuntimeWarning);comp=np.nanmedian(np.stack([date_features(d,policy) for d in BASE]),axis=0)
    tr=z['train_mask'];lab=z['train_labels'];rr,cc=find_colmajor(tr&np.isfinite(comp).all(2)&np.isin(lab,[2,3]));S=comp[rr,cc];Y=lab[rr,cc];r0=io.load_glint_stack().ref
    records={};polys=[];keys=(train['R'].astype(np.int64)*r0.width+train['C'])
    for i,f in enumerate(io.load_geojson(str(paths.STAGE1_LABELS['train']))):
        if f['properties'].get('Class',f['properties'].get('class'))!='HC':continue
        r,c=feature_pixels(f,r0,rule='stage1');pix=np.unique(r*r0.width+c);ix=(train['Y']==1)&np.isin(keys,pix);ok=np.isfinite(comp[r,c]).all(1);V=comp[r[ok],c[ok]]
        pid=str(f['properties'].get('id',i));records[f'held_{pid}']=ix;records[f'V_{pid}']=V;polys.append(dict(id=pid,n_test=len(V),held_training_n=int(ix.sum()),geometry_sha=hashlib.sha256(json.dumps(f['geometry'],sort_keys=True).encode()).hexdigest()))
    np.savez_compressed(OUT/'loo/inputs.npz',X=X,Y=train['Y'],S=S,SY=Y,**records);save(OUT/'loo/design.json',dict(policy=policy,polygons=polys,trees=50,cost=[[0,3,3],[1,0,1],[1,1,0]],baseline_dates=BASE,composite='pixel-wise nanmedian of coefficient-adjusted feature arrays; explicit new design',noncoral_cap='min(3 * retained coral training rows, available baseline noncoral rows) per class',calibration_exclusion='External sand polygons exclude final-test tiles. Sand reference is independent of HC holdout; no held habitat labels used for calibration.',scope='Polygon holdout across all training dates; not leave-date-out or external-site validation.'))
def loo_job(task):
    pid,seed=task;p=OUT/f'loo/polygon{pid}_seed{seed}.json'
    if p.exists():return f'LOO {pid}/{seed} reused'
    t=time.time();z=np.load(OUT/'loo/inputs.npz');keep=~z[f'held_{pid}'];coral=keep&(z['Y']==1);Xc=z['X'][coral];n=len(Xc);V=z[f'V_{pid}']
    if not len(V) or not n:save(p,dict(polygon_id=pid,seed=seed,evaluable=False,reason='No held evaluation or remaining coral rows'));return f'LOO {pid} unevaluable'
    rng=np.random.default_rng(seed);xx=[Xc];yy=[np.ones(n,np.uint8)]
    for k in (2,3):
        pool=z['S'][z['SY']==k];take=rng.choice(len(pool),min(n*3,len(pool)),replace=False);xx.append(pool[take]);yy.append(np.full(len(take),k,np.uint8))
    f=grow_forest(np.concatenate(xx),np.concatenate(yy),n_trees=50,seed=seed,classes=(1,2,3),cost=np.array([[0,3,3],[1,0,1],[1,1,0]]))
    pr=f.predict(V)[0];save_forest(f,OUT/f'loo/polygon{pid}_seed{seed}_forest.npz');np.savez_compressed(OUT/f'loo/polygon{pid}_seed{seed}_predictions.npz',pred=pr,y=np.ones(len(V),np.uint8));save(p,dict(polygon_id=pid,seed=seed,evaluable=True,n_train=sum(map(len,xx)),n_test=len(V),coral_recall=float(np.mean(pr==1)),seconds=time.time()-t))
    return f'LOO {pid}/{seed}: {time.time()-t:.1f}s recall={np.mean(pr==1):.3f}'
def run_loo():
    prepare_loo();polys=load(OUT/'loo/design.json')['polygons'];parallel(loo_job,[(p['id'],s) for s in SEEDS for p in polys],WORKERS)

def preflight(config,cache=None):
    # All identity checks precede creation of scientific outputs.
    for rel,expected in config.get('code_hashes',{}).items():
        # vendor/coralremote/* entries of the run configuration name this repository's coralremote/* (byte-identical).
        # The hashes of the two entry-point scripts are those of the released (path-templated) files.
        # Code hashes use LF line endings, so Windows (CRLF) and Unix checkouts give the same identity.
        f=REPO/rel[len('vendor/'):] if rel.startswith('vendor/') else HERE/rel
        if hashlib.sha256(f.read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=expected:raise ValueError(f'Release code identity mismatch: {rel}')
    for rec in config['input_hashes']:
        p=Path(expand(rec['path']))
        if not p.exists():raise FileNotFoundError(f'Missing required licensed/prepared input: {p}')
        if sha(p)!=rec['sha256']:raise ValueError(f'Input identity mismatch: {p}')
    if cache:
        os.environ['CORALREMOTE_CACHE']=str(cache)
        for path,var in ((paths.GLINT_STACK,'stack_glint'),(paths.RAW_STACK,'stack_raw')):
            st=io.StackFile(path,var,on_disk_scale=10000 if var=='stack_glint' else 1);cp=st._cache_path();meta=io.cache_meta(st.path,st.var,st._expected_shape(),st.dates)
            if (cp.exists() or cp.with_suffix('.meta.json').exists()) and not io.cache_is_valid(cp,meta):raise ValueError(f'Inconsistent cache metadata/content: {cp}; remove or regenerate in a separate clean cache')
    log('Release input/cache preflight passed before scientific outputs')

def primary(config_path,cache):
    config=load(config_path);preflight(config,cache);freeze();calibrate();prepare_stage1();policy=config['selected_policy'];seed=config['seed']
    parallel(stage1_job,[(policy,'mask',seed),(policy,'accuracy',seed)],2)
    tag=f'{policy}_seed{seed}';support=np.load(OUT/f'stage1/{policy}_mask_seed{seed}_masks.npz')['maskFinal'];z=prepare_support(support);z.pop('view');z.pop('training_polygons');z.pop('test_polygons');dump(OUT/f'stage2/{tag}_inputs.pkl',z);log(fit2(tag,OUT/'stage2',seed))
    save(OUT/'PRIMARY_COMPLETE.json',dict(policy=policy,seed=seed,config_sha256=sha(config_path),runner_sha256=sha(__file__),mask_n=int(support.sum()),entry_point='primary',scientific_output_root=str(OUT)))

def main():
    global OUT,WORKERS
    ap=argparse.ArgumentParser();ap.add_argument('command');ap.add_argument('--out');ap.add_argument('--workers',type=int,default=WORKERS);ap.add_argument('--config');ap.add_argument('--cache');a=ap.parse_args()
    if a.out:OUT=Path(a.out);os.environ['REVISION_OUT']=str(OUT)
    WORKERS=a.workers;os.environ['REVISION_WORKERS']=str(WORKERS)
    if a.command=='freeze':freeze()
    elif a.command=='calibrate':calibrate()
    elif a.command=='prepare':prepare_stage1()
    elif a.command=='stage1':parallel(stage1_job,[(p,k,s) for p in ('pooled','per_date') for s in SEEDS for k in ('mask','accuracy','cv0','cv1','cv2')],WORKERS)
    elif a.command=='select':select()
    elif a.command=='stage2':run_stage2()
    elif a.command=='omissions':run_omissions()
    elif a.command=='loo':run_loo()
    elif a.command=='primary':primary(a.config,a.cache)
    else:raise ValueError(a.command)
if __name__=='__main__':
    try:main()
    except Exception:
        save(OUT/f'failures/command_{time.time_ns()}.json',dict(command=sys.argv,traceback=traceback.format_exc()));raise
