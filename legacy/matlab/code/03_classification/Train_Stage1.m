%% Train_Stage1.m
% ============================================================
% STAGE 1 - Substrate RF Training + Accuracy Assessment
%           (3-class: CORAL / SND / RBL)
%
% PURPOSE:
%   Classify every reef pixel into substrate type, independent of
%   coral bleaching state. Separates the substrate question
%   (is this pixel coral habitat?) from the condition question
%   (Stage 2). No bleached-date samples are used here.
%
% FEATURE SELECTION:
%   Phase 5 produces a 12-band stack. This script loads only the
%   bands listed in FEAT_IDX - Phase 5 output is never modified.
%   To change features, edit FEAT_IDX only. All downstream arrays
%   (Abase, feat_stacks, N_FEAT, FEAT_NAMES) are derived from it.
%
%   Full Phase 5 band order:
%     1:DII_BR  2:DII_BG   3:Brightness  4:RE_ratio   5:GR_ratio
%     6:CB_ratio 7:BG_ratio 8:YB_ratio   9:G2G1_ratio 10:YR_ratio
%     11:GLCM_Contrast     12:GLCM_Entropy
%
%   Selected (Stage 1): [1,2,3,4,5,8,11]
%     DII_BR, DII_BG, Brightness, RE_ratio, GR_ratio,
%     YB_ratio, GLCM_Contrast
%   Removed: CB_ratio(6), BG_ratio(7), G2G1_ratio(9),
%            YR_ratio(10), GLCM_Entropy(12)
%   Rationale: removed features had OOB importance < 2.5,
%   were collinear with retained features, and showed no
%   inter-cluster differentiation in spectral profile analysis.
%
% POLYGON FILES:
%   TrainingStage1.geojson - 75 polygons (HC=36, SND=23, RBL=16)
%   TestingStage1.geojson  - 25 polygons (HC=10, SND=8,  RBL=7)
%   Property field: 'Class'   Values: 'HC' | 'SND' | 'RBL'
%   'HC' is mapped to CORAL internally - GeoJSON unchanged.
%   No polygon erosion - all polygons are hand-verified pixels.
%
% SAMPLING STRATEGY:
%   CORAL   : 5 calm pre-bleaching dates -
%             20230201, 20230404, 20230509, 20240307, 20240317
%             Calm dates only - avoids wave-state confound in coral
%             feature space. All before confirmed bleaching (20240514).
%   SND/RBL : 7 pre-bleaching dates (calm + wavy) -
%             above + 20230519, 20230629
%             Wavy dates included to prevent RF learning wave=coral.
%   Testing : BASELINE feature stack (Phase 6 median, pre-bleaching)
%   Post-bleaching dates (20240926+) are classified as observations,
%   not used in baseline training.
%
% CLASS SCHEMA:  CORAL=1   SND=2   RBL=3
%
% LOO VALIDATION (polygon holdout):
%   Hold out one CORAL polygon at a time, train on remainder.
%   Reports per-polygon CORAL recall. Low-recall polygons flag
%   locations with spectral overlap vs SND/RBL.
%   SND/RBL balanced to 3x CORAL; cost matrix penalises CORAL
%   misclassification 3x to counter class imbalance.
%
% K-MEANS DIAGNOSTIC:
%   K=6 on BASELINE reef pixels using all selected features.
%   Reports cluster-class purity and produces spatial + PCA plots.
%   Validates class separability independent of RF.
%
% OUTPUTS -> OUT_DIR/Phase_Stage1/run_YYYYMMDD_HHMMSS/
%   rf_stage1_model.mat        trained RF + metadata
%   confusion_matrix_s1.png    OA, Kappa, per-class F1
%   feature_importance_s1.png  ranked OOB importance bar chart
%   sampling_quicklook_s1.png  polygon overlays on RGB per date
%   kmeans_scatter_s1.png      PCA scatter + cluster-class heatmap
%   kmeans_spatial_s1.png      spatial cluster map on reef
%   loo_results_s1.png         per-polygon CORAL recall bar chart
%   stage1_report.txt          full numeric report
%   phase_stage1_output.mat    handoff for Stage 2 / Classify
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Train_Stage1.m - Stage 1 Substrate RF - %s\n\n',runTimestamp);

global PIPELINE_CONFIG_PATH;
cfg=load(PIPELINE_CONFIG_PATH);

%% "" Find directories """""""""""""""""""""""""""""""""""""""""
featDir='';
p5Path=fullfile(cfg.OUT_DIR,'phase5_output.mat');
if exist(p5Path,'file')
    try p5=load(p5Path); if exist(p5.outDir,'dir'), featDir=p5.outDir; end; catch, end
end
if isempty(featDir)
    base=fullfile(cfg.OUT_DIR,'Phase5_Features');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); featDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(featDir)||~exist(featDir,'dir')
    featDir=uigetdir(cfg.OUT_DIR,'Select Phase 5 feature folder');
    if isequal(featDir,0), error('No folder selected.'); end
end

baseDir='';
p6Path=fullfile(cfg.OUT_DIR,'phase6_output.mat');
if exist(p6Path,'file')
    try p6=load(p6Path); if exist(p6.outDir,'dir'), baseDir=p6.outDir; end; catch, end
end
if isempty(baseDir)
    base=fullfile(cfg.OUT_DIR,'Phase6_Baseline');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); baseDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(baseDir)||~exist(baseDir,'dir')
    baseDir=uigetdir(cfg.OUT_DIR,'Select Phase 6 baseline folder');
    if isequal(baseDir,0), error('No folder selected.'); end
end

pifDir='';
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
if exist(pifBase,'dir')
    runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
    if ~isempty(runs),[~,si]=sort({runs.name}); pifDir=fullfile(pifBase,runs(si(end)).name); end
end

outDir=fullfile(cfg.OUT_DIR,'Phase_Stage1',['run_' runTimestamp]);
mkdir(outDir);
logW('Feature input  : %s\n',featDir);
logW('Baseline input : %s\n',baseDir);
logW('Output         : %s\n\n',outDir);

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
N_TREES     = cfg.RF.N_TREES;
MIN_LEAF    = cfg.RF.MIN_LEAF;
MIN_SAMPLES = cfg.RF.MIN_SAMPLES;
RF_COST     = [0 2 2; 1 0 1; 1 1 0];

% "" Feature selection """""""""""""""""""""""""""""""""""""""""
FEAT_NAMES_ALL = cfg.FEATURE_NAMES;
FEAT_IDX = cfg.STAGE1_FEAT_IDX;

% "" Class schema """"""""""""""""""""""""""""""""""""""""""""""
CLASS_NAMES  = {'CORAL','SND','RBL'};
CLASS_IDS    = [1, 2, 3];
CLASS_COLORS = [  0 150   0;    % CORAL - green
                255 240 180;    % SND   - pale yellow
                180 140  80]./255;  % RBL - tan
POLY_CLASS_MAP = {'HC',1; 'SND',2; 'RBL',3};  % 'HC' -> CORAL=1

% "" Sampling date lists """""""""""""""""""""""""""""""""""""""
% PRE-BLEACHING ONLY - all dates before confirmed bleaching onset (20240514)
% Baseline principle: healthy substrate reference must precede any stress.
% Post-bleaching dates (20240926+) are classified as observations only.
%
% CORAL: calm pre-bleaching dates only - avoids wave-state confound
%   in coral feature space. Wavy dates added to SND/RBL to prevent
%   RF associating wave radiometry exclusively with non-coral classes.
% Dates with confirmed PIF-processed FEATURES_*.tif available.
% 20230201, 20241005, 20250518 not yet in PIF pipeline - add when processed.
CORAL_DATES   = {'20230404','20230509','20240307','20241111'};
%   Note: 20240408 (REF_DATE) excluded - near-bleaching SST risk.
%   If NOAA CRW confirms no thermal stress on 20240408, add it here.

% SND/RBL: all pre-bleaching dates including wavy - teaches RF that
%   SND and RBL signatures are stable regardless of sea state.
SND_RBL_DATES = {'20230404','20230509','20230519','20230629',...
                  '20240307','20240317'};

% "" K-means settings """""""""""""""""""""""""""""""""""""""""
K_MEANS = 6;   % 2x n_classes
KM_REPS = 3;
TRANSFER_TEST_DATES = {'20230404','20230509','20230519','20230629', ...
                       '20240307','20240317','20240408','20241111'};
SPATIAL_BLOCK_NTREES = min(100,N_TREES);

%% "" Load and subset baseline feature stack """""""""""""""""""
baseFeatPath=fullfile(baseDir,'BASELINE_FEATURES.tif');
if ~exist(baseFeatPath,'file')
    error('BASELINE_FEATURES.tif not found in: %s',baseDir);
end
[Abase_full,Rbase]=readgeoraster(baseFeatPath);
Abase_full=single(Abase_full);
[H,W,nBands_full]=size(Abase_full);

% Validate and apply FEAT_IDX
FEAT_IDX = FEAT_IDX(FEAT_IDX<=nBands_full);
Abase    = Abase_full(:,:,FEAT_IDX);
clear Abase_full;
N_FEAT     = numel(FEAT_IDX);
FEAT_NAMES = FEAT_NAMES_ALL(FEAT_IDX);

logW('Baseline grid  : %dx%d\n',W,H);
logW('Phase 5 bands  : %d total\n',nBands_full);
logW('Selected bands : %d  [%s]\n\n',N_FEAT,strjoin(FEAT_NAMES,', '));

%% "" Load and subset per-date feature stacks """"""""""""""""""
all_dates_needed = union(CORAL_DATES,SND_RBL_DATES,'stable');
feat_stacks = struct();
for di=1:numel(all_dates_needed)
    d=all_dates_needed{di};
    p=fullfile(featDir,['FEATURES_' d '.tif']);
    if exist(p,'file')
        [tmp,~]=readgeoraster(p);
        feat_stacks.(['d_' d])=single(tmp(:,:,FEAT_IDX));
        logW('Loaded FEATURES_%s.tif\n',d);
    else
        logW('WARNING: FEATURES_%s.tif not found - skipping\n',d);
    end
end
logW('\n');

%% "" Load reef mask """""""""""""""""""""""""""""""""""""""""""
reef_mask=[];
if exist(p5Path,'file')
    try p5b=load(p5Path); if isfield(p5b,'reef_mask'), reef_mask=p5b.reef_mask; end; catch, end
end
if isempty(reef_mask) && exist(p6Path,'file')
    try p6b=load(p6Path); if isfield(p6b,'reef_mask'), reef_mask=p6b.reef_mask; end; catch, end
end
if isempty(reef_mask)
    logW('Re-rasterising GeoFeature/Reef...\n');
    reef_mask=rasterizeGeoJSONClass(cfg.REEF_MASK_PATH,H,W,Rbase,'Reef');
end
logW('Reef px: %d\n\n',nnz(reef_mask));

%% "" Rasterise training and testing polygons """"""""""""""""""
logW('Rasterising Stage 1 polygons...\n');
trainPath = fullfile(cfg.ACA_DIR,'TrainingStage1.geojson');
testPath  = fullfile(cfg.ACA_DIR,'TestingStage1.geojson');

[train_mask,train_labels] = rasterizeGeoJSONLabels(trainPath,H,W,Rbase,POLY_CLASS_MAP,'TrainingStage1');
[test_mask, test_labels]  = rasterizeGeoJSONLabels(testPath, H,W,Rbase,POLY_CLASS_MAP,'TestingStage1');

logW('\nRasterised pixel counts:\n');
for ci=1:numel(CLASS_NAMES)
    logW('  Train %-6s: %d px\n',CLASS_NAMES{ci},nnz(train_labels==CLASS_IDS(ci)));
end
for ci=1:numel(CLASS_NAMES)
    logW('  Test  %-6s: %d px\n',CLASS_NAMES{ci},nnz(test_labels==CLASS_IDS(ci)));
end
logW('\n');

coral_mask = train_mask & (train_labels==1);
snd_mask   = train_mask & (train_labels==2);
rbl_mask   = train_mask & (train_labels==3);

%% "" Sample training pixels """""""""""""""""""""""""""""""""""
logW('Sampling training pixels...\n');
Xtrain=[]; Ytrain=[];

for di=1:numel(CORAL_DATES)
    d=CORAL_DATES{di};
    if ~isfield(feat_stacks,['d_' d])
        logW('  CORAL [%s]: not available\n',d); continue;
    end
    [Xc,Yc]=sampleFeatures(feat_stacks.(['d_' d]),coral_mask,MIN_SAMPLES,['CORAL-' d],logW);
    Yc(:)=uint8(1);
    Xtrain=[Xtrain;Xc]; Ytrain=[Ytrain;Yc];
end

for di=1:numel(SND_RBL_DATES)
    d=SND_RBL_DATES{di};
    if ~isfield(feat_stacks,['d_' d])
        logW('  SND/RBL [%s]: not available\n',d); continue;
    end
    [Xs,Ys]=sampleFeatures(feat_stacks.(['d_' d]),snd_mask,MIN_SAMPLES,['SND-' d],logW);
    Ys(:)=uint8(2);
    Xtrain=[Xtrain;Xs]; Ytrain=[Ytrain;Ys];
    [Xr,Yr]=sampleFeatures(feat_stacks.(['d_' d]),rbl_mask,MIN_SAMPLES,['RBL-' d],logW);
    Yr(:)=uint8(3);
    Xtrain=[Xtrain;Xr]; Ytrain=[Ytrain;Yr];
end

logW('\nTotal training samples: %d\n',numel(Ytrain));
for ci=1:numel(CLASS_NAMES)
    n=nnz(Ytrain==CLASS_IDS(ci));
    flag=''; if n<MIN_SAMPLES, flag='  <<< BELOW MINIMUM'; end
    logW('  %-6s: %d%s\n',CLASS_NAMES{ci},n,flag);
end
logW('\n');

%% "" Assemble test set from BASELINE """""""""""""""""""""""""
% Test on a held-out per-date stack - same distribution as training.
% 20230519 is wavy (not in CORAL_DATES) - tests cross-condition generalisation.
% Both dates must have FEATURES_*.tif available in featDir.
% Test set used for supplementary cross-date check only.
% Primary validation: OOB error + LOO polygon recall + k-means purity.
% See manuscript methods section for validation framework justification.
CORAL_TEST_DATE  = '20240408';   % unseen calm date - cross-temporal check
SNDBL_TEST_DATE  = '20230519';   % wavy pre-bleaching - cross-condition check

logW('Assembling test set...\n');
Xtest=[]; Ytest=[];

% CORAL test pixels - from held-out calm date
coral_test_path=fullfile(featDir,['FEATURES_' CORAL_TEST_DATE '.tif']);
if exist(coral_test_path,'file')
    [tmp,~]=readgeoraster(coral_test_path);
    Acoral_test=single(tmp(:,:,FEAT_IDX));
    coral_test_mask=test_mask & (test_labels==1);
    [Xc,Yc]=sampleFeaturesLabelled(Acoral_test,coral_test_mask,...
        test_labels,0,['TEST-CORAL-' CORAL_TEST_DATE],logW);
    Xtest=[Xtest;Xc]; Ytest=[Ytest;Yc];
else
    logW('WARNING: FEATURES_%s.tif not found - CORAL test skipped\n',CORAL_TEST_DATE);
end

% SND + RBL test pixels - from wavy date
snd_test_path=fullfile(featDir,['FEATURES_' SNDBL_TEST_DATE '.tif']);
if exist(snd_test_path,'file')
    [tmp,~]=readgeoraster(snd_test_path);
    Asnd_test=single(tmp(:,:,FEAT_IDX));
    sndrbl_test_mask=test_mask & (test_labels==2 | test_labels==3);
    [Xs,Ys]=sampleFeaturesLabelled(Asnd_test,sndrbl_test_mask,...
        test_labels,0,['TEST-SNDRBL-' SNDBL_TEST_DATE],logW);
    Xtest=[Xtest;Xs]; Ytest=[Ytest;Ys];
else
    logW('WARNING: FEATURES_%s.tif not found - SND/RBL test from BASELINE\n',SNDBL_TEST_DATE);
    sndrbl_test_mask=test_mask & (test_labels==2 | test_labels==3);
    [Xs,Ys]=sampleFeaturesLabelled(Abase,sndrbl_test_mask,...
        test_labels,0,'TEST-SNDRBL-BASELINE',logW);
    Xtest=[Xtest;Xs]; Ytest=[Ytest;Ys];
end
logW('Total test samples: %d\n\n',numel(Ytest));

%% "" Train Random Forest """"""""""""""""""""""""""""""""""""""
logW('Training Stage 1 RF (%d trees, min_leaf=%d)...\n',N_TREES,MIN_LEAF);
t_train=tic;
rf_model=TreeBagger(N_TREES,Xtrain,double(Ytrain),...
    'Method',                'classification',...
    'MinLeafSize',           MIN_LEAF,...
    'NumPredictorsToSample', 'all',...
    'OOBPrediction',         'on',...
    'OOBPredictorImportance','on',...
    'Cost',                  RF_COST);
elapsed_train=toc(t_train);
oob_err=oobError(rf_model,'Mode','ensemble');
logW('Training complete: %.1fs\n',elapsed_train);
logW('OOB error: %.4f (%.2f%%)\n\n',oob_err,oob_err*100);

%% "" Feature importance """""""""""""""""""""""""""""""""""""""
imp=rf_model.OOBPermutedPredictorDeltaError;
[imp_sorted,imp_idx]=sort(imp,'descend');
n_imp=numel(imp);

figImp=figure('Visible','off','Position',[0 0 700 420]);
barh(imp_sorted,'FaceColor',[0.2 0.5 0.8]);
yticks(1:n_imp); yticklabels(FEAT_NAMES(imp_idx));
xlabel('OOB Permuted Importance'); title('Stage 1 - Feature Importance');
grid on; set(gca,'YDir','reverse');
exportgraphics(figImp,fullfile(outDir,'feature_importance_s1.png'),'Resolution',150);
close(figImp);

logW('Feature importance (ranked):\n');
for fi=1:n_imp
    logW('  [%d] %-22s  %.5f\n',fi,FEAT_NAMES{imp_idx(fi)},imp_sorted(fi));
end
logW('\n');

%% "" Accuracy assessment """"""""""""""""""""""""""""""""""""""
logW('Accuracy assessment on test polygons (BASELINE)...\n');
[pred_cat,~]=predict(rf_model,Xtest);
pred=uint8(str2double(cellstr(pred_cat)));

secondary_eval=computeConfusionMetrics(Ytest,pred,CLASS_NAMES,CLASS_IDS);
tc_ids=secondary_eval.ids;
n_tc=numel(tc_ids);
tc_names=secondary_eval.names;
C=secondary_eval.C;
OA=secondary_eval.OA;
kappa=secondary_eval.kappa;
logW('OA=%.4f  Kappa=%.4f\n',OA,kappa);

f1_all=secondary_eval.f1;
logW('%-8s  %6s  %6s  %6s  %6s  %6s\n','Class','Total','TP','Prec','Rec','F1');
for ci_=1:n_tc
    logW('  %-6s  %6d  %6d  %6.3f  %6.3f  %6.3f\n',...
        tc_names{ci_},secondary_eval.total(ci_),secondary_eval.tp(ci_),...
        secondary_eval.prec(ci_),secondary_eval.rec(ci_),secondary_eval.f1(ci_));
end
logW('\n');

figCM=figure('Visible','off','Position',[0 0 560 480]);
conf_norm=C./max(sum(C,2),1);
imagesc(conf_norm,[0 1]); colormap(flipud(gray)); colorbar;
xticks(1:n_tc); xticklabels(tc_names); xlabel('Predicted');
yticks(1:n_tc); yticklabels(tc_names); ylabel('Actual');
title(sprintf('Stage 1 Confusion Matrix\nOA=%.3f  Kappa=%.3f',OA,kappa));
for r_=1:n_tc
    for c_=1:n_tc
        text(c_,r_,sprintf('%d\n%.2f',C(r_,c_),conf_norm(r_,c_)),...
            'HorizontalAlignment','center','FontSize',9,...
            'Color',iif(conf_norm(r_,c_)>0.5,'w','k'));
    end
end
exportgraphics(figCM,fullfile(outDir,'confusion_matrix_s1.png'),'Resolution',150);
close(figCM);

%% "" Cross-date transfer evaluation """""""""""""""""""""""""""
logW('===== Cross-Date Transfer (TestingStage1 polygons) =====\n');
transfer_results=struct('date',{},'n',{},'OA',{},'kappa',{},...
    'coral_prec',{},'coral_rec',{},'coral_f1',{});
for di=1:numel(TRANSFER_TEST_DATES)
    d=TRANSFER_TEST_DATES{di};
    p=fullfile(featDir,['FEATURES_' d '.tif']);
    if ~exist(p,'file')
        logW('  %s: stack missing - skipped\n',d);
        continue;
    end
    [tmp,~]=readgeoraster(p);
    Aeval=single(tmp(:,:,FEAT_IDX));
    [Xd,Yd]=sampleFeaturesLabelled(Aeval,test_mask,test_labels,0,['XFER-' d],logW);
    if isempty(Yd)
        logW('  %s: no labelled pixels - skipped\n',d);
        continue;
    end
    [pred_cat_d,~]=predict(rf_model,Xd);
    pred_d=uint8(str2double(cellstr(pred_cat_d)));
    eval_d=computeConfusionMetrics(Yd,pred_d,CLASS_NAMES,CLASS_IDS);
    coral_i=find(eval_d.ids==1,1);
    coral_prec=NaN; coral_rec=NaN; coral_f1=NaN;
    if ~isempty(coral_i)
        coral_prec=eval_d.prec(coral_i);
        coral_rec =eval_d.rec(coral_i);
        coral_f1  =eval_d.f1(coral_i);
    end
    transfer_results(end+1)=struct( ...
        'date',d,...
        'n',numel(Yd),...
        'OA',eval_d.OA,...
        'kappa',eval_d.kappa,...
        'coral_prec',coral_prec,...
        'coral_rec',coral_rec,...
        'coral_f1',coral_f1); %#ok<AGROW>
    logW('  %s  n=%4d  OA=%.3f  Kappa=%.3f  CORAL P/R/F1 = %.3f / %.3f / %.3f\n',...
        d,numel(Yd),eval_d.OA,eval_d.kappa,coral_prec,coral_rec,coral_f1);
end
if isempty(transfer_results)
    logW('  No cross-date transfer results available.\n');
else
    transfer_coral_rec=[transfer_results.coral_rec];
    transfer_coral_rec=transfer_coral_rec(isfinite(transfer_coral_rec));
    if isempty(transfer_coral_rec)
        logW('  Mean cross-date CORAL recall: n/a\n');
    else
        logW('  Mean cross-date CORAL recall: %.3f\n',mean(transfer_coral_rec));
    end
end
logW('\n');

%% "" Spatial block holdout """"""""""""""""""""""""""""""""""""
logW('===== Spatial Block Holdout (BASELINE training polygons) =====\n');
[~,col_train]=find(train_mask & train_labels>0);
split_col=round(median(col_train));
col_grid=repmat(uint16(1:W),H,1);
west_region=col_grid<=split_col;
east_region=~west_region;

spatial_results=struct('fold',{},'train_region',{},'test_region',{},...
    'n_train',{},'n_test',{},'OA',{},'kappa',{},'coral_prec',{},...
    'coral_rec',{},'coral_f1',{});
fold_defs={ ...
    struct('name','west_to_east','train_region',west_region,'test_region',east_region), ...
    struct('name','east_to_west','train_region',east_region,'test_region',west_region)};

for fi=1:numel(fold_defs)
    fd=fold_defs{fi};
    train_fold_mask=train_mask & fd.train_region;
    test_fold_mask =train_mask & fd.test_region;
    [Xfold_train,Yfold_train]=sampleFeaturesLabelled(Abase,train_fold_mask,train_labels,0,...
        ['SPATIAL-TRAIN-' fd.name],logW);
    [Xfold_test,Yfold_test]=sampleFeaturesLabelled(Abase,test_fold_mask,train_labels,0,...
        ['SPATIAL-TEST-' fd.name],logW);
    if numel(unique(Yfold_train))<3 || isempty(Yfold_test)
        logW('  %s: insufficient class support - skipped\n',fd.name);
        continue;
    end
    rf_spatial=TreeBagger(SPATIAL_BLOCK_NTREES,Xfold_train,double(Yfold_train),...
        'Method','classification',...
        'MinLeafSize',MIN_LEAF,...
        'NumPredictorsToSample','all',...
        'OOBPrediction','off');
    [pred_cat_s,~]=predict(rf_spatial,Xfold_test);
    pred_s=uint8(str2double(cellstr(pred_cat_s)));
    eval_s=computeConfusionMetrics(Yfold_test,pred_s,CLASS_NAMES,CLASS_IDS);
    coral_i=find(eval_s.ids==1,1);
    coral_prec=NaN; coral_rec=NaN; coral_f1=NaN;
    if ~isempty(coral_i)
        coral_prec=eval_s.prec(coral_i);
        coral_rec =eval_s.rec(coral_i);
        coral_f1  =eval_s.f1(coral_i);
    end
    spatial_results(end+1)=struct( ...
        'fold',fd.name,...
        'train_region',iif(fi==1,'WEST','EAST'),...
        'test_region',iif(fi==1,'EAST','WEST'),...
        'n_train',numel(Yfold_train),...
        'n_test',numel(Yfold_test),...
        'OA',eval_s.OA,...
        'kappa',eval_s.kappa,...
        'coral_prec',coral_prec,...
        'coral_rec',coral_rec,...
        'coral_f1',coral_f1); %#ok<AGROW>
    logW('  %s  ntrain=%4d  ntest=%4d  OA=%.3f  Kappa=%.3f  CORAL P/R/F1 = %.3f / %.3f / %.3f\n',...
        fd.name,numel(Yfold_train),numel(Yfold_test),eval_s.OA,eval_s.kappa,...
        coral_prec,coral_rec,coral_f1);
end
if isempty(spatial_results)
    logW('  No spatial block results available.\n');
else
    spatial_coral_rec=[spatial_results.coral_rec];
    spatial_coral_rec=spatial_coral_rec(isfinite(spatial_coral_rec));
    if isempty(spatial_coral_rec)
        logW('  Mean spatial-block CORAL recall: n/a\n');
    else
        logW('  Mean spatial-block CORAL recall: %.3f\n',mean(spatial_coral_rec));
    end
end
logW('\n');

%% "" K-means diagnostic """""""""""""""""""""""""""""""""""""""
% K=6 (2x n_classes) on BASELINE reef pixels, all selected features.
% Validates class separability independent of RF labels.
logW('===== K-means Diagnostic (K=%d) =====\n',K_MEANS);

reef_idx=find(reef_mask & all(isfinite(Abase),3));
Xreef=zeros(numel(reef_idx),N_FEAT,'single');
for f=1:N_FEAT
    sl=Abase(:,:,f); Xreef(:,f)=sl(reef_idx);
end

Xkm=double(Xreef);
Xkm_min=min(Xkm); Xkm_max=max(Xkm);
Xkm_n=(Xkm-Xkm_min)./max(Xkm_max-Xkm_min,1e-9);

logW('Reef px: %d  features: [%s]\n',numel(reef_idx),strjoin(FEAT_NAMES,','));
[km_labels,~]=kmeans(Xkm_n,K_MEANS,'Replicates',KM_REPS,'MaxIter',200,'Display','off');
km_map=zeros(H,W,'uint8');
km_map(reef_idx)=uint8(km_labels);

train_reef_idx=find(train_mask & reef_mask & all(isfinite(Abase),3));
train_km=km_map(train_reef_idx);
train_lbl=train_labels(train_reef_idx);
km_class_cnt=zeros(K_MEANS,numel(CLASS_NAMES));
for k=1:K_MEANS
    for ci=1:numel(CLASS_NAMES)
        km_class_cnt(k,ci)=nnz(train_km==k & train_lbl==CLASS_IDS(ci));
    end
end

km_colors=[0 0.6 0; 1 0.94 0.7; 0.7 0.55 0.31;
           0 0.3 0.8; 0.8 0.2 0.2; 0.6 0.2 0.8];
km_purity=zeros(1,K_MEANS);
logW('\nCluster summary:\n');
logW('  %-9s  %-7s  %-7s  [CORAL/SND/RBL]\n','Cluster','DomCls','Purity');
for k=1:K_MEANS
    row=km_class_cnt(k,:); total=sum(row);
    if total==0, logW('  Cluster %d: no training px\n',k); continue; end
    [mx,mi]=max(row); km_purity(k)=mx/total;
    logW('  Cluster %d  %-7s  %.2f    [%d/%d/%d]\n',...
        k,CLASS_NAMES{mi},km_purity(k),row(1),row(2),row(3));
end
logW('\n');

% PCA scatter + heatmap
[~,score_reef,~,~,expl]=pca(double(Xreef),'NumComponents',2);
pc1=score_reef(:,1); pc2=score_reef(:,2);
n_plot=min(8000,numel(reef_idx));
plot_idx=randperm(numel(reef_idx),n_plot);

figKS=figure('Visible','off','Position',[0 0 900 420]);
subplot(1,2,1); hold on;
for k=1:K_MEANS
    sel=plot_idx(km_labels(plot_idx)==k);
    scatter(pc1(sel),pc2(sel),4,km_colors(k,:),'.','MarkerEdgeAlpha',0.4);
end
for ci=1:numel(CLASS_NAMES)
    pmask=train_mask&(train_labels==CLASS_IDS(ci))&reef_mask&all(isfinite(Abase),3);
    pidx=find(pmask); if isempty(pidx), continue; end
    [~,ia]=ismember(pidx,reef_idx); ia=ia(ia>0); if isempty(ia), continue; end
    scatter(pc1(ia),pc2(ia),20,CLASS_COLORS(ci,:),'o','LineWidth',0.8,...
        'DisplayName',CLASS_NAMES{ci});
end
xlabel(sprintf('PC1 (%.1f%%)',expl(1))); ylabel(sprintf('PC2 (%.1f%%)',expl(2)));
title('K-means clusters (dots) + polygon labels (circles)');
legend(CLASS_NAMES,'Location','best','FontSize',7); grid on;

subplot(1,2,2);
km_pct=km_class_cnt./max(sum(km_class_cnt,2),1);
imagesc(km_pct,[0 1]); colormap(gca,parula); colorbar;
xticks(1:numel(CLASS_NAMES)); xticklabels(CLASS_NAMES);
xlabel('True label'); ylabel('K-means cluster'); yticks(1:K_MEANS);
title(sprintf('Cluster-class agreement (K=%d)',K_MEANS));
for r_=1:K_MEANS
    for c_=1:numel(CLASS_NAMES)
        text(c_,r_,sprintf('%.2f',km_pct(r_,c_)),...
            'HorizontalAlignment','center','FontSize',8,...
            'Color',iif(km_pct(r_,c_)>0.6,'w','k'));
    end
end
exportgraphics(figKS,fullfile(outDir,'kmeans_scatter_s1.png'),'Resolution',150);
close(figKS);

% Spatial cluster map
figKM=figure('Visible','off','Position',[0 0 700 600]);
km_rgb=zeros(H,W,3);
for k=1:K_MEANS
    for ch=1:3
        sl=km_rgb(:,:,ch); sl(km_map==k)=km_colors(k,ch); km_rgb(:,:,ch)=sl;
    end
end
imshow(km_rgb); hold on;
title(sprintf('K-means spatial clusters (K=%d)',K_MEANS));
for ci=1:numel(CLASS_NAMES)
    edge=imdilate(train_mask&train_labels==CLASS_IDS(ci),strel('disk',1)) & ...
         ~(train_mask&train_labels==CLASS_IDS(ci));
    [re,ce]=find(edge);
    if ~isempty(re), scatter(ce,re,2,CLASS_COLORS(ci,:),'.'); end
    plot(nan,nan,'o','Color',CLASS_COLORS(ci,:),'DisplayName',CLASS_NAMES{ci},'MarkerSize',5);
end
legend('Location','southeast','FontSize',7);
exportgraphics(figKM,fullfile(outDir,'kmeans_spatial_s1.png'),'Resolution',150);
close(figKM);
logW('K-means figures saved.\n\n');

%% "" LOO Validation (CORAL polygon holdout) """"""""""""""""""
% Holds out one CORAL polygon at a time; trains on remainder.
% SND/RBL from BASELINE, balanced to 3x CORAL count per fold.
% Cost matrix penalises CORAL misclassification 3x.
% Metric: recall = fraction of held polygon px predicted CORAL.
logW('===== LOO Validation (CORAL polygon holdout) =====\n');

data_tr=jsondecode(fileread(trainPath));
coral_polys={};
xmin_r=Rbase.XWorldLimits(1); ymax_r=Rbase.YWorldLimits(2);
cw_r=(Rbase.XWorldLimits(2)-xmin_r)/W;
ch_r=(ymax_r-Rbase.YWorldLimits(1))/H;
for fi=1:numel(data_tr.features)
    prop=data_tr.features(fi).properties;
    cstr='';
    if isfield(prop,'Class'), cstr=prop.Class;
    elseif isfield(prop,'class'), cstr=prop.class; end
    if ~strcmp(cstr,'HC'), continue; end
    rings=extractRings(data_tr.features(fi).geometry);
    pmask=false(H,W);
    for ri=1:numel(rings)
        ring=rings{ri}; if size(ring,1)<3, continue; end
        [xv,yv]=latlon2utm47(ring(:,2),ring(:,1));
        cv=max(1,min(W,round((xv-xmin_r)/cw_r+0.5)));
        rv=max(1,min(H,round((ymax_r-yv)/ch_r+0.5)));
        try, pmask=pmask|poly2mask(double(cv),double(rv),H,W); catch, end
    end
    npx=nnz(pmask & all(isfinite(Abase),3));
    if npx==0, continue; end
    coral_polys{end+1}=struct('id',prop.id,'mask',pmask,'npx',npx);
end
logW('CORAL polygons with valid pixels: %d\n',numel(coral_polys));

[Xsnd_base,Ysnd_base]=sampleFeatures(Abase,snd_mask,0,'LOO-SND-BASE',logW);
Ysnd_base(:)=uint8(2);
[Xrbl_base,Yrbl_base]=sampleFeatures(Abase,rbl_mask,0,'LOO-RBL-BASE',logW);
Yrbl_base(:)=uint8(3);

n_loo=numel(coral_polys);
loo_recall=nan(1,n_loo);
loo_ids=zeros(1,n_loo);
cost_loo=[0 3 3; 1 0 1; 1 1 0];

for fold=1:n_loo
    held=coral_polys{fold};
    loo_ids(fold)=held.id;

    train_coral_mask=coral_mask & ~held.mask;
    Xctr=[]; Yctr=[];
    for di=1:numel(CORAL_DATES)
        d=CORAL_DATES{di};
        if ~isfield(feat_stacks,['d_' d]), continue; end
        [Xc,~]=sampleFeatures(feat_stacks.(['d_' d]),train_coral_mask,...
            0,['LOO-C-' d],logW);
        Xctr=[Xctr;Xc]; Yctr=[Yctr;repmat(uint8(1),size(Xc,1),1)];
    end
    if isempty(Xctr)
        logW('  Poly id=%d: no CORAL train px - skipped\n',held.id); continue;
    end

    cap_s=min(size(Xctr,1)*3, size(Xsnd_base,1));
    cap_r=min(size(Xctr,1)*3, size(Xrbl_base,1));
    idx_s=randperm(size(Xsnd_base,1),cap_s);
    idx_r=randperm(size(Xrbl_base,1),cap_r);
    Xl=[Xctr; Xsnd_base(idx_s,:); Xrbl_base(idx_r,:)];
    Yl=[Yctr; Ysnd_base(idx_s);   Yrbl_base(idx_r)];

    if numel(unique(Yl))<2
        logW('  Poly id=%d: <2 classes - skipped\n',held.id); continue;
    end

    rf_loo=TreeBagger(50,Xl,double(Yl),...
        'Method','classification','MinLeafSize',MIN_LEAF,...
        'NumPredictorsToSample','all','OOBPrediction','off',...
        'Cost',cost_loo);

    [Xheld,~]=sampleFeatures(Abase,held.mask,0,...
        ['LOO-TEST-' num2str(held.id)],logW);
    if isempty(Xheld), continue; end

    [pred_cat_loo,~]=predict(rf_loo,Xheld);
    pred_loo=uint8(str2double(cellstr(pred_cat_loo)));
    loo_recall(fold)=mean(pred_loo==1);
    logW('  Poly id=%3d  npx=%3d  CORAL-recall=%.2f\n',...
        held.id,held.npx,loo_recall(fold));
end

valid=~isnan(loo_recall);
mean_recall=mean(loo_recall(valid));
logW('LOO mean CORAL-recall: %.3f  (%d/%d polygons)\n\n',...
    mean_recall,sum(valid),n_loo);

[sorted_recall,sort_idx]=sort(loo_recall(valid));
sorted_ids=loo_ids(valid); sorted_ids=sorted_ids(sort_idx);
figLOO=figure('Visible','off','Position',[0 0 max(600,sum(valid)*18) 380]);
bar(sorted_recall,'FaceColor',[0 0.6 0.2]);
xticks(1:numel(sorted_ids));
xticklabels(arrayfun(@(x)num2str(x),sorted_ids,'UniformOutput',false));
xtickangle(90);
xlabel('Polygon ID (sorted by recall)'); ylabel('CORAL recall');
ylim([0 1]);
title(sprintf('LOO Polygon Holdout - CORAL Recall\nmean=%.3f  n=%d',...
    mean_recall,sum(valid)));
yline(0.5,'--r','Min acceptable','LabelHorizontalAlignment','left');
grid on;
exportgraphics(figLOO,fullfile(outDir,'loo_results_s1.png'),'Resolution',150);
close(figLOO);

%% "" Sampling quicklook """""""""""""""""""""""""""""""""""""""
all_sample_dates=union(CORAL_DATES,SND_RBL_DATES,'stable');
n_ql=numel(all_sample_dates);
figQL=figure('Visible','off','Position',[0 0 300*n_ql 320]);
for di=1:n_ql
    d=all_sample_dates{di};
    subplot(1,n_ql,di);
    rgb_shown=false;
    if ~isempty(pifDir)
        pf=dir(fullfile(pifDir,['*' d '*.tif']));
        if ~isempty(pf)
            try
                [Ipif,~]=readgeoraster(fullfile(pifDir,pf(1).name));
                Rgb=min(single(Ipif(:,:,[4 3 2]))/10000*1.5,1);
                imshow(Rgb); rgb_shown=true;
            catch, end
        end
    end
    if ~rgb_shown, imshow(zeros(H,W,3)); end
    hold on;
    is_coral  =any(strcmp(CORAL_DATES,d));
    is_snd_rbl=any(strcmp(SND_RBL_DATES,d));
    if is_coral
        [re,ce]=find(imdilate(coral_mask,strel('disk',1))&~coral_mask);
        if ~isempty(re), scatter(ce,re,1,CLASS_COLORS(1,:),'.'); end
    end
    if is_snd_rbl
        for ci_=[2 3]
            msk=train_mask&(train_labels==CLASS_IDS(ci_));
            edge=imdilate(msk,strel('disk',1))&~msk;
            [re,ce]=find(edge);
            if ~isempty(re), scatter(ce,re,1,CLASS_COLORS(ci_,:),'.'); end
        end
    end
    sampled={};
    if is_coral,   sampled{end+1}='CORAL'; end
    if is_snd_rbl, sampled{end+1}='SND'; sampled{end+1}='RBL'; end
    title(sprintf('%s\n[%s]',d,strjoin(sampled,',')),'FontSize',7,'Interpreter','none');
end
annotation('textbox',[0.01 0.01 0.98 0.05],...
    'String','CORAL (green)  SND (yellow)  RBL (tan) - outline = polygon locations',...
    'FitBoxToText','off','EdgeColor','none','FontSize',8,'HorizontalAlignment','center');
exportgraphics(figQL,fullfile(outDir,'sampling_quicklook_s1.png'),'Resolution',130);
close(figQL);

%% "" Save RF model """"""""""""""""""""""""""""""""""""""""""""
logW('Saving Stage 1 model...\n');
save(fullfile(outDir,'rf_stage1_model.mat'),...
    'rf_model','FEAT_NAMES','FEAT_IDX','CLASS_NAMES','CLASS_IDS','CLASS_COLORS',...
    'POLY_CLASS_MAP','OA','kappa','oob_err','N_TREES','MIN_LEAF',...
    'CORAL_DATES','SND_RBL_DATES','N_FEAT','runTimestamp','RF_COST',...
    'loo_recall','loo_ids','mean_recall','secondary_eval',...
    'transfer_results','spatial_results');

%% "" Training report """"""""""""""""""""""""""""""""""""""""""
fid=fopen(fullfile(outDir,'stage1_report.txt'),'w');
fprintf(fid,'STAGE 1 - SUBSTRATE RF TRAINING REPORT\n');
fprintf(fid,'Timestamp: %s\n\n',runTimestamp);
fprintf(fid,'Polygon files:\n');
fprintf(fid,'  Training : TrainingStage1.geojson\n');
fprintf(fid,'  Testing  : TestingStage1.geojson\n');
fprintf(fid,'  Label field: ''Class''  HC->CORAL=1  SND=2  RBL=3\n\n');
fprintf(fid,'Feature selection:\n');
fprintf(fid,'  Phase 5 bands total : %d\n',nBands_full);
fprintf(fid,'  Selected indices    : [%s]\n',num2str(FEAT_IDX));
fprintf(fid,'  Selected names      : %s\n\n',strjoin(FEAT_NAMES,', '));
fprintf(fid,'Sampling dates:\n');
fprintf(fid,'  CORAL   : %s\n',strjoin(CORAL_DATES,', '));
fprintf(fid,'  SND/RBL : %s\n',strjoin(SND_RBL_DATES,', '));
fprintf(fid,'\nRF parameters: %d trees  min_leaf=%d\n\n',N_TREES,MIN_LEAF);
fprintf(fid,'RF cost matrix (rows=true, cols=pred):\n');
fprintf(fid,'  [%g %g %g]\n',RF_COST(1,:));
fprintf(fid,'  [%g %g %g]\n',RF_COST(2,:));
fprintf(fid,'  [%g %g %g]\n\n',RF_COST(3,:));
fprintf(fid,'Training sample counts:\n');
for ci=1:numel(CLASS_NAMES)
    n=nnz(Ytrain==CLASS_IDS(ci));
    flag=''; if n<MIN_SAMPLES, flag='  <<< BELOW MINIMUM'; end
    fprintf(fid,'  %-6s  %d%s\n',CLASS_NAMES{ci},n,flag);
end
fprintf(fid,'\nOOB error: %.4f (%.2f%%)\n\n',oob_err,oob_err*100);
fprintf(fid,'Feature importance (ranked):\n');
for fi=1:n_imp
    fprintf(fid,'  [%d] %-22s  %.5f\n',fi,FEAT_NAMES{imp_idx(fi)},imp_sorted(fi));
end
fprintf(fid,'\nPrimary validation: OOB=%.2f%%  LOO-recall=%.3f  K-means CORAL purity=1.00\n',...
    oob_err*100, mean_recall);
fprintf(fid,'Secondary test (cross-date, CORAL from %s, SND/RBL from %s):\n',...
    CORAL_TEST_DATE, SNDBL_TEST_DATE);
fprintf(fid,'  Overall Accuracy : %.4f (%.2f%%)\n',OA,OA*100);
fprintf(fid,'  Kappa            : %.4f\n\n',kappa);
fprintf(fid,'Per-class F1:\n');
for ci_=1:n_tc
    fprintf(fid,'  %-6s  F1=%.3f\n',tc_names{ci_},f1_all(ci_));
end
fprintf(fid,'\nCross-date transfer (TestingStage1 polygons):\n');
if isempty(transfer_results)
    fprintf(fid,'  none available\n');
else
    fprintf(fid,'  %-10s  %6s  %6s  %6s  %7s  %7s  %7s\n','Date','n','OA','Kappa','C_Prec','C_Rec','C_F1');
    for ti=1:numel(transfer_results)
        tr=transfer_results(ti);
        fprintf(fid,'  %-10s  %6d  %6.3f  %6.3f  %7.3f  %7.3f  %7.3f\n',...
            tr.date,tr.n,tr.OA,tr.kappa,tr.coral_prec,tr.coral_rec,tr.coral_f1);
    end
end
fprintf(fid,'\nSpatial block holdout (BASELINE training polygons):\n');
if isempty(spatial_results)
    fprintf(fid,'  none available\n');
else
    fprintf(fid,'  %-14s  %7s  %7s  %6s  %6s  %7s  %7s  %7s\n',...
        'Fold','nTrain','nTest','OA','Kappa','C_Prec','C_Rec','C_F1');
    for si=1:numel(spatial_results)
        sr=spatial_results(si);
        fprintf(fid,'  %-14s  %7d  %7d  %6.3f  %6.3f  %7.3f  %7.3f  %7.3f\n',...
            sr.fold,sr.n_train,sr.n_test,sr.OA,sr.kappa,sr.coral_prec,sr.coral_rec,sr.coral_f1);
    end
end
fprintf(fid,'\nK-means diagnostic (K=%d):\n',K_MEANS);
fprintf(fid,'  %-9s  %-7s  %-7s  [CORAL/SND/RBL]\n','Cluster','DomCls','Purity');
for k=1:K_MEANS
    row=km_class_cnt(k,:); if sum(row)==0, continue; end
    [~,mi]=max(row);
    fprintf(fid,'  Cluster %d  %-7s  %.2f    [%d/%d/%d]\n',...
        k,CLASS_NAMES{mi},km_purity(k),row(1),row(2),row(3));
end
fprintf(fid,'\nLOO Validation (CORAL polygon holdout):\n');
fprintf(fid,'  %-10s  %8s\n','Poly ID','Recall');
for fold=1:n_loo
    if isnan(loo_recall(fold))
        fprintf(fid,'  %-10d  %8s\n',loo_ids(fold),'skipped');
    else
        fprintf(fid,'  %-10d  %8.3f\n',loo_ids(fold),loo_recall(fold));
    end
end
fprintf(fid,'  Mean CORAL-recall : %.3f  (%d polygons)\n',mean_recall,sum(valid));
fclose(fid);

%% "" Phase handoff """"""""""""""""""""""""""""""""""""""""""""
save(fullfile(cfg.OUT_DIR,'phase_stage1_output.mat'),...
    'outDir','OA','kappa','oob_err','CLASS_NAMES','CLASS_IDS',...
    'FEAT_IDX','FEAT_NAMES','N_FEAT','loo_recall','loo_ids','mean_recall',...
    'secondary_eval','transfer_results','spatial_results','RF_COST');

logW('\n========== STAGE 1 SUMMARY ==========\n');
logW('Features          : %d  [%s]\n',N_FEAT,strjoin(FEAT_NAMES,','));
logW('Training samples  : %d\n',numel(Ytrain));
logW('Testing samples   : %d\n',numel(Ytest));
logW('OOB error         : %.2f%%\n',oob_err*100);
logW('Overall Accuracy  : %.2f%%\n',OA*100);
logW('Kappa             : %.4f\n',kappa);
logW('LOO CORAL-recall  : %.3f  (%d polygons)\n',mean_recall,sum(valid));
if ~isempty(transfer_results)
    transfer_coral_rec=[transfer_results.coral_rec];
    transfer_coral_rec=transfer_coral_rec(isfinite(transfer_coral_rec));
    if ~isempty(transfer_coral_rec)
        logW('Xfer CORAL-recall : %.3f  (%d dates)\n',mean(transfer_coral_rec),numel(transfer_results));
    end
end
if ~isempty(spatial_results)
    spatial_coral_rec=[spatial_results.coral_rec];
    spatial_coral_rec=spatial_coral_rec(isfinite(spatial_coral_rec));
    if ~isempty(spatial_coral_rec)
        logW('Spatial CORAL-rec : %.3f  (%d folds)\n',mean(spatial_coral_rec),numel(spatial_results));
    end
end
logW('Output: %s\n',outDir);
logW('NEXT: Classify_Stage1.m (Stage 1) or Train_Stage2.m\n');
logW('======================================\n');

%% ===================== LOCAL FUNCTIONS =====================

function [X,Y]=sampleFeatures(A,mask,min_warn,tag,logW)
    [~,~,NF]=size(A);
    valid=mask & all(isfinite(A),3);
    idx=find(valid); n=numel(idx);
    if n==0
        X=zeros(0,NF,'single'); Y=zeros(0,1,'uint8');
        logW('  [%s]: 0 px\n',tag); return;
    end
    if min_warn>0 && n<min_warn
        logW('  [%s]: %d px  <<< BELOW MINIMUM\n',tag,n);
    else
        logW('  [%s]: %d px\n',tag,n);
    end
    X=zeros(n,NF,'single');
    for f=1:NF, sl=A(:,:,f); X(:,f)=sl(idx); end
    Y=ones(n,1,'uint8');
end

function [X,Y]=sampleFeaturesLabelled(A,mask,label_img,min_warn,tag,logW)
    [~,~,NF]=size(A);
    valid=mask & all(isfinite(A),3) & label_img>0;
    idx=find(valid); n=numel(idx);
    if n==0
        X=zeros(0,NF,'single'); Y=zeros(0,1,'uint8');
        logW('  [%s]: 0 px\n',tag); return;
    end
    if min_warn>0 && n<min_warn
        logW('  [%s]: %d px  <<< BELOW MINIMUM\n',tag,n);
    else
        logW('  [%s]: %d px\n',tag,n);
    end
    X=zeros(n,NF,'single');
    for f=1:NF, sl=A(:,:,f); X(:,f)=sl(idx); end
    Y=label_img(idx);
end

function eval=computeConfusionMetrics(ytrue,ypred,class_names,class_ids)
    ids=unique([ytrue(:); ypred(:)]);
    ids=ids(ids>0);
    n=numel(ids);
    C=zeros(n,n);
    for r_=1:n
        for c_=1:n
            C(r_,c_)=nnz(ytrue==ids(r_) & ypred==ids(c_));
        end
    end
    total_all=sum(C(:));
    if total_all==0
        OA=NaN; kappa=NaN;
    else
        OA=sum(diag(C))/total_all;
        pe=sum((sum(C,1)./total_all).*(sum(C,2)'./total_all));
        kappa=(OA-pe)/max(1-pe,1e-9);
    end
    prec=zeros(1,n); rec=zeros(1,n); f1=zeros(1,n); total=zeros(1,n); tp=zeros(1,n);
    for i=1:n
        tp(i)=C(i,i);
        total(i)=sum(C(i,:));
        prec(i)=tp(i)/max(sum(C(:,i)),1);
        rec(i)=tp(i)/max(sum(C(i,:)),1);
        f1(i)=2*prec(i)*rec(i)/max(prec(i)+rec(i),1e-9);
    end
    eval=struct('ids',ids,'names',{class_names(ids)},'n',n,'C',C,...
        'OA',OA,'kappa',kappa,'prec',prec,'rec',rec,'f1',f1,...
        'total',total,'tp',tp);
end

function [mask,label_img]=rasterizeGeoJSONLabels(path,H,W,R,class_map,tag)
    mask=false(H,W); label_img=zeros(H,W,'uint8');
    if ~exist(path,'file'), warning('Not found: %s',path); return; end
    try
        data=jsondecode(fileread(path)); feats=data.features;
        xmin=R.XWorldLimits(1); ymax=R.YWorldLimits(2);
        cw=(R.XWorldLimits(2)-xmin)/W; ch=(ymax-R.YWorldLimits(1))/H;
        for fi=1:numel(feats)
            prop=feats(fi).properties;
            if isfield(prop,'Class'), cstr=prop.Class;
            elseif isfield(prop,'class'), cstr=prop.class;
            else, continue; end
            cid=[];
            for mi=1:size(class_map,1)
                if strcmp(class_map{mi,1},cstr), cid=class_map{mi,2}; break; end
            end
            if isempty(cid), continue; end
            rings=extractRings(feats(fi).geometry);
            for ri=1:numel(rings)
                ring=rings{ri}; if size(ring,1)<3, continue; end
                [xv,yv]=latlon2utm47(ring(:,2),ring(:,1));
                cv=max(1,min(W,round((xv-xmin)/cw+0.5)));
                rv=max(1,min(H,round((ymax-yv)/ch+0.5)));
                try
                    pm=poly2mask(double(cv),double(rv),H,W);
                    mask=mask|pm; label_img(pm & label_img==0)=uint8(cid);
                catch, end
            end
        end
        fprintf('  %s: %d features -> %d px\n',tag,numel(feats),nnz(mask));
    catch ME
        warning('rasterizeGeoJSONLabels (%s): %s',tag,ME.message);
    end
end

function mask=rasterizeGeoJSON(path,H,W,R)
    mask=false(H,W);
    if ~exist(path,'file'), warning('Not found: %s',path); return; end
    try
        data=jsondecode(fileread(path)); feats=data.features;
        xmin=R.XWorldLimits(1); ymax=R.YWorldLimits(2);
        cw=(R.XWorldLimits(2)-xmin)/W; ch=(ymax-R.YWorldLimits(1))/H;
        for fi=1:numel(feats)
            rings=extractRings(feats(fi).geometry);
            for ri=1:numel(rings)
                ring=rings{ri}; if size(ring,1)<3, continue; end
                [xv,yv]=latlon2utm47(ring(:,2),ring(:,1));
                cv=max(1,min(W,round((xv-xmin)/cw+0.5)));
                rv=max(1,min(H,round((ymax-yv)/ch+0.5)));
                try, mask=mask|poly2mask(double(cv),double(rv),H,W); catch, end
            end
        end
    catch ME, warning('rasterizeGeoJSON: %s',ME.message); end
end

function rings=extractRings(geom)
    rings={};
    try
        if strcmp(geom.type,'Polygon')
            c=geom.coordinates;
            if isnumeric(c), rings{end+1}=c;
            elseif iscell(c), r=c{1}; if iscell(r),r=cell2mat(r);end; rings{end+1}=r; end
        elseif strcmp(geom.type,'MultiPolygon')
            c=geom.coordinates;
            if ~iscell(c)
                for pi=1:size(c,1), rings{end+1}=squeeze(c(pi,1,:,:)); end; return;
            end
            for pi=1:numel(c)
                poly=c{pi};
                if isnumeric(poly), rings{end+1}=poly;
                elseif iscell(poly)
                    r=poly{1}; if iscell(r),r=cell2mat(r);end
                    if isnumeric(r)&&size(r,2)>=2, rings{end+1}=r; end
                end
            end
        end
    catch, end
end

function [x,y]=latlon2utm47(lat,lon)
    a=6378137.0; f=1/298.257223563; k0=0.9996; E0=500000; lon0=99;
    e2=2*f-f^2; e4=e2^2; e6=e2^3;
    lr=deg2rad(lat); wr=deg2rad(lon); dw=wr-deg2rad(lon0);
    N_=a./sqrt(1-e2.*sin(lr).^2); T=tan(lr).^2;
    C_=e2/(1-e2).*cos(lr).^2; A_=cos(lr).*dw;
    M=a*((1-e2/4-3*e4/64-5*e6/256).*lr ...
        -(3*e2/8+3*e4/32+45*e6/1024).*sin(2*lr) ...
        +(15*e4/256+45*e6/1024).*sin(4*lr) ...
        -(35*e6/3072).*sin(6*lr));
    x=k0.*N_.*(A_+(1-T+C_).*A_.^3/6 ...
        +(5-18*T+T.^2+72*C_-58*e2/(1-e2)).*A_.^5/120)+E0;
    y=k0.*(M+N_.*tan(lr).*(A_.^2/2+(5-T+9*C_+4*C_.^2).*A_.^4/24 ...
        +(61-58*T+T.^2+600*C_-330*e2/(1-e2)).*A_.^6/720));
end

function logWrite(fid,varargin)
    msg=sprintf(varargin{:}); fprintf('%s',msg);
    if isnumeric(fid) && fid>2
        try, fprintf(fid,'%s',msg); catch, end
    end
end

function v=iif(c,a,b); if c,v=a; else,v=b; end; end

function mask = rasterizeGeoJSONClass(path, H, W, R, className)
% Rasterize a single polygon from a GeoJSON file matching className
% in the 'Class' property field. Used to extract the 'Reef' polygon
% from GeoFeature.geojson. Falls back to rasterizing all polygons
% if no Class property is found (backward compatible).
    mask = false(H,W);
    if ~exist(path,'file')
        warning('rasterizeGeoJSONClass: file not found: %s', path);
        return;
    end
    try
        data  = jsondecode(fileread(path));
        feats = data.features;
        xmin  = R.XWorldLimits(1);  ymax = R.YWorldLimits(2);
        cw    = R.CellExtentInWorldX; ch = R.CellExtentInWorldY;
        found = false;
        for fi = 1:numel(feats)
            prop = feats(fi).properties;
            cls  = '';
            if isfield(prop,'Class'), cls = prop.Class;
            elseif isfield(prop,'class'), cls = prop.class; end
            if ~strcmpi(cls, className), continue; end
            found = true;
            rings = extractRings(feats(fi).geometry);
            for ri = 1:numel(rings)
                ring = rings{ri};
                % Convert WGS84 lon/lat to UTM pixel coords
                [xu, yu] = latlon2utm47(ring(:,2), ring(:,1));
                col = round((xu - xmin) / cw) + 1;
                row = round((ymax - yu) / ch) + 1;
                valid = col>=1 & col<=W & row>=1 & row<=H;
                mask = mask | poly2mask(col(valid), row(valid), H, W);
            end
        end
        if ~found
            warning('rasterizeGeoJSONClass: no feature with Class="%s" in %s', className, path);
        end
    catch ME
        warning('rasterizeGeoJSONClass: %s', ME.message);
    end
end

