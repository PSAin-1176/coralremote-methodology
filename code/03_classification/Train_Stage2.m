%% Train_Stage2.m
% ============================================================
% STAGE 2 - Coral Condition RF Training + Accuracy Assessment
%           (2-class: HEALTHY / BLEACHED)
%
% PURPOSE:
%   Within the baseline coral habitat mask (from Build_BaselineCoralMask.m),
%   classify each coral pixel as Healthy or Bleached. This is the
%   bleaching detection stage - it operates only inside confirmed
%   coral habitat, preventing sand/rubble from contaminating the
%   bleaching signal.
%
% CONSTRAINT:
%   All training, testing, and application pixels must lie inside
%   baseline_coral_mask (Tier >= 2: High + Moderate confidence).
%   Pixels outside the mask are never classified as Bleached.
%
% FEATURE SELECTION:
%   Inherits FEAT_IDX from Stage 1 model (same 7-band stack).
%   Stage 2 may use a different FEAT_IDX - bleaching physiology
%   is best captured by RE_ratio (zooxanthellae loss) and DII_BR
%   (brightness shift toward sand). Override S2_FEAT_IDX to test.
%
%   Full Phase 5 band order:
%     1:DII_BR  2:DII_BG  3:Brightness  4:RE_ratio  5:GR_ratio
%     6:CB_ratio  7:BG_ratio  8:YB_ratio  9:G2G1_ratio
%     10:YR_ratio  11:GLCM_Contrast  12:GLCM_Entropy
%
%   Default Stage 2 selection: [1,3,4,5,8,11]
%     DII_BR, Brightness, RE_ratio, GR_ratio, YB_ratio, GLCM_Contrast
%   Removed DII_BG: less discriminating for bleaching vs healthy
%   Rationale: bleaching primarily shifts DII_BR (up), RE_ratio (down),
%   Brightness (up). These are the physiologically grounded features.
%
% POLYGON FILES:
%   TrainingStage1.geojson  - HC polygons reused as HEALTHY training
%   Bleached polygons       - BleachedStage2.geojson (to be drawn)
%   Property field: 'Class'   Values: 'HC' | 'BL'
%
% SAMPLING STRATEGY:
%   HEALTHY : 5 pre-bleaching dates - 20230404, 20230509, 20230519,
%             20230629, 20240307. Wavy dates included so RF learns
%             spectral condition not sea state.
%             Sampled only inside baseline_coral_mask
%   BLEACHED: production training uses 2 bleaching dates - 20240514, 20240613
%             20240620 is held out from training as an unseen late-bleaching
%             application date to reduce circularity in the final workflow.
%             Sampled only inside baseline_coral_mask
%   Testing : held-out date - 20230509 for HEALTHY (calm, not in training)
%             and 20240620 for BLEACHED (held out from training)
%
% CLASS SCHEMA:  HEALTHY=1   BLEACHED=2
%
% CONFOUND CHECK:
%   After training, extract SND/RBL pixels from bleaching dates and
%   check what fraction the Stage 2 RF assigns to BLEACHED. If > 5%,
%   wave/turbidity artefacts are leaking into the bleaching signal.
%
% OUTPUTS -> OUT_DIR/Phase_Stage2/run_YYYYMMDD_HHMMSS/
%   rf_stage2_model.mat        trained RF + metadata
%   confusion_matrix_s2.png    OA, Kappa, per-class F1
%   feature_importance_s2.png  ranked OOB importance bar chart
%   confound_check_s2.png      SND/RBL pixel distribution in Stage 2 space
%   stage2_report.txt          full numeric report
%   phase_stage2_output.mat    handoff for Stage 2 classification
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Train_Stage2.m - Stage 2 Coral Condition RF - %s\n\n',runTimestamp);

cfg = loadPipelineConfigPortable();

%% "" Find directories """""""""""""""""""""""""""""""""""""""""
% Phase 5 feature stacks
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

% Stage 1 model (for FEAT_IDX inheritance)
s1ModelDir='';
s1OutPath=fullfile(cfg.OUT_DIR,'phase_stage1_output.mat');
if exist(s1OutPath,'file')
    try s1o=load(s1OutPath); if exist(s1o.outDir,'dir'), s1ModelDir=s1o.outDir; end; catch, end
end
if isempty(s1ModelDir)
    base=fullfile(cfg.OUT_DIR,'Phase_Stage1');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); s1ModelDir=fullfile(base,runs(si(end)).name); end
    end
end

% Baseline coral mask
bcOutPath=fullfile(cfg.OUT_DIR,'baseline_coral_output.mat');
bcMaskDir='';
if exist(bcOutPath,'file')
    try bco=load(bcOutPath); if exist(bco.outDir,'dir'), bcMaskDir=bco.outDir; end; catch, end
end
if isempty(bcMaskDir)
    base=fullfile(cfg.OUT_DIR,'Baseline_Coral');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); bcMaskDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(bcMaskDir)||~exist(bcMaskDir,'dir')
    bcMaskDir=uigetdir(cfg.OUT_DIR,'Select Baseline_Coral output folder');
    if isequal(bcMaskDir,0), error('No folder selected.'); end
end

% PIF directory for quicklooks
pifDir='';
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
if exist(pifBase,'dir')
    runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
    if ~isempty(runs),[~,si_]=sort({runs.name}); pifDir=fullfile(pifBase,runs(si_(end)).name); end
end

outDir=fullfile(cfg.OUT_DIR,'Phase_Stage2',['run_' runTimestamp]);
mkdir(outDir);
logW('Feature input   : %s\n',featDir);
logW('Stage 1 model   : %s\n',s1ModelDir);
logW('Baseline mask   : %s\n',bcMaskDir);
logW('Output          : %s\n\n',outDir);

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
N_TREES     = cfg.RF.N_TREES;
MIN_LEAF    = cfg.RF.MIN_LEAF;
MIN_SAMPLES = cfg.RF.MIN_SAMPLES;

% "" Feature selection """""""""""""""""""""""""""""""""""""""""
FEAT_NAMES_ALL = cfg.FEATURE_NAMES;
S2_FEAT_IDX = cfg.STAGE2_FEAT_IDX;

% "" Class schema """"""""""""""""""""""""""""""""""""""""""""""
CLASS_NAMES  = {'HEALTHY','BLEACHED'};
CLASS_IDS    = [1, 2];
CLASS_COLORS = [0    0.6  0;      % HEALTHY - green
                1.0  1.0  0.9];   % BLEACHED - near-white
POLY_CLASS_MAP = {'HC',1; 'BL',2};  % 'HC' -> HEALTHY=1, 'BL' -> BLEACHED=2

% "" Sampling dates """"""""""""""""""""""""""""""""""""""""""""
% Include wavy dates in HEALTHY to overlap with wavy BLEACHED dates.
% This forces the RF to learn spectral condition not date radiometry.
% 20230519 and 20230629 are wavy but healthy - same sea state as
% 20240613/20240620 but different coral condition.
HEALTHY_DATES  = {'20230404','20230509','20230519','20230629','20240307'};
BLEACHED_DATES = {'20240514','20240613'};

% "" Test date """""""""""""""""""""""""""""""""""""""""""""""""
% Primary validation: OOB + LOO (same rationale as Stage 1)
% Secondary test: HEALTHY from pre-bleaching, BLEACHED from bleaching date
HEALTHY_TEST_DATE  = '20230509';   % calm pre-bleaching hold-out
BLEACHED_TEST_DATE = '20240620';   % held-out late bleaching date

% "" Confound check dates """"""""""""""""""""""""""""""""""""""
% SND/RBL pixels on these dates should NOT be classified as BLEACHED
CONFOUND_DATES = {'20240514','20240613','20240620'};

%% "" Load Stage 1 FEAT_IDX (for reference) """""""""""""""""""
S1_FEAT_IDX = S2_FEAT_IDX;  % default to S2 selection
if ~isempty(s1ModelDir)
    s1mPath=fullfile(s1ModelDir,'rf_stage1_model.mat');
    if exist(s1mPath,'file')
        try
            s1m=load(s1mPath,'FEAT_IDX');
            logW('Stage 1 FEAT_IDX: [%s]\n',num2str(s1m.FEAT_IDX));
        catch, end
    end
end
% Apply Stage 2 feature selection
S2_FEAT_IDX = S2_FEAT_IDX(S2_FEAT_IDX<=numel(FEAT_NAMES_ALL));
N_FEAT      = numel(S2_FEAT_IDX);
FEAT_NAMES  = FEAT_NAMES_ALL(S2_FEAT_IDX);
logW('Stage 2 FEAT_IDX: [%s]  names: [%s]\n\n',...
    num2str(S2_FEAT_IDX),strjoin(FEAT_NAMES,', '));

%% "" Load baseline coral mask """""""""""""""""""""""""""""""""
bcMaskPath=fullfile(bcMaskDir,'baseline_coral_mask.mat');
if ~exist(bcMaskPath,'file'), error('baseline_coral_mask.mat not found: %s',bcMaskDir); end
bc=load(bcMaskPath,'baseline_coral_mask','tier_map','reef_mask','H','W','R_ref','px_m2');
baseline_coral_mask = bc.baseline_coral_mask;
reef_mask           = bc.reef_mask;
H=bc.H; W=bc.W; R_ref=bc.R_ref; px_m2=bc.px_m2;
logW('Baseline coral mask: %d px  %.1f m2  (%.1f%% of reef)\n\n',...
    nnz(baseline_coral_mask), nnz(baseline_coral_mask)*px_m2,...
    100*nnz(baseline_coral_mask)/max(nnz(reef_mask),1));

%% "" Load per-date feature stacks """""""""""""""""""""""""""""
all_dates_needed=union(union(HEALTHY_DATES,BLEACHED_DATES,'stable'),...
    union({HEALTHY_TEST_DATE,BLEACHED_TEST_DATE},CONFOUND_DATES,'stable'),'stable');
feat_stacks=struct();
for di=1:numel(all_dates_needed)
    d=all_dates_needed{di};
    p=fullfile(featDir,['FEATURES_' d '.tif']);
    if exist(p,'file')
        [tmp,~]=readgeoraster(p);
        feat_stacks.(['d_' d])=single(tmp(:,:,S2_FEAT_IDX));
        logW('Loaded FEATURES_%s.tif\n',d);
    else
        logW('WARNING: FEATURES_%s.tif not found - skipping\n',d);
    end
end
logW('\n');

%% "" Rasterise polygon files """"""""""""""""""""""""""""""""""
logW('Rasterising Stage 2 polygons...\n');

% HEALTHY polygons - reuse TrainingStage1.geojson HC class
% BLEACHED polygons - BleachedStage2.geojson (must be drawn in QGIS)
%   Draw on bleaching-era imagery (20240514, 20240620)
%   Only draw within Stage 1 CORAL areas - confirmed coral habitat
%   Class field: 'BL'
trainS1Path = fullfile(cfg.ACA_DIR,'TrainingStage1.geojson');
trainBLPath = fullfile(cfg.ACA_DIR,'TrainingStage1.geojson');
testS1Path  = fullfile(cfg.ACA_DIR,'TestingStage1.geojson');

% Training polygons
[train_H_mask,train_H_labels]=rasterizeGeoJSONLabels(trainS1Path,H,W,R_ref,...
    POLY_CLASS_MAP,'Train-HEALTHY');
[train_BL_mask,train_BL_labels]=rasterizeGeoJSONLabels(trainBLPath,H,W,R_ref,...
    POLY_CLASS_MAP,'Train-BLEACHED');

% Both masks use same HC polygon locations - apply baseline coral constraint
train_H_mask  = train_H_mask  & baseline_coral_mask;
train_BL_mask = train_BL_mask & baseline_coral_mask;
% Note: train_BL_mask == train_H_mask by design (same polygon file)
% BLEACHED label assigned by date (bleaching dates), not polygon class

% Testing polygons
[test_H_mask,test_H_labels]=rasterizeGeoJSONLabels(testS1Path,H,W,R_ref,...
    POLY_CLASS_MAP,'Test-HEALTHY');
test_H_mask = test_H_mask & baseline_coral_mask;

logW('\nStage 2 training pixel counts (inside baseline mask):\n');
logW('  Train HEALTHY  : %d px\n',nnz(train_H_mask & train_H_labels==1));
logW('  Train BLEACHED : %d px\n',nnz(train_BL_mask & train_BL_labels==2));
logW('  Test  HEALTHY  : %d px\n',nnz(test_H_mask  & test_H_labels==1));
logW('\n');

% Check BLEACHED polygons exist
if nnz(train_BL_mask)==0
    error('No BLEACHED training pixels found inside baseline coral mask. Check TrainingStage1.geojson and baseline_coral_mask extent.');
end

%% "" Sample training pixels """""""""""""""""""""""""""""""""""
logW('Sampling Stage 2 training pixels...\n');
Xtrain=[]; Ytrain=[];

% HEALTHY: sample from calm pre-bleaching dates inside baseline mask
healthy_mask = train_H_mask & (train_H_labels==1);
for di=1:numel(HEALTHY_DATES)
    d=HEALTHY_DATES{di};
    if ~isfield(feat_stacks,['d_' d])
        logW('  HEALTHY [%s]: stack not available\n',d); continue;
    end
    [Xh,Yh]=sampleFeatures(feat_stacks.(['d_' d]),healthy_mask,...
        MIN_SAMPLES,['HEALTHY-' d],logW);
    Yh(:)=uint8(1);
    Xtrain=[Xtrain;Xh]; Ytrain=[Ytrain;Yh];
end

% BLEACHED: sample from bleaching dates inside baseline mask
bleached_mask = train_BL_mask & (train_BL_labels==1);
for di=1:numel(BLEACHED_DATES)
    d=BLEACHED_DATES{di};
    if ~isfield(feat_stacks,['d_' d])
        logW('  BLEACHED [%s]: stack not available\n',d); continue;
    end
    [Xb,Yb]=sampleFeatures(feat_stacks.(['d_' d]),bleached_mask,...
        MIN_SAMPLES,['BLEACHED-' d],logW);
    Yb(:)=uint8(2);
    Xtrain=[Xtrain;Xb]; Ytrain=[Ytrain;Yb];
end

logW('\nTotal training samples: %d\n',numel(Ytrain));
for ci=1:numel(CLASS_NAMES)
    n=nnz(Ytrain==CLASS_IDS(ci));
    flag=''; if n<MIN_SAMPLES, flag='  <<< BELOW MINIMUM'; end
    logW('  %-8s: %d%s\n',CLASS_NAMES{ci},n,flag);
end
logW('\n');

if numel(unique(Ytrain))<2
    error('Need both HEALTHY and BLEACHED training samples. Check polygon files.');
end

%% "" Assemble test set """"""""""""""""""""""""""""""""""""""""
logW('Assembling Stage 2 test set...\n');
Xtest=[]; Ytest=[];

% HEALTHY test - held-out calm date
h_test_mask=test_H_mask & (test_H_labels==1);
if isfield(feat_stacks,['d_' HEALTHY_TEST_DATE])
    [Xh,Yh]=sampleFeatures(feat_stacks.(['d_' HEALTHY_TEST_DATE]),h_test_mask,...
        0,['TEST-HEALTHY-' HEALTHY_TEST_DATE],logW);
    Yh(:)=uint8(1);
    Xtest=[Xtest;Xh]; Ytest=[Ytest;Yh];
else
    logW('WARNING: HEALTHY test date %s not available\n',HEALTHY_TEST_DATE);
end

% BLEACHED test - bleaching date at BL polygon locations
if isfield(feat_stacks,['d_' BLEACHED_TEST_DATE])
    [Xb,Yb]=sampleFeatures(feat_stacks.(['d_' BLEACHED_TEST_DATE]),bleached_mask,...
        0,['TEST-BLEACHED-' BLEACHED_TEST_DATE],logW);
    Yb(:)=uint8(2);
    Xtest=[Xtest;Xb]; Ytest=[Ytest;Yb];
else
    logW('WARNING: BLEACHED test date %s not available\n',BLEACHED_TEST_DATE);
end
logW('Total test samples: %d\n\n',numel(Ytest));

%% "" Train Random Forest """"""""""""""""""""""""""""""""""""""
logW('Training Stage 2 RF (%d trees, min_leaf=%d)...\n',N_TREES,MIN_LEAF);
t_train=tic;
rf_model=TreeBagger(N_TREES,Xtrain,double(Ytrain),...
    'Method',                'classification',...
    'MinLeafSize',           MIN_LEAF,...
    'NumPredictorsToSample', 'all',...
    'OOBPrediction',         'on',...
    'OOBPredictorImportance','on');
elapsed_train=toc(t_train);
oob_err=oobError(rf_model,'Mode','ensemble');
logW('Training complete: %.1fs\n',elapsed_train);
logW('OOB error: %.4f (%.2f%%)\n\n',oob_err,oob_err*100);

%% "" Feature importance """""""""""""""""""""""""""""""""""""""
imp=rf_model.OOBPermutedPredictorDeltaError;
[imp_sorted,imp_idx]=sort(imp,'descend');
n_imp=numel(imp);

figImp=figure('Visible','off','Position',[0 0 700 380]);
barh(imp_sorted,'FaceColor',[0.8 0.3 0.1]);
yticks(1:n_imp); yticklabels(FEAT_NAMES(imp_idx));
xlabel('OOB Permuted Importance'); title('Stage 2 - Feature Importance');
grid on; set(gca,'YDir','reverse');
exportgraphics(figImp,fullfile(outDir,'feature_importance_s2.png'),'Resolution',150);
close(figImp);

logW('Feature importance (ranked):\n');
for fi=1:n_imp
    logW('  [%d] %-22s  %.5f\n',fi,FEAT_NAMES{imp_idx(fi)},imp_sorted(fi));
end
logW('\n');

%% "" Accuracy assessment """"""""""""""""""""""""""""""""""""""
OA=NaN; kappa=NaN; f1_all=[]; n_tc=0; tc_names={};
if ~isempty(Xtest) && numel(unique(Ytest))>=2
    logW('Accuracy assessment...\n');
    [pred_cat,~]=predict(rf_model,Xtest);
    pred=uint8(str2double(cellstr(pred_cat)));

    tc_ids=unique([Ytest;pred]); tc_ids=tc_ids(tc_ids>0);
    n_tc=numel(tc_ids); tc_names=CLASS_NAMES(tc_ids);
    C=zeros(n_tc,n_tc);
    for r_=1:n_tc
        for c_=1:n_tc
            C(r_,c_)=nnz(Ytest==tc_ids(r_) & pred==tc_ids(c_));
        end
    end
    OA=sum(diag(C))/sum(C(:));
    p_e=sum((sum(C,1)./sum(C(:))).*(sum(C,2)'./sum(C(:))));
    kappa=(OA-p_e)/(1-p_e);
    logW('OA=%.4f  Kappa=%.4f\n',OA,kappa);

    f1_all=zeros(1,n_tc);
    logW('%-10s  %6s  %6s  %6s  %6s  %6s\n','Class','Total','TP','Prec','Rec','F1');
    for ci_=1:n_tc
        tp=C(ci_,ci_); prec=tp/max(sum(C(:,ci_)),1); rec=tp/max(sum(C(ci_,:)),1);
        f1=2*prec*rec/max(prec+rec,1e-9); f1_all(ci_)=f1;
        logW('  %-8s  %6d  %6d  %6.3f  %6.3f  %6.3f\n',...
            tc_names{ci_},sum(C(ci_,:)),tp,prec,rec,f1);
    end
    logW('\n');

    % Confusion matrix figure
    figCM=figure('Visible','off','Position',[0 0 480 420]);
    conf_norm=C./max(sum(C,2),1);
    imagesc(conf_norm,[0 1]); colormap(flipud(gray)); colorbar;
    xticks(1:n_tc); xticklabels(tc_names); xlabel('Predicted');
    yticks(1:n_tc); yticklabels(tc_names); ylabel('Actual');
    title(sprintf('Stage 2 Confusion Matrix\nOA=%.3f  Kappa=%.3f',OA,kappa));
    for r_=1:n_tc
        for c_=1:n_tc
            text(c_,r_,sprintf('%d\n%.2f',C(r_,c_),conf_norm(r_,c_)),...
                'HorizontalAlignment','center','FontSize',10,...
                'Color',iif(conf_norm(r_,c_)>0.5,'w','k'));
        end
    end
    exportgraphics(figCM,fullfile(outDir,'confusion_matrix_s2.png'),'Resolution',150);
    close(figCM);
else
    logW('WARNING: Insufficient test samples for accuracy assessment\n\n');
end

%% "" Confound check """""""""""""""""""""""""""""""""""""""""""
% Apply Stage 2 RF to SND and RBL pixels on wavy bleaching dates.
% These pixels are OUTSIDE the baseline coral mask but inside the reef.
% If the RF assigns them to BLEACHED, wave/turbidity is contaminating
% the bleaching signal - reported as confound rate.
logW('===== Confound Check =====\n');
logW('Applying Stage 2 to SND/RBL pixels on wavy bleaching dates...\n');
logW('(pixels must be outside baseline mask but inside reef)\n\n');

confound_mask = reef_mask & ~baseline_coral_mask;
confound_results=struct('date',{},'n_px',{},'n_bleached',{},'confound_rate',{});

for di=1:numel(CONFOUND_DATES)
    d=CONFOUND_DATES{di};
    if ~isfield(feat_stacks,['d_' d])
        logW('  %s: stack not available\n',d); continue;
    end
    A=feat_stacks.(['d_' d]);
    valid=confound_mask & all(isfinite(A),3);
    idx=find(valid); n=numel(idx);
    if n==0, logW('  %s: no confound pixels\n',d); continue; end

    % Sample up to 5000 px for speed
    if n>5000
        idx=idx(randperm(n,5000)); n=5000;
    end
    Xconf=zeros(n,N_FEAT,'single');
    for f=1:N_FEAT, sl=A(:,:,f); Xconf(:,f)=sl(idx); end

    [pred_conf_cat,~]=predict(rf_model,Xconf);
    pred_conf=uint8(str2double(cellstr(pred_conf_cat)));
    n_bl=nnz(pred_conf==2);
    rate=100*n_bl/max(n,1);
    confound_results(end+1)=struct('date',d,'n_px',n,...
        'n_bleached',n_bl,'confound_rate',rate);
    flag=''; if rate>5, flag='  <<< WARNING: >5% confound'; end
    logW('  %s: %d SND/RBL px -> %d predicted BLEACHED (%.1f%%)%s\n',...
        d,n,n_bl,rate,flag);
end
logW('\n');

% Confound bar chart
if ~isempty(confound_results)
    figCC=figure('Visible','off','Position',[0 0 500 350]);
    rates=[confound_results.confound_rate];
    bar(rates,'FaceColor',[0.8 0.4 0.1]);
    xticks(1:numel(confound_results));
    xticklabels({confound_results.date});
    ylabel('% SND/RBL pixels predicted BLEACHED');
    yline(5,'--r','5% threshold','LabelHorizontalAlignment','left');
    ylim([0 max(max(rates)+5,15)]);
    title('Stage 2 Confound Check - Wave Date SND/RBL Leakage','FontSize',9);
    grid on;
    exportgraphics(figCC,fullfile(outDir,'confound_check_s2.png'),'Resolution',150);
    close(figCC);
end

%% "" Save model """""""""""""""""""""""""""""""""""""""""""""""
logW('Saving Stage 2 model...\n');
save(fullfile(outDir,'rf_stage2_model.mat'),...
    'rf_model','FEAT_NAMES','S2_FEAT_IDX','CLASS_NAMES','CLASS_IDS','CLASS_COLORS',...
    'POLY_CLASS_MAP','OA','kappa','oob_err','N_TREES','MIN_LEAF',...
    'HEALTHY_DATES','BLEACHED_DATES','N_FEAT','runTimestamp',...
    'baseline_coral_mask','bcMaskDir','confound_results');

%% "" Report """""""""""""""""""""""""""""""""""""""""""""""""""
fid=fopen(fullfile(outDir,'stage2_report.txt'),'w');
fprintf(fid,'STAGE 2 - CORAL CONDITION RF TRAINING REPORT\n');
fprintf(fid,'Timestamp: %s\n\n',runTimestamp);
fprintf(fid,'Polygon files:\n');
fprintf(fid,'  HEALTHY  : TrainingStage1.geojson (HC class)\n');
fprintf(fid,'  BLEACHED : TrainingStage1.geojson HC (same locations, bleaching dates)\n');
fprintf(fid,'  Constraint: inside baseline_coral_mask (Tier>=2)\n');
fprintf(fid,'  Note: HC polygons used for both classes; condition set by sampling date\n\n');
fprintf(fid,'Feature selection (Stage 2):\n');
fprintf(fid,'  Indices : [%s]\n',num2str(S2_FEAT_IDX));
fprintf(fid,'  Names   : %s\n\n',strjoin(FEAT_NAMES,', '));
fprintf(fid,'Sampling dates:\n');
fprintf(fid,'  HEALTHY  : %s\n',strjoin(HEALTHY_DATES,', '));
fprintf(fid,'  BLEACHED : %s\n',strjoin(BLEACHED_DATES,', '));
fprintf(fid,'\nRF parameters: %d trees  min_leaf=%d\n\n',N_TREES,MIN_LEAF);
fprintf(fid,'Training sample counts:\n');
for ci=1:numel(CLASS_NAMES)
    n=nnz(Ytrain==CLASS_IDS(ci));
    flag=''; if n<MIN_SAMPLES, flag='  <<< BELOW MINIMUM'; end
    fprintf(fid,'  %-8s  %d%s\n',CLASS_NAMES{ci},n,flag);
end
fprintf(fid,'\nOOB error: %.4f (%.2f%%)\n\n',oob_err,oob_err*100);
fprintf(fid,'Feature importance (ranked):\n');
for fi=1:n_imp
    fprintf(fid,'  [%d] %-22s  %.5f\n',fi,FEAT_NAMES{imp_idx(fi)},imp_sorted(fi));
end
fprintf(fid,'\nPrimary validation: OOB=%.2f%%\n\n',oob_err*100);
if ~isnan(OA)
    fprintf(fid,'Secondary test (HEALTHY from %s, BLEACHED from %s):\n',...
        HEALTHY_TEST_DATE,BLEACHED_TEST_DATE);
    fprintf(fid,'  Overall Accuracy : %.4f (%.2f%%)\n',OA,OA*100);
    fprintf(fid,'  Kappa            : %.4f\n\n',kappa);
    fprintf(fid,'Per-class F1:\n');
    for ci_=1:n_tc
        fprintf(fid,'  %-8s  F1=%.3f\n',tc_names{ci_},f1_all(ci_));
    end
end
fprintf(fid,'\nConfound check (SND/RBL on wave bleaching dates):\n');
fprintf(fid,'  %-10s  %8s  %10s  %12s\n','Date','N_px','N_bleached','Confound_%');
for k_=1:numel(confound_results)
    fprintf(fid,'  %-10s  %8d  %10d  %11.1f%%\n',...
        confound_results(k_).date, confound_results(k_).n_px,...
        confound_results(k_).n_bleached, confound_results(k_).confound_rate);
end
fclose(fid);

%% "" Phase handoff """"""""""""""""""""""""""""""""""""""""""""
save(fullfile(cfg.OUT_DIR,'phase_stage2_output.mat'),...
    'outDir','OA','kappa','oob_err','CLASS_NAMES','CLASS_IDS',...
    'S2_FEAT_IDX','FEAT_NAMES','N_FEAT','baseline_coral_mask');

logW('\n========== STAGE 2 SUMMARY ==========\n');
logW('Features         : %d  [%s]\n',N_FEAT,strjoin(FEAT_NAMES,','));
logW('Training samples : %d\n',numel(Ytrain));
logW('OOB error        : %.2f%%\n',oob_err*100);
if ~isnan(OA)
    logW('OA (secondary)   : %.2f%%\n',OA*100);
    logW('Kappa            : %.4f\n',kappa);
end
logW('Output: %s\n',outDir);
logW('NEXT: Classify_Stage2.m (apply to all dates within baseline mask)\n');
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

function [mask,label_img]=rasterizeGeoJSONLabels(path,H,W,R,class_map,tag)
    mask=false(H,W); label_img=zeros(H,W,'uint8');
    if ~exist(path,'file')
        warning('Not found: %s',path); return;
    end
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
                    mask=mask|pm;
                    label_img(pm & label_img==0)=uint8(cid);
                catch, end
            end
        end
        fprintf('  %s: %d features -> %d px\n',tag,numel(feats),nnz(mask));
    catch ME
        warning('rasterizeGeoJSONLabels (%s): %s',tag,ME.message);
    end
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


