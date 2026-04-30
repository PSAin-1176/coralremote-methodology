%% Classify_Stage2.m
% ============================================================
% STAGE 2 CLASSIFICATION - Coral Condition (HEALTHY / BLEACHED)
%
% Applies the Stage 2 RF model to all dates within the baseline
% coral habitat mask. Outside the mask: pixels are not classified
% (Stage 1 substrate class is retained). Inside the mask: each
% coral pixel is assigned HEALTHY, BLEACHED, or UNK.
%
% SPATIAL CONSTRAINT:
%   Only pixels inside baseline_coral_mask (Tier >= 2) are
%   classified by Stage 2. All other reef pixels retain their
%   Stage 1 substrate label (SND, RBL, UNK).
%
% UNK THRESHOLD:
%   Coral pixels where max Stage 2 RF probability < S2_UNK_THRESHOLD
%   are labelled UNK_S2 rather than forced into HEALTHY/BLEACHED.
%
% WAVE FLAG:
%   Dates with confirmed high confound rates are flagged in the
%   report as high-uncertainty. Bleaching estimates on these dates
%   should be interpreted as upper bounds.
%   High confound: 20240613 (95.9%), 20240620 (99.1%)
%   Lower confound: 20240514 (58.5%)
%
% OUTPUT CLASS SCHEMA (combined Stage 1 + Stage 2):
%   1 = HEALTHY  (coral, Stage 2 confirmed healthy)
%   2 = BLEACHED (coral, Stage 2 confirmed bleached)
%   3 = SND      (sand, Stage 1)
%   4 = RBL      (rubble, Stage 1)
%   5 = UNK      (uncertain - low confidence in either stage)
%
% DATES TO CLASSIFY (all 10):
%   Pre-bleaching: 20230404, 20230509, 20230519, 20230629
%   Bleaching    : 20240514, 20240613*, 20240620*
%   Recovery     : 20240926, 20241111
%   * high wave confound - flagged as upper bound estimates
%
% OUTPUTS -> OUT_DIR/Phase_Stage2_Classify/run_YYYYMMDD_HHMMSS/
%   S2CLASS_YYYYMMDD.tif         uint8 combined class map per date
%   S2CONF_YYYYMMDD.tif          float32 Stage 2 confidence per date
%   VIZ_S2_YYYYMMDD.png          RGB + S2 class overlay per date
%   s2_mosaic_all.png            all dates in one figure
%   s2_bleaching_timeseries.png  BLEACHED % of coral mask over time
%   s2_area_chart.png            stacked area chart per date
%   s2_classification_report.txt full numeric report with wave flags
%   phase_stage2_classify.mat    handoff for change analysis
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Classify_Stage2.m - Stage 2 Coral Condition - %s\n\n',runTimestamp);

cfg = loadPipelineConfigPortable();
FAST_MODE = isfield(cfg,'FAST_MODE') && cfg.FAST_MODE;

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
S2_UNK_THRESHOLD = 0.55;   % Stage 2 confidence below this -> UNK

% All 11 scientific observation dates
DATES_TO_CLASSIFY = {...
    '20230404','20230509','20230519','20230629','20240307'...  % pre-bleaching
    '20240514','20240613','20240620',...             % bleaching
    '20240926','20241111','20250518'};                          % recovery

% Wave-contaminated dates - flagged as upper bound in report
WAVE_FLAG_DATES = {'20240613','20240620'};
CONFOUND_RATES  = containers.Map(...
    {'20240514','20240613','20240620'},...
    [58.5, 95.9, 99.1]);

% Combined output class schema
OUT_CLASS.HEALTHY  = 1;
OUT_CLASS.BLEACHED = 2;
OUT_CLASS.SND      = 3;
OUT_CLASS.RBL      = 4;
OUT_CLASS.UNK      = 5;

OUT_CLASS_NAMES  = {'HEALTHY','BLEACHED','SND','RBL','UNK'};
OUT_CLASS_IDS    = [1, 2, 3, 4, 5];
OUT_CLASS_COLORS = [0.0  0.6  0.0;    % HEALTHY  - green
                    1.0  1.0  0.85;   % BLEACHED - near-white/yellow
                    1.0  0.94 0.70;   % SND      - pale yellow
                    0.7  0.55 0.31;   % RBL      - tan
                    0.65 0.65 0.65];  % UNK      - grey

%% "" Find directories """""""""""""""""""""""""""""""""""""""""
% Stage 2 model
s2ModelDir='';
s2OutPath=fullfile(cfg.OUT_DIR,'phase_stage2_output.mat');
if exist(s2OutPath,'file')
    try s2o=load(s2OutPath); if exist(s2o.outDir,'dir'), s2ModelDir=s2o.outDir; end; catch, end
end
if isempty(s2ModelDir)
    base=fullfile(cfg.OUT_DIR,'Phase_Stage2');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); s2ModelDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(s2ModelDir)||~exist(s2ModelDir,'dir')
    s2ModelDir=uigetdir(cfg.OUT_DIR,'Select Stage 2 model folder');
    if isequal(s2ModelDir,0), error('No folder selected.'); end
end

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
    bcMaskDir=uigetdir(cfg.OUT_DIR,'Select Baseline_Coral folder');
    if isequal(bcMaskDir,0), error('No folder selected.'); end
end

% Stage 1 classification output (for SND/RBL labels outside coral mask)
s1ClassDir='';
p9Path=fullfile(cfg.OUT_DIR,'phase9_output.mat');
if exist(p9Path,'file')
    try p9=load(p9Path); if exist(p9.outDir,'dir'), s1ClassDir=p9.outDir; end; catch, end
end
if isempty(s1ClassDir)
    base=fullfile(cfg.OUT_DIR,'Phase9_Classify');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); s1ClassDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(s1ClassDir)||~exist(s1ClassDir,'dir')
    s1ClassDir=uigetdir(cfg.OUT_DIR,'Select Stage 1 classification folder');
    if isequal(s1ClassDir,0), error('No folder selected.'); end
end

% Display RGB source (prefer current glint-corrected imagery)
rgbDir='';
g3Path=fullfile(cfg.OUT_DIR,'phase3_output.mat');
if exist(g3Path,'file')
    try
        g3=load(g3Path);
        if exist(g3.outDir,'dir'), rgbDir=g3.outDir; end
    catch
    end
end
if isempty(rgbDir)
    rgbBase=fullfile(cfg.OUT_DIR,'Phase3_LogRatio');
    if exist(rgbBase,'dir')
        runs=dir(fullfile(rgbBase,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si_]=sort({runs.name}); rgbDir=fullfile(rgbBase,runs(si_(end)).name); end
    end
end

outDir=fullfile(cfg.OUT_DIR,'Phase_Stage2_Classify',['run_' runTimestamp]);
mkdir(outDir);
logW('Stage 2 model    : %s\n',s2ModelDir);
logW('Feature input    : %s\n',featDir);
logW('Baseline mask    : %s\n',bcMaskDir);
logW('Stage 1 classes  : %s\n',s1ClassDir);
logW('Display RGB      : %s\n',rgbDir);
logW('Output           : %s\n\n',outDir);
logW('UNK threshold    : %.2f\n',S2_UNK_THRESHOLD);
logW('Wave flag dates  : %s\n',strjoin(WAVE_FLAG_DATES,', '));
logW('Fast mode        : %s\n\n', iif(FAST_MODE,'on','off'));

%% "" Load Stage 2 model """""""""""""""""""""""""""""""""""""""
s2mPath=fullfile(s2ModelDir,'rf_stage2_model.mat');
if ~exist(s2mPath,'file'), error('rf_stage2_model.mat not found: %s',s2ModelDir); end
s2m=load(s2mPath);
rf_s2      = s2m.rf_model;
S2_FEAT_IDX= s2m.S2_FEAT_IDX;
N_FEAT_S2  = s2m.N_FEAT;
FEAT_NAMES_S2 = s2m.FEAT_NAMES;
logW('Stage 2 RF loaded: OOB=%.2f%%  features=[%s]\n\n',...
    s2m.oob_err*100, strjoin(FEAT_NAMES_S2,','));

%% "" Load baseline coral mask """""""""""""""""""""""""""""""""
bcMaskPath=fullfile(bcMaskDir,'baseline_coral_mask.mat');
if ~exist(bcMaskPath,'file'), error('baseline_coral_mask.mat not found.'); end
bc=load(bcMaskPath,'baseline_coral_mask','reef_mask','H','W','R_ref','px_m2');
baseline_coral_mask = bc.baseline_coral_mask;
reef_mask           = bc.reef_mask;
H=bc.H; W=bc.W; R_ref=bc.R_ref; px_m2=bc.px_m2;
n_coral_mask = nnz(baseline_coral_mask);
logW('Baseline coral mask: %d px  %.1f m2  (%.1f%% reef)\n\n',...
    n_coral_mask, n_coral_mask*px_m2,...
    100*n_coral_mask/max(nnz(reef_mask),1));

%% "" Classify each date """""""""""""""""""""""""""""""""""""""
nDates   = numel(DATES_TO_CLASSIFY);
area_stats = struct('date',{},'class',{},'px',{},'m2',{},'pct_reef',{},...
                    'pct_coral',{},'wave_flag',{});
s2_class_maps = cell(nDates,1);
s2_conf_maps  = cell(nDates,1);
rgb_cache     = cell(nDates,1);

% Per-date summary: % BLEACHED of coral mask
bleach_pct = zeros(1,nDates);
healthy_pct= zeros(1,nDates);
reef_pct_s2= zeros(nDates, numel(OUT_CLASS_IDS));

for di=1:nDates
    d=DATES_TO_CLASSIFY{di};
    is_wave_flag = ismember(d,WAVE_FLAG_DATES);
    logW('============================================================\n');
    logW('[%d/%d] %s%s\n',di,nDates,d,...
        iif(is_wave_flag,'  *** WAVE FLAG - upper bound ***',''));

    % "" Load Stage 2 feature stack """"""""""""""""""""""""""""
    featPath=fullfile(featDir,['FEATURES_' d '.tif']);
    if ~exist(featPath,'file')
        logW('  FEATURES_%s.tif not found - skipping\n\n',d); continue;
    end
    [Afull,Ri]=readgeoraster(featPath);
    Ai=single(Afull(:,:,S2_FEAT_IDX));

    px_ha_=abs(Ri.CellExtentInWorldX*Ri.CellExtentInWorldY)/10000;

    % "" Load Stage 1 classification for this date """""""""""""
    % Used to assign SND/RBL labels outside coral mask
    s1_class=zeros(H,W,'uint8');
    s1cPath=fullfile(s1ClassDir,['CLASS_' d '.tif']);
    if exist(s1cPath,'file')
        [tmp,~]=readgeoraster(s1cPath);
        s1_class=uint8(tmp);
    else
        logW('  WARNING: CLASS_%s.tif not found - SND/RBL outside mask will be UNK\n',d);
    end

    % "" Initialise output maps """""""""""""""""""""""""""""""""
    % Start with Stage 1 substrate classes remapped to output schema
    % Stage 1: CORAL=1->handled by S2, SND=2->SND=3, RBL=3->RBL=4, UNK=4->UNK=5
    s1_remap = [0, OUT_CLASS.SND, OUT_CLASS.RBL, OUT_CLASS.UNK];
    out_map  = zeros(H,W,'uint8');
    conf_map = zeros(H,W,'single');

    % Fill non-coral pixels from Stage 1
    for ci=2:4   % SND, RBL, UNK from Stage 1
        out_map(s1_class==ci & ~baseline_coral_mask) = s1_remap(ci);
    end
    % Reef pixels not in Stage 1 output and not in coral mask -> UNK
    out_map(reef_mask & ~baseline_coral_mask & out_map==0) = OUT_CLASS.UNK;

    % "" Apply Stage 2 inside coral mask """""""""""""""""""""""
    classify_mask = baseline_coral_mask & all(isfinite(Ai),3);
    nonclass_mask = baseline_coral_mask & ~classify_mask;
    if any(nonclass_mask(:))
        % Keep every baseline-coral pixel in a reported class so downstream
        % area tables reconcile exactly to the coral mask area.
        out_map(nonclass_mask) = OUT_CLASS.UNK;
    end
    idx_px = find(classify_mask);
    n_px   = numel(idx_px);
    logW('  Coral mask pixels classifiable: %d / %d\n',n_px,n_coral_mask);
    if any(nonclass_mask(:))
        logW('  Non-classifiable coral-mask pixels -> UNK: %d px\n',nnz(nonclass_mask));
    end

    if n_px>0
        BATCH=50000;
        Ypred_all = zeros(n_px,1,'uint8');
        prob_all  = zeros(n_px,1,'single');

        for b_start=1:BATCH:n_px
            b_end=min(b_start+BATCH-1,n_px);
            b_idx=idx_px(b_start:b_end);
            Xb=zeros(numel(b_idx),N_FEAT_S2,'single');
            for f=1:N_FEAT_S2
                sl=Ai(:,:,f); Xb(:,f)=sl(b_idx);
            end
            [Ypred_cat,prob_b]=predict(rf_s2,Xb);
            Ypred_all(b_start:b_end)=uint8(str2double(cellstr(Ypred_cat)));
            prob_all(b_start:b_end) =single(max(prob_b,[],2));
        end

        % Map Stage 2 classes to output schema
        % S2: HEALTHY=1 -> OUT=1, BLEACHED=2 -> OUT=2
        out_coral = Ypred_all;   % already 1 or 2

        % Apply UNK threshold
        unk_s2 = prob_all < S2_UNK_THRESHOLD;
        out_coral(unk_s2) = uint8(OUT_CLASS.UNK);
        n_unk = nnz(unk_s2);
        logW('  Stage 2 UNK (conf<%.2f): %d px (%.1f%% of coral mask)\n',...
            S2_UNK_THRESHOLD,n_unk,100*n_unk/max(n_px,1));

        % Write Stage 2 results into output map
        out_map(idx_px) = out_coral;
        conf_map(idx_px)= prob_all;
    end

    s2_class_maps{di} = out_map;
    s2_conf_maps{di}  = conf_map;

    % "" Save GeoTIFFs """"""""""""""""""""""""""""""""""""""""""
    saveUint8Tif(fullfile(outDir,['S2CLASS_' d '.tif']), out_map, Ri);
    saveMultiBandTif(fullfile(outDir,['S2CONF_' d '.tif']), conf_map, Ri);
    logW('  Saved S2CLASS_ and S2CONF_ tifs\n');

    % "" Area statistics """"""""""""""""""""""""""""""""""""""""
    reef_px = nnz(reef_mask);
    logW('  Area stats (inside reef):\n');
    for ci=1:numel(OUT_CLASS_IDS)
        cid   = OUT_CLASS_IDS(ci);
        cname = OUT_CLASS_NAMES{ci};
        npx   = nnz(out_map==cid & reef_mask);
        m2_   = npx*px_m2;
        pct_r = 100*npx/max(reef_px,1);
        pct_c = 100*npx/max(n_coral_mask,1);
        reef_pct_s2(di,ci) = pct_r;
        logW('    %-8s  %6d px  %10.1f m2  %5.1f%% reef\n',cname,npx,m2_,pct_r);
        area_stats(end+1)=struct('date',d,'class',cname,...
            'px',npx,'m2',m2_,'pct_reef',pct_r,'pct_coral',pct_c,...
            'wave_flag',is_wave_flag);
    end

    % Store bleaching and healthy % of coral mask
    n_bl  = nnz(out_map==OUT_CLASS.BLEACHED & baseline_coral_mask);
    n_hy  = nnz(out_map==OUT_CLASS.HEALTHY  & baseline_coral_mask);
    bleach_pct(di)  = 100*n_bl/max(n_coral_mask,1);
    healthy_pct(di) = 100*n_hy/max(n_coral_mask,1);
    logW('  Bleaching severity (of coral mask): %.1f%%  Healthy: %.1f%%\n\n',...
        bleach_pct(di),healthy_pct(di));

    % "" Load PIF RGB """""""""""""""""""""""""""""""""""""""""""
    rgb_img=zeros(H,W,3,'single');
    if ~isempty(rgbDir)
        pf=dir(fullfile(rgbDir,['GLINT_' d '.tif']));
        if ~isempty(pf)
            try
                [Ap,~]=readgeoraster(fullfile(rgbDir,pf(1).name));
                Ap=single(Ap);
                r_=Ap(:,:,cfg.BANDS.RED);
                g_=Ap(:,:,cfg.BANDS.GREEN1);
                b_=Ap(:,:,cfg.BANDS.BLUE);
                rgb_img=makeDisplayRGB(cat(3,r_,g_,b_), reef_mask);
            catch, end
        end
    end
    rgb_cache{di}=rgb_img;

    % "" Per-date VIZ (3 panels) """""""""""""""""""""""""""""""
    if ~FAST_MODE
    figV=figure('Visible','off','Position',[0 0 1200 430]);

    % "" Compute coral mask bounding box for zoom """""""""""""
    [rows_c,cols_c]=find(baseline_coral_mask);
    margin=30;   % px padding around coral mask
    zoom_r1=max(1,min(rows_c)-margin); zoom_r2=min(H,max(rows_c)+margin);
    zoom_c1=max(1,min(cols_c)-margin); zoom_c2=min(W,max(cols_c)+margin);

    subplot(1,3,1);
    imshow(rgb_img);
    xlim([zoom_c1 zoom_c2]); ylim([zoom_r1 zoom_r2]);
    title_str=sprintf('%s - RGB',d);
    if is_wave_flag, title_str=[title_str '  WAVE']; end
    title(title_str,'FontSize',9,'Interpreter','none');

    subplot(1,3,2);
    cmap_out=[0 0 0; OUT_CLASS_COLORS];
    display_map=zeros(H,W,'uint8');
    display_map(baseline_coral_mask)=out_map(baseline_coral_mask);
    cm_shift=display_map+1; cm_shift(~reef_mask)=1;
    class_rgb=ind2rgb(cm_shift,cmap_out);
    blend=rgb_img;
    coral_disp=baseline_coral_mask & isfinite(conf_map);
    for ch=1:3
        tmp=blend(:,:,ch);
        cls=class_rgb(:,:,ch);
        tmp(coral_disp)=0.70.*cls(coral_disp)+0.30.*tmp(coral_disp);
        blend(:,:,ch)=tmp;
    end
    blend=max(0,min(1,blend));
    imshow(blend); hold on;
    coral_edge=imdilate(baseline_coral_mask,strel('disk',1))&~baseline_coral_mask;
    [re,ce]=find(coral_edge);
    if ~isempty(re), scatter(ce,re,1,[1 1 1],'.'); end
    xlim([zoom_c1 zoom_c2]); ylim([zoom_r1 zoom_r2]);
    title(sprintf('Stage 2 Classification\n(coral mask overlay)  BLEACHED %.1f%%  HEALTHY %.1f%%',...
        bleach_pct(di),healthy_pct(di)),'FontSize',8);
    for ci=1:numel(OUT_CLASS_IDS)
        patch(nan,nan,OUT_CLASS_COLORS(ci,:),'EdgeColor','none',...
            'DisplayName',OUT_CLASS_NAMES{ci});
    end
    legend('Location','southoutside','Orientation','horizontal','FontSize',6,'Box','off');

    subplot(1,3,3);
    conf_disp=conf_map; conf_disp(~baseline_coral_mask)=NaN;
    imagesc(conf_disp,[0 1]); axis image off;
    xlim([zoom_c1 zoom_c2]); ylim([zoom_r1 zoom_r2]);
    colormap(gca,hot); cb=colorbar; cb.Label.String='Stage 2 confidence';
    title(sprintf('Stage 2 Confidence\n(coral mask only)'),'FontSize',9);

    wave_str='';
    if is_wave_flag
        cr=CONFOUND_RATES(d);
        wave_str=sprintf(' - WAVE FLAG: confound=%.1f%%',cr);
    end
    sgtitle(sprintf('Stage 2 - %s%s',d,wave_str),'FontSize',10,'FontWeight','bold');
    exportgraphics(figV,fullfile(outDir,['VIZ_S2_' d '.png']),'Resolution',150);
    close(figV);
    logW('  Saved VIZ_S2_%s.png\n\n',d);
    end
end

%% "" All-dates mosaic """""""""""""""""""""""""""""""""""""""""
if ~FAST_MODE
logW('Generating all-dates mosaic...\n');
ncols=min(nDates,5); nrows=ceil(nDates/ncols);
figM=figure('Visible','off','Position',[0 0 260*ncols 320*nrows]);
cmap_out=[0 0 0; OUT_CLASS_COLORS];
for di=1:nDates
    if isempty(s2_class_maps{di}), continue; end
    d=DATES_TO_CLASSIFY{di};
    subplot(nrows,ncols,di);
    cm_shift=s2_class_maps{di}+1; cm_shift(~reef_mask)=1;
    class_rgb=ind2rgb(cm_shift,cmap_out);
    blend=0.65.*class_rgb + 0.35.*rgb_cache{di};
    blend=max(0,min(1,blend));
    imshow(blend);
    wave_str=''; if ismember(d,WAVE_FLAG_DATES), wave_str=' '; end
    title(sprintf('%s%s\nBL:%.1f%% HY:%.1f%%',...
        d,wave_str,bleach_pct(di),healthy_pct(di)),...
        'FontSize',7,'Interpreter','none');
end
axL=axes(figM,'Visible','off','Position',[0 0 1 0.04]);
hold(axL,'on');
for ci=1:numel(OUT_CLASS_IDS)
    patch(axL,nan,nan,OUT_CLASS_COLORS(ci,:),'EdgeColor','none',...
        'DisplayName',OUT_CLASS_NAMES{ci});
end
legend(axL,'Location','north','Orientation','horizontal','FontSize',8,'Box','off');
sgtitle('Stage 2 Classification - All Dates','FontSize',12,'FontWeight','bold');
exportgraphics(figM,fullfile(outDir,'s2_mosaic_all.png'),'Resolution',130);
close(figM);
logW('Mosaic saved.\n\n');
end

%% "" Bleaching time series """"""""""""""""""""""""""""""""""""
if ~FAST_MODE
logW('Generating bleaching time series...\n');
date_nums=zeros(1,nDates);
for di=1:nDates
    d=DATES_TO_CLASSIFY{di};
    try
        date_nums(di)=datenum(datetime(d,'InputFormat','yyyyMMdd'));
    catch
        date_nums(di)=di;
    end
end

figT=figure('Visible','off','Position',[0 0 950 450]);
hold on;

% Shade bleaching period
bl_nums=date_nums(ismember(DATES_TO_CLASSIFY,...
    {'20240514','20240613','20240620'}));
if ~isempty(bl_nums)
    fill([min(bl_nums)-15 max(bl_nums)+15 max(bl_nums)+15 min(bl_nums)-15],...
        [0 0 100 100],[1 0.9 0.85],'EdgeColor','none','FaceAlpha',0.5,...
        'DisplayName','Bleaching period');
end

% Plot BLEACHED % - wave-flagged dates as dashed markers
for di=1:nDates
    if isempty(s2_class_maps{di}), continue; end
    d=DATES_TO_CLASSIFY{di};
    is_wf=ismember(d,WAVE_FLAG_DATES);
    if is_wf
        plot(date_nums(di),bleach_pct(di),'v','Color',[0.8 0.3 0.1],...
            'MarkerSize',9,'MarkerFaceColor',[1 0.7 0.6],'LineWidth',1.5);
    else
        plot(date_nums(di),bleach_pct(di),'o','Color',[0.8 0.2 0.1],...
            'MarkerSize',9,'MarkerFaceColor',[0.8 0.2 0.1]);
    end
    text(date_nums(di),bleach_pct(di)+2,sprintf('%.1f%%',bleach_pct(di)),...
        'HorizontalAlignment','center','FontSize',7,...
        'Color',iif(is_wf,[0.6 0.2 0],'k'));
end
% Connect with line
valid=~cellfun(@isempty,s2_class_maps);
plot(date_nums(valid),bleach_pct(valid),'r-','LineWidth',1.5,...
    'DisplayName','BLEACHED % of coral mask');

% Plot HEALTHY %
plot(date_nums(valid),healthy_pct(valid),'g-o','LineWidth',1.5,...
    'MarkerFaceColor',[0 0.6 0],'MarkerSize',6,...
    'DisplayName','HEALTHY % of coral mask');

% Wave flag legend entry
plot(nan,nan,'v','Color',[0.8 0.3 0.1],'MarkerSize',8,...
    'MarkerFaceColor',[1 0.7 0.6],'DisplayName','Wave-flagged (upper bound)');

datetick('x','mmm-yy','keepticks');
xlabel('Date','FontSize',10); ylabel('% of coral habitat mask','FontSize',10);
ylim([0 110]);
title('Stage 2 - Coral Condition Over Time','FontSize',11);
legend('Location','best','FontSize',8); grid on;
set(gca,'GridAlpha',0.3);
exportgraphics(figT,fullfile(outDir,'s2_bleaching_timeseries.png'),'Resolution',150);
close(figT);
logW('Time series saved.\n\n');
end

%% "" Stacked area chart """""""""""""""""""""""""""""""""""""""
if ~FAST_MODE
logW('Generating area chart...\n');
figA=figure('Visible','off','Position',[0 0 max(800,nDates*80) 480]);
valid_idx=find(valid);
bar_data=reef_pct_s2(valid_idx,:);
b=bar(bar_data,'stacked','BarWidth',0.7);
for ci=1:numel(OUT_CLASS_IDS)
    b(ci).FaceColor=OUT_CLASS_COLORS(ci,:);
    b(ci).DisplayName=OUT_CLASS_NAMES{ci};
end
xticks(1:sum(valid));
xticklabels(arrayfun(@(i) iif(ismember(DATES_TO_CLASSIFY{i},WAVE_FLAG_DATES),...
    [DATES_TO_CLASSIFY{i} ' '],DATES_TO_CLASSIFY{i}),...
    valid_idx,'UniformOutput',false));
xtickangle(45);
ylabel('% of reef area','FontSize',10);
title('Stage 2 - Substrate + Condition by Date','FontSize',11);
legend('Location','eastoutside','FontSize',9);
grid on; set(gca,'GridAlpha',0.3);
exportgraphics(figA,fullfile(outDir,'s2_area_chart.png'),'Resolution',150);
close(figA);
logW('Area chart saved.\n\n');
end

%% "" Classification report """"""""""""""""""""""""""""""""""""
fid=fopen(fullfile(outDir,'s2_classification_report.txt'),'w');
fprintf(fid,'STAGE 2 CLASSIFICATION REPORT\n');
fprintf(fid,'Timestamp  : %s\n',runTimestamp);
fprintf(fid,'S2 model   : %s\n',s2mPath);
fprintf(fid,'OOB error  : %.2f%%\n',s2m.oob_err*100);
fprintf(fid,'UNK thresh : %.2f\n',S2_UNK_THRESHOLD);
fprintf(fid,'Features   : [%s]\n\n',strjoin(FEAT_NAMES_S2,', '));
fprintf(fid,'Baseline coral mask: %d px  %.1f m2\n\n',n_coral_mask,n_coral_mask*px_m2);
fprintf(fid,'Wave-flagged dates (upper bound estimates): %s\n',strjoin(WAVE_FLAG_DATES,', '));
fprintf(fid,'Confound rates: 20240514=58.5%% 20240613=95.9%% 20240620=99.1%%\n\n');

fprintf(fid,'%-10s  %-9s  %-8s  %8s  %14s  %10s  %11s\n',...
    'Date','Class','WaveFlag','Pixels','Area_m2','Pct_reef','Pct_coral');
CORAL_CLASSES = {'HEALTHY','BLEACHED'};
for k_=1:numel(area_stats)
    as=area_stats(k_);
    wf=iif(as.wave_flag,'YES','no');
    if ismember(as.class,CORAL_CLASSES)
        fprintf(fid,'%-10s  %-9s  %-8s  %8d  %14.1f  %9.1f%%  %10.1f%%\n',...
            as.date,as.class,wf,as.px,as.m2,as.pct_reef,as.pct_coral);
    else
        fprintf(fid,'%-10s  %-9s  %-8s  %8d  %14.1f  %9.1f%%  %11s\n',...
            as.date,as.class,wf,as.px,as.m2,as.pct_reef,'-');
    end
end

fprintf(fid,'\nBleaching severity summary (within coral mask):\n');
fprintf(fid,'%-10s  %12s  %12s  %10s\n','Date','BLEACHED_%','HEALTHY_%','WaveFlag');
for di=1:nDates
    if isempty(s2_class_maps{di}), continue; end
    wf=iif(ismember(DATES_TO_CLASSIFY{di},WAVE_FLAG_DATES),' YES','no');
    fprintf(fid,'%-10s  %11.1f%%  %11.1f%%  %10s\n',...
        DATES_TO_CLASSIFY{di},bleach_pct(di),healthy_pct(di),wf);
end
fclose(fid);

%% "" Crop outputs to GeoFeature bounding box """""""""""""""""
% Crop all S2CLASS_*.tif to rectangular extent defined by
% GeoFeature.geojson. Runtime check: bbox must be inside image extent.
if ~FAST_MODE
logW('Cropping outputs to GeoFeature bounding box...\n');
gfPath=fullfile(cfg.ACA_DIR,'GeoFeature.geojson');
cropDir=fullfile(outDir,'Cropped');
do_crop=false;
if ~exist(gfPath,'file')
    logW('WARNING: GeoFeature.geojson not found - skipping crop\n');
else
    try
        gfData_=jsondecode(fileread(gfPath));
        gf_lons_=[]; gf_lats_=[];
        for fi_=1:numel(gfData_.features)
            % Use only the polygon named 'Reef' as crop boundary
            prop_=gfData_.features(fi_).properties;
            cls_='';
            if isfield(prop_,'Class'), cls_=prop_.Class;
            elseif isfield(prop_,'class'), cls_=prop_.class; end
            if ~strcmpi(cls_,'Reef'), continue; end
            rings__=extractRings(gfData_.features(fi_).geometry);
            for ri_=1:numel(rings__)
                ring__=rings__{ri_};
                gf_lons_=[gf_lons_;ring__(:,1)];
                gf_lats_=[gf_lats_;ring__(:,2)];
            end
        end
        if isempty(gf_lons_)
            error('No polygon with Class=Reef found in GeoFeature.geojson');
        end
        % Convert bbox corners to UTM47N
        [x_sw,y_sw]=latlon2utm47(min(gf_lats_),min(gf_lons_));
        [x_ne,y_ne]=latlon2utm47(max(gf_lats_),max(gf_lons_));
        gf_xmin=min(x_sw,x_ne); gf_xmax=max(x_sw,x_ne);
        gf_ymin=min(y_sw,y_ne); gf_ymax=max(y_sw,y_ne);
        logW('  GeoFeature UTM bbox: X[%.1f %.1f]  Y[%.1f %.1f]\n',...
            gf_xmin,gf_xmax,gf_ymin,gf_ymax);

        % Get image extent from R_ref
        img_xmin=R_ref.XWorldLimits(1); img_xmax=R_ref.XWorldLimits(2);
        img_ymin=R_ref.YWorldLimits(1); img_ymax=R_ref.YWorldLimits(2);
        logW('  Image UTM extent:    X[%.1f %.1f]  Y[%.1f %.1f]\n',...
            img_xmin,img_xmax,img_ymin,img_ymax);

        % Check bbox is inside image extent
        inside = gf_xmin>=img_xmin && gf_xmax<=img_xmax && ...
                 gf_ymin>=img_ymin && gf_ymax<=img_ymax;
        if ~inside
            logW('  WARNING: GeoFeature bbox partially outside image extent\n');
            logW('  Clamping to image extent - output may be smaller than full bbox\n');
            gf_xmin=max(gf_xmin,img_xmin); gf_xmax=min(gf_xmax,img_xmax);
            gf_ymin=max(gf_ymin,img_ymin); gf_ymax=min(gf_ymax,img_ymax);
        else
            logW('  Bbox is fully inside image extent "\n');
        end

        % Convert bbox to pixel row/col indices
        cw=R_ref.CellExtentInWorldX; ch_=R_ref.CellExtentInWorldY;
        col1=max(1, round((gf_xmin-img_xmin)/cw)+1);
        col2=min(W, round((gf_xmax-img_xmin)/cw));
        row1=max(1, round((img_ymax-gf_ymax)/ch_)+1);
        row2=min(H, round((img_ymax-gf_ymin)/ch_));
        logW('  Crop rows %d:%d  cols %d:%d  -> %dx%d px\n',...
            row1,row2,col1,col2,col2-col1+1,row2-row1+1);

        if col2>col1 && row2>row1
            mkdir(cropDir);
            do_crop=true;
        else
            logW('  WARNING: Invalid crop region - skipping\n');
        end
    catch ME
        logW('  ERROR reading GeoFeature.geojson: %s\n',ME.message);
    end
end

if do_crop
    % Build cropped spatial reference
    x_origin = img_xmin + (col1-1)*cw + cw/2;
    y_origin  = img_ymax - (row1-1)*ch_ - ch_/2;
    R_crop=R_ref;
    R_crop.XWorldLimits=[img_xmin+(col1-1)*cw, img_xmin+col2*cw];
    R_crop.YWorldLimits=[img_ymax-row2*ch_,     img_ymax-(row1-1)*ch_];
    R_crop.RasterSize=[row2-row1+1, col2-col1+1];

    n_cropped=0;
    for di=1:nDates
        if isempty(s2_class_maps{di}), continue; end
        d=DATES_TO_CLASSIFY{di};
        cropped=s2_class_maps{di}(row1:row2,col1:col2);
        saveUint8Tif(fullfile(cropDir,['S2CLASS_' d '_crop.tif']),cropped,R_crop);
        n_cropped=n_cropped+1;
    end
    logW('  Cropped %d S2CLASS tifs -> %s\n\n',n_cropped,cropDir);
else
    logW('  Crop skipped.\n\n');
end
end

%% "" Phase handoff """"""""""""""""""""""""""""""""""""""""""""
save(fullfile(cfg.OUT_DIR,'phase_stage2_classify.mat'),...
    'outDir','DATES_TO_CLASSIFY','OUT_CLASS_NAMES','OUT_CLASS_IDS',...
    'bleach_pct','healthy_pct','reef_pct_s2','baseline_coral_mask',...
    'WAVE_FLAG_DATES','n_coral_mask','px_m2');

classified_count = sum(~cellfun(@isempty, s2_class_maps));
logW('\n========== STAGE 2 SUMMARY ==========\n');
logW('Dates classified  : %d\n',classified_count);
logW('Coral mask        : %d px  %.1f m2\n',n_coral_mask,n_coral_mask*px_m2);
logW('\n%-12s  %10s  %10s  %10s\n','Date','BLEACHED%','HEALTHY%','WaveFlag');
for di=1:nDates
    if isempty(s2_class_maps{di}), continue; end
    wf=iif(ismember(DATES_TO_CLASSIFY{di},WAVE_FLAG_DATES),'','');
    logW('%-12s  %9.1f%%  %9.1f%%  %10s\n',...
        DATES_TO_CLASSIFY{di},bleach_pct(di),healthy_pct(di),wf);
end
logW('\nOutput: %s\n',outDir);
logW('NEXT: Diag_SpectralProfile.m or Analyze_Change.m\n');
logW('======================================\n');

%% ===================== LOCAL FUNCTIONS =====================

function saveUint8Tif(path,A,R)
    try geotiffwrite(path,uint8(A),R); return; catch, end
    try geotiffwrite(path,uint8(A),R,'CoordRefSysCode',32647); return; catch, end
    try
        t=Tiff(path,'w'); [H_,W_]=size(A);
        t.setTag('Photometric',1); t.setTag('ImageLength',H_);
        t.setTag('ImageWidth',W_); t.setTag('BitsPerSample',8);
        t.setTag('SampleFormat',1); t.setTag('SamplesPerPixel',1);
        t.setTag('Compression',5); t.setTag('PlanarConfig',1);
        t.write(uint8(A)); t.close();
    catch ME, warning('saveUint8Tif: %s',ME.message); end
end

function saveMultiBandTif(path,A,R)
    if ~isa(A,'single'), A=single(A); end
    try geotiffwrite(path,A,R); return; catch, end
    try geotiffwrite(path,A,R,'CoordRefSysCode',32647); return; catch, end
    try
        [H_,W_,nb]=size(A); t=Tiff(path,'w');
        t.setTag('Photometric',1); t.setTag('ImageLength',H_);
        t.setTag('ImageWidth',W_); t.setTag('BitsPerSample',32);
        t.setTag('SampleFormat',3); t.setTag('SamplesPerPixel',nb);
        t.setTag('Compression',5); t.setTag('PlanarConfig',1);
        t.write(A); t.close();
    catch ME, error('saveMultiBandTif: %s',ME.message); end
end

function logWrite(fid,varargin)
    msg=sprintf(varargin{:}); fprintf('%s',msg);
    if isnumeric(fid) && fid>2
        try, fprintf(fid,'%s',msg); catch, end
    end
end

function rgb_out=makeDisplayRGB(rgb_in, valid_mask)
    rgb_in=single(rgb_in);
    rgb_out=zeros(size(rgb_in),'single');
    if nargin<2 || isempty(valid_mask)
        valid_mask=all(isfinite(rgb_in),3);
    end
    valid_mask=logical(valid_mask);
    for ch=1:size(rgb_in,3)
        sl=rgb_in(:,:,ch);
        v=sl(valid_mask & isfinite(sl));
        v=v(v>0);
        if isempty(v)
            continue;
        end
        lo=prctile(v,1);
        hi=prctile(v,99.5);
        if ~isfinite(lo), lo=min(v); end
        if ~isfinite(hi), hi=max(v); end
        if hi<=lo
            hi=max(v);
            lo=min(v);
        end
        sl=(sl-lo)./max(hi-lo,1e-6);
        sl=max(0,min(1,sl));
        sl(~valid_mask)=0;
        rgb_out(:,:,ch)=sl;
    end
    gamma_=0.85;
    rgb_out=max(rgb_out,0).^gamma_;
end

function v=iif(c,a,b); if c,v=a; else,v=b; end; end


