%% Classify_Stage1.m
% ============================================================
% PHASE 9 - Apply RF Classifier to All Dates
%
% Loads the trained RF model (Stage 1 or legacy Phase 8) and
% classifies each date's Phase 5 feature stack. Feature
% subsetting is controlled by FEAT_IDX loaded from the model
% file - Phase 5 stacks are never modified.
%
% MODEL LOOKUP ORDER:
%   1. phase_stage1_output.mat  -> Phase_Stage1 folder (preferred)
%   2. phase8_output.mat        -> Phase8_TrainRF folder (legacy)
%   Model file: rf_stage1_model.mat  or  rf_model.mat
%
% FEATURE SUBSETTING:
%   FEAT_IDX loaded from model mat - applied to each date stack
%   immediately after readgeoraster. Legacy models without
%   FEAT_IDX default to using all bands (backward compatible).
%
% UNK CLASS:
%   Pixels where max RF class probability < UNK_THRESHOLD are
%   relabelled as UNK=4 (uncertain). Threshold default 0.55.
%   UNK pixels are excluded from area statistics and Stage 2.
%
% CLASS SCHEMA:
%   CORAL=1   SND=2   RBL=3   UNK=4 (uncertain)
%
% OUTPUTS -> OUT_DIR/Phase9_Classify/run_YYYYMMDD_HHMMSS/
%   CLASS_YYYYMMDD.tif              uint8 class map per date
%   CONF_YYYYMMDD.tif               float32 confidence map per date
%   VIZ_YYYYMMDD.png                RGB + class overlay + conf (per date)
%   classify_mosaic_all.png         all dates classified in one figure
%   classify_area_chart.png         stacked area bar chart across dates
%   classify_coral_timeseries.png   CORAL % of reef over time
%   classify_quicklook.png          legacy 2-row RGB+class mosaic
%   classification_report.txt       area statistics per date per class
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Classify_Stage1.m - Phase 9 - %s\n\n',runTimestamp);

cfg = loadPipelineConfigPortable();

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
UNK_THRESHOLD = 0.55;   % max RF probability below this -> UNK class
CLASS_UNK     = 4;      % UNK class ID (appended after RF classes)
UNK_COLOR     = [0.65 0.65 0.65];  % grey

% Dates to classify - confirmed PIF-processed dates only.
% Excludes: 20240408 (REF_DATE - anomalous Stage 1 output, not a scientific observation)
% Leave empty [] to classify all FEATURES_*.tif found in featDir.
DATES_TO_CLASSIFY = {...
    '20230404','20230509','20230519','20230629',...   % pre-bleaching healthy
    '20240307','20240317',...                         % calm pre-bleaching
    '20240514','20240613','20240620',...              % bleaching era
    '20240926','20241111','20250518'};                % recovery + post-bleaching

%% "" Find RF model """"""""""""""""""""""""""""""""""""""""""""
modelDir='';
s1Path=fullfile(cfg.OUT_DIR,'phase_stage1_output.mat');
p8Path=fullfile(cfg.OUT_DIR,'phase8_output.mat');
if exist(s1Path,'file')
    try s1=load(s1Path); if exist(s1.outDir,'dir'), modelDir=s1.outDir; end; catch, end
end
if isempty(modelDir) && exist(p8Path,'file')
    try p8=load(p8Path); if exist(p8.outDir,'dir'), modelDir=p8.outDir; end; catch, end
end
if isempty(modelDir)
    base=fullfile(cfg.OUT_DIR,'Phase_Stage1');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); modelDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(modelDir)
    base=fullfile(cfg.OUT_DIR,'Phase8_TrainRF');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); modelDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(modelDir)||~exist(modelDir,'dir')
    modelDir=uigetdir(cfg.OUT_DIR,'Select RF model folder (Stage1 or Phase8)');
    if isequal(modelDir,0), error('No folder selected.'); end
end

% Load model
modelPath=fullfile(modelDir,'rf_stage1_model.mat');
if ~exist(modelPath,'file'), modelPath=fullfile(modelDir,'rf_model.mat'); end
if ~exist(modelPath,'file'), error('No RF model mat found in: %s',modelDir); end
md=load(modelPath);
rf_model     = md.rf_model;
FEAT_NAMES   = md.FEAT_NAMES;
CLASS_NAMES  = md.CLASS_NAMES;   % e.g. {'CORAL','SND','RBL'}
CLASS_IDS    = md.CLASS_IDS;     % e.g. [1 2 3]
CLASS_COLORS = md.CLASS_COLORS;
N_FEAT       = md.N_FEAT;
if isfield(md,'FEAT_IDX')
    FEAT_IDX = md.FEAT_IDX;
    logW('FEAT_IDX: [%s]\n',num2str(FEAT_IDX));
else
    FEAT_IDX = 1:N_FEAT;
    logW('WARNING: FEAT_IDX not in model - using all %d bands\n',N_FEAT);
end
logW('RF model: OA=%.3f  Kappa=%.4f  features=[%s]\n\n',...
    md.OA,md.kappa,strjoin(FEAT_NAMES,','));

% Extended class schema including UNK
CLASS_NAMES_EXT  = [CLASS_NAMES, {'UNK'}];
CLASS_IDS_EXT    = [CLASS_IDS,    CLASS_UNK];
CLASS_COLORS_EXT = [CLASS_COLORS; UNK_COLOR];
n_cls            = numel(CLASS_NAMES);
n_cls_ext        = numel(CLASS_NAMES_EXT);

%% "" Find Phase 5 feature directory """"""""""""""""""""""""""
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

% Find PIF directory for RGB quicklooks
pif_dir='';
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
if exist(pifBase,'dir')
    runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
    if ~isempty(runs),[~,si_]=sort({runs.name}); pif_dir=fullfile(pifBase,runs(si_(end)).name); end
end

outDir=fullfile(cfg.OUT_DIR,'Phase9_Classify',['run_' runTimestamp]);
mkdir(outDir);
logW('Feature input: %s\n',featDir);
logW('Output:        %s\n\n',outDir);
logW('UNK threshold: %.2f\n\n',UNK_THRESHOLD);

%% "" Load masks """""""""""""""""""""""""""""""""""""""""""""""
reef_mask=[]; deep_mask=[];
if exist(p5Path,'file')
    try
        p5=load(p5Path);
        if isfield(p5,'reef_mask'), reef_mask=p5.reef_mask; end
        if isfield(p5,'deep_mask'), deep_mask=p5.deep_mask; end
    catch, end
end

%% "" Find all feature stacks """"""""""""""""""""""""""""""""""
ff=dir(fullfile(featDir,'FEATURES_*.tif'));
if isempty(ff), error('No FEATURES_*.tif in %s',featDir); end
dates_feat=cellfun(@(n) regexp(n,'FEATURES_(\d{8})\.tif','tokens','once'),...
    {ff.name},'UniformOutput',false);
dates_feat=cellfun(@(t) t{1},dates_feat,'UniformOutput',false);
[dates_feat,si]=sort(dates_feat); ff=ff(si);
% Filter to intended dates only
if ~isempty(DATES_TO_CLASSIFY)
    keep=ismember(dates_feat,DATES_TO_CLASSIFY);
    if ~any(keep)
        logW('WARNING: DATES_TO_CLASSIFY matched no files - classifying all dates\n');
    else
        dates_feat=dates_feat(keep); ff=ff(keep);
        logW('Date filter applied: %d of %d stacks selected\n',...
            sum(keep),numel(keep));
    end
end
nDates=numel(ff);
logW('Dates to classify: %s\n\n',strjoin(dates_feat,', '));

%% "" Classify each date """""""""""""""""""""""""""""""""""""""
area_stats  = struct('date',{},'class',{},'px',{},'ha',{},'pct',{});
class_maps  = cell(nDates,1);
conf_maps   = cell(nDates,1);
rgb_cache   = cell(nDates,1);   % PIF RGB for visualisation
reef_pct    = zeros(nDates,n_cls_ext);  % % of reef per class per date
px_ha_global= 0;

for di=1:nDates
    logW('============================================================\n');
    logW('[%d/%d] %s\n',di,nDates,dates_feat{di});

    [Ai,Ri]=readgeoraster(fullfile(ff(di).folder,ff(di).name));
    Ai=single(Ai(:,:,FEAT_IDX));
    [H,W,~]=size(Ai);

    if isempty(reef_mask)
        logW('  Re-rasterising GeoFeature/Reef...\n');
        reef_mask=rasterizeGeoJSONClass(cfg.REEF_MASK_PATH,H,W,Ri,'Reef');
    end

    px_ha=abs(Ri.CellExtentInWorldX*Ri.CellExtentInWorldY)/10000;
    if px_ha_global==0, px_ha_global=px_ha; end

    % "" RF classification """""""""""""""""""""""""""""""""""""
    classify_mask=reef_mask & all(isfinite(Ai),3);
    idx_px=find(classify_mask);
    n_px=numel(idx_px);
    logW('  Classifiable pixels: %d\n',n_px);

    class_map=zeros(H,W,'uint8');
    conf_map =zeros(H,W,'single');

    if n_px>0
        if isfield(cfg,'RF') && isfield(cfg.RF,'STAGE1_PREDICT_BATCH') ...
                && ~isempty(cfg.RF.STAGE1_PREDICT_BATCH)
            BATCH = cfg.RF.STAGE1_PREDICT_BATCH;
        else
            BATCH = 10000;
        end
        logW('  Predict batch size: %d\n', BATCH);
        Ypred_all=zeros(n_px,1,'uint8');
        prob_all =zeros(n_px,1,'single');
        for b_start=1:BATCH:n_px
            b_end=min(b_start+BATCH-1,n_px);
            b_idx=idx_px(b_start:b_end);
            Xb=zeros(numel(b_idx),N_FEAT,'single');
            for f=1:N_FEAT
                sl=Ai(:,:,f); Xb(:,f)=sl(b_idx);
            end
            [Ypred_cat,prob_b]=predict(rf_model,Xb);
            Ypred_all(b_start:b_end)=uint8(str2double(cellstr(Ypred_cat)));
            prob_all(b_start:b_end) =single(max(prob_b,[],2));
            clear Xb Ypred_cat prob_b;
        end
        class_map(idx_px)=Ypred_all;
        conf_map( idx_px)=prob_all;

        % "" UNK threshold """"""""""""""""""""""""""""""""""""""
        % Pixels below confidence threshold -> UNK=4
        unk_mask=classify_mask & (conf_map<UNK_THRESHOLD);
        class_map(unk_mask)=uint8(CLASS_UNK);
        n_unk=nnz(unk_mask);
        logW('  UNK pixels (conf<%.2f): %d (%.1f%% of reef)\n',...
            UNK_THRESHOLD,n_unk,100*n_unk/max(n_px,1));
    end

    class_maps{di}=class_map;
    conf_maps{di} =conf_map;

    % "" Save GeoTIFFs """"""""""""""""""""""""""""""""""""""""""
    saveUint8Tif(fullfile(outDir,['CLASS_' dates_feat{di} '.tif']),class_map,Ri);
    saveMultiBandTif(fullfile(outDir,['CONF_' dates_feat{di} '.tif']),conf_map,Ri);
    logW('  Saved CLASS_ and CONF_ tifs\n');

    % "" Area statistics """"""""""""""""""""""""""""""""""""""""
    reef_px=nnz(reef_mask);
    logW('  Area stats:\n');
    for ci=1:n_cls_ext
        cid=CLASS_IDS_EXT(ci); cname=CLASS_NAMES_EXT{ci};
        npx=nnz(class_map==cid);
        ha=npx*px_ha;
        pct=100*npx/max(reef_px,1);
        reef_pct(di,ci)=pct;
        logW('    %-6s  %6d px  %7.2f ha  %5.1f%%\n',cname,npx,ha,pct);
        area_stats(end+1)=struct('date',dates_feat{di},'class',cname,...
            'px',npx,'ha',ha,'pct',pct);
    end
    logW('\n');

    % "" Load PIF RGB for visualisation """""""""""""""""""""""""
    rgb_img=[];
    if ~isempty(pif_dir)
        pif_f=fullfile(pif_dir,[dates_feat{di} '_PIF.tif']);
        if ~exist(pif_f,'file')
            pf=dir(fullfile(pif_dir,['*' dates_feat{di} '*.tif']));
            if ~isempty(pf), pif_f=fullfile(pif_dir,pf(1).name); end
        end
        if exist(pif_f,'file')
            try
                [Ap,~]=readgeoraster(pif_f); Ap=single(Ap)/10000;
                r_=Ap(:,:,cfg.BANDS.RED);
                g_=Ap(:,:,cfg.BANDS.GREEN1);
                b_=Ap(:,:,cfg.BANDS.BLUE);
                rgb_img=cat(3,r_,g_,b_);
                for ch=1:3   % min/max stretch per channel
                    sl=rgb_img(:,:,ch); v=sl(isfinite(sl)&sl>0);
                    if ~isempty(v)
                        lo=min(v); hi=max(v);
                        sl=max(0,min(1,(sl-lo)/(max(hi-lo,1e-4))));
                        rgb_img(:,:,ch)=sl;
                    end
                end
            catch, end
        end
    end
    if isempty(rgb_img), rgb_img=zeros(H,W,3,'single'); end
    rgb_cache{di}=rgb_img;

    % "" Per-date visualisation (3 panels) """""""""""""""""""""
    % Panel 1: RGB  Panel 2: Class overlay on RGB  Panel 3: Confidence
    figV=figure('Visible','off','Position',[0 0 1200 450]);

    % Panel 1 - RGB
    subplot(1,3,1);
    imshow(rgb_img); title(sprintf('%s - RGB',dates_feat{di}),'FontSize',9);

    % Panel 2 - Class overlay (70% class, 30% RGB)
    subplot(1,3,2);
    cmap_ext=[0 0 0; CLASS_COLORS_EXT];   % index 0 = black (outside reef)
    cm_shift=class_map+1;
    cm_shift(~reef_mask)=1;               % outside reef -> black slot
    class_rgb=ind2rgb(cm_shift,cmap_ext);
    % Blend: 70% class colour, 30% RGB for texture visibility
    blend=0.70.*class_rgb + 0.30.*rgb_img;
    blend=max(0,min(1,blend));
    imshow(blend);
    title('Classification (70% overlay)','FontSize',9);
    % Legend patches
    hold on;
    for ci=1:n_cls_ext
        patch(nan,nan,CLASS_COLORS_EXT(ci,:),...
            'DisplayName',CLASS_NAMES_EXT{ci},'EdgeColor','none');
    end
    legend('Location','southoutside','Orientation','horizontal',...
        'FontSize',7,'Box','off');

    % Panel 3 - Confidence map
    subplot(1,3,3);
    conf_disp=conf_map; conf_disp(~reef_mask)=NaN;
    imagesc(conf_disp,[0 1]); axis image off;
    colormap(gca,hot); cb=colorbar; cb.Label.String='Max class probability';
    hold on;
    % Mark UNK pixels in cyan
    unk_ovl=class_map==CLASS_UNK & reef_mask;
    if nnz(unk_ovl)>0
        [ru,cu]=find(unk_ovl);
        scatter(cu,ru,1,[0 0.9 0.9],'.','MarkerEdgeAlpha',0.3);
    end
    title(sprintf('Confidence  UNK=%d px',nnz(unk_ovl)),'FontSize',9);

    sgtitle(sprintf('Stage 1 Classification - %s',dates_feat{di}),...
        'FontSize',11,'FontWeight','bold');
    exportgraphics(figV,fullfile(outDir,['VIZ_' dates_feat{di} '.png']),'Resolution',150);
    close(figV);
    logW('  Saved VIZ_%s.png\n\n',dates_feat{di});
end

%% "" Mosaic: all dates classified """"""""""""""""""""""""""""
logW('Generating classification mosaic (all dates)...\n');
ncols=min(nDates,6); nrows=ceil(nDates/ncols);
figM=figure('Visible','off','Position',[0 0 260*ncols 320*nrows]);
cmap_ext=[0 0 0; CLASS_COLORS_EXT];
for di=1:nDates
    subplot(nrows,ncols,di);
    cm_shift=class_maps{di}+1; cm_shift(~reef_mask)=1;
    class_rgb=ind2rgb(cm_shift,cmap_ext);
    % Light RGB blend for context
    blend=0.65.*class_rgb + 0.35.*rgb_cache{di};
    blend=max(0,min(1,blend));
    imshow(blend);
    % UNK count annotation
    n_unk=nnz(class_maps{di}==CLASS_UNK & reef_mask);
    coral_pct=reef_pct(di,CLASS_IDS==1);
    if isempty(coral_pct), coral_pct=0; end
    title(sprintf('%s\nCORAL %.1f%%  UNK %d px',...
        dates_feat{di},coral_pct,n_unk),'FontSize',7,'Interpreter','none');
end
% Shared legend
axL=axes(figM,'Visible','off','Position',[0 0 1 0.05]);
hold(axL,'on');
for ci=1:n_cls_ext
    patch(axL,nan,nan,CLASS_COLORS_EXT(ci,:),...
        'DisplayName',CLASS_NAMES_EXT{ci},'EdgeColor','none');
end
legend(axL,'Location','north','Orientation','horizontal','FontSize',8,'Box','off');
sgtitle('Stage 1 Classification - All Dates','FontSize',12,'FontWeight','bold');
exportgraphics(figM,fullfile(outDir,'classify_mosaic_all.png'),'Resolution',130);
close(figM);
logW('Mosaic saved.\n\n');

%% "" Area bar chart (stacked, all dates) """""""""""""""""""""
logW('Generating area chart...\n');
figA=figure('Visible','off','Position',[0 0 max(800,nDates*70) 480]);

% reef_pct: nDates x n_cls_ext  (% of reef per class)
% Stacked bar - classes in order CORAL, SND, RBL, UNK
bar_data=reef_pct;   % already in CLASS_IDS_EXT order
b=bar(bar_data,'stacked','BarWidth',0.7);
for ci=1:n_cls_ext
    b(ci).FaceColor=CLASS_COLORS_EXT(ci,:);
    b(ci).DisplayName=CLASS_NAMES_EXT{ci};
end
xticks(1:nDates); xticklabels(dates_feat); xtickangle(45);
ylabel('% of reef area'); ylim([0 100]);
title('Stage 1 - Substrate Area by Date','FontSize',11);
legend('Location','eastoutside','FontSize',9);
grid on; set(gca,'GridAlpha',0.3);
exportgraphics(figA,fullfile(outDir,'classify_area_chart.png'),'Resolution',150);
close(figA);
logW('Area chart saved.\n\n');

%% "" CORAL time series """"""""""""""""""""""""""""""""""""""""
logW('Generating CORAL time series...\n');

% Extract CORAL % - find column index matching CLASS_IDS==1
coral_col=find(CLASS_IDS_EXT==1,1);
coral_ts=reef_pct(:,coral_col);

% Parse dates for x-axis
date_nums=zeros(1,nDates);
for di=1:nDates
    d=dates_feat{di};
    try
        date_nums(di)=datenum(datetime(d,'InputFormat','yyyyMMdd'));
    catch
        date_nums(di)=di;
    end
end

figT=figure('Visible','off','Position',[0 0 900 420]);
hold on;

% Shade bleaching period - dates with year 2024 between May and Aug
is_bleach=false(1,nDates);
for di=1:nDates
    d=dates_feat{di};
    yr=str2double(d(1:4)); mo=str2double(d(5:6));
    if yr==2024 && mo>=5 && mo<=8, is_bleach(di)=true; end
end
if any(is_bleach)
    bl_dates=date_nums(is_bleach);
    fill([min(bl_dates)-15 max(bl_dates)+15 max(bl_dates)+15 min(bl_dates)-15],...
        [0 0 100 100],[1 0.9 0.8],'EdgeColor','none','FaceAlpha',0.5,...
        'DisplayName','Bleaching period');
end

% Plot CORAL %
plot(date_nums,coral_ts,'o-','Color',[0 0.6 0],'LineWidth',2,...
    'MarkerSize',8,'MarkerFaceColor',[0 0.6 0],'DisplayName','CORAL %');

% Annotate each point
for di=1:nDates
    text(date_nums(di),coral_ts(di)+1.5,sprintf('%.1f%%',coral_ts(di)),...
        'HorizontalAlignment','center','FontSize',8,'Color',[0 0.4 0]);
end

datetick('x','mmm-yy','keepticks');
xlabel('Date'); ylabel('% of reef area');
ylim([0 max(coral_ts)*1.3+5]);
title('CORAL Coverage Over Time - Stage 1 Classification','FontSize',11);
legend('Location','best','FontSize',9); grid on;
set(gca,'GridAlpha',0.3);
exportgraphics(figT,fullfile(outDir,'classify_coral_timeseries.png'),'Resolution',150);
close(figT);
logW('CORAL time series saved.\n\n');

%% "" Legacy quicklook mosaic (2-row RGB + class) """"""""""""""
if ~isfield(cfg,'SAVE_LEGACY_QUICKLOOK') || cfg.SAVE_LEGACY_QUICKLOOK
    logW('Generating legacy quicklook mosaic...\n');
    figQ=figure('Visible','off','Position',[0 0 200*nDates 420]);
    for di=1:nDates
        subplot(2,nDates,di);
        imshow(rgb_cache{di});
        title(dates_feat{di},'FontSize',7,'Interpreter','none');

        subplot(2,nDates,nDates+di);
        cm_shift=class_maps{di}+1; cm_shift(~reef_mask)=1;
        imshow(ind2rgb(cm_shift,cmap_ext));
        coral_pct=reef_pct(di,coral_col);
        title(sprintf('CORAL %.0f%%',coral_pct),'FontSize',7);
    end
    annotation('textbox',[0.01 0.01 0.98 0.04],...
        'String',strjoin(arrayfun(@(i)sprintf('%s=%d',...
        CLASS_NAMES_EXT{i},CLASS_IDS_EXT(i)),1:n_cls_ext,'UniformOutput',false),'  '),...
        'FitBoxToText','off','EdgeColor','none','FontSize',8,'HorizontalAlignment','center');
    exportgraphics(figQ,fullfile(outDir,'classify_quicklook.png'),'Resolution',120);
    close(figQ);
    logW('Quicklook saved.\n\n');
else
    logW('Skipping legacy quicklook mosaic (disabled in config).\n\n');
end

%% "" Classification report """"""""""""""""""""""""""""""""""""
fid=fopen(fullfile(outDir,'classification_report.txt'),'w');
fprintf(fid,'PHASE 9 - CLASSIFICATION REPORT\n');
fprintf(fid,'Timestamp : %s\n',runTimestamp);
fprintf(fid,'RF model  : %s\n',modelPath);
fprintf(fid,'OA=%.4f  Kappa=%.4f  OOB=%.4f\n',md.OA,md.kappa,md.oob_err);
fprintf(fid,'UNK threshold : %.2f\n',UNK_THRESHOLD);
fprintf(fid,'Features  : [%s]\n\n',strjoin(FEAT_NAMES,', '));
fprintf(fid,'%-10s  %-6s  %8s  %10s  %7s\n','Date','Class','Pixels','Area_ha','Pct');
for k_=1:numel(area_stats)
    fprintf(fid,'%-10s  %-6s  %8d  %10.2f  %6.1f%%\n',...
        area_stats(k_).date, area_stats(k_).class,...
        area_stats(k_).px,   area_stats(k_).ha, area_stats(k_).pct);
end
fclose(fid);

%% "" Crop outputs to GeoFeature bounding box """""""""""""""""
% Crop all CLASS_*.tif to rectangular extent from GeoFeature.geojson.
% Runtime check confirms bbox is inside image working space.
logW('Cropping CLASS outputs to GeoFeature bounding box...\n');
gfPath_=fullfile(cfg.ACA_DIR,'GeoFeature.geojson');
cropDir_=fullfile(outDir,'Cropped');
do_crop_=false; col1_=1; col2_=W; row1_=1; row2_=H; R_crop_=[];
if ~exist(gfPath_,'file')
    logW('  GeoFeature.geojson not found - skipping crop\n');
else
    try
        gfData_=jsondecode(fileread(gfPath_));
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
        [x_sw_,y_sw_]=latlon2utm47(min(gf_lats_),min(gf_lons_));
        [x_ne_,y_ne_]=latlon2utm47(max(gf_lats_),max(gf_lons_));
        gf_xmin_=min(x_sw_,x_ne_); gf_xmax_=max(x_sw_,x_ne_);
        gf_ymin_=min(y_sw_,y_ne_); gf_ymax_=max(y_sw_,y_ne_);

        % Use spatial ref from first loaded date
        if isempty(R_crop_) && ~isempty(ff)
            [~,Rtmp_]=readgeoraster(fullfile(ff(1).folder,ff(1).name));
        else
            Rtmp_=R_crop_;
        end
        img_xmin_=Rtmp_.XWorldLimits(1); img_xmax_=Rtmp_.XWorldLimits(2);
        img_ymin_=Rtmp_.YWorldLimits(1); img_ymax_=Rtmp_.YWorldLimits(2);
        logW('  GeoFeature UTM bbox: X[%.1f %.1f]  Y[%.1f %.1f]\n',...
            gf_xmin_,gf_xmax_,gf_ymin_,gf_ymax_);
        logW('  Image extent:        X[%.1f %.1f]  Y[%.1f %.1f]\n',...
            img_xmin_,img_xmax_,img_ymin_,img_ymax_);

        % Clamp to image extent if needed
        inside_=(gf_xmin_>=img_xmin_)&&(gf_xmax_<=img_xmax_)&&...
                 (gf_ymin_>=img_ymin_)&&(gf_ymax_<=img_ymax_);
        if ~inside_
            logW('  Bbox partially outside - clamping to image extent\n');
            gf_xmin_=max(gf_xmin_,img_xmin_); gf_xmax_=min(gf_xmax_,img_xmax_);
            gf_ymin_=max(gf_ymin_,img_ymin_); gf_ymax_=min(gf_ymax_,img_ymax_);
        else
            logW('  Bbox fully inside image extent "\n');
        end

        cw_=Rtmp_.CellExtentInWorldX; ch__=Rtmp_.CellExtentInWorldY;
        col1_=max(1,round((gf_xmin_-img_xmin_)/cw_)+1);
        col2_=min(W,round((gf_xmax_-img_xmin_)/cw_));
        row1_=max(1,round((img_ymax_-gf_ymax_)/ch__)+1);
        row2_=min(H,round((img_ymax_-gf_ymin_)/ch__));
        logW('  Crop rows %d:%d  cols %d:%d  -> %dx%d px\n',...
            row1_,row2_,col1_,col2_,col2_-col1_+1,row2_-row1_+1);

        if col2_>col1_ && row2_>row1_
            R_crop_=Rtmp_;
            R_crop_.XWorldLimits=[img_xmin_+(col1_-1)*cw_, img_xmin_+col2_*cw_];
            R_crop_.YWorldLimits=[img_ymax_-row2_*ch__, img_ymax_-(row1_-1)*ch__];
            R_crop_.RasterSize=[row2_-row1_+1, col2_-col1_+1];
            mkdir(cropDir_); do_crop_=true;
        else
            logW('  Invalid crop region - skipping\n');
        end
    catch ME
        logW('  ERROR: %s\n',ME.message);
    end
end

if do_crop_
    n_cr=0;
    for di=1:nDates
        cpath_=fullfile(outDir,['CLASS_' dates_feat{di} '.tif']);
        if ~exist(cpath_,'file'), continue; end
        [tmp_,~]=readgeoraster(cpath_);
        saveUint8Tif(fullfile(cropDir_,['CLASS_' dates_feat{di} '_crop.tif']),...
            uint8(tmp_(row1_:row2_,col1_:col2_)),R_crop_);
        n_cr=n_cr+1;
    end
    logW('  Cropped %d CLASS tifs -> %s\n\n',n_cr,cropDir_);
else
    logW('  Crop skipped.\n\n');
end

%% "" Phase handoff """"""""""""""""""""""""""""""""""""""""""""
save(fullfile(cfg.OUT_DIR,'phase9_output.mat'),...
    'outDir','dates_feat','CLASS_NAMES','CLASS_IDS',...
    'CLASS_NAMES_EXT','CLASS_IDS_EXT','CLASS_UNK','UNK_THRESHOLD',...
    'reef_pct','coral_ts');

logW('\n========== PHASE 9 SUMMARY ==========\n');
logW('Classified dates  : %d\n',nDates);
logW('UNK threshold     : %.2f\n',UNK_THRESHOLD);
logW('CORAL range       : %.1f%% - %.1f%%\n',min(coral_ts),max(coral_ts));
logW('Output: %s\n',outDir);
logW('NEXT: Analyze_Change.m (Phase 10)\n');
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



