%% Analyze_Change.m
% ============================================================
% CHANGE ANALYSIS - Stage 2 Bleaching Spatial & Temporal Analysis
%
% Operates on Stage 2 combined classification maps (S2CLASS_*.tif)
% and the baseline coral mask to produce spatial change products
% and temporal statistics for the manuscript.
%
% INPUT CLASS SCHEMA (from Classify_Stage2.m):
%   1=HEALTHY  2=BLEACHED  3=SND  4=RBL  5=UNK
%
% BASELINE FOR CHANGE:
%   Pre-bleaching consensus from 4 dates (20230404, 20230509,
%   20230519, 20230629). A pixel is "baseline HEALTHY" if it was
%   HEALTHY on >= 3/4 pre-bleaching dates (Moderate+ confidence).
%   Post-bleaching reference: 20250518 (healthy, no thermal stress).
%   This is more robust than using a single reference date.
%
% ANALYSES:
%   1. Bleaching frequency map - how many bleaching dates each
%      coral pixel was classified BLEACHED (0-4 dates)
%   2. Spatial overlap - pixel-level agreement between bleaching dates
%   3. Transition analysis - pre-bleaching HEALTHY -> each date class
%   4. Recovery trajectory - BLEACHED% decline from peak to recovery
%   5. Persistent bleaching map - pixels BLEACHED on ALL bleaching dates
%
% WAVE FLAG HANDLING:
%   20240613, 20240620 flagged as upper bound - reported separately
%   Primary bleaching estimate: 20240514
%
% OUTPUTS -> OUT_DIR/Change_Analysis/run_YYYYMMDD_HHMMSS/
%   BLEACH_FREQ.tif              uint8 bleaching frequency (0-4)
%   BLEACH_PERSISTENT.tif        uint8 pixels bleached on all primary dates
%   RECOVERY_MAP.tif             uint8 recovery state at 20240926
%   bleach_frequency_map.png     spatial map of bleaching frequency
%   bleach_spatial_overlap.png   Venn-style overlap between dates
%   recovery_trajectory.png      time series BLEACHED% primary dates
%   transition_heatmap.png       pre-bleach HEALTHY -> each date
%   spatial_comparison.png       side-by-side maps per bleaching date
%   change_report.txt            full numeric report (m)
%   change_output.mat            handoff
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Analyze_Change.m - %s\n\n',runTimestamp);

cfg = loadPipelineConfigPortable();

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
% Class IDs from Stage 2 combined schema
ID_HEALTHY  = 1;
ID_BLEACHED = 2;
ID_SND      = 3;
ID_RBL      = 4;
ID_UNK      = 5;
CLASS_NAMES = {'HEALTHY','BLEACHED','SND','RBL','UNK'};
CLASS_COLORS= [0.0  0.6  0.0;
               1.0  1.0  0.85;
               1.0  0.94 0.70;
               0.7  0.55 0.31;
               0.65 0.65 0.65];

% Date groups
PRE_BLEACH_DATES = {'20230404','20230509','20230519','20230629', '20240307'};
BLEACH_DATES_ALL = {'20240514','20240613','20240620'};
BLEACH_PRIMARY   = {'20240514'};
BLEACH_WAVE      = {'20240613','20240620'};   % upper bound
RECOVERY_DATES   = {'20240926','20241111'};
POST_BLEACH_DATES = {'20250518'};
ALL_DATES        = [PRE_BLEACH_DATES, BLEACH_DATES_ALL, RECOVERY_DATES, POST_BLEACH_DATES];

%% "" Find Stage 2 classification output """"""""""""""""""""""
s2ClassDir='';
s2Path=fullfile(cfg.OUT_DIR,'phase_stage2_classify.mat');
if exist(s2Path,'file')
    try s2c=load(s2Path); if exist(s2c.outDir,'dir'), s2ClassDir=s2c.outDir; end; catch, end
end
if isempty(s2ClassDir)
    base=fullfile(cfg.OUT_DIR,'Phase_Stage2_Classify');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); s2ClassDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(s2ClassDir)||~exist(s2ClassDir,'dir')
    s2ClassDir=uigetdir(cfg.OUT_DIR,'Select Stage 2 classification folder');
    if isequal(s2ClassDir,0), error('No folder selected.'); end
end

% Baseline coral mask
bcPath=fullfile(cfg.OUT_DIR,'baseline_coral_output.mat');
bcMaskDir='';
if exist(bcPath,'file')
    try bco=load(bcPath); if exist(bco.outDir,'dir'), bcMaskDir=bco.outDir; end; catch, end
end
if isempty(bcMaskDir)
    base=fullfile(cfg.OUT_DIR,'Baseline_Coral');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); bcMaskDir=fullfile(base,runs(si(end)).name); end
    end
end

% PIF for RGB
pifDir='';
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
if exist(pifBase,'dir')
    runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
    if ~isempty(runs),[~,si_]=sort({runs.name}); pifDir=fullfile(pifBase,runs(si_(end)).name); end
end

outDir=fullfile(cfg.OUT_DIR,'Change_Analysis',['run_' runTimestamp]);
mkdir(outDir);
logW('S2 classification : %s\n',s2ClassDir);
logW('Baseline mask     : %s\n',bcMaskDir);
logW('Output            : %s\n\n',outDir);

%% "" Load baseline coral mask """""""""""""""""""""""""""""""""
bc=load(fullfile(bcMaskDir,'baseline_coral_mask.mat'),...
    'baseline_coral_mask','tier_map','vote_map','reef_mask','H','W','R_ref','px_m2');
baseline_coral_mask = bc.baseline_coral_mask;
reef_mask = bc.reef_mask;
H=bc.H; W=bc.W; R_ref=bc.R_ref; px_m2=bc.px_m2;
px_ha = px_m2 / 10000;
n_coral = nnz(baseline_coral_mask);
logW('Baseline coral mask: %d px  %.1f m2  (%.1f%% reef)\n\n',...
    n_coral,n_coral*px_m2,100*n_coral/max(nnz(reef_mask),1));

%% "" Load S2CLASS maps """"""""""""""""""""""""""""""""""""""""
logW('Loading Stage 2 classification maps...\n');
maps=struct();
for di=1:numel(ALL_DATES)
    d=ALL_DATES{di};
    p=fullfile(s2ClassDir,['S2CLASS_' d '.tif']);
    if exist(p,'file')
        [tmp,~]=readgeoraster(p);
        maps.(['d_' d])=uint8(tmp);
        logW('  Loaded S2CLASS_%s.tif\n',d);
    else
        logW('  MISSING: S2CLASS_%s.tif\n',d);
    end
end
logW('\n');

%% "" Pre-bleaching HEALTHY consensus """""""""""""""""""""""""
% A pixel is "pre-bleaching HEALTHY" if it was HEALTHY on >=3/4
% pre-bleaching dates AND is inside baseline_coral_mask
logW('Computing pre-bleaching HEALTHY consensus...\n');
map_names=fieldnames(maps);
for mi=1:numel(map_names)
    m=maps.(map_names{mi});
    zero_in_coral = baseline_coral_mask & (m==0);
    if any(zero_in_coral(:))
        m(zero_in_coral)=ID_UNK;
        maps.(map_names{mi})=m;
    end
end
prebl_healthy_votes = zeros(H,W,'uint8');
for di=1:numel(PRE_BLEACH_DATES)
    d=PRE_BLEACH_DATES{di};
    if isfield(maps,['d_' d])
        prebl_healthy_votes = prebl_healthy_votes + uint8(maps.(['d_' d])==ID_HEALTHY);
    end
end
n_prebl = sum(cellfun(@(d) isfield(maps,['d_' d]), PRE_BLEACH_DATES));
prebl_threshold = max(1, ceil(n_prebl*0.75));   % >=75% of available dates
prebl_healthy = baseline_coral_mask & (prebl_healthy_votes >= prebl_threshold);
logW('Pre-bleaching HEALTHY pixels (>=%d/%d votes): %d  (%.1f%% of coral mask)\n\n',...
    prebl_threshold,n_prebl,nnz(prebl_healthy),100*nnz(prebl_healthy)/max(n_coral,1));

%% "" Bleaching frequency map """"""""""""""""""""""""""""""""""
% Count how many bleaching dates each coral pixel was BLEACHED
logW('Computing bleaching frequency map (primary dates only)...\n');
bl_freq_primary = zeros(H,W,'uint8');
bl_freq_all     = zeros(H,W,'uint8');
for di=1:numel(BLEACH_DATES_ALL)
    d=BLEACH_DATES_ALL{di};
    if isfield(maps,['d_' d])
        bl_map = maps.(['d_' d])==ID_BLEACHED & baseline_coral_mask;
        bl_freq_all = bl_freq_all + uint8(bl_map);
        if ismember(d,BLEACH_PRIMARY)
            bl_freq_primary = bl_freq_primary + uint8(bl_map);
        end
    end
end
bl_freq_primary(~baseline_coral_mask)=0;
bl_freq_all(~baseline_coral_mask)=0;

saveUint8Tif(fullfile(outDir,'BLEACH_FREQ_PRIMARY.tif'),bl_freq_primary,R_ref);
saveUint8Tif(fullfile(outDir,'BLEACH_FREQ_ALL.tif'),bl_freq_all,R_ref);

logW('Bleaching frequency (primary dates, inside coral mask):\n');
n_primary = numel(BLEACH_PRIMARY);
for v=0:n_primary
    npx=nnz(bl_freq_primary==v & baseline_coral_mask);
    logW('  %d/%d primary dates bleached: %6d px  %10.1f m2  (%.1f%%)\n',...
        v,n_primary,npx,npx*px_m2,100*npx/max(n_coral,1));
end
logW('\n');

%% "" Persistent bleaching map """""""""""""""""""""""""""""""""
% Pixels BLEACHED on ALL primary bleaching dates
bl_persistent = bl_freq_primary == n_primary & baseline_coral_mask;
saveUint8Tif(fullfile(outDir,'BLEACH_PERSISTENT.tif'),uint8(bl_persistent),R_ref);
logW('Persistent BLEACHED (all %d primary dates): %d px  %.1f m2  (%.1f%% coral mask)\n\n',...
    n_primary,nnz(bl_persistent),nnz(bl_persistent)*px_m2,...
    100*nnz(bl_persistent)/max(n_coral,1));

%% "" Recovery map at 20240926 """"""""""""""""""""""""""""""""
% Recovery state: was BLEACHED in primary dates, now HEALTHY at 20240926
d_rec='20240926';
if isfield(maps,['d_' d_rec])
    recovered  = bl_persistent & maps.(['d_' d_rec])==ID_HEALTHY;
    still_bl   = bl_persistent & maps.(['d_' d_rec])==ID_BLEACHED;
    uncertain  = bl_persistent & maps.(['d_' d_rec])==ID_UNK;
    saveUint8Tif(fullfile(outDir,'RECOVERY_MAP.tif'),...
        uint8(recovered)*2 + uint8(still_bl)*1 + uint8(uncertain)*3, R_ref);
    logW('Recovery at %s (of persistent bleached %d px):\n',d_rec,nnz(bl_persistent));
    logW('  Recovered  (HEALTHY): %d px  %.1f m2  (%.1f%%)\n',...
        nnz(recovered),nnz(recovered)*px_m2,100*nnz(recovered)/max(nnz(bl_persistent),1));
    logW('  Still BL  (BLEACHED): %d px  %.1f m2  (%.1f%%)\n',...
        nnz(still_bl),nnz(still_bl)*px_m2,100*nnz(still_bl)/max(nnz(bl_persistent),1));
    logW('  Uncertain      (UNK): %d px  %.1f m2  (%.1f%%)\n\n',...
        nnz(uncertain),nnz(uncertain)*px_m2,100*nnz(uncertain)/max(nnz(bl_persistent),1));
end

%% "" Per-date statistics inside coral mask """""""""""""""""""
logW('Per-date statistics (inside coral mask):\n');
logW('%-10s  %10s  %10s  %8s  %8s  %8s\n',...
    'Date','HEALTHY_m2','BLEACH_m2','BL_pct','HY_pct','UNK_pct');

date_stats=struct('date',{},...
    'healthy_px',{},'bleached_px',{},'unk_px',{},...
    'healthy_m2',{},'bleached_m2',{},'unk_m2',{},...
    'healthy_ha',{},'bleached_ha',{},'unk_ha',{},...
    'bleach_pct',{},'healthy_pct',{},'unk_pct',{},...
    'wave_flag',{});

for di=1:numel(ALL_DATES)
    d=ALL_DATES{di};
    if ~isfield(maps,['d_' d]), continue; end
    m=maps.(['d_' d]);
    n_hy = nnz(m==ID_HEALTHY  & baseline_coral_mask);
    n_bl = nnz(m==ID_BLEACHED & baseline_coral_mask);
    n_uk = nnz(m==ID_UNK      & baseline_coral_mask);
    hy_m2=n_hy*px_m2; bl_m2=n_bl*px_m2; uk_m2=n_uk*px_m2;
    hy_ha=n_hy*px_ha; bl_ha=n_bl*px_ha; uk_ha=n_uk*px_ha;
    bl_pct=100*n_bl/max(n_coral,1);
    hy_pct=100*n_hy/max(n_coral,1);
    uk_pct=100*n_uk/max(n_coral,1);
    wf=ismember(d,BLEACH_WAVE);
    logW('  %-10s  %10.1f  %10.1f  %7.1f%%  %7.1f%%  %7.1f%%%s\n',...
        d,hy_m2,bl_m2,bl_pct,hy_pct,uk_pct,iif(wf,'  ',''));
    date_stats(end+1)=struct(...
        'date',d,...
        'healthy_px',n_hy,'bleached_px',n_bl,'unk_px',n_uk,...
        'healthy_m2',hy_m2,'bleached_m2',bl_m2,'unk_m2',uk_m2,...
        'healthy_ha',hy_ha,'bleached_ha',bl_ha,'unk_ha',uk_ha,...
        'bleach_pct',bl_pct,'healthy_pct',hy_pct,'unk_pct',uk_pct,...
        'wave_flag',wf);
end
logW('\n');

%% "" Figure 1: Bleaching frequency spatial map """"""""""""""""
logW('Generating bleaching frequency map...\n');
figF=figure('Visible','off','Position',[0 0 1000 480]);

% Load reference RGB
rgb_ref=loadRGB(pifDir,cfg.REF_DATE,cfg.BANDS,H,W);

subplot(1,3,1);
imshow(rgb_ref); title(sprintf('RGB - %s',cfg.REF_DATE),'FontSize',9);

subplot(1,3,2);
freq_disp=single(bl_freq_primary); freq_disp(~baseline_coral_mask)=NaN;
imagesc(freq_disp,[0 n_primary]); axis image off;
colormap(gca,flipud(hot));
cb=colorbar; cb.Ticks=0:n_primary;
cb.TickLabels=arrayfun(@(v)sprintf('%d/%d',v,n_primary),0:n_primary,'UniformOutput',false);
cb.Label.String='Primary dates bleached';
title(sprintf('Bleaching frequency\n(primary dates: %s)',strjoin(BLEACH_PRIMARY,' & ')),'FontSize',8);

subplot(1,3,3);
% Persistent bleaching overlay on RGB
overlay=rgb_ref;
for ch=1:3
    sl=overlay(:,:,ch);
    sl(bl_persistent)=CLASS_COLORS(2,ch);   % BLEACHED colour
    if exist('recovered','var') && any(recovered(:))
        sl(recovered)=CLASS_COLORS(1,ch);   % recovered - green
    end
    overlay(:,:,ch)=sl;
end
imshow(overlay);
title(sprintf('Persistent bleached (red)\nRecovered at %s (green)',d_rec),'FontSize',8);

sgtitle('Bleaching Frequency & Persistence - Koh Man Nai','FontSize',11,'FontWeight','bold');
exportgraphics(figF,fullfile(outDir,'bleach_frequency_map.png'),'Resolution',150);
close(figF);

%% "" Figure 2: Recovery trajectory """""""""""""""""""""""""""
logW('Generating recovery trajectory...\n');

all_dates_sorted = ALL_DATES;
date_nums=zeros(1,numel(all_dates_sorted));
for di=1:numel(all_dates_sorted)
    try
        date_nums(di)=datenum(datetime(all_dates_sorted{di},'InputFormat','yyyyMMdd'));
    catch
        date_nums(di)=di;
    end
end

bl_ts=zeros(1,numel(all_dates_sorted));
hy_ts=zeros(1,numel(all_dates_sorted));
has_data=false(1,numel(all_dates_sorted));
for di=1:numel(all_dates_sorted)
    d=all_dates_sorted{di};
    idx=find(strcmp({date_stats.date},d),1);
    if ~isempty(idx)
        bl_ts(di)=date_stats(idx).bleach_pct;
        hy_ts(di)=date_stats(idx).healthy_pct;
        has_data(di)=true;
    end
end

figR=figure('Visible','off','Position',[0 0 960 460]);
hold on;

% Shade bleaching period
bl_nums=date_nums(ismember(all_dates_sorted,BLEACH_DATES_ALL));
if ~isempty(bl_nums)
    fill([min(bl_nums)-20 max(bl_nums)+20 max(bl_nums)+20 min(bl_nums)-20],...
        [0 0 105 105],[1 0.93 0.88],'EdgeColor','none','FaceAlpha',0.5,...
        'DisplayName','Bleaching period');
end

% Plot lines for all valid dates
for di=1:numel(all_dates_sorted)
    if ~has_data(di), continue; end
    d=all_dates_sorted{di};
    is_wf=ismember(d,BLEACH_WAVE);
    mk=iif(is_wf,'v','o');
    plot(date_nums(di),bl_ts(di),mk,'Color',[0.8 0.1 0.1],...
        'MarkerSize',iif(is_wf,9,8),'MarkerFaceColor',...
        iif(is_wf,[1 0.7 0.6],[0.8 0.1 0.1]),'LineWidth',1.5);
    plot(date_nums(di),hy_ts(di),'o','Color',[0 0.55 0],...
        'MarkerSize',8,'MarkerFaceColor',[0 0.55 0]);
    text(date_nums(di),bl_ts(di)+3,sprintf('%.1f%%',bl_ts(di)),...
        'HorizontalAlignment','center','FontSize',7,...
        'Color',iif(is_wf,[0.6 0.2 0],'k'));
end

% Connect primary dates with solid line, wave dates with dashed
valid_primary_idx = find(has_data & ~ismember(all_dates_sorted,BLEACH_WAVE));
valid_wave_idx    = find(has_data &  ismember(all_dates_sorted,BLEACH_WAVE));
if ~isempty(valid_primary_idx)
    plot(date_nums(valid_primary_idx),bl_ts(valid_primary_idx),'r-',...
        'LineWidth',2,'DisplayName','BLEACHED % (primary)');
    plot(date_nums(valid_primary_idx),hy_ts(valid_primary_idx),'g-',...
        'LineWidth',2,'DisplayName','HEALTHY %');
end
for di=valid_wave_idx
    plot(date_nums([di-1 di]),bl_ts([di-1 di]),'r--','LineWidth',1.2,...
        'HandleVisibility','off');
end

% Wave flag and recovery annotations
for d_={'20240613','20240620'}
    d=d_{1};
    di=find(strcmp(all_dates_sorted,d));
    if ~isempty(di) && has_data(di)
        text(date_nums(di),bl_ts(di)-6,' WAVE','HorizontalAlignment','center',...
            'FontSize',7,'Color',[0.6 0.2 0]);
    end
end

plot(nan,nan,'v','Color',[0.8 0.3 0.1],'MarkerSize',8,...
    'MarkerFaceColor',[1 0.7 0.6],'DisplayName','Wave-flagged (upper bound)');

datetick('x','mmm-yy','keepticks');
xlabel('Date','FontSize',10); ylabel('% of coral habitat mask','FontSize',10);
ylim([0 108]);
title(sprintf('Bleaching & Recovery Trajectory - Koh Man Nai\nBaseline coral mask: %d px  %.0f m',...
    n_coral,n_coral*px_m2),'FontSize',11);
legend('Location','best','FontSize',8); grid on; set(gca,'GridAlpha',0.3);
exportgraphics(figR,fullfile(outDir,'recovery_trajectory.png'),'Resolution',150);
close(figR);

%% "" Figure 3: Side-by-side bleaching date maps """""""""""""""
logW('Generating spatial comparison of bleaching dates...\n');
n_bl_dates=numel(BLEACH_DATES_ALL);
figS=figure('Visible','off','Position',[0 0 300*n_bl_dates 560]);
cmap_full=[0 0 0; CLASS_COLORS];

for di=1:n_bl_dates
    d=BLEACH_DATES_ALL{di};
    is_wf=ismember(d,BLEACH_WAVE);
    subplot(2,n_bl_dates,di);
    rgb_d=loadRGB(pifDir,d,cfg.BANDS,H,W);
    imshow(rgb_d);
    title(sprintf('%s%s',d,iif(is_wf,' ','')),'FontSize',8,'Interpreter','none');

    subplot(2,n_bl_dates,n_bl_dates+di);
    if isfield(maps,['d_' d])
        cm_shift=maps.(['d_' d])+1; cm_shift(~reef_mask)=1;
        blend=0.70.*ind2rgb(cm_shift,cmap_full) + 0.30.*rgb_d;
        blend=max(0,min(1,blend));
        imshow(blend);
        n_bl_=nnz(maps.(['d_' d])==ID_BLEACHED & baseline_coral_mask);
        title(sprintf('BL=%.1f%%',100*n_bl_/max(n_coral,1)),'FontSize',8);
    else
        text(0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        axis off;
    end
end

% Common legend
axL=axes(figS,'Visible','off','Position',[0 0 1 0.05]);
hold(axL,'on');
leg_names={'HEALTHY','BLEACHED','SND','RBL','UNK'};
for ci=1:size(CLASS_COLORS,1)
    patch(axL,[nan nan],[nan nan],CLASS_COLORS(ci,:),...
        'EdgeColor','none','DisplayName',leg_names{ci});
end
legend(axL,'Location','north','Orientation','horizontal','FontSize',8,'Box','off');
sgtitle('Bleaching Dates - Stage 2 Classification','FontSize',11,'FontWeight','bold');
exportgraphics(figS,fullfile(outDir,'spatial_comparison.png'),'Resolution',130);
close(figS);

%% "" Figure 4: Transition heatmap """""""""""""""""""""""""""""
% Pre-bleaching HEALTHY -> each bleaching date class
logW('Generating transition heatmap...\n');
n_trans_cls=5;
trans_names={'HEALTHY','BLEACHED','SND','RBL','UNK'};
trans_ids=[1,2,3,4,5];

figTR=figure('Visible','off','Position',[0 0 280*numel(BLEACH_DATES_ALL) 420]);
for di=1:numel(BLEACH_DATES_ALL)
    d=BLEACH_DATES_ALL{di};
    subplot(1,numel(BLEACH_DATES_ALL),di);
    if ~isfield(maps,['d_' d]), axis off; continue; end

    % Transition: prebl_healthy -> this date (inside coral mask)
    T=zeros(1,n_trans_cls);
    for ci=1:n_trans_cls
        T(ci)=nnz(prebl_healthy & maps.(['d_' d])==trans_ids(ci));
    end
    T_pct=100*T/max(sum(T),1);
    b_=bar(T_pct,'FaceColor','flat');
    b_.CData=CLASS_COLORS(1:n_trans_cls,:);
    set(gca,'XTickLabel',trans_names,'FontSize',7);
    xtickangle(30);
    ylabel(iif(di==1,'% of pre-bleach HEALTHY',''),'FontSize',8);
    title(sprintf('%s%s',d,iif(ismember(d,BLEACH_WAVE),' ','')),'FontSize',8,'Interpreter','none');
    ylim([0 105]); grid on;
    for ci=1:n_trans_cls
        if T_pct(ci)>2
            text(ci,T_pct(ci)+2,sprintf('%.1f%%',T_pct(ci)),...
                'HorizontalAlignment','center','FontSize',6);
        end
    end
end
sgtitle(sprintf('Transition from Pre-bleaching HEALTHY\n(%d/%d vote consensus) -> Bleaching Dates',...
    prebl_threshold,n_prebl),'FontSize',10,'FontWeight','bold');
exportgraphics(figTR,fullfile(outDir,'transition_heatmap.png'),'Resolution',150);
close(figTR);

%% "" Change report """"""""""""""""""""""""""""""""""""""""""""
logW('Writing change report...\n');
fid=fopen(fullfile(outDir,'change_report.txt'),'w');
fprintf(fid,'CHANGE ANALYSIS REPORT - Koh Man Nai\n');
fprintf(fid,'Timestamp : %s\n\n',runTimestamp);
fprintf(fid,'Baseline coral mask: %d px  %.1f m2  (%.1f%% reef)\n',...
    n_coral,n_coral*px_m2,100*n_coral/max(nnz(reef_mask),1));
fprintf(fid,'Pre-bleaching HEALTHY consensus (>=%d/%d dates): %d px  %.1f m2\n\n',...
    prebl_threshold,n_prebl,nnz(prebl_healthy),nnz(prebl_healthy)*px_m2);

fprintf(fid,'Wave-flagged dates (upper bound): %s\n',strjoin(BLEACH_WAVE,', '));
fprintf(fid,'Primary bleaching dates          : %s\n\n',strjoin(BLEACH_PRIMARY,', '));

fprintf(fid,'Per-date statistics (inside coral mask):\n');
fprintf(fid,'%-10s  %12s  %12s  %8s  %8s  %8s  %8s\n',...
    'Date','HEALTHY_m2','BLEACH_m2','BL_%','HY_%','UNK_%','WaveFlag');
for k_=1:numel(date_stats)
    ds=date_stats(k_);
    fprintf(fid,'%-10s  %12.1f  %12.1f  %7.1f%%  %7.1f%%  %7.1f%%  %8s\n',...
        ds.date,ds.healthy_m2,ds.bleached_m2,...
        ds.bleach_pct,ds.healthy_pct,...
        ds.unk_pct,...
        iif(ds.wave_flag,' YES','no'));
end

fprintf(fid,'\nBleaching frequency (primary dates, coral mask):\n');
for v=0:n_primary
    npx=nnz(bl_freq_primary==v & baseline_coral_mask);
    fprintf(fid,'  %d/%d dates bleached: %6d px  %10.1f m2  (%.1f%%)\n',...
        v,n_primary,npx,npx*px_m2,100*npx/max(n_coral,1));
end

fprintf(fid,'\nPersistent bleaching (all %d primary dates): %d px  %.1f m2  (%.1f%% coral mask)\n',...
    n_primary,nnz(bl_persistent),nnz(bl_persistent)*px_m2,...
    100*nnz(bl_persistent)/max(n_coral,1));

if isfield(maps,['d_' d_rec])
    fprintf(fid,'\nRecovery at %s (of %d persistent bleached px):\n',d_rec,nnz(bl_persistent));
    fprintf(fid,'  Recovered  (HEALTHY): %d px  %.1f m2  (%.1f%%)\n',...
        nnz(recovered),nnz(recovered)*px_m2,...
        100*nnz(recovered)/max(nnz(bl_persistent),1));
    fprintf(fid,'  Still BL  (BLEACHED): %d px  %.1f m2  (%.1f%%)\n',...
        nnz(still_bl),nnz(still_bl)*px_m2,...
        100*nnz(still_bl)/max(nnz(bl_persistent),1));
    fprintf(fid,'  Uncertain      (UNK): %d px  %.1f m2  (%.1f%%)\n',...
        nnz(uncertain),nnz(uncertain)*px_m2,...
        100*nnz(uncertain)/max(nnz(bl_persistent),1));
end
fclose(fid);

table6Path=fullfile(outDir,'table6_bleaching_summary.csv');
fid=fopen(table6Path,'w');
fprintf(fid,['Date,WaveFlag,Healthy_px,Bleached_px,Unk_px,' ...
    'Healthy_m2,Bleached_m2,Unk_m2,Healthy_ha,Bleached_ha,Unk_ha,' ...
    'Healthy_pct_mask,Bleached_pct_mask,Unk_pct_mask,' ...
    'CoralMask_px,CoralMask_m2,CoralMask_ha\n']);
for k_=1:numel(date_stats)
    ds=date_stats(k_);
    fprintf(fid,'%s,%s,%d,%d,%d,%.1f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.1f,%.4f\n',...
        ds.date,tf2csv(ds.wave_flag),...
        ds.healthy_px,ds.bleached_px,ds.unk_px,...
        ds.healthy_m2,ds.bleached_m2,ds.unk_m2,...
        ds.healthy_ha,ds.bleached_ha,ds.unk_ha,...
        ds.healthy_pct,ds.bleach_pct,ds.unk_pct,...
        n_coral,n_coral*px_m2,n_coral*px_ha);
end
fclose(fid);

%% "" Phase handoff """"""""""""""""""""""""""""""""""""""""""""
save(fullfile(cfg.OUT_DIR,'change_output.mat'),...
    'outDir','date_stats','bl_freq_primary','bl_freq_all',...
    'bl_persistent','baseline_coral_mask','n_coral','px_m2','px_ha',...
    'BLEACH_PRIMARY','BLEACH_WAVE','PRE_BLEACH_DATES','RECOVERY_DATES');

logW('\n========== CHANGE ANALYSIS SUMMARY ==========\n');
logW('Coral mask        : %d px  %.0f m2\n',n_coral,n_coral*px_m2);
logW('Pre-bleach HEALTHY: %d px  (%.1f%% of mask)\n',...
    nnz(prebl_healthy),100*nnz(prebl_healthy)/max(n_coral,1));
logW('Persistent BL     : %d px  %.0f m2  (%.1f%% mask)\n',...
    nnz(bl_persistent),nnz(bl_persistent)*px_m2,...
    100*nnz(bl_persistent)/max(n_coral,1));
if isfield(maps,['d_' d_rec])
    logW('Recovered at %s  : %d px  (%.1f%% of persistent)\n',...
        d_rec,nnz(recovered),100*nnz(recovered)/max(nnz(bl_persistent),1));
end
logW('\nKey bleaching numbers:\n');
for k_=1:numel(date_stats)
    ds=date_stats(k_);
    logW('  %-10s  BL=%.1f%%  HY=%.1f%%%s\n',...
        ds.date,ds.bleach_pct,ds.healthy_pct,...
        iif(ds.wave_flag,'  ',''));
end
logW('\nOutput: %s\n',outDir);
logW('Table 6 CSV: %s\n',table6Path);
logW('NEXT: Diag_SpectralProfile.m\n');
logW('==============================================\n');

%% ===================== LOCAL FUNCTIONS =====================

function rgb=loadRGB(pifDir,date,BANDS,H,W)
    rgb=zeros(H,W,3,'single');
    if isempty(pifDir), return; end
    pf=dir(fullfile(pifDir,['*' date '*.tif']));
    if isempty(pf), return; end
    try
        [Ap,~]=readgeoraster(fullfile(pifDir,pf(1).name));
        Ap=single(Ap)/10000;
        r_=Ap(:,:,BANDS.RED); g_=Ap(:,:,BANDS.GREEN1); b_=Ap(:,:,BANDS.BLUE);
        rgb_=cat(3,r_,g_,b_);
        for ch=1:3
            sl=rgb_(:,:,ch); v=sl(isfinite(sl)&sl>0);
            if ~isempty(v)
                sl=max(0,min(1,(sl-min(v))/max(max(v)-min(v),1e-4)));
                rgb_(:,:,ch)=sl;
            end
        end
        rgb=rgb_;
    catch, end
end

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

function logWrite(fid,varargin)
    msg=sprintf(varargin{:}); fprintf('%s',msg);
    if isnumeric(fid) && fid>2
        try, fprintf(fid,'%s',msg); catch, end
    end
end

function v=iif(c,a,b); if c,v=a; else,v=b; end; end
function s=tf2csv(tf); if tf,s='YES'; else,s='NO'; end; end



