%% Build_BaselineCoralMask.m
% ============================================================
% BASELINE CORAL MASK - Consensus from Stage 1 Pre-bleaching Maps
%
% PURPOSE:
%   Build a persistent coral habitat mask from Stage 1 substrate
%   classification maps on pre-bleaching healthy dates. This mask
%   defines WHERE coral habitat exists regardless of condition,
%   and serves as the spatial prior for Stage 2 (Healthy/Bleached).
%
% DESIGN PRINCIPLE:
%   Only pre-bleaching dates contribute to the baseline.
%   Post-bleaching dates are excluded - spectral recovery is incomplete
%   and Stage 1 cannot reliably detect coral on those dates (CORAL F1=0
%   on 20241111 confirmed). Bleached dates excluded because Stage 1
%   misclassifies bleached coral as SND/UNK.
%
% BASELINE DATES (pre-bleaching Stage 1 maps):
%   Confirmed available: 20230404, 20230509, 20230519
%   Optional (if classified): 20230629
%   Excluded: all 2024+ dates
%
% CONSENSUS VOTING:
%   Per-pixel vote: how many baseline dates classified this pixel as CORAL?
%   Confidence tiers:
%     High     - voted CORAL on >= HIGH_THRESH dates (default: all N dates)
%     Moderate - voted CORAL on >= MOD_THRESH dates (default: N-1 dates)
%     Low      - voted CORAL on >= LOW_THRESH dates (default: 1 date)
%     Exclude  - never voted CORAL
%
% CONFIDENCE ALSO WEIGHTED BY RF CONFIDENCE:
%   Mean confidence score at CORAL pixels across baseline dates added
%   as a supplementary layer - useful for Stage 2 pixel weighting.
%
% INPUTS:
%   CLASS_YYYYMMDD.tif   - Stage 1 uint8 class maps (Phase 9 output)
%   CONF_YYYYMMDD.tif    - Stage 1 confidence maps (Phase 9 output)
%
% OUTPUTS -> OUT_DIR/Baseline_Coral/run_YYYYMMDD_HHMMSS/
%   baseline_coral_mask.mat     consensus mask + confidence tier + vote count
%   BASELINE_CORAL_VOTE.tif     uint8 vote count map (0..N dates)
%   BASELINE_CORAL_CONF.tif     float32 mean RF confidence at CORAL pixels
%   BASELINE_CORAL_TIER.tif     uint8 confidence tier (3=High 2=Mod 1=Low 0=Excl)
%   baseline_coral_quicklook.png spatial map of tiers
%   baseline_coral_report.txt   pixel counts per tier, area in m
%   baseline_coral_output.mat   handoff for Stage 2
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if ~isnumeric(PIPELINE_LOG_FID)||isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<3
    PIPELINE_LOG_FID=-1;
end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Build_BaselineCoralMask.m - %s\n\n',runTimestamp);

cfg = loadPipelineConfigPortable();

%% "" Settings """""""""""""""""""""""""""""""""""""""""""""""""
% Pre-bleaching dates to use for consensus - must have CLASS_*.tif available
% All confirmed before bleaching onset (20240514).
% Calm pre-bleaching dates only - wavy dates excluded (wave confound
% Confirmed clean dates: DHW<1 on all four dates (NOAA CRW verified).
% 20230404: DHW~0, 20230509: DHW~0, 20240307: DHW~0, 20240408: DHW=0.5
% 20240317 excluded - noisy acquisition (CORAL=452px vs ~18700px on clean dates)
% 20230519, 20230629 excluded - wavy dates (near-zero CORAL from wave confound)
% 20240408 excluded: REF_DATE, anomalous spectral profile (CORAL=419px,
%   SND=73.1% vs ~45% on comparable dates) despite DHW=0.5.
%   Likely tide/atmospheric effect on REF_DATE acquisition.
% 20240317 excluded: noisy acquisition (CORAL=452px).
% Three confirmed clean dates: 20230404, 20230509, 20240307
%   all show consistent CORAL=8.1-8.2% and DHW~0.
BASELINE_DATES = {'20230404','20230509','20240307'};

% Confidence tier thresholds (votes out of N available dates)
% Expressed as fraction - applied after determining N loaded dates.
% High:     >= HIGH_FRAC of dates voted CORAL
% Moderate: >= MOD_FRAC  of dates voted CORAL
% Low:      >= 1 date voted CORAL
HIGH_FRAC = 1.0;   % all dates agree
MOD_FRAC  = 0.6;   % majority agree (e.g. 3/4 or 2/3)

% Class IDs from Stage 1 model
CORAL_ID = 1;
UNK_ID   = 4;

% Colours for tier map
TIER_COLORS = [0.85 0.85 0.85;   % Tier 0 - Excluded (light grey)
               0.6  0.9  0.6;    % Tier 1 - Low (pale green)
               0.1  0.7  0.1;    % Tier 2 - Moderate (green)
               0.0  0.4  0.0];   % Tier 3 - High (dark green)

%% "" Find Phase 9 classification output """"""""""""""""""""""
classDir='';
p9Path=fullfile(cfg.OUT_DIR,'phase9_output.mat');
if exist(p9Path,'file')
    try p9=load(p9Path); if exist(p9.outDir,'dir'), classDir=p9.outDir; end; catch, end
end
if isempty(classDir)
    base=fullfile(cfg.OUT_DIR,'Phase9_Classify');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs),[~,si]=sort({runs.name}); classDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(classDir)||~exist(classDir,'dir')
    classDir=uigetdir(cfg.OUT_DIR,'Select Phase 9 classification folder');
    if isequal(classDir,0), error('No folder selected.'); end
end

outDir=fullfile(cfg.OUT_DIR,'Baseline_Coral',['run_' runTimestamp]);
mkdir(outDir);
logW('Classification input: %s\n',classDir);
logW('Output:              %s\n\n',outDir);

%% "" Load available baseline CLASS and CONF maps """""""""""""
logW('Loading baseline classification maps...\n');
loaded_dates={};
class_stack=[];  % H x W x nDates uint8
conf_stack =[];  % H x W x nDates single
H=0; W=0; R_ref=[];

for di=1:numel(BASELINE_DATES)
    d=BASELINE_DATES{di};
    cpath=fullfile(classDir,['CLASS_' d '.tif']);
    fpath=fullfile(classDir,['CONF_'  d '.tif']);
    if ~exist(cpath,'file')
        logW('  MISSING: CLASS_%s.tif - skipping\n',d); continue;
    end
    [Ac,Rc]=readgeoraster(cpath);
    if H==0
        H=size(Ac,1); W=size(Ac,2); R_ref=Rc;
        class_stack=zeros(H,W,0,'uint8');
        conf_stack =zeros(H,W,0,'single');
    end
    class_stack(:,:,end+1)=uint8(Ac);
    if exist(fpath,'file')
        [Af,~]=readgeoraster(fpath);
        conf_stack(:,:,end+1)=single(Af);
    else
        conf_stack(:,:,end+1)=zeros(H,W,'single');
        logW('  WARNING: CONF_%s.tif not found - using zero confidence\n',d);
    end
    loaded_dates{end+1}=d;
    logW('  Loaded CLASS_%s.tif\n',d);
end

nDates=numel(loaded_dates);
if nDates==0, error('No baseline CLASS_*.tif found in %s',classDir); end
logW('\n%d baseline dates loaded: %s\n\n',nDates,strjoin(loaded_dates,', '));

%% "" Load reef mask """""""""""""""""""""""""""""""""""""""""""
reef_mask=[];
p5Path=fullfile(cfg.OUT_DIR,'phase5_output.mat');
p6Path=fullfile(cfg.OUT_DIR,'phase6_output.mat');
if exist(p5Path,'file')
    try p5b=load(p5Path); if isfield(p5b,'reef_mask'), reef_mask=p5b.reef_mask; end; catch, end
end
if isempty(reef_mask) && exist(p6Path,'file')
    try p6b=load(p6Path); if isfield(p6b,'reef_mask'), reef_mask=p6b.reef_mask; end; catch, end
end
if isempty(reef_mask)
    logW('Re-rasterising GeoFeature/Reef...\n');
    reef_mask=rasterizeGeoJSONClass(cfg.REEF_MASK_PATH,H,W,R_ref,'Reef');
end
logW('Reef px: %d\n\n',nnz(reef_mask));

%% "" Compute vote count and mean confidence """""""""""""""""""
logW('Computing consensus vote map...\n');

% vote_map: per-pixel count of dates that classified as CORAL
vote_map=zeros(H,W,'uint8');
conf_sum=zeros(H,W,'single');
conf_n  =zeros(H,W,'uint8');

for di=1:nDates
    coral_px = class_stack(:,:,di)==CORAL_ID;
    vote_map = vote_map + uint8(coral_px);
    % Extract this date's confidence slice first, then index with 2D mask
    conf_slice = conf_stack(:,:,di);
    conf_sum(coral_px) = conf_sum(coral_px) + conf_slice(coral_px);
    conf_n(coral_px)   = conf_n(coral_px) + uint8(1);
end

% Mean RF confidence at CORAL votes
conf_mean=zeros(H,W,'single');
has_vote=conf_n>0;
conf_mean(has_vote)=conf_sum(has_vote)./single(conf_n(has_vote));

% Apply reef mask
vote_map(~reef_mask)=0;
conf_mean(~reef_mask)=0;

logW('Vote distribution (reef pixels):\n');
for v=0:nDates
    n=nnz(vote_map==v & reef_mask);
    logW('  %d/%d votes: %6d px  (%.1f%%)\n',v,nDates,n,100*n/max(nnz(reef_mask),1));
end
logW('\n');

%% "" Assign confidence tiers """"""""""""""""""""""""""""""""""
% Tier 3 - High:     voted CORAL on ALL dates
% Tier 2 - Moderate: voted CORAL on >= MOD_FRAC dates
% Tier 1 - Low:      voted CORAL on >= 1 date
% Tier 0 - Excluded: never voted CORAL
high_thresh = nDates;
mod_thresh  = max(1, ceil(nDates * MOD_FRAC));

tier_map=zeros(H,W,'uint8');
tier_map(reef_mask & vote_map>=1)            = 1;   % Low
tier_map(reef_mask & vote_map>=mod_thresh)   = 2;   % Moderate
tier_map(reef_mask & vote_map>=high_thresh)  = 3;   % High

logW('Confidence tiers  (High threshold=%d/%d  Mod threshold=%d/%d):\n',...
    high_thresh,nDates,mod_thresh,nDates);
tier_names={'Excluded','Low','Moderate','High'};
tier_px=zeros(1,4); tier_m2=zeros(1,4);
px_m2=abs(R_ref.CellExtentInWorldX*R_ref.CellExtentInWorldY);
for t=0:3
    n=nnz(tier_map==t & reef_mask);
    tier_px(t+1)=n; tier_m2(t+1)=n*px_m2;
    logW('  Tier %d %-9s: %6d px  %10.1f m2  (%.1f%%)\n',...
        t,tier_names{t+1},n,n*px_m2,100*n/max(nnz(reef_mask),1));
end
logW('\n');

% Baseline coral mask - High + Moderate confidence (min Stage 2 input)
baseline_coral_mask = tier_map >= 2;
logW('Baseline coral mask (Tier>=2): %d px  %.1f m2\n',...
    nnz(baseline_coral_mask), nnz(baseline_coral_mask)*px_m2);
logW('  As %% of reef: %.1f%%\n\n',100*nnz(baseline_coral_mask)/max(nnz(reef_mask),1));

%% "" Save GeoTIFFs """"""""""""""""""""""""""""""""""""""""""""
logW('Saving GeoTIFFs...\n');
saveUint8Tif(fullfile(outDir,'BASELINE_CORAL_VOTE.tif'), vote_map, R_ref);
saveUint8Tif(fullfile(outDir,'BASELINE_CORAL_TIER.tif'), tier_map, R_ref);
saveMultiBandTif(fullfile(outDir,'BASELINE_CORAL_CONF.tif'), conf_mean, R_ref);
logW('Saved VOTE, TIER, CONF tifs.\n\n');

%% "" Quicklook """"""""""""""""""""""""""""""""""""""""""""""""
logW('Generating quicklook...\n');

% Build tier RGB image
tier_rgb=zeros(H,W,3,'single');
for t=0:3
    for ch=1:3
        sl=tier_rgb(:,:,ch);
        sl(tier_map==t & reef_mask)=TIER_COLORS(t+1,ch);
        tier_rgb(:,:,ch)=sl;
    end
end

% Load PIF RGB for reference
rgb_ref=zeros(H,W,3,'single');
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
if exist(pifBase,'dir')
    runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
    if ~isempty(runs)
        [~,si_]=sort({runs.name});
        pif_f=fullfile(pifBase,runs(si_(end)).name,[cfg.REF_DATE '_PIF.tif']);
        if exist(pif_f,'file')
            try
                [Ap,~]=readgeoraster(pif_f); Ap=single(Ap)/10000;
                r_=Ap(:,:,cfg.BANDS.RED);
                g_=Ap(:,:,cfg.BANDS.GREEN1);
                b_=Ap(:,:,cfg.BANDS.BLUE);
                rgb_ref=cat(3,r_,g_,b_);
                for ch=1:3
                    sl=rgb_ref(:,:,ch); v=sl(isfinite(sl)&sl>0);
                    if ~isempty(v)
                        lo=quantile(v,0.02); hi=quantile(v,0.98);
                        sl=max(0,min(1,(sl-lo)/max(hi-lo,1e-4)));
                        rgb_ref(:,:,ch)=sl;
                    end
                end
            catch, end
        end
    end
end

figBC=figure('Visible','off','Position',[0 0 1100 480]);

% Panel 1: RGB reference
subplot(1,3,1);
imshow(rgb_ref); title(sprintf('RGB - %s',cfg.REF_DATE),'FontSize',9);

% Panel 2: Vote count map
subplot(1,3,2);
vote_disp=single(vote_map); vote_disp(~reef_mask)=NaN;
imagesc(vote_disp,[0 nDates]); axis image off;
colormap(gca,flipud(hot)); cb=colorbar;
cb.Label.String=sprintf('CORAL votes (max %d)',nDates);
cb.Ticks=0:nDates; cb.TickLabels=arrayfun(@num2str,0:nDates,'UniformOutput',false);
title(sprintf('Consensus vote  N=%d dates',nDates),'FontSize',9);

% Panel 3: Tier map with blend
subplot(1,3,3);
blend=0.65.*tier_rgb + 0.35.*rgb_ref;
blend=max(0,min(1,blend));
imshow(blend); hold on;
title('Coral habitat tiers','FontSize',9);
tier_vote_thresh=[1, mod_thresh, high_thresh];   % votes needed per tier 1/2/3
for t=1:3
    patch(nan,nan,TIER_COLORS(t+1,:),'EdgeColor','none',...
        'DisplayName',sprintf('Tier %d - %s (>=%d/%d votes)',...
        t,tier_names{t+1},tier_vote_thresh(t),nDates));
end
legend('Location','southoutside','Orientation','horizontal','FontSize',7,'Box','off');

sgtitle(sprintf('Baseline Coral Habitat Mask - %s',strjoin(loaded_dates,', ')),...
    'FontSize',11,'FontWeight','bold');
exportgraphics(figBC,fullfile(outDir,'baseline_coral_quicklook.png'),'Resolution',150);
close(figBC);
logW('Quicklook saved.\n\n');

%% "" Report """""""""""""""""""""""""""""""""""""""""""""""""""
fid=fopen(fullfile(outDir,'baseline_coral_report.txt'),'w');
fprintf(fid,'BASELINE CORAL HABITAT MASK REPORT\n');
fprintf(fid,'Timestamp : %s\n\n',runTimestamp);
fprintf(fid,'Baseline dates (%d): %s\n',nDates,strjoin(loaded_dates,', '));
fprintf(fid,'Classification source: %s\n\n',classDir);
fprintf(fid,'Confidence thresholds:\n');
fprintf(fid,'  High     : %d/%d dates voted CORAL\n',high_thresh,nDates);
fprintf(fid,'  Moderate : >= %d/%d dates voted CORAL (frac=%.1f)\n',mod_thresh,nDates,MOD_FRAC);
fprintf(fid,'  Low      : >= 1 date voted CORAL\n\n');
fprintf(fid,'Vote distribution (reef = %d px  %.1f m2):\n',nnz(reef_mask),nnz(reef_mask)*px_m2);
for v=0:nDates
    n=nnz(vote_map==v & reef_mask);
    fprintf(fid,'  %d/%d votes: %6d px  %12.1f m2  (%.1f%%)\n',...
        v,nDates,n,n*px_m2,100*n/max(nnz(reef_mask),1));
end
fprintf(fid,'\nConfidence tier areas:\n');
fprintf(fid,'  %-9s  %8s  %14s  %7s\n','Tier','Pixels','Area_m2','Pct');
for t=0:3
    fprintf(fid,'  Tier %d %-9s  %8d  %14.1f  %6.1f%%\n',...
        t,tier_names{t+1},tier_px(t+1),tier_m2(t+1),...
        100*tier_px(t+1)/max(nnz(reef_mask),1));
end
fprintf(fid,'\nBaseline coral mask (Tier >= 2 = High + Moderate):\n');
fprintf(fid,'  Pixels : %d\n',nnz(baseline_coral_mask));
fprintf(fid,'  Area   : %.1f m2\n',nnz(baseline_coral_mask)*px_m2);
fprintf(fid,'  %% reef : %.1f%%\n',100*nnz(baseline_coral_mask)/max(nnz(reef_mask),1));
fclose(fid);

%% "" Save handoff """""""""""""""""""""""""""""""""""""""""""""
save(fullfile(outDir,'baseline_coral_mask.mat'),...
    'baseline_coral_mask','tier_map','vote_map','conf_mean',...
    'loaded_dates','nDates','high_thresh','mod_thresh',...
    'H','W','R_ref','reef_mask','px_m2','runTimestamp');
save(fullfile(cfg.OUT_DIR,'baseline_coral_output.mat'),...
    'outDir','baseline_coral_mask','tier_map','loaded_dates');

logW('========== BASELINE CORAL SUMMARY ==========\n');
logW('Baseline dates    : %s\n',strjoin(loaded_dates,', '));
logW('Reef px           : %d\n',nnz(reef_mask));
logW('High confidence   : %d px  (%.1f%%)\n',...
    tier_px(4),100*tier_px(4)/max(nnz(reef_mask),1));
logW('Moderate conf     : %d px  (%.1f%%)\n',...
    tier_px(3),100*tier_px(3)/max(nnz(reef_mask),1));
logW('Baseline mask     : %d px  %.1f m2  (%.1f%% reef)\n',...
    nnz(baseline_coral_mask),nnz(baseline_coral_mask)*px_m2,...
    100*nnz(baseline_coral_mask)/max(nnz(reef_mask),1));
logW('Output: %s\n',outDir);
logW('NEXT: Train_Stage2.m (Healthy/Bleached within baseline mask)\n');
logW('=============================================\n');

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


