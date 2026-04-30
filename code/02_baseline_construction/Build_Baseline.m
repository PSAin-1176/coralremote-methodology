%% Build_Baseline.m
% ============================================================
% PHASE 6 - 2023 Median Baseline + Stability Map
%
% Pixel-wise temporal median of the 3 baseline date feature
% stacks (20230404, 20230509, 20240307) inside GeoFeature "Reef" polygon.
% CV stability map per feature identifies temporally unstable
% pixels that should get lower confidence in training.
%
% Masks are loaded from phase5_output.mat (no re-rasterisation).
% Sand stability check runs on MSandMask pixels.
% ACA benthic is NOT used.
%
% OUTPUTS -> OUT_DIR/Phase6_Baseline/run_YYYYMMDD_HHMMSS/
%   BASELINE_FEATURES.tif    9-band median feature stack
%   STABILITY_CV.tif         9-band CV map (lower = more stable)
%   BASELINE_quicklook.png
%   baseline_report.txt
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<0, PIPELINE_LOG_FID=-1; end
logW=@(varargin) logWrite(PIPELINE_LOG_FID,varargin{:});

runTime=datetime('now','TimeZone','local');
runTimestamp=char(runTime,'yyyyMMdd_HHmmss');
logW('Build_Baseline.m - Phase 6 - %s\n\n',runTimestamp);

global PIPELINE_CONFIG_PATH;
cfg=load(PIPELINE_CONFIG_PATH);

%% "" Find Phase 5 output """"""""""""""""""""""""""""""""""""""
featDir='';
p5Path=fullfile(cfg.OUT_DIR,'phase5_output.mat');
if exist(p5Path,'file')
    try p5=load(p5Path); if exist(p5.outDir,'dir'), featDir=p5.outDir; end; catch, end
end
if isempty(featDir)
    base=fullfile(cfg.OUT_DIR,'Phase5_Features');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs), [~,si]=sort({runs.name}); featDir=fullfile(base,runs(si(end)).name); end
    end
end
if isempty(featDir)||~exist(featDir,'dir')
    featDir=uigetdir(cfg.OUT_DIR,'Select Phase 5 feature folder');
    if isequal(featDir,0), error('No folder selected.'); end
end

outDir=fullfile(cfg.OUT_DIR,'Phase6_Baseline',['run_' runTimestamp]);
mkdir(outDir);
logW('Feature input: %s\n',featDir);
logW('Output:        %s\n\n',outDir);

%% "" Load masks from Phase 5 handoff """""""""""""""""""""""""
reef_mask=[]; sand_mask=[];
if exist(p5Path,'file')
    try
        p5=load(p5Path);
        if isfield(p5,'reef_mask'), reef_mask=p5.reef_mask; end
        if isfield(p5,'sand_mask'), sand_mask=p5.sand_mask; end
        logW('Masks loaded from phase5_output.mat\n');
    catch ME
        logW('Could not load masks from phase5_output: %s\n',ME.message);
    end
end

%% "" Baseline dates """""""""""""""""""""""""""""""""""""""""""
DATE_CFG=cfg.DATE_CFG;
baseline_dates={};
for k=1:numel(DATE_CFG)
    if DATE_CFG(k).baseline_coral, baseline_dates{end+1}=DATE_CFG(k).date; end
end
logW('Baseline dates: %s\n\n',strjoin(baseline_dates,', '));

%% "" Find feature stacks """"""""""""""""""""""""""""""""""""""
feat_files=struct('date',{},'path',{});
for k=1:numel(baseline_dates)
    p=fullfile(featDir,['FEATURES_' baseline_dates{k} '.tif']);
    if exist(p,'file')
        feat_files(end+1)=struct('date',baseline_dates{k},'path',p);
        logW('  Found: FEATURES_%s.tif\n',baseline_dates{k});
    else
        logW('  MISSING: FEATURES_%s.tif\n',baseline_dates{k});
    end
end
if numel(feat_files)<2, error('Need >=2 baseline stacks. Found %d.',numel(feat_files)); end
logW('\n%d baseline stacks ready.\n\n',numel(feat_files));

%% "" Load into 4D stack """""""""""""""""""""""""""""""""""""""
[Afirst,Rfirst]=readgeoraster(feat_files(1).path);
[H,W,N_FEAT]=size(Afirst); nBase=numel(feat_files);

% Fall back: re-rasterise masks if Phase 5 handoff unavailable
if isempty(reef_mask)
    logW('Re-rasterising GeoFeature/Reef...\n');
    reef_mask=rasterizeGeoJSONClass(cfg.REEF_MASK_PATH,H,W,Rfirst,'Reef');
end
if isempty(sand_mask)
    logW('Re-rasterising MSandMask.geojson...\n');
    sand_mask=rasterizeGeoJSON(fullfile(cfg.ACA_DIR,'MSandMask.geojson'),H,W,Rfirst);
end

% Road mask - spectrally stable on-land target for radiometric QC
road_mask=false(H,W);
roadPath='';
if isfield(cfg,'ROAD_MASK'), roadPath=cfg.ROAD_MASK; end
if isempty(roadPath), roadPath=fullfile(cfg.ACA_DIR,'RoadMark.geojson'); end
if exist(roadPath,'file')
    [~,~,roadExt]=fileparts(roadPath);
    if any(strcmpi(roadExt,{'.tif','.tiff'}))
        [roadTmp,~]=readgeoraster(roadPath);
        if ndims(roadTmp) > 2
            roadTmp = roadTmp(:,:,1);
        end
        road_mask = logical(roadTmp);
        if ~isequal(size(road_mask), [H W])
            warning('Road mask raster size mismatch: %s', roadPath);
            road_mask = false(H,W);
        end
        logW('Road mask raster (%s): %d px\n',roadPath,nnz(road_mask));
    else
        road_mask=rasterizeGeoJSON(roadPath,H,W,Rfirst);
        logW('Road mask vector (%s): %d px\n',roadPath,nnz(road_mask));
    end
else
    logW('Road mask: RoadMark.geojson not found - stability check will skip road\n');
end
logW('Reef px: %d  Sand px: %d  Road px: %d\n\n',nnz(reef_mask),nnz(sand_mask),nnz(road_mask));

stack4D=nan(H,W,N_FEAT,nBase,'single');
stack4D(:,:,:,1)=single(Afirst);
for k=2:nBase
    [Ak,~]=readgeoraster(feat_files(k).path);
    stack4D(:,:,:,k)=single(Ak);
    logW('Loaded [%d/%d]: %s\n',k,nBase,feat_files(k).date);
end

%% "" Median baseline + CV """""""""""""""""""""""""""""""""""""
logW('\nComputing pixel-wise median and CV...\n');
baseline_median=nan(H,W,N_FEAT,'single');
baseline_cv    =nan(H,W,N_FEAT,'single');
for f=1:N_FEAT
    sl=squeeze(stack4D(:,:,f,:));         % H x W x nBase
    baseline_median(:,:,f)=median(sl,3,'omitnan');
    mu=mean(sl,3,'omitnan'); sg=std(sl,0,3,'omitnan');
    baseline_cv(:,:,f)=sg./max(abs(mu),1e-6);
end

%% "" Sand + Road stability check """"""""""""""""""""""""""""""
% DII_BR (feature 1) CV should be low at both sand and road pixels.
% Sand: reef substrate stability. Road: on-land radiometric anchor.
% High CV at either suggests baseline date inconsistency.

% Load master sand mask for QC (must be defined before use below)
master_sand_qc = false(H,W);
masterSandPath = '';
if isfield(cfg,'MASTER_SAND_MASK'), masterSandPath = cfg.MASTER_SAND_MASK; end
if ~isempty(masterSandPath) && exist(masterSandPath,'file')
    master_sand_qc = rasterizeGeoJSON(masterSandPath,H,W,Rfirst);
    logW('Sand stability QC using master sand mask: %d px\n', nnz(master_sand_qc));
end
if nnz(master_sand_qc) == 0
    master_sand_qc = sand_mask;
    logW('Sand stability QC fallback to MSandMask: %d px\n', nnz(master_sand_qc));
end

% Collect per-pixel DII_BR time series at sand and road locations
dii_br_4d = squeeze(stack4D(:,:,1,:));   % H x W x nBase

% Sand stability
if nnz(master_sand_qc) > 0
    cv_sand     = baseline_cv(:,:,1);
    cv_sand_vals= cv_sand(master_sand_qc & isfinite(cv_sand));
    sand_mean_ts= squeeze(mean(mean(dii_br_4d .* repmat(master_sand_qc,[1 1 nBase]),1,'omitnan'),2,'omitnan'));
    logW('\nSand DII_BR stability (CV target <0.05):\n');
    logW('  Median CV=%.4f  P75=%.4f  P90=%.4f\n',...
        median(cv_sand_vals),prctile(cv_sand_vals,75),prctile(cv_sand_vals,90));
    if median(cv_sand_vals)>0.05
        logW('  WARNING: sand CV > 0.05 - check 20230629 quality or PIF calibration\n');
    end
else
    cv_sand_vals=[]; sand_mean_ts=nan(nBase,1);
    logW('\nSand stability: no sand pixels available\n');
end

% Road stability
if nnz(road_mask) > 0
    % Road pixels are on-land - outside reef_mask, so baseline_cv is NaN there.
    % Compute CV directly from the 4D stack at road locations.
    road_ts = reshape(dii_br_4d(repmat(road_mask,[1 1 nBase])),[nnz(road_mask),nBase]);
    road_mu  = mean(road_ts,2,'omitnan');
    road_sg  = std(road_ts,0,2,'omitnan');
    cv_road_vals = road_sg ./ max(abs(road_mu),1e-6);
    cv_road_vals = cv_road_vals(isfinite(cv_road_vals));
    road_mean_ts= squeeze(mean(mean(dii_br_4d .* repmat(road_mask,[1 1 nBase]),1,'omitnan'),2,'omitnan'));
    logW('Road DII_BR stability (CV target <0.05):\n');
    logW('  Median CV=%.4f  P75=%.4f  P90=%.4f\n',...
        median(cv_road_vals),prctile(cv_road_vals,75),prctile(cv_road_vals,90));
    if median(cv_road_vals)>0.05
        logW('  WARNING: road CV > 0.05 - check PIF calibration\n');
    end
else
    cv_road_vals=[]; road_mean_ts=nan(nBase,1);
    logW('Road stability: no road pixels available\n');
end

% "" 2-panel stability quicklook """""""""""""""""""""""""""""""
figStab=figure('Visible','off','Position',[0 0 900 380]);
bdate_labels={feat_files.date};
bdate_labels=strrep(bdate_labels,'2023','23');
bdate_labels=strrep(bdate_labels,'2024','24');

subplot(1,2,1);
if ~isempty(cv_sand_vals)
    histogram(cv_sand_vals,30,'FaceColor',[1 0.8 0],'EdgeColor','none'); hold on;
    xline(median(cv_sand_vals),'k-','LineWidth',2,...
        'Label',sprintf('Med=%.3f',median(cv_sand_vals)));
    xline(0.05,'r--','LineWidth',1.5,'Label','target');
    xlabel('DII\_BR CV'); ylabel('Pixel count'); grid on;
    title(sprintf('Sand pixels  (n=%d)',nnz(master_sand_qc)),'FontSize',9);
else
    text(0.5,0.5,'No sand data','HorizontalAlignment','center','Units','normalized');
    axis off;
end

subplot(1,2,2);
if ~isempty(cv_road_vals)
    histogram(cv_road_vals,20,'FaceColor',[1 0.4 0.4],'EdgeColor','none'); hold on;
    xline(median(cv_road_vals),'k-','LineWidth',2,...
        'Label',sprintf('Med=%.3f',median(cv_road_vals)));
    xline(0.05,'r--','LineWidth',1.5,'Label','target');
    xlabel('DII\_BR CV'); ylabel('Pixel count'); grid on;
    title(sprintf('Road pixels  (n=%d)',nnz(road_mask)),'FontSize',9);
else
    text(0.5,0.5,'No road data','HorizontalAlignment','center','Units','normalized');
    axis off;
end

sgtitle('Baseline DII\_BR Stability CV - Sand vs Road','FontSize',11,'FontWeight','bold');
exportgraphics(figStab,fullfile(outDir,'BASELINE_stability_QC.png'),'Resolution',150);
close(figStab);

%% "" Save """""""""""""""""""""""""""""""""""""""""""""""""""""
logW('\nSaving outputs...\n');
saveMultiBandTif(fullfile(outDir,'BASELINE_FEATURES.tif'),baseline_median,Rfirst);
saveMultiBandTif(fullfile(outDir,'STABILITY_CV.tif'),     baseline_cv,    Rfirst);

%% "" Quicklook """"""""""""""""""""""""""""""""""""""""""""""""
figB=figure('Visible','off','Position',[0 0 1500 420]);

% Reference RGB (from Phase 3 PIF)
rgb_img=[];
pifBase=fullfile(cfg.OUT_DIR,'Phase3_PIF');
runs=dir(fullfile(pifBase,'run_*')); runs=runs([runs.isdir]);
if ~isempty(runs)
    [~,si]=sort({runs.name});
    pifRef=fullfile(pifBase,runs(si(end)).name,[cfg.REF_DATE '_PIF.tif']);
    if exist(pifRef,'file')
        try
            [Ap,~]=readgeoraster(pifRef); Ap=single(Ap);
            rgb_img=imadjust(cat(3,Ap(:,:,cfg.BANDS.RED),...
                Ap(:,:,cfg.BANDS.GREEN1),Ap(:,:,cfg.BANDS.BLUE)),[0 0.12],[0 1]);
        catch, end
    end
end

subplot(1,5,1);
if ~isempty(rgb_img), imshow(rgb_img); else, imagesc(zeros(H,W)); axis image off; end
title(sprintf('RGB %s',cfg.REF_DATE),'FontSize',9);

panelNames = cfg.LR.FEATURE_NAMES;
for pi_=1:4
    subplot(1,5,pi_+1);
    if pi_ <= min(3, N_FEAT)
        d = baseline_median(:,:,pi_);
        cmap_ = parula;
        title_ = strrep(panelNames{pi_}, '_', '\_');
    else
        d = baseline_cv(:,:,1);
        cmap_ = hot;
        title_ = sprintf('%s CV', strrep(panelNames{1}, '_', '\_'));
    end
    d(~reef_mask)=NaN;
    dv=d(reef_mask&isfinite(d));
    if ~isempty(dv)
        imagesc(d,[prctile(dv,2) prctile(dv,98)]);
    else
        imagesc(d);
    end
    colormap(gca,cmap_); axis image off; colorbar('FontSize',7);
    title(title_,'FontSize',9);
end

% Overlay sand mask on DII_BR panel
subplot(1,5,2); hold on;
if nnz(sand_mask)>0
    ov=zeros(H,W,3); ov(:,:,3)=1;
    h=imshow(ov); set(h,'AlphaData',double(sand_mask)*0.5);
end

sgtitle(sprintf('2023 Baseline (%s)',strjoin(baseline_dates,', ')),'FontSize',11);
exportgraphics(figB,fullfile(outDir,'BASELINE_quicklook.png'),'Resolution',150);
close(figB);

%% "" Report """""""""""""""""""""""""""""""""""""""""""""""""""
FEAT_NAMES = cfg.LR.FEATURE_NAMES;
fid=fopen(fullfile(outDir,'baseline_report.txt'),'w');
fprintf(fid,'PHASE 6 - BASELINE REPORT\n');
fprintf(fid,'Timestamp: %s\n',runTimestamp);
fprintf(fid,'Baseline dates: %s\n',strjoin(baseline_dates,', '));
fprintf(fid,'Grid: %dx%d  Reef px: %d  Sand px: %d\n\n',H,W,nnz(reef_mask),nnz(sand_mask));
fprintf(fid,'Masks: GeoFeature.geojson "Reef" polygon (reef)  MSandMask.geojson (sand)\n');
fprintf(fid,'ACA benthic: NOT used\n\n');
fprintf(fid,'%-22s  %8s  %8s  %8s\n','Feature','Med CV','P25 CV','P75 CV');
for f=1:N_FEAT
    cv_f=baseline_cv(:,:,f); cv_f=cv_f(reef_mask&isfinite(cv_f));
    fname=iif(f<=numel(FEAT_NAMES),FEAT_NAMES{f},sprintf('F%d',f));
    if isempty(cv_f)
        fprintf(fid,'%-22s  %8s  %8s  %8s\n',fname,'--','--','--');
    else
        fprintf(fid,'%-22s  %8.4f  %8.4f  %8.4f\n',...
            fname,median(cv_f),prctile(cv_f,25),prctile(cv_f,75));
    end
end
fclose(fid);

save(fullfile(cfg.OUT_DIR,'phase6_output.mat'),'outDir','baseline_dates','N_FEAT',...
    'reef_mask','sand_mask');

logW('\n========== PHASE 6 SUMMARY ==========\n');
logW('Baseline dates  : %s\n',strjoin(baseline_dates,', '));
logW('Grid            : %dx%d  Features: %d\n',H,W,N_FEAT);
logW('Reef px         : %d\n',nnz(reef_mask));
logW('Sand px         : %d\n',nnz(sand_mask));
logW('Output: %s\n',outDir);
logW('NEXT: Explore_Spectral.m (Phase 7)\n');
logW('======================================\n');

%% ===================== LOCAL FUNCTIONS =====================

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
                ring=rings{ri};
                if size(ring,1)<3, continue; end
                [xv,yv]=latlon2utm47(ring(:,2),ring(:,1));
                cv=max(1,min(W,round((xv-xmin)/cw+0.5)));
                rv=max(1,min(H,round((ymax-yv)/ch+0.5)));
                try, mask=mask|poly2mask(double(cv),double(rv),H,W); catch, end
            end
        end
    catch ME, warning('rasterizeGeoJSON(%s): %s',path,ME.message); end
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
    a=6378137.0;f=1/298.257223563;k0=0.9996;E0=500000;lon0=99;
    e2=2*f-f^2;e4=e2^2;e6=e2^3;
    lr=deg2rad(lat);wr=deg2rad(lon);dw=wr-deg2rad(lon0);
    N_=a./sqrt(1-e2.*sin(lr).^2);T=tan(lr).^2;
    C=e2/(1-e2).*cos(lr).^2;A_=cos(lr).*dw;
    M=a*((1-e2/4-3*e4/64-5*e6/256).*lr...
        -(3*e2/8+3*e4/32+45*e6/1024).*sin(2*lr)...
        +(15*e4/256+45*e6/1024).*sin(4*lr)...
        -(35*e6/3072).*sin(6*lr));
    x=k0.*N_.*(A_+(1-T+C).*A_.^3/6+(5-18*T+T.^2+72*C-58*e2/(1-e2)).*A_.^5/120)+E0;
    y=k0.*(M+N_.*tan(lr).*(A_.^2/2+(5-T+9*C+4*C.^2).*A_.^4/24+...
        (61-58*T+T.^2+600*C-330*e2/(1-e2)).*A_.^6/720));
end

function saveMultiBandTif(path,A,R)
    if ~isa(A,'single'), A=single(A); end
    try geotiffwrite(path,A,R); return; catch, end
    try geotiffwrite(path,A,R,'CoordRefSysCode',32647); return; catch, end
    try
        [H,W,nb]=size(A); t=Tiff(path,'w');
        t.setTag('Photometric',1);t.setTag('ImageLength',H);
        t.setTag('ImageWidth',W);t.setTag('BitsPerSample',32);
        t.setTag('SampleFormat',3);t.setTag('SamplesPerPixel',nb);
        t.setTag('Compression',5);t.setTag('PlanarConfig',1);
        t.write(A);t.close();
    catch ME, error('saveMultiBandTif: %s',ME.message); end
end

function logWrite(fid,varargin)
    msg=sprintf(varargin{:}); fprintf('%s',msg);
    if fid>0, fprintf(fid,'%s',msg); end
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


