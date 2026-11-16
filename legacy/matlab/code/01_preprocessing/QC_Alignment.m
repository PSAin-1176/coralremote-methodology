%% QC_Alignment.m
% ============================================================
% PHASE 0 - Co-registration quality check on RAW SR images
%
% Runs before any correction to establish baseline offset measurements.
% Results determine whether AROSICS correction is needed and
% which images are most problematic.
%
% OUTPUTS -> OUT_DIR/Phase0_AlignCheck_raw/run_YYYYMMDD_HHMMSS/
%   alignment_report.txt     full text report + recommendation
%   offset_summary.mat       all metrics
%   offset_map.png           scatter + bar chart of offsets
%   coastline_overlay.png    NDWI coastline stability
%   nir_land_profile.png     NIR temporal stability at land pixels
%   xcorr_<date>.png         cross-correlation surface per image
% ============================================================

clear; clc;
runTime      = datetime('now','TimeZone','local');
runTimestamp = char(runTime,'yyyyMMdd_HHmmss');
fprintf('QC_Alignment.m - Phase 0 - %s\n\n', runTimestamp);

%% Load shared config
cfg = loadPipelineConfig();
if isempty(cfg)
    % Standalone fallback - ask user for config file
    [cfgFile, cfgFolder] = uigetfile('pipeline_config.mat','Select pipeline_config.mat');
    if isequal(cfgFile,0)
        % Manual defaults
        cfg.RAW_DIR        = uigetdir(pwd,'Select RAW SR image folder');
        cfg.OUT_DIR        = uigetdir(pwd,'Select output folder');
        cfg.RAW_PATTERN    = '*3B_AnalyticMS_SR_8b*.tif';
        cfg.QC.MAX_OFFSET_PX = 1.0;
        cfg.BANDS.NIR      = 8;
        cfg.BANDS.GREEN1   = 3;
        cfg.BANDS.BLUE     = 2;
        cfg.SR_SCALE       = 10000;
        cfg.REF_DATE       = '20240408';
    else
        cfg = load(fullfile(cfgFolder, cfgFile));
    end
end

%% Settings
inDir          = cfg.RAW_DIR;
pat            = cfg.RAW_PATTERN;
MAX_OFFSET_PX  = cfg.QC.MAX_OFFSET_PX;
bNIR           = cfg.BANDS.NIR;
bGREEN         = cfg.BANDS.GREEN1;
bBLUE          = cfg.BANDS.BLUE;
SR_SCALE       = cfg.SR_SCALE;
REF_DATE       = cfg.REF_DATE;

XCORR_ROI_FRAC    = 0.60;   % central crop fraction for cross-correlation
WINDOW_SIZE       = [256 256];
MAX_SHIFT_SEARCH  = 20;     % max lag to search in xcorr (pixels)
COAST_DILATION    = 3;      % px dilation for coastline comparison
LAND_SAMPLE_N     = 500;    % land pixels for NIR stability check
NDWI_THRESHOLD    = 0.0;

%% Output folder
outDir = fullfile(cfg.OUT_DIR, 'Phase0_AlignCheck_raw', ['run_' runTimestamp]);
mkdir(outDir);
fprintf('Output: %s\n\n', outDir);

%% Find files
ff = dir(fullfile(inDir, pat));
if isempty(ff)
    error('No files matched: %s in %s', pat, inDir);
end
fprintf('Found %d SR files.\n', numel(ff));

% Sort chronologically by filename date
dates = cellfun(@(n) extractDateFromName(n), {ff.name}, 'UniformOutput', false);
[dates_sorted, sortIdx] = sort(dates);
ff = ff(sortIdx);
fprintf('Dates found: %s\n\n', strjoin(dates_sorted, ', '));

% Find reference image index
refIdx = find(strcmp(dates_sorted, REF_DATE));
if isempty(refIdx)
    warning('Reference date %s not found - using first image.', REF_DATE);
    refIdx = 1;
end
fprintf('Reference image: %s (%s)\n\n', ff(refIdx).name, REF_DATE);

%% ===================== PASS 1: METADATA =====================
fprintf('============================================================\n');
fprintf('PASS 1: Georeference Metadata\n');
fprintf('============================================================\n');

metaLog = struct('file',{},'H',{},'W',{},'cellX',{},'cellY',{},...
    'originX',{},'originY',{},'status',{},'notes',{});
refMeta = [];

for i = 1:numel(ff)
    srPath = fullfile(ff(i).folder, ff(i).name);
    fprintf('[%d/%d] %s\n', i, numel(ff), ff(i).name);
    try
        [A, R] = readgeoraster(srPath);
    catch ME
        fprintf('  SKIP: %s\n', ME.message);
        metaLog(end+1) = makeMetaEntry(ff(i).name,0,0,NaN,NaN,NaN,NaN,'SKIP',ME.message);
        continue;
    end

    [H,W,~] = size(A);
    [cellX, cellY, originX, originY] = extractSpatialRef(R);

    if isempty(refMeta)
        refMeta = struct('H',H,'W',W,'cellX',cellX,'cellY',cellY,...
            'originX',originX,'originY',originY);
        fprintf('  REFERENCE: %d%d px | cell=(%.6f,%.6f) | origin=(%.4f,%.4f)\n',...
            H,W,cellX,cellY,originX,originY);
        metaLog(end+1) = makeMetaEntry(ff(i).name,H,W,cellX,cellY,originX,originY,'REF','');
        continue;
    end

    notes = {}; status = 'OK';
    if H~=refMeta.H || W~=refMeta.W
        notes{end+1} = sprintf('SIZE: %dx%d vs %dx%d',H,W,refMeta.H,refMeta.W);
        status = 'WARN';
    end
    dOX = (originX - refMeta.originX);
    dOY = (originY - refMeta.originY);
    dOX_px = dOX / cellX;   dOY_px = dOY / cellY;
    offPx = sqrt(dOX_px^2 + dOY_px^2);
    notes{end+1} = sprintf('origin offset=%.3f px (%.2f m E, %.2f m N)',...
        offPx, dOX*111320*cosd(originY), dOY*111320);
    if offPx > MAX_OFFSET_PX, status = 'WARN'; end

    fprintf('  %s - %s\n', status, strjoin(notes,' | '));
    metaLog(end+1) = makeMetaEntry(ff(i).name,H,W,cellX,cellY,originX,originY,status,strjoin(notes,' | '));
end

%% ===================== PASS 2: CROSS-CORRELATION =====================
fprintf('\n============================================================\n');
fprintf('PASS 2: Image Cross-Correlation Offsets (NIR band)\n');
fprintf('============================================================\n');

% Load reference image
[Aref, Rref] = readgeoraster(fullfile(ff(refIdx).folder, ff(refIdx).name));
[H,W,~] = size(Aref);
cellSizeM = extractCellSizeM(Rref);

% Central ROI
r1=round(H*(1-XCORR_ROI_FRAC)/2)+1; r2=round(H*(1+XCORR_ROI_FRAC)/2);
c1=round(W*(1-XCORR_ROI_FRAC)/2)+1; c2=round(W*(1+XCORR_ROI_FRAC)/2);
refROI = normalizeROI(double(Aref(r1:r2,c1:c2,bNIR)));

xcorrLog = struct('file',{},'offsetRow',{},'offsetCol',{},...
    'offsetPx',{},'offsetM',{},'peakCorr',{},'status',{});

for i = 1:numel(ff)
    srPath = fullfile(ff(i).folder, ff(i).name);
    fprintf('[%d/%d] %s\n', i, numel(ff), ff(i).name);

    if i == refIdx
        fprintf('  REFERENCE\n');
        xcorrLog(end+1) = struct('file',ff(i).name,'offsetRow',0,'offsetCol',0,...
            'offsetPx',0,'offsetM',0,'peakCorr',1,'status','REF');
        continue;
    end

    try
        [Ai,~] = readgeoraster(srPath);
    catch ME
        fprintf('  SKIP: %s\n', ME.message);
        xcorrLog(end+1) = struct('file',ff(i).name,'offsetRow',NaN,'offsetCol',NaN,...
            'offsetPx',NaN,'offsetM',NaN,'peakCorr',NaN,'status','SKIP');
        continue;
    end

    if size(Ai,1)~=H || size(Ai,2)~=W
        fprintf('  SKIP: size mismatch\n');
        xcorrLog(end+1) = struct('file',ff(i).name,'offsetRow',NaN,'offsetCol',NaN,...
            'offsetPx',NaN,'offsetM',NaN,'peakCorr',NaN,'status','SKIP');
        continue;
    end

    tgtROI = normalizeROI(double(Ai(r1:r2,c1:c2,bNIR)));
    xc     = normxcorr2(refROI, tgtROI);

    % Find peak within MAX_SHIFT_SEARCH lag window
    [roiH, roiW] = size(refROI);
    midR = roiH; midC = roiW;   % zero-lag location in normxcorr2 output
    sr1=max(1,midR-MAX_SHIFT_SEARCH); sr2=min(size(xc,1),midR+MAX_SHIFT_SEARCH);
    sc1=max(1,midC-MAX_SHIFT_SEARCH); sc2=min(size(xc,2),midC+MAX_SHIFT_SEARCH);
    xcSub = xc(sr1:sr2, sc1:sc2);
    [peakVal, peakIdx] = max(xcSub(:));
    [peakR, peakC] = ind2sub(size(xcSub), peakIdx);
    peakR = peakR + sr1 - 1;
    peakC = peakC + sc1 - 1;

    [subR, subC, peakFit] = subpixelPeak(xc, peakR, peakC);
    rowOff = subR - roiH;
    colOff = subC - roiW;
    offPx  = sqrt(rowOff^2 + colOff^2);
    offM   = offPx * cellSizeM;

    if offPx > MAX_OFFSET_PX
        st = 'WARN'; flag = '!';
    else
        st = 'OK'; flag = ' ';
    end
    fprintf('  [%s] row=%+.2fpx  col=%+.2fpx  |%.2f|px  (%.1fm)  corr=%.4f\n',...
        flag, rowOff, colOff, offPx, offM, peakFit);

    xcorrLog(end+1) = struct('file',ff(i).name,'offsetRow',rowOff,'offsetCol',colOff,...
        'offsetPx',offPx,'offsetM',offM,'peakCorr',peakFit,'status',st);

    % Save xcorr quicklook
    fig=figure('Visible','off');
    imagesc(xc(sr1:sr2,sc1:sc2)); colormap(hot); colorbar; axis image;
    hold on;
    plot(colOff+MAX_SHIFT_SEARCH+1, rowOff+MAX_SHIFT_SEARCH+1,'c+','MarkerSize',14,'LineWidth',2);
    title(sprintf('%s | row=%+.2f col=%+.2f px (%.1fm)', ...
        dates_sorted{i}, rowOff, colOff, offM),'Interpreter','none');
    [~,fn]=fileparts(ff(i).name);
    exportgraphics(fig, fullfile(outDir,['xcorr_' fn '.png']),'Resolution',120);
    close(fig);
end

%% ===================== PASS 3: NDWI COASTLINE STABILITY =====================
fprintf('\n============================================================\n');
fprintf('PASS 3: NDWI Coastline Stability\n');
fprintf('============================================================\n');

coastLines = cell(numel(ff),1);
cmap = lines(numel(ff));
for i = 1:numel(ff)
    try
        [Ai,~] = readgeoraster(fullfile(ff(i).folder, ff(i).name));
        G = double(Ai(:,:,bGREEN))/SR_SCALE;
        N = double(Ai(:,:,bNIR))/SR_SCALE;
        ndwi = (G-N)./(G+N+eps);
        land = ndwi < NDWI_THRESHOLD & isfinite(ndwi);
        edge = imdilate(bwperim(land), strel('disk',COAST_DILATION));
        coastLines{i} = edge;
        fprintf('[%d/%d] Edge pixels: %d\n', i, numel(ff), nnz(edge));
    catch
        coastLines{i} = [];
        fprintf('[%d/%d] SKIP\n', i, numel(ff));
    end
end

% Agreement map
coastAgree = zeros(H,W,'uint8');
figC = figure('Visible','off');
refNIR = double(Aref(:,:,bNIR))/SR_SCALE;
imagesc(refNIR); colormap(gray); axis image; hold on;
title('NDWI Coastline Overlay - all dates');
for i=1:numel(coastLines)
    if isempty(coastLines{i}), continue; end
    [ry,cx]=find(coastLines{i});
    scatter(cx,ry,1,repmat(cmap(i,:),numel(ry),1),'.');
    coastAgree = coastAgree + uint8(coastLines{i});
end
nCoastImgs = sum(~cellfun(@isempty, coastLines));
pctAgree = 100*nnz(coastAgree >= round(0.7*nCoastImgs))/max(nnz(coastAgree>0),1);
fprintf('High-agreement coastline: %.1f%%\n', pctAgree);
exportgraphics(figC, fullfile(outDir,'coastline_overlay.png'),'Resolution',150);
close(figC);

%% ===================== PASS 4: NIR LAND STABILITY =====================
fprintf('\n============================================================\n');
fprintf('PASS 4: NIR Temporal Stability at Land Pixels\n');
fprintf('============================================================\n');

G=double(Aref(:,:,bGREEN))/SR_SCALE;
N=double(Aref(:,:,bNIR))/SR_SCALE;
ndwiRef=(G-N)./(G+N+eps);
landRef = ndwiRef < NDWI_THRESHOLD & isfinite(ndwiRef);
landInt = imerode(landRef, strel('disk',5));
fprintf('Interior land pixels: %d\n', nnz(landInt));

nirProfile = nan(min(LAND_SAMPLE_N,nnz(landInt)), numel(ff));
medTemporalStd = NaN;

if nnz(landInt) >= 10
    lx = find(landInt);
    rng(42);
    idx = lx(randperm(numel(lx), min(LAND_SAMPLE_N,numel(lx))));
    [sampR,sampC] = ind2sub([H,W], idx);

    for i=1:numel(ff)
        try
            [Ai,~]=readgeoraster(fullfile(ff(i).folder,ff(i).name));
            Ni=double(Ai(:,:,bNIR))/SR_SCALE;
            for p=1:numel(idx)
                nirProfile(p,i)=Ni(sampR(p),sampC(p));
            end
        catch, end
    end

    pixStd = nanstd(nirProfile,[],2);
    medTemporalStd = median(pixStd(isfinite(pixStd)));
    fprintf('Median temporal NIR std at land: %.4f  [threshold: 0.02]\n', medTemporalStd);
    if medTemporalStd > 0.02
        fprintf('WARNING: High temporal variance - possible misalignment\n');
    else
        fprintf('OK: Land NIR temporally stable\n');
    end

    % Plot
    figN=figure('Visible','off');
    mu=nanmean(nirProfile,1); sig=nanstd(nirProfile,[],1);
    dIdx=1:numel(ff);
    fill([dIdx fliplr(dIdx)],[mu+sig fliplr(mu-sig)],[0.2 0.6 1],'FaceAlpha',0.25,'EdgeColor','none');
    hold on;
    plot(dIdx,mu,'b-o','LineWidth',2,'MarkerFaceColor','b');
    xticks(dIdx);
    xticklabels(dates_sorted); xtickangle(35);
    ylabel('NIR reflectance'); title('NIR Land Stability (interior pixels)');
    grid on;
    exportgraphics(figN, fullfile(outDir,'nir_land_profile.png'),'Resolution',150);
    close(figN);
end

%% ===================== OFFSET SUMMARY PLOT =====================
figO=figure('Visible','off','Position',[0 0 1000 450]);
rows=[xcorrLog.offsetRow]; cols=[xcorrLog.offsetCol]; mags=[xcorrLog.offsetPx];

subplot(1,2,1);
scatter(cols,rows,80,mags,'filled'); colormap(gca,'hot'); cb=colorbar;
cb.Label.String='Offset (px)'; clim([0 max(MAX_OFFSET_PX*2,max(mags(isfinite(mags))))]);
hold on;
theta=linspace(0,2*pi,200);
plot(MAX_OFFSET_PX*cos(theta),MAX_OFFSET_PX*sin(theta),'g--','LineWidth',1.5);
plot(0,0,'g+','MarkerSize',12,'LineWidth',2);
for k=1:numel(xcorrLog)
    if isfinite(rows(k))
        text(cols(k),rows(k),sprintf('  %s',dates_sorted{k}),'FontSize',7);
    end
end
axis equal; grid on;
xlabel('Col offset (px)'); ylabel('Row offset (px)');
title(sprintf('Cross-correlation offsets\n(green=%.1fpx tolerance)',MAX_OFFSET_PX));

subplot(1,2,2);
mags_plot = mags; mags_plot(isnan(mags_plot))=0;
b=bar(mags_plot,'FaceColor',[0.3 0.6 0.9]); hold on;
yline(MAX_OFFSET_PX,'r--','LineWidth',1.5);
% colour bars exceeding threshold red
for k=1:numel(mags)
    if isfinite(mags(k)) && mags(k) > MAX_OFFSET_PX
        b.FaceColor='flat'; b.CData(k,:)=[0.9 0.2 0.2];
    end
end
xticks(1:numel(ff)); xticklabels(dates_sorted); xtickangle(35);
ylabel('Offset magnitude (px)'); title('Per-image offset magnitude'); grid on;

exportgraphics(figO, fullfile(outDir,'offset_map.png'),'Resolution',150);
close(figO);

%% ===================== WRITE REPORT =====================
nWarn = sum(strcmp({xcorrLog.status},'WARN'));
reportPath = fullfile(outDir,'alignment_report.txt');
fid = fopen(reportPath,'w');
fprintf(fid,'================================================\n');
fprintf(fid,'  PHASE 0 - ALIGNMENT CHECK (RAW IMAGES)\n');
fprintf(fid,'================================================\n');
fprintf(fid,'Timestamp      : %s\n', runTimestamp);
fprintf(fid,'Input folder   : %s\n', inDir);
fprintf(fid,'Reference date : %s\n', REF_DATE);
fprintf(fid,'Files checked  : %d\n', numel(ff));
fprintf(fid,'Tolerance      : %.1f px\n', MAX_OFFSET_PX);
fprintf(fid,'Cell size      : %.2f m\n', cellSizeM);
fprintf(fid,'\n');
fprintf(fid,'PASS 1 - Metadata\n');
fprintf(fid,'%-60s  %-6s  %s\n','File','Status','Notes');
for k=1:numel(metaLog)
    fprintf(fid,'%-60s  %-6s  %s\n',metaLog(k).file,metaLog(k).status,metaLog(k).notes);
end
fprintf(fid,'\nPASS 2 - Cross-Correlation Offsets\n');
fprintf(fid,'%-60s  %-6s  %8s  %8s  %8s  %8s\n','File','Status','Row(px)','Col(px)','Mag(px)','Mag(m)');
for k=1:numel(xcorrLog)
    fprintf(fid,'%-60s  %-6s  %8.3f  %8.3f  %8.3f  %8.2f\n',...
        xcorrLog(k).file, xcorrLog(k).status,...
        xcorrLog(k).offsetRow, xcorrLog(k).offsetCol,...
        xcorrLog(k).offsetPx, xcorrLog(k).offsetM);
end
fprintf(fid,'\nPASS 3 - Coastline Agreement: %.1f%%\n', pctAgree);
fprintf(fid,'PASS 4 - Land NIR temporal std: %.4f\n', medTemporalStd);
fprintf(fid,'\n================================================\n');
fprintf(fid,'  RECOMMENDATION\n');
fprintf(fid,'================================================\n');
if nWarn == 0
    fprintf(fid,'All images within %.1f px tolerance.\n', MAX_OFFSET_PX);
    fprintf(fid,'PROCEED to Phase 1 (AROSICS) if you want sub-pixel\n');
    fprintf(fid,'refinement, otherwise skip to Phase 3 (PIF).\n');
elseif nWarn <= 2
    fprintf(fid,'%d image(s) exceed tolerance - AROSICS recommended.\n', nWarn);
    fprintf(fid,'PROCEED to Phase 1 (AROSICS co-registration).\n');
else
    fprintf(fid,'%d images exceed tolerance - AROSICS required.\n', nWarn);
    fprintf(fid,'Do NOT proceed to PIF until offsets are corrected.\n');
end
fprintf(fid,'================================================\n');
fclose(fid);

% Save mat
save(fullfile(outDir,'alignment_summary_raw.mat'),...
    'metaLog','xcorrLog','coastAgree','nirProfile','medTemporalStd',...
    'nWarn','pctAgree','MAX_OFFSET_PX','cellSizeM','dates_sorted');

fprintf('\n========== PHASE 0 SUMMARY ==========\n');
fprintf('Images flagged: %d / %d\n', nWarn, numel(ff));
fprintf('Coastline agreement: %.1f%%\n', pctAgree);
fprintf('Land NIR std: %.4f\n', medTemporalStd);
fprintf('Report: %s\n', reportPath);
if nWarn > 0
    fprintf('ACTION REQUIRED: Run Phase 1 (AROSICS) before PIF.\n');
else
    fprintf('OK: Proceed to Phase 1 (AROSICS) for sub-pixel refinement.\n');
end
fprintf('======================================\n');

%% ===================== LOCAL FUNCTIONS =====================

function cfg = loadPipelineConfig()
    cfg = [];
    global PIPELINE_CONFIG_PATH;
    % 1. Try global path set by RunAll
    if ~isempty(PIPELINE_CONFIG_PATH) && exist(PIPELINE_CONFIG_PATH,'file')
        cfg = load(PIPELINE_CONFIG_PATH);
        return;
    end
    % 2. Try current working directory
    p = fullfile(pwd, 'pipeline_config.mat');
    if exist(p,'file'), cfg = load(p); return; end
    % 3. Try same folder as this script (works when run standalone)
    p = fullfile(fileparts(mfilename('fullpath')), 'pipeline_config.mat');
    if exist(p,'file'), cfg = load(p); return; end
end

function d = extractDateFromName(fname)
    tok = regexp(fname, '(\d{8})', 'tokens');
    if ~isempty(tok), d = tok{1}{1}; else, d = fname(1:8); end
end

function entry = makeMetaEntry(file,H,W,cellX,cellY,originX,originY,status,notes)
    entry = struct('file',file,'H',H,'W',W,'cellX',cellX,'cellY',cellY,...
        'originX',originX,'originY',originY,'status',status,'notes',notes);
end

function [cellX,cellY,originX,originY] = extractSpatialRef(R)
    try
        if isprop(R,'CellExtentInWorldX')
            cellX=R.CellExtentInWorldX; cellY=R.CellExtentInWorldY;
            originX=R.XWorldLimits(1);  originY=R.YWorldLimits(2);
        elseif isprop(R,'CellExtentInLongitude')
            cellX=R.CellExtentInLongitude; cellY=R.CellExtentInLatitude;
            originX=R.LongitudeLimits(1);  originY=R.LatitudeLimits(2);
        else
            cellX=R.DeltaX; cellY=abs(R.DeltaY);
            originX=R.XWorldLimits(1); originY=R.YWorldLimits(2);
        end
    catch
        cellX=NaN; cellY=NaN; originX=NaN; originY=NaN;
    end
end

function m = extractCellSizeM(R)
    try
        if isprop(R,'CellExtentInWorldX')
            m = R.CellExtentInWorldX;
        elseif isprop(R,'CellExtentInLongitude')
            m = R.CellExtentInLongitude * 111320;
        else
            m = R.DeltaX;
        end
    catch
        m = 3.0;
    end
end

function roi = normalizeROI(roi)
    roi = roi - mean(roi(:));
    s   = std(roi(:));
    if s > 0, roi = roi / s; end
end

function [subR,subC,peakVal] = subpixelPeak(xc, peakR, peakC)
    [nr,nc] = size(xc);
    rSlice = double(xc(max(peakR-1,1):min(peakR+1,nr), peakC));
    cSlice = double(xc(peakR, max(peakC-1,1):min(peakC+1,nc)));
    subR=peakR; subC=peakC;
    if numel(rSlice)==3 && rSlice(2)>rSlice(1) && rSlice(2)>rSlice(3)
        d=2*(rSlice(1)-2*rSlice(2)+rSlice(3));
        if abs(d)>1e-10, subR=peakR+(rSlice(1)-rSlice(3))/d; end
    end
    if numel(cSlice)==3 && cSlice(2)>cSlice(1) && cSlice(2)>cSlice(3)
        d=2*(cSlice(1)-2*cSlice(2)+cSlice(3));
        if abs(d)>1e-10, subC=peakC+(cSlice(1)-cSlice(3))/d; end
    end
    peakVal=xc(peakR,peakC);
end
