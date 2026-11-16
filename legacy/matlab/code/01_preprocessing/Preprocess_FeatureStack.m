%% Preprocess_FeatureStack.m
clear; clc;
global PIPELINE_LOG_FID PIPELINE_CONFIG_PATH;
if isempty(PIPELINE_LOG_FID), PIPELINE_LOG_FID = -1; end
cfg = load(PIPELINE_CONFIG_PATH);

runTimestamp = char(datetime('now','TimeZone','local'),'yyyyMMdd_HHmmss');
outDir = fullfile(cfg.OUT_DIR, 'Phase5_Features', ['run_' runTimestamp]);
if ~exist(outDir, 'dir'), mkdir(outDir); end
logWrite(PIPELINE_LOG_FID, '  Feature mode: %s\n', cfg.FEATURE_MODE);
logWrite(PIPELINE_LOG_FID, 'Preprocess_FeatureStack.m - %s\n', runTimestamp);

p2 = load(fullfile(cfg.OUT_DIR, 'phase2_masking_output.mat'));
p4 = load(fullfile(cfg.OUT_DIR, 'phase4_wave_output.mat'));
cleanDir = p4.outDir;

ff = dir(fullfile(cleanDir, 'CLEAN_NORM_*.tif'));
if isempty(ff)
    error('No CLEAN_NORM rasters found in %s', cleanDir);
end
featNames = cfg.LR.FEATURE_NAMES;

msgs = cell(numel(ff),1);
usePar = isfield(cfg, 'PAR') && isfield(cfg.PAR, 'ENABLE') && cfg.PAR.ENABLE;
pool = gcp('nocreate');
if usePar && ~isempty(pool)
    try
        isThreadPool = isa(pool, 'parallel.ThreadPool');
        if ~isThreadPool
            try
                isThreadPool = contains(class(pool), 'ThreadPool', 'IgnoreCase', true);
            catch
            end
        end
        if ~isThreadPool
            try
                isThreadPool = strcmpi(pool.Cluster.Type, 'threads');
            catch
            end
        end
        if isThreadPool
            usePar = false;
            logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase5 (thread pool incompatible with geospatial IO)\n');
        end
    catch
        usePar = false;
        logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase5 (could not verify pool type safely)\n');
    end
end
if usePar
    parfor i = 1:numel(ff)
        msgs{i} = processFeatureDate(i, ff, cleanDir, outDir, p2.reef_mask, cfg);
    end
else
    for i = 1:numel(ff)
        msgs{i} = processFeatureDate(i, ff, cleanDir, outDir, p2.reef_mask, cfg);
    end
end
for i = 1:numel(msgs)
    if ~isempty(msgs{i})
        logWrite(PIPELINE_LOG_FID, '%s', msgs{i});
    end
end

reef_mask = p2.reef_mask;
deep_mask = p2.deep_mask;
sand_mask = p2.sand_mask;
save(fullfile(outDir, 'phase5_feature_output.mat'), ...
    'outDir','reef_mask','deep_mask','sand_mask','featNames');
save(fullfile(cfg.OUT_DIR, 'phase5_output.mat'), ...
    'outDir','reef_mask','deep_mask','sand_mask','featNames');

fid = fopen(fullfile(outDir, 'feature_report.txt'), 'w');
fprintf(fid, 'Feature mode: %s\n', cfg.FEATURE_MODE);
fprintf(fid, 'Feature order:\n');
for k = 1:numel(featNames)
    fprintf(fid, '  %d. %s\n', k, featNames{k});
end
fprintf(fid, '\nDII constants:\n');
fprintf(fid, '  k_BR = %.3f\n', cfg.LR.K_BR);
fprintf(fid, '  k_GR = %.3f\n', cfg.LR.K_GR);
if isfield(cfg.LR, 'GLCM_WINDOW')
    fprintf(fid, '  glcm_window = %d\n', cfg.LR.GLCM_WINDOW);
end
fclose(fid);

logWrite(PIPELINE_LOG_FID, '  Output   : %s\n\n', outDir);
return;

function msg = processFeatureDate(i, ff, cleanDir, outDir, reefMask, cfg)
msg = '';
tok = regexp(ff(i).name, 'CLEAN_NORM_(\d{8})\.tif', 'tokens', 'once');
if isempty(tok), return; end
d = tok{1};
[A, Rref] = readgeoraster(fullfile(cleanDir, ff(i).name));
A = single(A);
[feat, ~] = buildFeatureStack_LogRatio(A, cfg);
feat(~repmat(reefMask, [1 1 size(feat, 3)])) = nan;
saveMultiBandTif(fullfile(outDir, ['FEATURES_' d '.tif']), feat, Rref);
msg = sprintf('  Features %s (%d bands)\n', d, size(feat, 3));
end

%% Preprocess_DII.m
% ============================================================
% PHASE 5 - Sunglint Correction + DII + Feature Stack
%
% SUNGLINT CORRECTION - STANDARD HEDLEY (all dates):
%   Hedley regression on raw NIR over BM deep water mask.
%   Mask : BM union (bathy>300cm UNION MDeepWater.geojson)
%   ROI  : P95 outlier filter within BM (no global NIR ceiling)
%   Regression: polyfit(NIR_raw, band_raw) over reg_mask
%   Subtraction: band_corr = band - slope*(NIR - NIR_min)
%   REF_DATE (20230509): skip (no correction)
%
%   NOTE: Track B (temporal Hedley / low_nir) removed after diagnostic
%   showed T_scale collapse on heavy-wave dates. All dates standardised.
%
% MASK STRATEGY:
%   Deep water  - BM union:
%                   bathy_mask : bathymetry_0.tif >300cm (Method C)
%                   mask_M     : MDeepWater.geojson (~26180 px)
%                   mask_BM    = bathy_mask | mask_M
%                 NO global NIR ceiling applied to mask construction.
%                 P95 outlier filter applied inside ROI only.
%   Sand calib  - MSandMask.geojson   (manual, hand-verified)
%   Reef extent - GeoFeature.geojson "Reef" polygon
% NEW (line 47)
% Features cover the full image extent; reef mask applied at classification.
% NEW (line 637)
% fprintf(fid,'  Reef (GeoFeature/Reef): %d px\n\n', nnz(reef_mask));
%   Study area  - full Planet image extent (no additional clip)
%   ACA benthic - NOT used
%
% DII FORMULA:
%   DII_BR = log(max(Blue,  1e-4)) - k_BR * log(max(Red, 1e-4))
%   DII_GR = log(max(Green, 1e-4)) - k_GR * log(max(Red, 1e-4))
%   Log floor raised from 1e-6 to 1e-4 to avoid floor domination.
%
% K-RATIO STRATEGY:
%   REF_DATE uses shared k (it IS the shared source).
%   All other dates attempt per-date k from sand mask.
%   Falls back to shared if delta > K_MAX_DELTA.
%
% 12-band feature stack per image:
%   [1]  DII_BR          [2]  DII_GR         [3]  Brightness
%   [4]  RE_ratio        [5]  GR_ratio       [6]  CB_ratio
%   [7]  BG_ratio        [8]  YB_ratio       [9]  G2G1_ratio
%   [10] YR_ratio        [11] GLCM_Contrast (5x5)
%   [12] GLCM_Entropy (5x5)
%   NOTE: NIR_corr removed - encodes date identity (sea state) not
%         benthic class. RF exploits it as a date proxy.
%
% OUTPUTS -> OUT_DIR/Phase5_Features/run_YYYYMMDD_HHMMSS/
%   FEATURES_YYYYMMDD.tif     12-band float32 feature stack
%   GLINT_YYYYMMDD.tif        glint-corrected 8-band image
%   deep_mask.tif             BM deep water mask (uint8, shared)
%   sand_mask.tif             sand calibration mask (uint8)
%   reef_mask.tif             reef analysis extent (uint8)
%   DII_params.mat            all correction + QC parameters
%   phase5_viz.mat            full reconstruction data for figures
%   preprocessing_report.txt  human-readable QC report
%   quicklook_DII_YYYYMMDD.png 6-panel QC figure per date
% ============================================================

clear; clc;
global PIPELINE_LOG_FID;
if isempty(PIPELINE_LOG_FID)||PIPELINE_LOG_FID<0, PIPELINE_LOG_FID=-1; end
logW = @(varargin) logWrite(PIPELINE_LOG_FID, varargin{:});

runTime      = datetime('now','TimeZone','local');
runTimestamp = char(runTime,'yyyyMMdd_HHmmss');
logW('============================================================\n');
logW('Preprocess_DII.m - Phase 5 - %s\n', runTimestamp);
logW('Two-track BM Hedley sunglint correction\n');
logW('============================================================\n\n');

global PIPELINE_CONFIG_PATH;
cfg = load(PIPELINE_CONFIG_PATH);

%% "" Constants """"""""""""""""""""""""""""""""""""""""""""""""
BANDS            = cfg.BANDS;
DATE_CFG         = cfg.DATE_CFG;
REF_DATE         = cfg.REF_DATE;
FEAT_NAMES       = {'DII_BR','DII_GR','Brightness','RedEdge_ratio',...
                    'Green_Red_ratio','Coastal_Blue_ratio','Blue_Green_ratio',...
                    'Yellow_Blue_ratio','Green2_Green1_ratio','Yellow_Red_ratio',...
                    'GLCM_Contrast','GLCM_Entropy'};
N_FEATURES       = 12;
GLCM_WIN         = 5;
DEEP_THRESH      = 300;    % cm - bathy threshold for deep water
NIR_P95_FRAC     = 0.95;   % outlier filter within ROI (soft, not global ceiling)
K_MAX_DELTA      = 0.25;   % max allowed deviation of per-date k from shared k
K_CLAMP_LO       = 0.3;
K_CLAMP_HI       = 3.0;
LOG_FLOOR        = 1e-4;   % DII log floor (raised from 1e-6)
MIN_REG_PX       = 50;     % minimum pixels for valid regression
SAND_NIR_LO      = 0.005;   % k-ratio sand filter: exclude deeply submerged (NIR too low)
SAND_NIR_HI      = 0.12;   % k-ratio sand filter: exclude exposed/dry sand (NIR too high)

logW('CONSTANTS:\n');
logW('  DEEP_THRESH      = %d cm\n',   DEEP_THRESH);
logW('  NIR_P95_FRAC     = %.2f (outlier filter within ROI)\n', NIR_P95_FRAC);
logW('  K_MAX_DELTA      = %.2f (per-date vs shared fallback)\n', K_MAX_DELTA);
logW('  K_CLAMP          = [%.2f, %.2f]\n', K_CLAMP_LO, K_CLAMP_HI);
logW('  LOG_FLOOR        = %.0e (DII formula)\n', LOG_FLOOR);
logW('  MIN_REG_PX       = %d\n', MIN_REG_PX);
logW('  SAND_NIR_LO      = %.2f (k-ratio: exclude deeply submerged sand)\n', SAND_NIR_LO);
logW('  SAND_NIR_HI      = %.2f (k-ratio: exclude exposed/dry sand)\n\n', SAND_NIR_HI);

%% "" Find PIF input """""""""""""""""""""""""""""""""""""""""""
pifDir='';
p4Path=fullfile(cfg.OUT_DIR,'phase4_output.mat');
if exist(p4Path,'file')
    try p4=load(p4Path); if exist(p4.pifDir,'dir'), pifDir=p4.pifDir; end; catch, end
end
if isempty(pifDir)
    base=fullfile(cfg.OUT_DIR,'Phase3_PIF');
    if exist(base,'dir')
        runs=dir(fullfile(base,'run_*')); runs=runs([runs.isdir]);
        if ~isempty(runs)
            [~,si]=sort({runs.name});
            pifDir=fullfile(base,runs(si(end)).name);
        end
    end
end
if isempty(pifDir)||~exist(pifDir,'dir')
    pifDir=uigetdir(cfg.OUT_DIR,'Select Phase 3 PIF output folder');
    if isequal(pifDir,0), error('No folder selected.'); end
end

outDir=fullfile(cfg.OUT_DIR,'Phase5_Features',['run_' runTimestamp]);
mkdir(outDir);
logW('PIF input : %s\n', pifDir);
logW('Output    : %s\n\n', outDir);

%% "" Find + sort PIF files """"""""""""""""""""""""""""""""""""
ff=dir(fullfile(pifDir,'*_PIF.tif'));
if isempty(ff), error('No *_PIF.tif in %s',pifDir); end
dates=cellfun(@(n) regexp(n,'(\d{8})','tokens','once'),{ff.name},'UniformOutput',false);
dates=cellfun(@(t) t{1},dates,'UniformOutput',false);
[dates,si]=sort(dates); ff=ff(si); nDates=numel(ff);
logW('PIF images found: %d\n',nDates);
for i=1:nDates
    dcfg_tmp=getDatConfig(DATE_CFG,dates{i});
    logW('  [%d] %s  mode=%-8s\n',i,dates{i},dcfg_tmp.hedley);
end
logW('\n');

%% "" Load reference image """""""""""""""""""""""""""""""""""""
refIdx=find(strcmp(dates,REF_DATE));
if isempty(refIdx), error('Reference date %s not found in PIF folder',REF_DATE); end
[Aref,Rref]=readgeoraster(fullfile(ff(refIdx).folder,ff(refIdx).name));
Aref=single(Aref); [H,W,nBands]=size(Aref);
logW('Reference image : %s  (%d x %d px, %d bands)\n\n',REF_DATE,W,H,nBands);

%% "" BUILD SHARED MASKS """""""""""""""""""""""""""""""""""""""
logW('------------------------------------------------------------\n');
logW('BUILDING MASKS\n');
logW('------------------------------------------------------------\n');

% 1a) Bathy deep water - Method C scatteredInterpolant WGS84->UTM47N
acaPath=fullfile(cfg.ACA_DIR,'bathymetry_0.tif');
if ~exist(acaPath,'file'), error('bathymetry_0.tif not found: %s',acaPath); end
[Ab,Rab]=readgeoraster(acaPath);
depth_grid=reprojectDepth(single(Ab),Rab,Rref,H,W);
bathy_mask = isfinite(depth_grid) & (depth_grid > DEEP_THRESH);
logW('  Bathy deep (>%dcm)    : %d px  (%.1f%% coverage in depth grid)\n',...
    DEEP_THRESH, nnz(bathy_mask), ...
    100*nnz(bathy_mask)/max(1,nnz(isfinite(depth_grid))));

% 1b) Manual deep water polygon - MDeepWater.geojson
deepWaterPath=fullfile(cfg.ACA_DIR,'MDeepWater.geojson');
if ~exist(deepWaterPath,'file')
    logW('  WARNING: MDeepWater.geojson not found - BM mask = bathy only\n');
    mask_M=false(H,W);
else
    mask_M=rasterizeGeoJSON(deepWaterPath,H,W,Rref,'MDeepWater');
end
logW('  MDeepWater polygon    : %d px\n', nnz(mask_M));

% 1c) BM union - no NIR ceiling
mask_BM = bathy_mask | mask_M;
logW('  BM union (bathy|M)    : %d px  <-- primary regression ROI\n', nnz(mask_BM));
if nnz(mask_BM) < MIN_REG_PX
    error('FATAL: BM mask too small (%d px). Check bathymetry and MDeepWater paths.',nnz(mask_BM));
end

% 1d) Reference NIR statistics for QC
nir_ref_full = Aref(:,:,BANDS.NIR);
nir_ref_BM   = nir_ref_full(mask_BM);
NIR_ref_BM_mean = mean(nir_ref_BM,'omitnan');
NIR_ref_BM_std  = std(nir_ref_BM,'omitnan');
logW('  NIR_ref over BM mask  : mean=%.5f  std=%.5f  range=[%.5f, %.5f]\n',...
    NIR_ref_BM_mean, NIR_ref_BM_std, ...
    min(nir_ref_BM), max(nir_ref_BM));

% 2) Sand mask
sandPath=fullfile(cfg.ACA_DIR,'MSandMask.geojson');
sand_mask=rasterizeGeoJSON(sandPath,H,W,Rref,'MSandMask');
logW('  Sand (manual)         : %d px\n', nnz(sand_mask));
if nnz(sand_mask)<20, logW('  WARNING: very few sand pixels\n'); end

% 3) Reef extent
reefPath=cfg.REEF_MASK_PATH;   % GeoFeature.geojson set in RunAll_New
reef_mask=rasterizeGeoJSONClass(reefPath,H,W,Rref,'Reef');
logW('  Reef extent (manual)  : %d px\n', nnz(reef_mask));

% Reference NIR statistics for QC
valid_ref = isfinite(nir_ref_full) & nir_ref_full > 0;
nir_sand  = nir_ref_full(sand_mask & valid_ref);
if ~isempty(nir_sand)
    logW('  Sand NIR (REF, expect<0.05): mean=%.5f  std=%.5f\n',...
        mean(nir_sand), std(nir_sand));
end
nir_reef = nir_ref_full(reef_mask & valid_ref);
if ~isempty(nir_reef)
    logW('  Reef NIR (REF)     : mean=%.5f  P10=%.5f  P90=%.5f\n',...
        mean(nir_reef), prctile(nir_reef,10), prctile(nir_reef,90));
end
logW('\n');

% Save shared masks
saveMask(fullfile(outDir,'deep_mask.tif'), mask_BM,   Rref);
saveMask(fullfile(outDir,'sand_mask.tif'), sand_mask, Rref);
saveMask(fullfile(outDir,'reef_mask.tif'), reef_mask, Rref);

%% "" Master sand mask for k-ratio calibration """""""""""""""""
% Sand pixels at multiple depths give reef-specific water column attenuation.
% Falls back to BM deep water mask if file not found.
logW('------------------------------------------------------------\n');
logW('MASTER SAND MASK (for k-ratio calibration)\n');
logW('------------------------------------------------------------\n');

master_sand_mask = [];
% Try StableSand first (near-island, spectrally stable, better depth range)
if isfield(cfg,'STABLE_SAND_MASK') && exist(cfg.STABLE_SAND_MASK,'file')
    master_sand_mask = rasterizeGeoJSON(cfg.STABLE_SAND_MASK,H,W,Rref,'StableSand');
    logW('  StableSand mask: %d px\n', nnz(master_sand_mask));
    if nnz(master_sand_mask) < 30
        logW('  WARNING: too few StableSand pixels - trying MasterSand\n');
        master_sand_mask = [];
    end
end
% Fallback to MasterSand
if isempty(master_sand_mask) && isfield(cfg,'MASTER_SAND_MASK') && exist(cfg.MASTER_SAND_MASK,'file')
    master_sand_mask = rasterizeGeoJSON(cfg.MASTER_SAND_MASK,H,W,Rref,'MasterSand');

    logW('  Master sand mask: %d px\n', nnz(master_sand_mask));
    if nnz(master_sand_mask) < 30
        logW('  WARNING: too few pixels - falling back to BM mask for k-ratios\n');
        master_sand_mask = [];
    end
end
if isempty(master_sand_mask)
    logW('  MasterSand.geojson not found - using BM mask for k-ratios\n');
    master_sand_mask = mask_BM;
end
logW('\n');

%% "" Shared k-ratios from master sand mask on REF_DATE """"""""
logW('------------------------------------------------------------\n');
logW('SHARED K-RATIOS (from %s, master sand mask + NIR filter)\n', REF_DATE);
logW('------------------------------------------------------------\n');
% Per-date NIR filter: exclude exposed dry sand (NIR>SAND_NIR_HI)
% and deeply submerged sand (NIR<SAND_NIR_LO) - intertidal transition pixels
nir_ref_img = Aref(:,:,BANDS.NIR);
sand_k_ref  = master_sand_mask & (nir_ref_img > SAND_NIR_LO) & (nir_ref_img < SAND_NIR_HI);
logW('  Master sand: %d px  after NIR filter [%.2f,%.2f]: %d px\n',...
    nnz(master_sand_mask), SAND_NIR_LO, SAND_NIR_HI, nnz(sand_k_ref));
if nnz(sand_k_ref) < MIN_REG_PX
    logW('  WARNING: too few filtered sand px (%d) - using unfiltered mask\n', nnz(sand_k_ref));
    sand_k_ref = master_sand_mask;
end
[k_BR_shared, k_GR_shared, k_diag_shared] = computeKratios(Aref, sand_k_ref, BANDS, ...
    K_CLAMP_LO, K_CLAMP_HI, LOG_FLOOR);
logW('  k_BR_shared = %.4f  (R2=%.4f, n=%d)\n', ...
    k_BR_shared, k_diag_shared.r2_BR, k_diag_shared.n_px);
logW('  k_GR_shared = %.4f  (R2=%.4f)\n', k_GR_shared, k_diag_shared.r2_GR);
logW('  lR range    = [%.4f, %.4f]\n\n', k_diag_shared.lR_min, k_diag_shared.lR_max);

%% "" Main processing loop """""""""""""""""""""""""""""""""""""
% Per-date result structs - preallocate
dii_params = struct('date',{},'track',{},'hedley_mode',{},...
    'NIR_BM_mean',{},'NIR_BM_std',{},'NIR_BM_range',{},...
    'reg_px',{},...
    'slopes',{},'slope_nir',{},...
    'k_BR',{},'k_GR',{},'k_source',{},...
    'k_BR_raw',{},'k_GR_raw',{},'k_BR_delta',{},'k_GR_delta',{},...
    'DII_BR_mean',{},'DII_BR_std',{},'DII_BR_p5',{},'DII_BR_p95',{},...
    'DII_GR_mean',{},'DII_GR_std',{},'DII_GR_p5',{},'DII_GR_p95',{},...
    'NIR_raw_BM_mean',{},'NIR_corr_BM_mean',{},'NIR_reduction',{},...
    'hedley_ok',{},'hedley_skip_reason',{});

% Accumulate per-date data for viz .mat
viz_data = struct();
viz_data.dates      = dates;
viz_data.REF_DATE   = REF_DATE;
viz_data.mask_BM    = mask_BM;
viz_data.mask_M     = mask_M;
viz_data.bathy_mask = bathy_mask;
viz_data.sand_mask  = sand_mask;
viz_data.reef_mask  = reef_mask;
viz_data.depth_grid = depth_grid;
viz_data.BANDS      = BANDS;
viz_data.FEAT_NAMES = FEAT_NAMES;
viz_data.k_BR_shared     = k_BR_shared;
viz_data.k_GR_shared     = k_GR_shared;
viz_data.k_diag_shared   = k_diag_shared;
viz_data.NIR_ref_BM_mean = NIR_ref_BM_mean;
viz_data.per_date        = struct([]);   % filled in loop

report = {};   % text report lines

for i = 1:nDates
    logW('============================================================\n');
    logW('[%d/%d] %s\n', i, nDates, dates{i});
    logW('============================================================\n');

    dcfg  = getDatConfig(DATE_CFG, dates{i});
    track = ternary(strcmp(dcfg.hedley,'skip'),'C','A');
    logW('  Mode : %s  (hedley=%s)\n', track, dcfg.hedley);

    [Ai,Ri] = readgeoraster(fullfile(ff(i).folder, ff(i).name));
    Ai = single(Ai);
    nir_i = Ai(:,:,BANDS.NIR);

    % "" BM mask statistics for this date """"""""""""""""""""""
    nir_BM_i    = nir_i(mask_BM);
    NIR_BM_mean = mean(nir_BM_i, 'omitnan');
    NIR_BM_std  = std(nir_BM_i,  'omitnan');
    NIR_BM_min  = min(nir_BM_i);
    NIR_BM_max  = max(nir_BM_i);
    NIR_BM_range= NIR_BM_max - NIR_BM_min;
    logW('  NIR over BM mask : mean=%.5f  std=%.5f  range=%.5f  [%.5f, %.5f]\n',...
        NIR_BM_mean, NIR_BM_std, NIR_BM_range, NIR_BM_min, NIR_BM_max);

    % P95 outlier filter - applied within BM, not as a global ceiling
    nir_p95_i  = prctile(nir_BM_i, NIR_P95_FRAC * 100);
    reg_mask_i = mask_BM & (nir_i < nir_p95_i);
    n_reg      = nnz(reg_mask_i);
    logW('  P95 NIR cutoff   : %.5f  ->  reg_mask = %d px\n', nir_p95_i, n_reg);

    % "" Sunglint correction (unified standard Hedley) """""""""
    hedley_ok       = false;
    skip_reason     = '';
    slopes_struct   = struct();
    slope_nir       = 0;
    NIR_resid       = [];   % unused - kept for quicklook compatibility

    if strcmp(dcfg.hedley,'skip')
        logW('  --- SKIP (REF_DATE) ---\n');
        Ai_corr     = Ai;
        hedley_ok   = true;
        skip_reason = 'REF_DATE - no correction applied';
        logW('  %s\n', skip_reason);
    elseif n_reg < MIN_REG_PX
        skip_reason = sprintf('reg_mask too small (%d px)',n_reg);
        logW('  SKIP: %s\n', skip_reason);
        Ai_corr = Ai;
    else
        logW('  --- Standard Hedley ---\n');
        [Ai_corr, slope_nir, slopes_struct, hedley_ok, skip_reason] = ...
            hedleyCorrectA(Ai, reg_mask_i, BANDS, MIN_REG_PX);
        logW('  slope_NIR (RED band) : %.5f\n', slope_nir);
        logW('  Hedley OK            : %d\n',   hedley_ok);
        if ~isempty(skip_reason)
            logW('  Skip reason          : %s\n', skip_reason);
        end
        logW('  Band slopes:\n');
        fns = fieldnames(slopes_struct);
        for fi=1:numel(fns)
            logW('    %-10s: %.5f\n', fns{fi}, slopes_struct.(fns{fi}));
        end
    end

    % "" NIR reduction QC """"""""""""""""""""""""""""""""""""""
    NIR_raw_BM_mean  = NIR_BM_mean;
    NIR_corr_BM_mean = mean(Ai_corr(:,:,BANDS.NIR), 'all','omitnan');  % scene-wide
    nir_corr_flat    = Ai_corr(:,:,BANDS.NIR);
    NIR_corr_BM_deep = mean(nir_corr_flat(mask_BM), 'omitnan');        % over BM only
    NIR_reduction    = NIR_raw_BM_mean - NIR_corr_BM_deep;
    logW('  NIR BM mean  raw : %.5f\n', NIR_raw_BM_mean);
    logW('  NIR BM mean corr : %.5f\n', NIR_corr_BM_deep);
    logW('  NIR reduction    : %.5f  (%.1f%%)\n',...
        NIR_reduction, 100*NIR_reduction/max(NIR_raw_BM_mean,1e-9));

    % "" DII k-ratios """"""""""""""""""""""""""""""""""""""""""
    logW('  --- K-ratios ---\n');
    use_shared_k = strcmp(dcfg.hedley,'skip');  % only REF uses shared

    if use_shared_k
        k_BR     = k_BR_shared;
        k_GR     = k_GR_shared;
        k_BR_raw = k_BR_shared;
        k_GR_raw = k_GR_shared;
        k_src    = 'shared';
        k_BR_delta = 0;
        k_GR_delta = 0;
        logW('  k-ratios forced shared (REF_DATE)\n');
    else
        % Track A: compute per-date, fall back if deviation too large
        % Apply per-date NIR filter to sand mask - excludes intertidal
        % transition pixels (exposed dry sand NIR>HI, deeply submerged NIR<LO)
        sand_k_i = master_sand_mask & (nir_i > SAND_NIR_LO) & (nir_i < SAND_NIR_HI);
        logW('  Sand k-mask after NIR filter [%.2f,%.2f]: %d px\n',...
            SAND_NIR_LO, SAND_NIR_HI, nnz(sand_k_i));
        if nnz(sand_k_i) < MIN_REG_PX
            logW('  Too few filtered sand px - using unfiltered mask\n');
            sand_k_i = master_sand_mask;
        end
        [k_BR_raw, k_GR_raw, k_diag_i] = computeKratios(Ai_corr, sand_k_i, BANDS,...
            K_CLAMP_LO, K_CLAMP_HI, LOG_FLOOR);
        k_BR_delta = abs(k_BR_raw - k_BR_shared);
        k_GR_delta = abs(k_GR_raw - k_GR_shared);
        logW('  Per-date k_BR=%.4f  k_GR=%.4f  (R2_BR=%.4f, n=%d)\n',...
            k_BR_raw, k_GR_raw, k_diag_i.r2_BR, k_diag_i.n_px);
        logW('  Delta from shared: k_BR=%.4f  k_GR=%.4f  (threshold=%.2f)\n',...
            k_BR_delta, k_GR_delta, K_MAX_DELTA);

        if k_BR_delta > K_MAX_DELTA || k_GR_delta > K_MAX_DELTA
            k_BR = k_BR_shared;
            k_GR = k_GR_shared;
            k_src = 'shared_fallback';
            logW('  FALLBACK to shared (delta exceeded threshold)\n');
        else
            k_BR = k_BR_raw;
            k_GR = k_GR_raw;
            k_src = 'per_date';
        end
    end
    logW('  k_BR=%.4f  k_GR=%.4f  [%s]\n', k_BR, k_GR, k_src);

    % "" DII computation """""""""""""""""""""""""""""""""""""""
    logW('  --- DII ---\n');
    B_log = log(max(Ai_corr(:,:,BANDS.BLUE),  LOG_FLOOR));
    G_log = log(max(Ai_corr(:,:,BANDS.GREEN1),LOG_FLOOR));
    R_log = log(max(Ai_corr(:,:,BANDS.RED),   LOG_FLOOR));
    DII_BR = B_log - k_BR .* R_log;
    DII_GR = G_log - k_GR .* R_log;

    % DII statistics over reef only
    dii_br_reef = DII_BR(reef_mask);
    dii_gr_reef = DII_GR(reef_mask);
    DII_BR_mean = mean(dii_br_reef,'omitnan');
    DII_BR_std  = std( dii_br_reef,'omitnan');
    DII_BR_p5   = prctile(dii_br_reef, 5);
    DII_BR_p95  = prctile(dii_br_reef,95);
    DII_GR_mean = mean(dii_gr_reef,'omitnan');
    DII_GR_std  = std( dii_gr_reef,'omitnan');
    DII_GR_p5   = prctile(dii_gr_reef, 5);
    DII_GR_p95  = prctile(dii_gr_reef,95);
    logW('  DII_BR (reef): mean=%.4f  std=%.4f  [P5=%.4f, P95=%.4f]\n',...
        DII_BR_mean, DII_BR_std, DII_BR_p5, DII_BR_p95);
    logW('  DII_GR (reef): mean=%.4f  std=%.4f  [P5=%.4f, P95=%.4f]\n',...
        DII_GR_mean, DII_GR_std, DII_GR_p5, DII_GR_p95);

    % Sanity check: DII should be negative over water
    pct_positive_BR = 100 * mean(dii_br_reef > 0, 'omitnan');
    pct_positive_GR = 100 * mean(dii_gr_reef > 0, 'omitnan');
    logW('  DII_BR >0 in reef: %.1f%%  DII_GR >0 in reef: %.1f%%\n',...
        pct_positive_BR, pct_positive_GR);
    if pct_positive_BR > 20
        logW('  WARNING: >20%% positive DII_BR - check glint correction\n');
    end
    if DII_BR_p95 > 1.0
        logW('  WARNING: DII_BR P95=%.2f > 1.0 - possible glint residual or k issue\n',...
            DII_BR_p95);
    end

    % "" Spectral features """""""""""""""""""""""""""""""""""""
    Brightness = (Ai_corr(:,:,BANDS.BLUE) + Ai_corr(:,:,BANDS.GREEN1) + ...
                  Ai_corr(:,:,BANDS.RED)) / 3;
    RE_ratio   = Ai_corr(:,:,BANDS.REDEDGE) ./ max(Ai_corr(:,:,BANDS.RED),   1e-6);
    GR_ratio   = Ai_corr(:,:,BANDS.GREEN1)  ./ max(Ai_corr(:,:,BANDS.RED),   1e-6);
    CB_ratio   = Ai_corr(:,:,BANDS.COASTAL) ./ max(Ai_corr(:,:,BANDS.BLUE),  1e-6);
    BG_ratio   = Ai_corr(:,:,BANDS.BLUE)    ./ max(Ai_corr(:,:,BANDS.GREEN1),1e-6);
    YB_ratio   = Ai_corr(:,:,BANDS.YELLOW)  ./ max(Ai_corr(:,:,BANDS.BLUE),  1e-6);
    G2G1_ratio = Ai_corr(:,:,BANDS.GREEN2)  ./ max(Ai_corr(:,:,BANDS.GREEN1),1e-6);
    YR_ratio   = Ai_corr(:,:,BANDS.YELLOW)  ./ max(Ai_corr(:,:,BANDS.RED),   1e-6);

    % Brightness statistics over reef
    bright_reef = Brightness(reef_mask);
    logW('  Brightness (reef): mean=%.5f  std=%.5f\n',...
        mean(bright_reef,'omitnan'), std(bright_reef,'omitnan'));

    % "" GLCM texture """"""""""""""""""""""""""""""""""""""""""
    logW('  GLCM (%dx%d window)...\n', GLCM_WIN, GLCM_WIN);
    [glcm_contrast, glcm_entropy] = computeGLCM(Brightness, GLCM_WIN);

    % "" Assemble 12-band feature stack - NaN outside reef """""
    feat = cat(3, DII_BR, DII_GR, Brightness, RE_ratio, GR_ratio,...
                  CB_ratio, BG_ratio, YB_ratio, G2G1_ratio, YR_ratio,...
                  glcm_contrast, glcm_entropy);

    % "" Save TIFs """""""""""""""""""""""""""""""""""""""""""""
    saveMultiBandTif(fullfile(outDir,['GLINT_'    dates{i} '.tif']), Ai_corr, Ri);
    saveMultiBandTif(fullfile(outDir,['FEATURES_' dates{i} '.tif']), feat,    Ri);
    logW('  Saved: GLINT_%s.tif  FEATURES_%s.tif\n', dates{i}, dates{i});

    % "" 6-panel quicklook """""""""""""""""""""""""""""""""""""
    makeQuicklook(outDir, dates{i}, Ai, Ai_corr, nir_i, NIR_resid, ...
        DII_BR, DII_GR, depth_grid, reef_mask, mask_BM, reg_mask_i, ...
        track, NaN, BANDS);

    % "" Store per-date params """"""""""""""""""""""""""""""""""
    dii_params(end+1) = struct(...
        'date',          dates{i},...
        'track',         track,...
        'hedley_mode',   dcfg.hedley,...
        'NIR_BM_mean',   NIR_BM_mean,...
        'NIR_BM_std',    NIR_BM_std,...
        'NIR_BM_range',  NIR_BM_range,...
        'reg_px',        n_reg,...
        'slopes',        slopes_struct,...
        'slope_nir',     slope_nir,...
        'k_BR',          k_BR,...
        'k_GR',          k_GR,...
        'k_source',      k_src,...
        'k_BR_raw',      k_BR_raw,...
        'k_GR_raw',      k_GR_raw,...
        'k_BR_delta',    k_BR_delta,...
        'k_GR_delta',    k_GR_delta,...
        'DII_BR_mean',   DII_BR_mean,...
        'DII_BR_std',    DII_BR_std,...
        'DII_BR_p5',     DII_BR_p5,...
        'DII_BR_p95',    DII_BR_p95,...
        'DII_GR_mean',   DII_GR_mean,...
        'DII_GR_std',    DII_GR_std,...
        'DII_GR_p5',     DII_GR_p5,...
        'DII_GR_p95',    DII_GR_p95,...
        'NIR_raw_BM_mean',  NIR_raw_BM_mean,...
        'NIR_corr_BM_mean', NIR_corr_BM_deep,...
        'NIR_reduction',    NIR_reduction,...
        'hedley_ok',     hedley_ok,...
        'hedley_skip_reason', skip_reason);

    % Store in viz struct (lightweight - indices only, no full images)
    viz_data.per_date(i).date          = dates{i};
    viz_data.per_date(i).track         = track;
    viz_data.per_date(i).k_BR          = k_BR;
    viz_data.per_date(i).k_GR          = k_GR;
    viz_data.per_date(i).k_source      = k_src;
    viz_data.per_date(i).k_BR_raw      = k_BR_raw;
    viz_data.per_date(i).k_GR_raw      = k_GR_raw;
    viz_data.per_date(i).slopes        = slopes_struct;
    viz_data.per_date(i).slope_nir     = slope_nir;
    viz_data.per_date(i).NIR_BM_mean   = NIR_BM_mean;
    viz_data.per_date(i).NIR_BM_std    = NIR_BM_std;
    viz_data.per_date(i).NIR_BM_range  = NIR_BM_range;
    viz_data.per_date(i).reg_px        = n_reg;
    viz_data.per_date(i).DII_BR_mean   = DII_BR_mean;
    viz_data.per_date(i).DII_BR_std    = DII_BR_std;
    viz_data.per_date(i).DII_BR_p5     = DII_BR_p5;
    viz_data.per_date(i).DII_BR_p95    = DII_BR_p95;
    viz_data.per_date(i).DII_GR_mean   = DII_GR_mean;
    viz_data.per_date(i).DII_GR_std    = DII_GR_std;
    viz_data.per_date(i).DII_GR_p5     = DII_GR_p5;
    viz_data.per_date(i).DII_GR_p95    = DII_GR_p95;
    viz_data.per_date(i).NIR_raw_BM_mean  = NIR_raw_BM_mean;
    viz_data.per_date(i).NIR_corr_BM_mean = NIR_corr_BM_deep;
    viz_data.per_date(i).NIR_reduction    = NIR_reduction;
    viz_data.per_date(i).hedley_ok        = hedley_ok;
    viz_data.per_date(i).hedley_skip_reason = skip_reason;
    viz_data.per_date(i).pct_positive_BR  = pct_positive_BR;
    viz_data.per_date(i).pct_positive_GR  = pct_positive_GR;

    % Report line
    report{end+1} = sprintf(...
        '%s | track=%-1s mode=%-8s ok=%-1d | reg_px=%6d | slope_nir=%7.4f | k_BR=%6.4f k_GR=%6.4f [%-16s] | DII_BR=%6.3f%5.3f [%6.3f,%6.3f] | DII_GR=%6.3f%5.3f',...
        dates{i}, track, dcfg.hedley, hedley_ok,...
        n_reg, slope_nir,...
        k_BR, k_GR, k_src,...
        DII_BR_mean, DII_BR_std, DII_BR_p5, DII_BR_p95,...
        DII_GR_mean, DII_GR_std);

    logW('\n');
end  % date loop

%% "" Save DII_params.mat """"""""""""""""""""""""""""""""""""""
save(fullfile(outDir,'DII_params.mat'),...
    'dii_params',...
    'k_BR_shared','k_GR_shared','k_diag_shared',...
    'FEAT_NAMES','GLCM_WIN','DEEP_THRESH',...
    'NIR_P95_FRAC',...
    'K_MAX_DELTA','K_CLAMP_LO','K_CLAMP_HI','LOG_FLOOR',...
    'MIN_REG_PX','dates','REF_DATE',...
    'NIR_ref_BM_mean','NIR_ref_BM_std');
logW('Saved: DII_params.mat\n');

%% "" Save phase5_viz.mat (full reconstruction) """""""""""""""
% Add run metadata
viz_data.runTimestamp    = runTimestamp;
viz_data.outDir          = outDir;
viz_data.pifDir          = pifDir;
viz_data.DEEP_THRESH     = DEEP_THRESH;
viz_data.NIR_P95_FRAC    = NIR_P95_FRAC;
viz_data.K_MAX_DELTA     = K_MAX_DELTA;
viz_data.K_CLAMP_LO      = K_CLAMP_LO;
viz_data.K_CLAMP_HI      = K_CLAMP_HI;
viz_data.LOG_FLOOR       = LOG_FLOOR;
viz_data.MIN_REG_PX      = MIN_REG_PX;
viz_data.GLCM_WIN        = GLCM_WIN;
viz_data.N_FEATURES      = N_FEATURES;
viz_data.image_size      = [H, W];
viz_data.Rref            = Rref;
viz_data.k_BR_shared     = k_BR_shared;
viz_data.k_GR_shared     = k_GR_shared;
viz_data.k_diag_shared   = k_diag_shared;
viz_data.NIR_ref_BM_mean = NIR_ref_BM_mean;
viz_data.NIR_ref_BM_std  = NIR_ref_BM_std;
viz_data.report          = report;

save(fullfile(outDir,'phase5_viz.mat'), 'viz_data', '-v7.3');
logW('Saved: phase5_viz.mat\n');

%% "" Save phase5_output.mat (handoff to Phase 6/7) """""""""""
save(fullfile(cfg.OUT_DIR,'phase5_output.mat'),...
    'outDir','dates','FEAT_NAMES',...
    'mask_BM','mask_M','bathy_mask','sand_mask','reef_mask',...
    'depth_grid','k_BR_shared','k_GR_shared','dii_params');
logW('Saved: phase5_output.mat (handoff)\n\n');

%% "" Write preprocessing_report.txt """"""""""""""""""""""""""
fid=fopen(fullfile(outDir,'preprocessing_report.txt'),'w');
fprintf(fid,'PHASE 5 - PREPROCESSING REPORT\n');
fprintf(fid,'Timestamp : %s\n', runTimestamp);
fprintf(fid,'PIF input : %s\n', pifDir);
fprintf(fid,'Output    : %s\n\n', outDir);

fprintf(fid,'MASK SOURCES\n');
fprintf(fid,'  Bathy deep (>%dcm)    : %d px\n', DEEP_THRESH, nnz(bathy_mask));
fprintf(fid,'  MDeepWater polygon    : %d px\n', nnz(mask_M));
fprintf(fid,'  BM union              : %d px\n', nnz(mask_BM));
fprintf(fid,'  Sand (MSandMask)      : %d px\n', nnz(sand_mask));
fprintf(fid,'  Reef (GeoFeature/Reef): %d px\n\n', nnz(reef_mask));

fprintf(fid,'CONSTANTS\n');
fprintf(fid,'  NIR_P95_FRAC     = %.2f\n', NIR_P95_FRAC);
fprintf(fid,'  K_MAX_DELTA      = %.2f\n', K_MAX_DELTA);
fprintf(fid,'  K_CLAMP          = [%.2f, %.2f]\n', K_CLAMP_LO, K_CLAMP_HI);
fprintf(fid,'  LOG_FLOOR        = %.0e\n', LOG_FLOOR);
fprintf(fid,'  MIN_REG_PX       = %d\n\n', MIN_REG_PX);

fprintf(fid,'SHARED K-RATIOS (from %s, BM mask, n=%d px)\n', REF_DATE, k_diag_shared.n_px);
fprintf(fid,'  k_BR = %.4f  (R2=%.4f)\n', k_BR_shared, k_diag_shared.r2_BR);
fprintf(fid,'  k_GR = %.4f  (R2=%.4f)\n\n', k_GR_shared, k_diag_shared.r2_GR);

fprintf(fid,'PER-DATE RESULTS\n');
fprintf(fid,'%s\n', repmat('-',1,140));
hdr = sprintf('%-8s | Trk mode     ok | reg_px | slope_nir | k_BR   k_GR  [source          ] | DII_BR meanstd [P5,P95]       | DII_GR meanstd',...
    'date');
fprintf(fid,'%s\n',hdr);
fprintf(fid,'%s\n', repmat('-',1,140));
for k_=1:numel(report), fprintf(fid,'%s\n',report{k_}); end
fprintf(fid,'%s\n\n', repmat('-',1,140));

fprintf(fid,'TRACK LEGEND\n');
fprintf(fid,'  A = Standard Hedley on raw NIR\n');
fprintf(fid,'  C = Skip (REF_DATE)\n\n');
fprintf(fid,'FEATURE STACK BANDS\n');
for k_=1:N_FEATURES, fprintf(fid,'  [%d] %s\n',k_,FEAT_NAMES{k_}); end
fclose(fid);
logW('Saved: preprocessing_report.txt\n\n');

%% "" Summary """""""""""""""""""""""""""""""""""""""""""""""""
logW('============================================================\n');
logW('PHASE 5 SUMMARY\n');
logW('============================================================\n');
logW('Dates processed : %d\n', nDates);
logW('BM mask px      : %d\n', nnz(mask_BM));
logW('Sand mask px    : %d\n', nnz(sand_mask));
logW('Reef mask px    : %d\n', nnz(reef_mask));
logW('k_BR_shared     : %.4f\n', k_BR_shared);
logW('k_GR_shared     : %.4f\n', k_GR_shared);
logW('\nPer-date summary:\n');
logW('%-10s  Trk  ok  reg_px  slope_nir  k_BR    k_GR    DII_BR_mean  DII_GR_mean\n','Date');
for i=1:numel(dii_params)
    p=dii_params(i);
    logW('%-10s   %s   %d  %6d   %7.4f  %6.4f  %6.4f   %8.4f    %8.4f\n',...
        p.date, p.track, p.hedley_ok,...
        p.reg_px, p.slope_nir, p.k_BR, p.k_GR,...
        p.DII_BR_mean, p.DII_GR_mean);
end
logW('\nOutput: %s\n', outDir);
logW('NEXT  : Build_Baseline.m (Phase 6)\n');
logW('============================================================\n');


%% =====================================================================
%% LOCAL FUNCTIONS
%% =====================================================================

% ------------------------------------------------------------------
function [Ac, slope_nir, slopes, ok, reason] = hedleyCorrectA(A, reg_mask, BANDS, min_px)
% Track A - Standard Hedley on raw NIR
% Regresses each visible band against raw NIR over reg_mask.
% Subtraction: band_corr = band - slope*(NIR - NIR_min)
    Ac = A; slope_nir = 0; slopes = struct(); ok = false; reason = '';

    nir     = double(A(:,:,BANDS.NIR));
    nir_vec = nir(reg_mask);
    nir_min = min(nir_vec);
    nir_range = max(nir_vec) - nir_min;

    if numel(nir_vec) < min_px
        reason = sprintf('reg_mask %d px < min %d', numel(nir_vec), min_px);
        return;
    end
    if nir_range < 1e-5
        reason = sprintf('NIR range %.2e < 1e-5 (degenerate)', nir_range);
        return;
    end

    vis_bands = [BANDS.COASTAL, BANDS.BLUE, BANDS.GREEN1, BANDS.GREEN2,...
                 BANDS.YELLOW, BANDS.RED, BANDS.REDEDGE];
    band_names = {'Coastal','Blue','Green1','Green2','Yellow','Red','RedEdge'};

    for bi = 1:numel(vis_bands)
        b = vis_bands(bi);
        bv = double(A(:,:,b));
        band_vec = bv(reg_mask);
        if std(band_vec) < 1e-7, continue; end
        [p,~,mu] = polyfit(nir_vec, band_vec, 1);
        slope    = p(1) / mu(2);
        % Reject physically impossible slopes
        if slope < 0 || slope > 8
            slopes.(band_names{bi}) = NaN;
            continue;
        end
        slopes.(band_names{bi}) = slope;
        if b == BANDS.RED, slope_nir = slope; end
        corr = single(A(:,:,b)) - single(slope) .* (single(nir) - single(nir_min));
        Ac(:,:,b) = max(0, corr);
    end
    ok = true;
end

% ------------------------------------------------------------------
function [k_BR, k_GR, diag] = computeKratios(A, mask, BANDS, clo, chi, logfloor)
% Log-log regression: log(Blue) ~ k_BR*log(Red), log(Green) ~ k_GR*log(Red)
% Returns k values clamped to [clo, chi] and diagnostics.
    diag = struct('n_px',0,'r2_BR',NaN,'r2_GR',NaN,'lR_min',NaN,'lR_max',NaN);
    k_BR = 1.2; k_GR = 1.1;
    dm   = mask(:);
    if nnz(dm) < 30
        warning('computeKratios: too few px (%d)', nnz(dm));
        return;
    end
    Bv = double(A(:,:,BANDS.BLUE));   Bv = Bv(:);
    Gv = double(A(:,:,BANDS.GREEN1)); Gv = Gv(:);
    Rv = double(A(:,:,BANDS.RED));    Rv = Rv(:);
    lB = log(max(Bv, logfloor));
    lG = log(max(Gv, logfloor));
    lR = log(max(Rv, logfloor));
    ok = isfinite(lB) & isfinite(lG) & isfinite(lR) & dm;
    n  = nnz(ok);
    diag.n_px   = n;
    diag.lR_min = min(lR(ok));
    diag.lR_max = max(lR(ok));
    if n > 10
        p = polyfit(lR(ok), lB(ok), 1);
        k_BR = max(clo, min(chi, p(1)));
        yhat = p(1)*lR(ok) + p(2);
        ss_res = sum((lB(ok)-yhat).^2);
        ss_tot = sum((lB(ok)-mean(lB(ok))).^2);
        diag.r2_BR = 1 - ss_res/max(ss_tot,1e-12);

        p = polyfit(lR(ok), lG(ok), 1);
        k_GR = max(clo, min(chi, p(1)));
        yhat = p(1)*lR(ok) + p(2);
        ss_res = sum((lG(ok)-yhat).^2);
        ss_tot = sum((lG(ok)-mean(lG(ok))).^2);
        diag.r2_GR = 1 - ss_res/max(ss_tot,1e-12);
    end
end

% ------------------------------------------------------------------
function makeQuicklook(outDir, dateStr, A_raw, A_corr, nir_raw, NIR_resid, ...
        DII_BR, DII_GR, depth_grid, reef_mask, mask_BM, reg_mask, ...
        track, T_scale, BANDS)
% 6-panel QC quicklook per date:
%   [1] RGB raw (P2-P98 stretch)
%   [2] RGB corrected (same stretch)
%   [3] NIR raw + BM mask edge
%   [4] NIR histogram over BM
%   [5] DII_BR over reef
%   [6] DII_GR over reef

    fig = figure('Visible','off','Position',[0 0 1800 620]);

    % Shared RGB stretch - per-image P2/P98
    rgb_raw  = cat(3, A_raw(:,:,BANDS.RED),  A_raw(:,:,BANDS.GREEN1),  A_raw(:,:,BANDS.BLUE));
    rgb_corr = cat(3, A_corr(:,:,BANDS.RED), A_corr(:,:,BANDS.GREEN1), A_corr(:,:,BANDS.BLUE));
    lo = double(prctile(rgb_raw(:), 2));
    hi = double(prctile(rgb_raw(:),98));
    if hi <= lo, hi = lo + 0.01; end
    stretch = @(x) min(max((double(x)-lo)./(hi-lo),0),1);

    % Panel 1 - RGB raw
    subplot(1,6,1);
    imshow(stretch(rgb_raw));
    title(sprintf('RGB Raw\n%s', dateStr),'FontSize',8,'Interpreter','none');

    % Panel 2 - RGB corrected
    subplot(1,6,2);
    imshow(stretch(rgb_corr));
    title(sprintf('RGB Corrected\n%s',track),'FontSize',8);

    % Panel 3 - NIR raw + mask overlay
    subplot(1,6,3);
    nir_lo = prctile(nir_raw(:),1); nir_hi = prctile(nir_raw(:),99);
    if nir_hi <= nir_lo, nir_hi = nir_lo+0.01; end
    imagesc(nir_raw,[nir_lo nir_hi]); colormap(gca,gray); axis image off;
    hold on;
    bm_edge = edge(mask_BM,'Canny');
    [re,ce]  = find(bm_edge);
    plot(ce, re, '.r', 'MarkerSize',1);
    rm_edge = edge(reg_mask,'Canny');
    [re2,ce2] = find(rm_edge);
    plot(ce2, re2, '.y','MarkerSize',1);
    hold off;
    title(sprintf('NIR raw\nRed=BM  Yel=reg'),'FontSize',8);
    colorbar('FontSize',7);

    % Panel 4 - NIR histogram over BM mask
    subplot(1,6,4);
    nir_bm_vals = nir_raw(mask_BM);
    histogram(nir_bm_vals, 40, 'FaceColor',[0.2 0.5 0.8]);
    xlabel('NIR','FontSize',7); ylabel('px','FontSize',7);
    title('NIR hist (BM)','FontSize',8);

    % Panel 5 - DII_BR
    subplot(1,6,5);
    d = DII_BR; d(~reef_mask) = NaN;
    dv = d(reef_mask);
    if any(isfinite(dv))
        clim_lo = prctile(dv,2); clim_hi = prctile(dv,98);
        if clim_hi <= clim_lo, clim_hi = clim_lo+0.1; end
    else
        clim_lo=-2; clim_hi=0;
    end
    imagesc(d,[clim_lo clim_hi]); colormap(gca,parula); axis image off; colorbar('FontSize',7);
    title(sprintf('DII\_BR\nmean=%.3f',mean(dv,'omitnan')),'FontSize',8);

    % Panel 6 - DII_GR
    subplot(1,6,6);
    d = DII_GR; d(~reef_mask) = NaN;
    dv = d(reef_mask);
    if any(isfinite(dv))
        clim_lo = prctile(dv,2); clim_hi = prctile(dv,98);
        if clim_hi <= clim_lo, clim_hi = clim_lo+0.1; end
    else
        clim_lo=-2; clim_hi=0;
    end
    imagesc(d,[clim_lo clim_hi]); colormap(gca,parula); axis image off; colorbar('FontSize',7);
    title(sprintf('DII\_GR\nmean=%.3f',mean(dv,'omitnan')),'FontSize',8);

    sgtitle(sprintf('Phase 5 QC - %s  (Track %s)', dateStr, track),...
        'FontSize',10,'FontWeight','bold','Interpreter','none');

    exportgraphics(fig, fullfile(outDir,['quicklook_DII_' dateStr '.png']),'Resolution',120);
    close(fig);
end

% ------------------------------------------------------------------
function cmap = bwr_colormap()
% Blue-White-Red diverging colormap (64 levels)
    n=32;
    r=[linspace(0,1,n), ones(1,n)];
    g=[linspace(0,1,n), linspace(1,0,n)];
    b=[ones(1,n), linspace(1,0,n)];
    cmap=[r',g',b'];
end

% ------------------------------------------------------------------
function [contrast,entropy]=computeGLCM(img,win)
    img_norm = mat2gray(img);
    img_u8   = im2uint8(img_norm);
    se       = true(win,win);
    contrast = single(stdfilt(img_norm,se).^2);
    entropy  = single(entropyfilt(img_u8,se));
end

% ------------------------------------------------------------------
function dcfg=getDatConfig(DATE_CFG,dateStr)
    dcfg=struct('hedley','standard','baseline_coral',false,...
                'standalone',true,'dii_k','per_date');
    for k=1:numel(DATE_CFG)
        if strcmp(DATE_CFG(k).date,dateStr), dcfg=DATE_CFG(k); return; end
    end
    warning('Date %s not in DATE_CFG - using defaults.',dateStr);
end

% ------------------------------------------------------------------
function depth_out=reprojectDepth(Ab,Rab,Rimg,H,W)
    depth_out=nan(H,W,'single');
    try
        [Hb,Wb]=size(Ab);
        [cb,rb]=meshgrid(1:Wb,1:Hb);
        [lat_b,lon_b]=intrinsicToGeographic(Rab,cb,rb);
        [x_b,y_b]=latlon2utm47(lat_b(:),lon_b(:));
        [ci,ri]=meshgrid(1:W,1:H);
        [x_q,y_q]=intrinsicToWorld(Rimg,ci,ri);
        F=scatteredInterpolant(double(x_b),double(y_b),double(Ab(:)),'linear','none');
        depth_out=single(F(double(x_q),double(y_q)));
        cov=100*nnz(isfinite(depth_out))/(H*W);
        fprintf('  Depth reprojected: %.1f%% coverage\n',cov);
        if cov < 50
            fprintf('  WARNING: depth coverage < 50%% - check bathymetry_0.tif\n');
        end
    catch ME
        warning('reprojectDepth (Method C): %s',ME.message);
    end
end

% ------------------------------------------------------------------
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

% ------------------------------------------------------------------
function mask=rasterizeGeoJSON(path,H,W,R,label)
    mask=false(H,W);
    if ~exist(path,'file'), warning('Not found: %s',path); return; end
    try
        data=jsondecode(fileread(path)); feats=data.features;
        xmin=R.XWorldLimits(1); ymax=R.YWorldLimits(2);
        cw=(R.XWorldLimits(2)-xmin)/W; ch=(ymax-R.YWorldLimits(1))/H;
        n_rings=0;
        for fi=1:numel(feats)
            rings=extractRings(feats(fi).geometry);
            for ri=1:numel(rings)
                ring=rings{ri};
                if size(ring,1)<3, continue; end
                [xv,yv]=latlon2utm47(ring(:,2),ring(:,1));
                cv=max(1,min(W,round((xv-xmin)/cw+0.5)));
                rv=max(1,min(H,round((ymax-yv)/ch+0.5)));
                try
                    mask=mask|poly2mask(double(cv),double(rv),H,W);
                    n_rings=n_rings+1;
                catch
                end
            end
        end
        fprintf('  %s: %d features / %d rings -> %d px\n',...
            label,numel(feats),n_rings,nnz(mask));
    catch ME
        warning('rasterizeGeoJSON (%s): %s',label,ME.message);
    end
end

% ------------------------------------------------------------------
function rings=extractRings(geom)
    rings={};
    try
        if strcmp(geom.type,'Polygon')
            c=geom.coordinates;
            if isnumeric(c), rings{end+1}=c;
            elseif iscell(c)
                r=c{1}; if iscell(r),r=cell2mat(r);end; rings{end+1}=r;
            end
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
    catch
    end
end

% ------------------------------------------------------------------
function saveMask(path,mask,R)
    try geotiffwrite(path,uint8(mask),R); return; catch, end
    try geotiffwrite(path,uint8(mask),R,'CoordRefSysCode',32647); return; catch, end
    try
        t=Tiff(path,'w'); [H,W]=size(mask);
        t.setTag('Photometric',1);
        t.setTag('PlanarConfig',1);
        t.setTag('ImageLength',H);
        t.setTag('ImageWidth',W);
        t.setTag('BitsPerSample',8);
        t.setTag('SampleFormat',1);
        t.setTag('SamplesPerPixel',1);
        t.setTag('Compression',5);
        t.write(uint8(mask)); t.close();
    catch ME, warning('saveMask: %s',ME.message); end
end

% ------------------------------------------------------------------
function saveMultiBandTif(path,A,R)
    if ~isa(A,'single'), A=single(A); end
    try geotiffwrite(path,A,R); return; catch, end
    try geotiffwrite(path,A,R,'CoordRefSysCode',32647); return; catch, end
    try
        [H,W,nb]=size(A); t=Tiff(path,'w');
        t.setTag('Photometric',1);
        t.setTag('PlanarConfig',1);
        t.setTag('ImageLength',H);
        t.setTag('ImageWidth',W);
        t.setTag('BitsPerSample',32);
        t.setTag('SampleFormat',3);
        t.setTag('SamplesPerPixel',nb);
        t.setTag('Compression',5);
        t.write(A); t.close();
    catch ME, error('saveMultiBandTif: %s',ME.message); end
end

% ------------------------------------------------------------------
function logWrite(fid,varargin)
    msg=sprintf(varargin{:}); fprintf('%s',msg);
    if fid>0, fprintf(fid,'%s',msg); end
end

% ------------------------------------------------------------------
function out = ternary(cond, a, b)
    if cond, out=a; else, out=b; end
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

