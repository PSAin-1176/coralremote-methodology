%% Preprocess_Masking.m
clear; clc;
global PIPELINE_LOG_FID PIPELINE_CONFIG_PATH;
if isempty(PIPELINE_LOG_FID), PIPELINE_LOG_FID = -1; end
cfg = load(PIPELINE_CONFIG_PATH);

runTimestamp = char(datetime('now','TimeZone','local'),'yyyyMMdd_HHmmss');
outDir = fullfile(cfg.OUT_DIR, 'Phase2_Masking', ['run_' runTimestamp]);
if ~exist(outDir, 'dir'), mkdir(outDir); end
logWrite(PIPELINE_LOG_FID, 'Preprocess_Masking.m - %s\n', runTimestamp);

inputDir = cfg.RAW_DIR;
sourceTable = table();
prepared_source_dir = '';
prepared_mode = false;

prepOutPath = '';
if isfield(cfg, 'USE_PREPARED_SOURCE') && cfg.USE_PREPARED_SOURCE && isfield(cfg, 'PREP_SOURCE_ROOT')
    prepOutPath = fullfile(cfg.PREP_SOURCE_ROOT, 'prepared_source_output.mat');
end
if ~isempty(prepOutPath) && exist(prepOutPath, 'file')
    prep = load(prepOutPath);
    manifestPath = fullfile(prep.outDir, 'prepared_source_manifest.mat');
    if exist(manifestPath, 'file')
        man = load(manifestPath);
        sourceTable = man.preparedTable;
        prepared_source_dir = prep.outDir;
        prepared_mode = height(sourceTable) > 0;
        if prepared_mode
            inputDir = fileparts(sourceTable.prepared_sr_path{1});
        end
    end
end

coregBase = fullfile(cfg.OUT_DIR, 'Phase1_Coreg');
if ~prepared_mode && exist(coregBase, 'dir')
    runs = dir(fullfile(coregBase, 'run_*'));
    runs = runs([runs.isdir]);
    if ~isempty(runs)
        [~, idx] = sort({runs.name});
        inputDir = fullfile(coregBase, runs(idx(end)).name);
    end
end

ff = dir(fullfile(inputDir, cfg.COREG_PATTERN));
if isempty(ff), ff = dir(fullfile(inputDir, cfg.RAW_PATTERN)); end
if isempty(ff), error('No input imagery found in %s', inputDir); end

refPath = fullfile(ff(1).folder, ff(1).name);
for i = 1:numel(ff)
    tok = regexp(ff(i).name, '(\d{8})', 'tokens', 'once');
    if ~isempty(tok) && strcmp(tok{1}, cfg.GEOM_REF_DATE)
        refPath = fullfile(ff(i).folder, ff(i).name);
        break;
    end
end

[Aref, Rref] = readgeoraster(refPath); %#ok<ASGLU>
[H, W, ~] = size(Aref);

reef_mask = rasterizeGeoJSONClass(cfg.REEF_MASK_PATH, H, W, Rref, 'Reef');
mdeep_mask = rasterizeGeoJSON(cfg.DEEP_WATER_MASK, H, W, Rref, 'Deep water');
sand_mask = rasterizeGeoJSON(cfg.MASTER_SAND_MASK, H, W, Rref, 'Master sand');
stable_sand_mask = rasterizeGeoJSON(cfg.STABLE_SAND_MASK, H, W, Rref, 'Stable sand');
valid_mask = all(single(Aref) > 0, 3);
if prepared_mode
    refDateIdx = find(strcmp(sourceTable.date, cfg.GEOM_REF_DATE), 1);
    if ~isempty(refDateIdx) && exist(sourceTable.prepared_valid_mask_path{refDateIdx}, 'file')
        valid_mask = logical(readgeoraster(sourceTable.prepared_valid_mask_path{refDateIdx}));
    end
end

bathy_deep_mask = false(H, W);
depth_grid = nan(H, W, 'single');
if isfield(cfg, 'BATHY_PATH') && exist(cfg.BATHY_PATH, 'file')
    try
        [Ab, Rab] = readgeoraster(cfg.BATHY_PATH);
        depth_grid = reprojectDepth(single(Ab), Rab, Rref, H, W);
        bathy_deep_mask = isfinite(depth_grid) & (depth_grid > cfg.LR.SUNGLINT_DEEP_THRESH_CM);
    catch ME
        warning('Phase2 bathymetry reprojection failed: %s', ME.message);
    end
end
bm_deep_mask = bathy_deep_mask | mdeep_mask;
deep_mask = bm_deep_mask;

saveMask(fullfile(outDir, 'reef_mask.tif'), reef_mask, Rref);
saveMask(fullfile(outDir, 'deep_mask.tif'), deep_mask, Rref);
saveMask(fullfile(outDir, 'mdeep_mask.tif'), mdeep_mask, Rref);
saveMask(fullfile(outDir, 'bathy_deep_mask.tif'), bathy_deep_mask, Rref);
saveMask(fullfile(outDir, 'bm_deep_mask.tif'), bm_deep_mask, Rref);
saveMask(fullfile(outDir, 'sand_mask.tif'), sand_mask, Rref);
saveMask(fullfile(outDir, 'stable_sand_mask.tif'), stable_sand_mask, Rref);
saveMask(fullfile(outDir, 'valid_mask.tif'), valid_mask, Rref);

save(fullfile(outDir, 'mask_output.mat'), ...
    'outDir','inputDir','refPath','Rref','H','W', ...
    'reef_mask','deep_mask','mdeep_mask','bathy_deep_mask','bm_deep_mask','depth_grid', ...
    'sand_mask','stable_sand_mask','valid_mask', ...
    'sourceTable','prepared_source_dir','prepared_mode');
save(fullfile(cfg.OUT_DIR, 'phase2_masking_output.mat'), ...
    'outDir','inputDir','refPath','Rref','H','W', ...
    'reef_mask','deep_mask','mdeep_mask','bathy_deep_mask','bm_deep_mask','depth_grid', ...
    'sand_mask','stable_sand_mask','valid_mask', ...
    'sourceTable','prepared_source_dir','prepared_mode');

logWrite(PIPELINE_LOG_FID, '  Input dir : %s\n', inputDir);
if prepared_mode
    logWrite(PIPELINE_LOG_FID, '  Prepared source manifest : %s\n', prepared_source_dir);
end
logWrite(PIPELINE_LOG_FID, '  Reference : %s\n', refPath);
logWrite(PIPELINE_LOG_FID, '  Reef px   : %d\n', nnz(reef_mask));
logWrite(PIPELINE_LOG_FID, '  MDeep px  : %d\n', nnz(mdeep_mask));
logWrite(PIPELINE_LOG_FID, '  Bathy px  : %d\n', nnz(bathy_deep_mask));
logWrite(PIPELINE_LOG_FID, '  BM deep px: %d\n', nnz(bm_deep_mask));
logWrite(PIPELINE_LOG_FID, '  Sand px   : %d\n', nnz(sand_mask));
logWrite(PIPELINE_LOG_FID, '  Output    : %s\n\n', outDir);

function depth_out = reprojectDepth(Ab, Rab, Rimg, H, W)
depth_out = nan(H, W, 'single');
try
    [Hb, Wb] = size(Ab);
    [cb, rb] = meshgrid(1:Wb, 1:Hb);
    [lat_b, lon_b] = intrinsicToGeographic(Rab, cb, rb);
    [x_b, y_b] = latlon2utm47(lat_b(:), lon_b(:));
    [ci, ri] = meshgrid(1:W, 1:H);
    [x_q, y_q] = intrinsicToWorld(Rimg, ci, ri);
    F = scatteredInterpolant(double(x_b), double(y_b), double(Ab(:)), 'linear', 'none');
    depth_out = single(F(double(x_q), double(y_q)));
catch ME
    warning('reprojectDepth: %s', ME.message);
end
end
