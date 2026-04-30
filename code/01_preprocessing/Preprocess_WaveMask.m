%% Preprocess_WaveMask.m
clear; clc;
global PIPELINE_LOG_FID PIPELINE_CONFIG_PATH;
if isempty(PIPELINE_LOG_FID), PIPELINE_LOG_FID = -1; end
cfg = load(PIPELINE_CONFIG_PATH);

runTimestamp = char(datetime('now','TimeZone','local'),'yyyyMMdd_HHmmss');
outDir = fullfile(cfg.OUT_DIR, 'Phase4_Wave', ['run_' runTimestamp]);
if ~exist(outDir, 'dir'), mkdir(outDir); end
logWrite(PIPELINE_LOG_FID, 'Preprocess_WaveMask.m - %s\n', runTimestamp);

p2 = load(fullfile(cfg.OUT_DIR, 'phase2_masking_output.mat'));
p3 = load(fullfile(cfg.OUT_DIR, 'phase3_logratio_output.mat'));
inputDir = p2.inputDir;
normDir = p3.outDir;
sourceTable = table();
if isfield(p2, 'sourceTable')
    sourceTable = p2.sourceTable;
end

ff = dir(fullfile(normDir, 'NORM_*.tif'));
if isempty(ff), error('No normalized rasters found in %s', normDir); end
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
            logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase4 (thread pool incompatible with geospatial IO)\n');
        end
    catch
        usePar = false;
        logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase4 (could not verify pool type safely)\n');
    end
end
if usePar
    parfor i = 1:numel(ff)
        msgs{i} = processWaveDate(i, ff, normDir, inputDir, sourceTable, outDir, p2.reef_mask, cfg);
    end
else
    for i = 1:numel(ff)
        msgs{i} = processWaveDate(i, ff, normDir, inputDir, sourceTable, outDir, p2.reef_mask, cfg);
    end
end
for i = 1:numel(msgs)
    if ~isempty(msgs{i})
        logWrite(PIPELINE_LOG_FID, '%s', msgs{i});
    end
end

save(fullfile(outDir, 'wave_output.mat'), 'outDir','inputDir','normDir');
save(fullfile(cfg.OUT_DIR, 'phase4_wave_output.mat'), 'outDir','inputDir','normDir');
logWrite(PIPELINE_LOG_FID, '  Output       : %s\n\n', outDir);

function msg = processWaveDate(i, ff, normDir, inputDir, sourceTable, outDir, reefMask, cfg)
msg = '';
tok = regexp(ff(i).name, 'NORM_(\d{8})\.tif', 'tokens', 'once');
if isempty(tok), return; end
d = tok{1};
[Anorm, Rref] = readgeoraster(fullfile(normDir, ff(i).name));
Anorm = single(Anorm);
rawPath = '';
validMask = [];
if ~isempty(sourceTable)
    idx = find(strcmp(sourceTable.date, d), 1);
    if ~isempty(idx)
        rawPath = sourceTable.prepared_sr_path{idx};
        if exist(sourceTable.prepared_valid_mask_path{idx}, 'file')
            validMask = logical(readgeoraster(sourceTable.prepared_valid_mask_path{idx}));
        end
    end
else
    hits = dir(fullfile(inputDir, ['*' d '*.tif']));
    if ~isempty(hits)
        rawPath = fullfile(hits(1).folder, hits(1).name);
    end
end
if isempty(rawPath)
    error('No source raster found for %s in %s', d, inputDir);
end
[Araw, ~] = readgeoraster(rawPath);
Araw = single(Araw) ./ cfg.SR_SCALE;
if isempty(validMask)
    validMask = all(isfinite(Araw) & Araw > 0, 3);
end
waveThr = cfg.LR.NIR_WAVE_THRESHOLD;
if isfield(cfg.LR, 'WAVE_THRESHOLD_OVERRIDE_DATES') && isfield(cfg.LR, 'WAVE_THRESHOLD_OVERRIDE_VALUES')
    idxThr = find(strcmp(cfg.LR.WAVE_THRESHOLD_OVERRIDE_DATES, d), 1);
    if ~isempty(idxThr) && idxThr <= numel(cfg.LR.WAVE_THRESHOLD_OVERRIDE_VALUES)
        waveThr = cfg.LR.WAVE_THRESHOLD_OVERRIDE_VALUES(idxThr);
    end
end
disableWave = isfield(cfg.LR, 'WAVE_DISABLE_GLOBAL') && cfg.LR.WAVE_DISABLE_GLOBAL;
if ~disableWave && isfield(cfg.LR, 'WAVE_DISABLE_DATES')
    disableWave = any(strcmp(cfg.LR.WAVE_DISABLE_DATES, d));
end
flagOnly = false;
if isfield(cfg.LR, 'WAVE_FLAG_ONLY_DATES')
    flagOnly = any(strcmp(cfg.LR.WAVE_FLAG_ONLY_DATES, d));
end
wave_mask = Araw(:,:,cfg.BANDS.NIR) > waveThr;
wave_mask = wave_mask & reefMask & validMask;
if disableWave
    wave_mask(:) = false;
end
Aclean = Anorm;
Aclean(repmat(~validMask, [1 1 size(Aclean, 3)])) = nan;
if ~flagOnly
    Aclean(repmat(wave_mask, [1 1 size(Aclean, 3)])) = nan;
end
saveMask(fullfile(outDir, ['WAVE_MASK_' d '.tif']), wave_mask, Rref);
saveMultiBandTif(fullfile(outDir, ['CLEAN_NORM_' d '.tif']), Aclean, Rref);
modeTag = '';
if disableWave
    modeTag = ' [disabled]';
elseif flagOnly
    modeTag = ' [flag-only]';
elseif abs(waveThr - cfg.LR.NIR_WAVE_THRESHOLD) > eps('single')
    modeTag = sprintf(' [thr=%.3f]', waveThr);
end
msg = sprintf('  Wave mask %s : %d px%s\n', d, nnz(wave_mask), modeTag);
end
