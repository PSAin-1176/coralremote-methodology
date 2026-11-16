%% Preprocess_Normalize.m
clear; clc;
global PIPELINE_LOG_FID PIPELINE_CONFIG_PATH;
if isempty(PIPELINE_LOG_FID), PIPELINE_LOG_FID = -1; end
cfg = load(PIPELINE_CONFIG_PATH);

runTimestamp = char(datetime('now','TimeZone','local'),'yyyyMMdd_HHmmss');
outDir = fullfile(cfg.OUT_DIR, 'Phase3_LogRatio', ['run_' runTimestamp]);
if ~exist(outDir, 'dir'), mkdir(outDir); end
logWrite(PIPELINE_LOG_FID, 'Preprocess_Normalize.m - %s\n', runTimestamp);

m2 = load(fullfile(cfg.OUT_DIR, 'phase2_masking_output.mat'));
inputDir = m2.inputDir;
reef_mask = m2.reef_mask;
Aref = [];
Rref = [];
sourceTable = table();
if isfield(m2, 'sourceTable')
    sourceTable = m2.sourceTable;
end
if ~isempty(sourceTable)
    sourceTable = sortrows(sourceTable, 'date');
    refIdx = find(strcmp(sourceTable.date, cfg.RAD_REF_DATE), 1);
    if isempty(refIdx)
        error('Radiometric reference %s not found in prepared manifest.', cfg.RAD_REF_DATE);
    end
    [Aref, Rref] = readgeoraster(sourceTable.prepared_sr_path{refIdx});
    dateList = sourceTable.date;
    inputPaths = sourceTable.prepared_sr_path;
    validMaskPaths = sourceTable.prepared_valid_mask_path;
else
    ff = dir(fullfile(inputDir, cfg.COREG_PATTERN));
    if isempty(ff), ff = dir(fullfile(inputDir, cfg.RAW_PATTERN)); end
    if isempty(ff), error('No imagery available for normalization.'); end
    dateList = cell(numel(ff), 1);
    for i = 1:numel(ff)
        tok = regexp(ff(i).name, '(\d{8})', 'tokens', 'once');
        if isempty(tok), dateList{i} = ''; else, dateList{i} = tok{1}; end
    end
    [~, ord] = sort(dateList);
    ff = ff(ord);
    dateList = dateList(ord);
    refIdx = find(strcmp(dateList, cfg.RAD_REF_DATE), 1);
    if isempty(refIdx)
        error('Radiometric reference %s not found in %s', cfg.RAD_REF_DATE, inputDir);
    end
    [Aref, Rref] = readgeoraster(fullfile(ff(refIdx).folder, ff(refIdx).name));
    inputPaths = cell(numel(dateList),1);
    validMaskPaths = cell(numel(dateList),1);
    for i = 1:numel(dateList)
        inputPaths{i} = fullfile(ff(i).folder, ff(i).name);
        validMaskPaths{i} = '';
    end
end
Aref = single(Aref) ./ cfg.SR_SCALE;
valid_ref = reef_mask & all(isfinite(Aref) & Aref > 0, 3);
if nnz(valid_ref) == 0
    error('No valid reef pixels found for reference statistics.');
end

nBands = size(Aref, 3);
applyGlint = isfield(cfg.LR, 'APPLY_SUNGLINT_CORRECTION') && cfg.LR.APPLY_SUNGLINT_CORRECTION;
applyNorm = ~(isfield(cfg.LR, 'APPLY_NORMALIZATION') && ~cfg.LR.APPLY_NORMALIZATION);
mdeep_mask = false(size(reef_mask));
if isfield(m2, 'mdeep_mask') && ~isempty(m2.mdeep_mask)
    mdeep_mask = m2.mdeep_mask;
end
bm_deep_mask = m2.deep_mask;
if isfield(m2, 'bm_deep_mask') && ~isempty(m2.bm_deep_mask)
    bm_deep_mask = m2.bm_deep_mask;
end
Aref_mu = Aref;
if applyGlint
    if ~isempty(sourceTable)
        refValidMask = valid_ref;
        if exist(sourceTable.prepared_valid_mask_path{refIdx}, 'file')
            refValidMask = logical(readgeoraster(sourceTable.prepared_valid_mask_path{refIdx}));
        end
    else
        refValidMask = all(isfinite(Aref) & Aref > 0, 3);
    end
    [Aref_mu, gsRef] = maybeApplySunglintCorrection(Aref, refValidMask, reef_mask, mdeep_mask, bm_deep_mask, cfg, cfg.RAD_REF_DATE); %#ok<NASGU>
end
mu_ref = nan(1, nBands, 'single');
for b = 1:nBands
    band = Aref_mu(:,:,b);
    mu_ref(b) = mean(band(valid_ref), 'omitnan');
end
glintStats = struct('date', {}, 'applied', {}, 'roi_name', {}, 'reg_px', {}, 'nir_min', {}, ...
    'nir_p95', {}, 'blue_slope', {}, 'green1_slope', {}, 'red_slope', {}, 'reason', {});
msgs = cell(numel(dateList),1);
glintStatsCell = cell(numel(dateList),1);
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
            logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase3 (thread pool incompatible with projcrs)\n');
        end
    catch
        usePar = false;
        logWrite(PIPELINE_LOG_FID, '  Parallel note : serial fallback in Phase3 (could not verify pool type safely)\n');
    end
end
if usePar
    parfor i = 1:numel(dateList)
        [msgs{i}, glintStatsCell{i}] = processOneDate(i, dateList, inputPaths, validMaskPaths, ...
            sourceTable, outDir, Rref, mu_ref, nBands, applyGlint, reef_mask, mdeep_mask, bm_deep_mask, cfg);
    end
else
    for i = 1:numel(dateList)
        [msgs{i}, glintStatsCell{i}] = processOneDate(i, dateList, inputPaths, validMaskPaths, ...
            sourceTable, outDir, Rref, mu_ref, nBands, applyGlint, reef_mask, mdeep_mask, bm_deep_mask, cfg);
    end
end
for i = 1:numel(msgs)
    if ~isempty(msgs{i})
        logWrite(PIPELINE_LOG_FID, '%s', msgs{i});
    end
    if ~isempty(glintStatsCell{i})
        glintStats(end+1) = glintStatsCell{i}; %#ok<SAGROW>
    end
end

save(fullfile(outDir, 'logratio_reference.mat'), ...
    'outDir','inputDir','mu_ref','Rref','reef_mask','sourceTable','glintStats','applyNorm');
save(fullfile(cfg.OUT_DIR, 'phase3_logratio_output.mat'), ...
    'outDir','inputDir','mu_ref','Rref','reef_mask','sourceTable','glintStats','applyNorm');

fid = fopen(fullfile(outDir, 'logratio_report.txt'), 'w');
fprintf(fid, 'Radiometric reference date: %s\n', cfg.RAD_REF_DATE);
fprintf(fid, 'Normalization: %s\n', ternary(applyNorm, 'enabled', 'disabled'));
if applyNorm
    fprintf(fid, 'Mandatory formula: log(R) - log(mu_ref)\n');
    fprintf(fid, 'Stored linear product: R ./ mu_ref\n');
else
    fprintf(fid, 'Stored linear product: glint-corrected reflectance (no cross-date normalization)\n');
    fprintf(fid, 'Stored logratio product: not used; written as NaN placeholder for interface consistency\n');
end
fprintf(fid, 'Sunglint correction: %s\n', ternary(applyGlint, 'enabled', 'disabled'));
if applyGlint
    fprintf(fid, '  Method: Hedley regression on raw NIR over constrained deep-water ROI\n');
    fprintf(fid, '  Primary ROI : %s\n', cfg.LR.SUNGLINT_PRIMARY_ROI);
    fprintf(fid, '  Fallback ROI: %s\n', cfg.LR.SUNGLINT_FALLBACK_ROI);
    fprintf(fid, '  BM mask: bathymetry > %d cm OR MDeepWater.geojson\n', cfg.LR.SUNGLINT_DEEP_THRESH_CM);
    fprintf(fid, '  Min regression pixels: %d\n', cfg.LR.SUNGLINT_MIN_REG_PX);
    fprintf(fid, '  NIR P95 fraction: %.2f\n', cfg.LR.SUNGLINT_NIR_P95_FRAC);
end
for b = 1:numel(mu_ref)
    fprintf(fid, 'Band %d mean: %.6f\n', b, mu_ref(b));
end
if applyGlint && ~isempty(glintStats)
    fprintf(fid, '\nPer-date sunglint summary:\n');
    fprintf(fid, '%-10s  %-6s  %-16s  %8s  %10s  %10s  %10s\n', ...
        'Date','Applied','ROI','RegPx','BlueSlope','G1Slope','RedSlope');
    for k = 1:numel(glintStats)
        gs = glintStats(k);
        fprintf(fid, '%-10s  %-6s  %-16s  %8d  %10.5f  %10.5f  %10.5f', ...
            gs.date, ternary(gs.applied, 'yes', 'no'), gs.roi_name, gs.reg_px, ...
            gs.blue_slope, gs.green1_slope, gs.red_slope);
        if ~isempty(gs.reason)
            fprintf(fid, '   %s', gs.reason);
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);

logWrite(PIPELINE_LOG_FID, '  Reference date : %s\n', cfg.RAD_REF_DATE);
logWrite(PIPELINE_LOG_FID, '  Normalize      : %s\n', ternary(applyNorm, 'on', 'off'));
if applyGlint
    logWrite(PIPELINE_LOG_FID, '  Sunglint      : Hedley on constrained deep-water ROI\n');
end
logWrite(PIPELINE_LOG_FID, '  Output         : %s\n\n', outDir);

function [A_corr, gs] = maybeApplySunglintCorrection(A, validMask, reefMask, mdeepMask, bmDeepMask, cfg, dateStr)
A_corr = A;
gs = struct('date', '', 'applied', false, 'roi_name', '', 'reg_px', 0, 'nir_min', nan, ...
    'nir_p95', nan, 'blue_slope', nan, 'green1_slope', nan, 'red_slope', nan, 'reason', '');
if ~(isfield(cfg.LR, 'APPLY_SUNGLINT_CORRECTION') && cfg.LR.APPLY_SUNGLINT_CORRECTION)
    gs.reason = 'disabled';
    return;
end
roiPrimary = logical(mdeepMask) & logical(validMask);
roiFallback = logical(bmDeepMask) & ~logical(reefMask) & logical(validMask);
regMask = roiPrimary;
gs.roi_name = 'MDEEP';
if nnz(regMask) < cfg.LR.SUNGLINT_MIN_REG_PX
    regMask = roiFallback;
    gs.roi_name = 'BM_OUTSIDE_REEF';
end
nir = double(A(:,:,cfg.BANDS.NIR));
nirVals = nir(regMask);
nirVals = nirVals(isfinite(nirVals) & nirVals > 0);
if numel(nirVals) < cfg.LR.SUNGLINT_MIN_REG_PX
    gs.reason = sprintf('few ROI px (%d)', numel(nirVals));
    return;
end
nirP95 = prctile(nirVals, cfg.LR.SUNGLINT_NIR_P95_FRAC * 100);
regMask = regMask & isfinite(nir) & nir < nirP95;
nirVec = nir(regMask);
if numel(nirVec) < cfg.LR.SUNGLINT_MIN_REG_PX
    gs.reason = sprintf('few ROI px after P95 (%d)', numel(nirVec));
    return;
end
nirMin = min(nirVec);
nirRange = max(nirVec) - nirMin;
if nirRange < 1e-5
    gs.reason = 'degenerate NIR range';
    return;
end
visBands = [cfg.BANDS.COASTAL cfg.BANDS.BLUE cfg.BANDS.GREEN1 cfg.BANDS.GREEN2 ...
    cfg.BANDS.YELLOW cfg.BANDS.RED cfg.BANDS.REDEDGE];
slopes = nan(1, numel(visBands));
for bi = 1:numel(visBands)
    b = visBands(bi);
    band = double(A(:,:,b));
    y = band(regMask);
    if numel(y) < cfg.LR.SUNGLINT_MIN_REG_PX || std(y, 'omitnan') < 1e-8
        continue;
    end
    p = polyfit(nirVec, y, 1);
    slope = p(1);
    if ~isfinite(slope) || slope < 0 || slope > 8
        continue;
    end
    slopes(bi) = slope;
    corr = single(A(:,:,b)) - single(slope) .* (single(A(:,:,cfg.BANDS.NIR)) - single(nirMin));
    A_corr(:,:,b) = max(0, corr);
end
gs.applied = any(isfinite(slopes));
gs.reg_px = numel(nirVec);
gs.nir_min = nirMin;
gs.nir_p95 = nirP95;
gs.blue_slope = getBandSlope(visBands, slopes, cfg.BANDS.BLUE);
gs.green1_slope = getBandSlope(visBands, slopes, cfg.BANDS.GREEN1);
gs.red_slope = getBandSlope(visBands, slopes, cfg.BANDS.RED);
if ~gs.applied
    gs.reason = 'no valid slopes';
end
end

function s = getBandSlope(visBands, slopes, bandId)
idx = find(visBands == bandId, 1);
if isempty(idx)
    s = nan;
else
    s = slopes(idx);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [msg, gs] = processOneDate(i, dateList, inputPaths, validMaskPaths, ...
    sourceTable, outDir, Rref, mu_ref, nBands, applyGlint, reefMask, mdeepMask, bmDeepMask, cfg)
msg = '';
gs = [];
if isempty(dateList{i})
    return;
end
inPath = inputPaths{i};
[A, ~] = readgeoraster(inPath);
A = single(A) ./ cfg.SR_SCALE;
valid = all(isfinite(A) & A > 0, 3);
if ~isempty(sourceTable) && i <= numel(validMaskPaths) && ~isempty(validMaskPaths{i}) && exist(validMaskPaths{i}, 'file')
    valid = valid & logical(readgeoraster(validMaskPaths{i}));
end
[A_corr, gs] = maybeApplySunglintCorrection(A, valid, reefMask, mdeepMask, bmDeepMask, cfg, dateList{i});
gs.date = dateList{i};
applyNorm = ~(isfield(cfg.LR, 'APPLY_NORMALIZATION') && ~cfg.LR.APPLY_NORMALIZATION);
baseA = A;
if applyGlint
    baseA = A_corr;
end
if applyNorm
    A_norm = baseA ./ reshape(mu_ref, 1, 1, []);
    A_logratio = log(baseA + cfg.LR.EPS) - log(reshape(mu_ref, 1, 1, []));
else
    A_norm = baseA;
    A_logratio = nan(size(baseA), 'single');
end
A_norm(~repmat(valid, [1 1 nBands])) = nan;
A_logratio(~repmat(valid, [1 1 nBands])) = nan;
saveMultiBandTif(fullfile(outDir, ['NORM_' dateList{i} '.tif']), A_norm, Rref);
saveMultiBandTif(fullfile(outDir, ['LOGRATIO_' dateList{i} '.tif']), A_logratio, Rref);
if applyGlint
    saveMultiBandTif(fullfile(outDir, ['GLINT_' dateList{i} '.tif']), A_corr, Rref);
    msg = sprintf('  %s %s  [glint %s, roi=%s, reg=%d]\n', ...
        ternary(applyNorm, 'Normalized', 'Prepared'), dateList{i}, ternary(gs.applied, 'on', 'skip'), gs.roi_name, gs.reg_px);
else
    msg = sprintf('  %s %s\n', ternary(applyNorm, 'Normalized', 'Prepared'), dateList{i});
end
end
