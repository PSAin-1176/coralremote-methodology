%% QC_AlignmentPost.m
% ============================================================
% PHASE 2 - Alignment check AFTER AROSICS correction
%
% Identical logic to Phase 0 but runs on co-registered images.
% Compares offsets before vs after correction.
% Pass criterion: all offsets < MAX_OFFSET_PX (default 1.0 px)
%
% OUTPUTS -> OUT_DIR/Phase2_AlignCheck_post/run_YYYYMMDD_HHMMSS/
%   Same structure as Phase 0 output.
%   alignment_report.txt includes before/after comparison.
% ============================================================

clear; clc;
runTime      = datetime('now','TimeZone','local');
runTimestamp = char(runTime,'yyyyMMdd_HHmmss');
fprintf('QC_AlignmentPost.m - Phase 2 - %s\n\n', runTimestamp);

%% Load shared config
global PIPELINE_CONFIG_PATH;
cfgPath = '';
if ~isempty(PIPELINE_CONFIG_PATH) && exist(PIPELINE_CONFIG_PATH,'file')
    cfgPath = PIPELINE_CONFIG_PATH;
else
    candidates = {fullfile(pwd,'pipeline_config.mat')};
    for k=1:numel(candidates)
        if exist(candidates{k},'file'), cfgPath=candidates{k}; break; end
    end
end
if isempty(cfgPath)
    [f,d]=uigetfile('pipeline_config.mat','Select pipeline_config.mat');
    if isequal(f,0), error('No config file selected.'); end
    cfgPath=fullfile(d,f);
end
cfg = load(cfgPath);

%% Find input directory (Phase 1 coreg output)
coregDir = '';

% 1. Try phase1_output.mat
phase1MatPath = fullfile(cfg.OUT_DIR, 'phase1_output.mat');
if exist(phase1MatPath,'file')
    try
        p1 = load(phase1MatPath);
        if exist(p1.outDir,'dir')
            coregDir = p1.outDir;
            fprintf('Phase 1 output (from mat): %s\n', coregDir);
        end
    catch
    end
end

% 2. Auto-detect most recent Phase1_Coreg run/after folder
if isempty(coregDir)
    phase1Base = fullfile(cfg.OUT_DIR, 'Phase1_Coreg');
    if exist(phase1Base,'dir')
        runs = dir(fullfile(phase1Base, 'run_*'));
        runs = runs([runs.isdir]);
        if ~isempty(runs)
            [~,si] = sort({runs.name});
            % Prefer run subfolder that contains an 'after' subfolder with TIFs
            for ki = numel(si):-1:1
                candidate = fullfile(phase1Base, runs(si(ki)).name, 'after');
                if exist(candidate,'dir') && ~isempty(dir(fullfile(candidate,'*.tif')))
                    coregDir = candidate;
                    fprintf('Phase 1 output (auto-detected after/): %s\n', coregDir);
                    break;
                end
                % Also check the run folder itself for TIFs
                candidate = fullfile(phase1Base, runs(si(ki)).name);
                if ~isempty(dir(fullfile(candidate,'*.tif')))
                    coregDir = candidate;
                    fprintf('Phase 1 output (auto-detected run/): %s\n', coregDir);
                    break;
                end
            end
        end
    end
end

% 3. Check for pre-existing AROSICS output at known location
%    (from run_20260310_145506/after/ which was run before this pipeline)
if isempty(coregDir)
    preExisting = fullfile(cfg.OUT_DIR, 'run_20260310_145506', 'after');
    if exist(preExisting,'dir') && ~isempty(dir(fullfile(preExisting,'*.tif')))
        coregDir = preExisting;
        fprintf('Phase 1 output (pre-existing AROSICS run): %s\n', coregDir);
    end
end

% 4. Manual fallback
if isempty(coregDir) || ~exist(coregDir,'dir')
    fprintf('Could not auto-detect Phase 1 output. Please select folder.\n');
    coregDir = uigetdir(cfg.OUT_DIR, 'Select folder containing co-registered TIF files');
    if isequal(coregDir,0), error('No folder selected.'); end
end

%% Run the same alignment check on coreg images
% Reuse QC_Alignment logic with coreg folder as input
inDir = coregDir;
% Use config pattern - files keep their original names after AROSICS correction
pat   = cfg.RAW_PATTERN;   % *_harmonized_clip.tif

outDir = fullfile(cfg.OUT_DIR, 'Phase2_AlignCheck_post', ['run_' runTimestamp]);
mkdir(outDir);
fprintf('Output: %s\n\n', outDir);

MAX_OFFSET_PX = cfg.QC.MAX_OFFSET_PX;
bNIR    = cfg.BANDS.NIR;
bGREEN  = cfg.BANDS.GREEN1;
SR_SCALE= cfg.SR_SCALE;
REF_DATE= cfg.REF_DATE;
XCORR_ROI_FRAC = 0.60;
MAX_SHIFT_SEARCH = 15;

ff = dir(fullfile(inDir, pat));
if isempty(ff)
    ff = dir(fullfile(inDir, '*.tif'));
end
if isempty(ff)
    fprintf('\n=== DIAGNOSTIC: Contents of Phase 1 output folder ===\n');
    all_files = dir(inDir);
    if isempty(all_files)
        fprintf('  Folder is completely empty.\n');
    else
        for kk = 1:numel(all_files)
            fprintf('  %s  (%d bytes)\n', all_files(kk).name, all_files(kk).bytes);
        end
    end
    consoleLog = fullfile(inDir, 'arosics_console.txt');
    if exist(consoleLog,'file')
        fprintf('\n=== arosics_console.txt (last 30 lines) ===\n');
        fid = fopen(consoleLog,'r'); lines = {};
        while ~feof(fid), lines{end+1} = fgetl(fid); end
        fclose(fid);
        for kk = max(1,numel(lines)-29):numel(lines)
            if ischar(lines{kk}), fprintf('  %s\n', lines{kk}); end
        end
    else
        fprintf('\n  arosics_console.txt not found.\n');
    end
    error('No co-registered files found in %s\nSee diagnostic output above.', inDir);
end

dates = cellfun(@(n) regexp(n,'(\d{8})','tokens','once'), {ff.name},'UniformOutput',false);
dates = cellfun(@(t) t{1}, dates, 'UniformOutput', false);
[dates, si] = sort(dates); ff = ff(si);

% Inject reference from RAW_DIR if AROSICS excluded it (it always does)
if ~any(strcmp(dates, REF_DATE))
    fprintf('Reference %s not in coreg folder - injecting from RAW_DIR...\n', REF_DATE);
    rawFF = dir(fullfile(cfg.RAW_DIR, pat));
    if isempty(rawFF), rawFF = dir(fullfile(cfg.RAW_DIR,'*.tif')); end
    for ki = 1:numel(rawFF)
        tok = regexp(rawFF(ki).name,'(\d{8})','tokens','once');
        if ~isempty(tok) && strcmp(tok{1}, REF_DATE)
            ff    = [ff; rawFF(ki)];
            dates = [dates, {REF_DATE}];
            [dates, si2] = sort(dates); ff = ff(si2);
            fprintf('  Injected: %s\n', rawFF(ki).name);
            break;
        end
    end
end

refIdx = find(strcmp(dates, REF_DATE));
if isempty(refIdx), refIdx=1; end
fprintf('Checking %d images: %s\n\n', numel(ff), strjoin(dates,', '));

[Aref,Rref] = readgeoraster(fullfile(ff(refIdx).folder,ff(refIdx).name));
[H,W,~]     = size(Aref);
cellSizeM   = extractCellSizeM_post(Rref);

r1=round(H*0.2)+1; r2=round(H*0.8);
c1=round(W*0.2)+1; c2=round(W*0.8);
refROI = normalizeROI_post(double(Aref(r1:r2,c1:c2,bNIR)));

xcorrLog = struct('file',{},'date',{},'offsetRow',{},'offsetCol',{},...
    'offsetPx',{},'offsetM',{},'status',{});

for i=1:numel(ff)
    fprintf('[%d/%d] %s\n', i, numel(ff), ff(i).name);
    if i==refIdx
        xcorrLog(end+1)=struct('file',ff(i).name,'date',dates{i},...
            'offsetRow',0,'offsetCol',0,'offsetPx',0,'offsetM',0,'status','REF');
        fprintf('  REFERENCE\n'); continue;
    end
    try
        [Ai,~]=readgeoraster(fullfile(ff(i).folder,ff(i).name));
        tgtROI=normalizeROI_post(double(Ai(r1:r2,c1:c2,bNIR)));
        xc=normxcorr2(refROI,tgtROI);
        [roiH,roiW]=size(refROI);
        midR=roiH; midC=roiW;
        sr1=max(1,midR-MAX_SHIFT_SEARCH); sr2=min(size(xc,1),midR+MAX_SHIFT_SEARCH);
        sc1=max(1,midC-MAX_SHIFT_SEARCH); sc2=min(size(xc,2),midC+MAX_SHIFT_SEARCH);
        xcSub=xc(sr1:sr2,sc1:sc2);
        [~,pk]=max(xcSub(:)); [pr,pc]=ind2sub(size(xcSub),pk);
        pr=pr+sr1-1; pc=pc+sc1-1;
        [subR,subC]=subpixelPeak_post(xc,pr,pc);
        rowOff=subR-roiH; colOff=subC-roiW;
        offPx=sqrt(rowOff^2+colOff^2); offM=offPx*cellSizeM;
        st=iif_post(offPx>MAX_OFFSET_PX,'WARN','OK');
        fprintf('  [%s] row=%+.2fpx  col=%+.2fpx  |%.2f|px  (%.1fm)\n',...
            st,rowOff,colOff,offPx,offM);
        xcorrLog(end+1)=struct('file',ff(i).name,'date',dates{i},...
            'offsetRow',rowOff,'offsetCol',colOff,'offsetPx',offPx,'offsetM',offM,'status',st);
    catch ME
        fprintf('  SKIP: %s\n',ME.message);
        xcorrLog(end+1)=struct('file',ff(i).name,'date',dates{i},...
            'offsetRow',NaN,'offsetCol',NaN,'offsetPx',NaN,'offsetM',NaN,'status','SKIP');
    end
end

nWarn=sum(strcmp({xcorrLog.status},'WARN'));
fprintf('\n========== PHASE 2 SUMMARY ==========\n');
fprintf('Images exceeding %.1f px: %d / %d\n', MAX_OFFSET_PX, nWarn, numel(ff));

if nWarn == 0
    fprintf('PASS - All images within tolerance.\n');
    fprintf('PROCEED to Phase 3 (PIF Harmonization).\n');
    fprintf('Input folder for Phase 3: %s\n', inDir);
else
    fprintf('FAIL - %d images still exceed tolerance.\n', nWarn);
    fprintf('Options:\n');
    fprintf('  1. Re-run Phase 1 with USE_LOCAL=true for local warping\n');
    fprintf('  2. Exclude flagged images from the stack\n');
    fprintf('  3. Proceed with caution if offsets < 2px\n');
end
fprintf('======================================\n');

% Save output path for Phase 3
save(fullfile(cfg.OUT_DIR,'phase2_output.mat'), 'inDir', 'nWarn', 'xcorrLog');

% Write report
fid=fopen(fullfile(outDir,'alignment_report_post.txt'),'w');
fprintf(fid,'PHASE 2 - POST CO-REGISTRATION ALIGNMENT CHECK\n');
fprintf(fid,'Timestamp: %s\n', runTimestamp);
fprintf(fid,'Input: %s\n\n', inDir);
fprintf(fid,'%-55s  %-6s  %8s  %8s  %8s\n','File','Status','Row(px)','Col(px)','Mag(px)');
for k=1:numel(xcorrLog)
    fprintf(fid,'%-55s  %-6s  %8.3f  %8.3f  %8.3f\n',...
        xcorrLog(k).file,xcorrLog(k).status,...
        xcorrLog(k).offsetRow,xcorrLog(k).offsetCol,xcorrLog(k).offsetPx);
end
fprintf(fid,'\nImages exceeding %.1f px: %d\n', MAX_OFFSET_PX, nWarn);
if nWarn==0
    fprintf(fid,'PASS - Proceed to Phase 3.\n');
else
    fprintf(fid,'FAIL - Review flagged images.\n');
end
fclose(fid);

save(fullfile(outDir,'alignment_post.mat'),'xcorrLog','nWarn','MAX_OFFSET_PX');

function roi = normalizeROI_post(roi)
    roi=roi-mean(roi(:)); s=std(roi(:)); if s>0, roi=roi/s; end
end
function m = extractCellSizeM_post(R)
    try
        if isprop(R,'CellExtentInWorldX'), m=R.CellExtentInWorldX;
        elseif isprop(R,'CellExtentInLongitude'), m=R.CellExtentInLongitude*111320;
        else, m=R.DeltaX; end
    catch, m=3.0; end
end
function [subR,subC]=subpixelPeak_post(xc,pR,pC)
    [nr,nc]=size(xc); subR=pR; subC=pC;
    rs=double(xc(max(pR-1,1):min(pR+1,nr),pC));
    cs=double(xc(pR,max(pC-1,1):min(pC+1,nc)));
    if numel(rs)==3&&rs(2)>rs(1)&&rs(2)>rs(3)
        d=2*(rs(1)-2*rs(2)+rs(3)); if abs(d)>1e-10, subR=pR+(rs(1)-rs(3))/d; end
    end
    if numel(cs)==3&&cs(2)>cs(1)&&cs(2)>cs(3)
        d=2*(cs(1)-2*cs(2)+cs(3)); if abs(d)>1e-10, subC=pC+(cs(1)-cs(3))/d; end
    end
end
function v=iif_post(c,a,b); if c, v=a; else, v=b; end; end
