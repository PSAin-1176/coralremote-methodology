function RunAll_LogRatio_Mode(featureMode, override)
% RunAll_LogRatio_Mode  Shared runner for mode-specific log-ratio tests.

if nargin < 1 || isempty(featureMode)
    featureMode = 'A3';
end
if nargin < 2 || isempty(override)
    override = struct();
end
featureMode = upper(char(string(featureMode)));

fprintf('============================================================\n');
fprintf('  CORAL REEF BLEACHING PIPELINE - Log-Ratio Branch\n');
fprintf('  Feature mode: %s\n', featureMode);
fprintf('  Started: %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('============================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
REPO_ROOT  = fileparts(SCRIPT_DIR);
addpath(genpath(SCRIPT_DIR));
P = loadLocalPaths(REPO_ROOT);

RAW_DIR  = getOverride(override, 'RAW_DIR', P.RAW_DIR);
ACA_DIR  = getOverride(override, 'ACA_DIR', P.ACA_DIR);
OUT_ROOT = getOverride(override, 'OUT_ROOT', P.OUT_ROOT);
runLabel = getOverride(override, 'RUN_LABEL', '');
if ~isempty(runLabel)
    OUT_DIR  = fullfile(OUT_ROOT, sprintf('%s_%s', featureMode, runLabel));
else
    OUT_DIR  = fullfile(OUT_ROOT, featureMode);
end
UDM_DIR  = '';
PREP_SOURCE_ROOT = getOverride(override, 'PREP_SOURCE_ROOT', P.PREP_SOURCE_ROOT);
USE_PREPARED_SOURCE = true;
SAVE_LEGACY_QUICKLOOK = getOverride(override, 'SAVE_LEGACY_QUICKLOOK', false);
FAST_MODE = getOverride(override, 'FAST_MODE', false);
RESUME_CACHE = getOverride(override, 'RESUME_CACHE', true);

RAW_PATTERN   = '*_cropped.tif';
COREG_PATTERN = '*_cropped_coreg.tif';
PIF_PATTERN   = '*_PIF.tif'; %#ok<NASGU>

REEF_MASK_PATH   = fullfile(ACA_DIR, 'GeoFeature.geojson');
DEEP_WATER_MASK  = fullfile(ACA_DIR, 'MDeepWater.geojson');
STABLE_SAND_MASK = fullfile(ACA_DIR, 'StableSand.geojson');
MASTER_SAND_MASK = fullfile(ACA_DIR, 'MasterSand.geojson');
ROAD_MASK        = fullfile(ACA_DIR, 'road_consensus_mask.tif');
WAVE_DEEP_MASK   = fullfile(ACA_DIR, 'WaveWDMask.geojson');
BATHY_PATH       = fullfile(ACA_DIR, 'bathymetry_0.tif');

GEOM_REF_DATE = '20240408';
RAD_REF_DATE  = '20241111';
REF_DATE      = GEOM_REF_DATE;

BANDS.COASTAL = 1;
BANDS.BLUE    = 2;
BANDS.GREEN1  = 3;
BANDS.GREEN2  = 4;
BANDS.YELLOW  = 5;
BANDS.RED     = 6;
BANDS.REDEDGE = 7;
BANDS.NIR     = 8;
SR_SCALE      = 10000;

LR.EPS                = 1e-6;
LR.NIR_WAVE_THRESHOLD = getOverride(override, 'NIR_WAVE_THRESHOLD', 0.07);
LR.APPLY_NORMALIZATION = getOverride(override, 'APPLY_NORMALIZATION', false);
LR.APPLY_MEDIAN       = false;
LR.MEDIAN_WINDOW      = 3;
LR.K_BR               = 0.82;
LR.K_GR               = 0.80;
LR.GLCM_WINDOW        = 5;
LR.APPLY_SUNGLINT_CORRECTION = getOverride(override, 'APPLY_SUNGLINT_CORRECTION', true);
LR.SUNGLINT_MIN_REG_PX       = getOverride(override, 'SUNGLINT_MIN_REG_PX', 50);
LR.SUNGLINT_NIR_P95_FRAC     = getOverride(override, 'SUNGLINT_NIR_P95_FRAC', 0.95);
LR.SUNGLINT_DEEP_THRESH_CM   = getOverride(override, 'SUNGLINT_DEEP_THRESH_CM', 300);
LR.SUNGLINT_PRIMARY_ROI      = getOverride(override, 'SUNGLINT_PRIMARY_ROI', 'MDEEP');
LR.SUNGLINT_FALLBACK_ROI     = getOverride(override, 'SUNGLINT_FALLBACK_ROI', 'BM_OUTSIDE_REEF');
LR.WAVE_TEST_NAME     = getOverride(override, 'WAVE_TEST_NAME', 'soft_wave_global_0p07');
LR.WAVE_DISABLE_GLOBAL = getOverride(override, 'WAVE_DISABLE_GLOBAL', false);
LR.WAVE_DISABLE_DATES  = getOverride(override, 'WAVE_DISABLE_DATES', {});
LR.WAVE_FLAG_ONLY_DATES = getOverride(override, 'WAVE_FLAG_ONLY_DATES', {});
LR.WAVE_THRESHOLD_OVERRIDE_DATES  = getOverride(override, 'WAVE_THRESHOLD_OVERRIDE_DATES', {});
LR.WAVE_THRESHOLD_OVERRIDE_VALUES = getOverride(override, 'WAVE_THRESHOLD_OVERRIDE_VALUES', []);

switch featureMode
    case 'A3'
        FEATURE_NAMES   = {'B_over_G1','G1_over_G2','B_over_R','DII_BR','DII_GR'};
        STAGE1_FEAT_IDX = 1:5;
        STAGE2_FEAT_IDX = 1:5;
    case 'A4'
        FEATURE_NAMES   = {'B_over_G1','G1_over_G2','B_over_R','DII_BR','DII_GR', ...
                           'VisMean','Whiteness','LocalVar_DII_GR'};
        STAGE1_FEAT_IDX = 1:8;
        STAGE2_FEAT_IDX = 1:8;
    case 'A5B'
        FEATURE_NAMES   = {'B_over_G1','G1_over_G2','B_over_R','DII_BR','DII_GR'};
        STAGE1_FEAT_IDX = 1:5;
        STAGE2_FEAT_IDX = 1:5;
    case 'C1'
        FEATURE_NAMES   = {'B_over_G1','G1_over_G2','B_over_R','DII_BR','DII_GR', ...
                           'Blue_GLCM_Contrast','Blue_GLCM_Entropy', ...
                           'NIR_GLCM_Contrast','NIR_GLCM_Entropy'};
        STAGE1_FEAT_IDX = 1:9;
        STAGE2_FEAT_IDX = 1:9;
    otherwise
        error('Unsupported feature mode: %s', featureMode);
end
FEATURE_MODE     = featureMode;
LR.FEATURE_NAMES = FEATURE_NAMES;

DATE_CFG = struct( ...
    'date', { ...
        '20230404','20230509','20230519','20230629', ...
        '20240307','20240408', ...
        '20240514','20240613','20240620', ...
        '20240926','20241111','20250518'}, ...
    'baseline_coral', { ...
        true, true, true, true, ...
        true, false, ...
        false, false, false, ...
        false, false, false}, ...
    'standalone', { ...
        true, true, true, true, ...
        true, false, ...
        true, true, true, ...
        true, true, true} ...
);

QC.MAX_OFFSET_PX = 1.0;
RF.N_TREES       = 200;
RF.MIN_LEAF      = 5;
RF.MIN_SAMPLES   = 300;
RF.SLIC_K        = 2000;
RF.SLIC_M        = 10;
RF.STAGE1_PREDICT_BATCH = getOverride(override, 'STAGE1_PREDICT_BATCH', 10000);
PAR.ENABLE       = getOverride(override, 'ENABLE_PARALLEL', true);
PAR.MODE         = getOverride(override, 'PARALLEL_MODE', 'processes');

ACA_CLASSES = {'Coral/Algae','Rubble','Sand','Seagrass','Rock'};
CLASS_SCHEMA = {
    'CA',  'Coral/Algae healthy',   [0   150  0  ];
    'BL',  'Coral/Algae bleached',  [255 255 255];
    'RBL', 'Rubble',                [180 140  80];
    'SND', 'Sand',                  [255 240 180];
};

UTIL_DIR   = fullfile(SCRIPT_DIR, 'utils');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
MASTERLOG_DIR = getOverride(override, 'MASTERLOG_DIR', P.MASTERLOG_DIR);
if ~exist(MASTERLOG_DIR, 'dir'), mkdir(MASTERLOG_DIR); end

configPath = fullfile(OUT_DIR, 'pipeline_logratio_config.mat');
save(configPath, ...
    'RAW_DIR','ACA_DIR','OUT_ROOT','OUT_DIR','UDM_DIR','PREP_SOURCE_ROOT','USE_PREPARED_SOURCE','SAVE_LEGACY_QUICKLOOK','SCRIPT_DIR','UTIL_DIR', ...
    'MASTERLOG_DIR', ...
    'RAW_PATTERN','COREG_PATTERN','PIF_PATTERN', ...
    'REEF_MASK_PATH','DEEP_WATER_MASK','STABLE_SAND_MASK', ...
    'MASTER_SAND_MASK','ROAD_MASK','WAVE_DEEP_MASK','BATHY_PATH', ...
    'GEOM_REF_DATE','RAD_REF_DATE','REF_DATE', ...
    'BANDS','SR_SCALE','DATE_CFG','LR','QC','RF','PAR', ...
    'FEATURE_MODE','FEATURE_NAMES','STAGE1_FEAT_IDX','STAGE2_FEAT_IDX', ...
    'ACA_CLASSES','CLASS_SCHEMA', ...
    'FAST_MODE','RESUME_CACHE');

global PIPELINE_CONFIG_PATH PIPELINE_SCRIPT_DIR PIPELINE_LOG_FID;
PIPELINE_CONFIG_PATH = configPath;
PIPELINE_SCRIPT_DIR  = SCRIPT_DIR;

logTimestamp = char(datetime('now','TimeZone','local'),'yyyyMMdd_HHmmss');
LOG_FILE = fullfile(MASTERLOG_DIR, sprintf('pipeline_logratio_%s_master_log_%s.txt', lower(FEATURE_MODE), logTimestamp));
PIPELINE_LOG_FID = fopen(LOG_FILE, 'a');
if PIPELINE_LOG_FID < 0
    warning('Could not open master log file: %s', LOG_FILE);
    PIPELINE_LOG_FID = -1;
end

logWrite(PIPELINE_LOG_FID, '============================================================\n');
logWrite(PIPELINE_LOG_FID, '  LOG-RATIO PIPELINE RUN - %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
logWrite(PIPELINE_LOG_FID, '  FEATURE_MODE   : %s\n', FEATURE_MODE);
logWrite(PIPELINE_LOG_FID, '  RAW_DIR        : %s\n', RAW_DIR);
logWrite(PIPELINE_LOG_FID, '  OUT_DIR        : %s\n', OUT_DIR);
logWrite(PIPELINE_LOG_FID, '  GEOM_REF_DATE  : %s\n', GEOM_REF_DATE);
logWrite(PIPELINE_LOG_FID, '  RAD_REF_DATE   : %s\n', RAD_REF_DATE);
logWrite(PIPELINE_LOG_FID, '  NIR wave thr   : %.3f\n', LR.NIR_WAVE_THRESHOLD);
logWrite(PIPELINE_LOG_FID, '  Sunglint       : %s\n', ternary(LR.APPLY_SUNGLINT_CORRECTION, 'Hedley', 'off'));
logWrite(PIPELINE_LOG_FID, '  Normalize      : %s\n', ternary(LR.APPLY_NORMALIZATION, 'on', 'off'));
logWrite(PIPELINE_LOG_FID, '  Parallel       : %s\n', ternary(PAR.ENABLE, PAR.MODE, 'off'));
logWrite(PIPELINE_LOG_FID, '  WAVE_TEST      : %s\n', LR.WAVE_TEST_NAME);
logWrite(PIPELINE_LOG_FID, '  FEATURES       : %s\n', strjoin(FEATURE_NAMES, ', '));
logWrite(PIPELINE_LOG_FID, '  RESUME_CACHE   : %s\n', ternary(RESUME_CACHE, 'on', 'off'));
logWrite(PIPELINE_LOG_FID, '  FAST_MODE      : %s\n', ternary(FAST_MODE, 'on', 'off'));
logWrite(PIPELINE_LOG_FID, '============================================================\n\n');

if PAR.ENABLE
    try
        p = gcp('nocreate');
        if isempty(p)
            if strcmpi(PAR.MODE, 'threads')
                parpool('threads');
            else
                parpool;
            end
            logWrite(PIPELINE_LOG_FID, '  Parallel pool  : started (%s)\n\n', PAR.MODE);
        else
            logWrite(PIPELINE_LOG_FID, '  Parallel pool  : reusing existing pool\n\n');
        end
    catch ME
        PAR.ENABLE = false;
        save(configPath, 'PAR', '-append');
        logWrite(PIPELINE_LOG_FID, '  Parallel pool  : unavailable (%s)\n\n', ME.message);
    end
end

phaseList = getOverride(override, 'PHASE_LIST', { ...
    '01_preprocessing\Preprocess_Run.m', ...
    '02_baseline_construction\Build_Baseline.m', ...
    '03_classification\Train_Stage1.m', ...
    '03_classification\Classify_Stage1.m', ...
    '02_baseline_construction\Build_BaselineCoralMask.m', ...
    '03_classification\Train_Stage2.m', ...
    '03_classification\Classify_Stage2.m', ...
    '03_classification\Analyze_Change.m'});

if RESUME_CACHE
    cachedSkipPreproc    = phaseArtifactExists(OUT_DIR, {'phase2_output.mat','phase3_output.mat','phase4_output.mat'});
    cachedSkipFeatures   = phaseArtifactExists(OUT_DIR, {'phase5_output.mat','phase6_output.mat'});
    cachedSkipStage1     = phaseArtifactExists(OUT_DIR, {'phase_stage1_output.mat'});
    cachedSkipCoralBase  = phaseArtifactExists(OUT_DIR, {'baseline_coral_output.mat'});
    cachedSkipStage2     = phaseArtifactExists(OUT_DIR, {'phase_stage2_output.mat','phase_stage2_classify.mat'});
    cachedSkipChange     = phaseArtifactExists(OUT_DIR, {'change_output.mat'});

    logWrite(PIPELINE_LOG_FID, '  Cache status   : preproc=%s features=%s stage1=%s coral=%s stage2=%s change=%s\n\n', ...
        ternary(cachedSkipPreproc,'hit','miss'), ternary(cachedSkipFeatures,'hit','miss'), ...
        ternary(cachedSkipStage1,'hit','miss'), ternary(cachedSkipCoralBase,'hit','miss'), ...
        ternary(cachedSkipStage2,'hit','miss'), ternary(cachedSkipChange,'hit','miss'));
else
    cachedSkipPreproc = false;
    cachedSkipFeatures = false;
    cachedSkipStage1   = false;
    cachedSkipCoralBase = false;
    cachedSkipStage2   = false;
    cachedSkipChange   = false;
end

if cachedSkipPreproc
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phases 0-4 (cached outputs found)\n');
end
if cachedSkipFeatures
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phases 5-7 (cached outputs found)\n');
end
if cachedSkipStage1
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phases 8-9 (cached outputs found)\n');
end
if cachedSkipCoralBase
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phase 10 (cached output found)\n');
end
if cachedSkipStage2
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phases 11-12 (cached outputs found)\n');
end
if cachedSkipChange
    logWrite(PIPELINE_LOG_FID, '  Resume action  : skipping Phase 13 (cached output found)\n');
end
if cachedSkipPreproc || cachedSkipFeatures || cachedSkipStage1 || cachedSkipCoralBase || cachedSkipStage2 || cachedSkipChange
    logWrite(PIPELINE_LOG_FID, '\n');
end

for iPhase = 1:numel(phaseList)
    scriptName = phaseList{iPhase};
    if RESUME_CACHE && shouldSkipPhase(scriptName, cachedSkipPreproc, cachedSkipFeatures, cachedSkipStage1, cachedSkipCoralBase, cachedSkipStage2, cachedSkipChange)
        logWrite(PIPELINE_LOG_FID, '[SKIP] %s  (cached)\n', scriptName);
        continue;
    end
    runPhase(SCRIPT_DIR, scriptName, PIPELINE_LOG_FID);
end

logWrite(PIPELINE_LOG_FID, '\n============================================================\n');
logWrite(PIPELINE_LOG_FID, '  LOG-RATIO PIPELINE COMPLETE - %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
logWrite(PIPELINE_LOG_FID, '  Master log: %s\n', LOG_FILE);
logWrite(PIPELINE_LOG_FID, '============================================================\n');
fprintf('\nMaster log saved: %s\n', LOG_FILE);
if PIPELINE_LOG_FID > 0, fclose(PIPELINE_LOG_FID); end
end

function v = getOverride(override, name, defaultValue)
if isstruct(override) && isfield(override, name)
    v = override.(name);
else
    v = defaultValue;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function tf = phaseArtifactExists(outDir, fileNames)
tf = false;
for i = 1:numel(fileNames)
    if exist(fullfile(outDir, fileNames{i}), 'file')
        tf = true;
        return;
    end
end
end

function tf = shouldSkipPhase(scriptName, skipPreproc, skipFeatures, skipStage1, skipCoralBase, skipStage2, skipChange)
scriptName = strrep(scriptName, '/', '\');
switch scriptName
    case {'01_preprocessing\QC_Alignment.m','01_preprocessing\Coreg_AROSICS_runner.m','01_preprocessing\QC_AlignmentPost.m','PIF_Harmonize.m','QC_Radiometric.m'}
        tf = skipPreproc;
    case {'Preprocess_DII.m','Explore_Spectral.m', ...
          '01_preprocessing\Preprocess_Run.m', ...
          '02_baseline_construction\Build_Baseline.m', ...
          '01_preprocessing\Preprocess_Masking.m', ...
          '01_preprocessing\Preprocess_Normalize.m', ...
          '01_preprocessing\Preprocess_WaveMask.m', ...
          '01_preprocessing\Preprocess_FeatureStack.m'}
        tf = skipFeatures;
    case {'Train_Stage1.m','Classify_AllDates.m', ...
          '03_classification\Train_Stage1.m', ...
          '03_classification\Classify_Stage1.m'}
        tf = skipStage1;
    case {'Baseline_Coral.m','02_baseline_construction\Build_BaselineCoralMask.m'}
        tf = skipCoralBase;
    case {'Train_Stage2.m','Classify_Stage2.m', ...
          '03_classification\Train_Stage2.m', ...
          '03_classification\Classify_Stage2.m'}
        tf = skipStage2;
    case {'Change_Analysis.m','03_classification\Analyze_Change.m'}
        tf = skipChange;
    otherwise
        tf = false;
end
end


