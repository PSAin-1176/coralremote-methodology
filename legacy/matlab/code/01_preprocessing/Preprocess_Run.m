%% Preprocess_Run.m
clear; clc;
global PIPELINE_LOG_FID PIPELINE_SCRIPT_DIR;
if isempty(PIPELINE_LOG_FID), PIPELINE_LOG_FID = -1; end
if isempty(PIPELINE_SCRIPT_DIR)
    PIPELINE_SCRIPT_DIR = fileparts(fileparts(mfilename('fullpath')));
end

preprocessSteps = { ...
    '01_preprocessing\Preprocess_Masking.m', ...
    '01_preprocessing\Preprocess_Normalize.m', ...
    '01_preprocessing\Preprocess_WaveMask.m', ...
    '01_preprocessing\Preprocess_FeatureStack.m'};

logWrite(PIPELINE_LOG_FID, 'Preprocess_Run.m - start\n');
for iStep = 1:numel(preprocessSteps)
    runPhase(PIPELINE_SCRIPT_DIR, preprocessSteps{iStep}, PIPELINE_LOG_FID);
end
logWrite(PIPELINE_LOG_FID, 'Preprocess_Run.m - complete\n');
