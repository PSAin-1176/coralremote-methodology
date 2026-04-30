%% RunAll_LogRatio.m
% ============================================================
% CORAL REEF BLEACHING PIPELINE - LOG-RATIO BRANCH
% ============================================================
clear; clc;
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(genpath(SCRIPT_DIR));
if ~exist('FEATURE_MODE_OVERRIDE', 'var') || isempty(FEATURE_MODE_OVERRIDE)
    FEATURE_MODE_OVERRIDE = 'A3';
end
RunAll_LogRatio_Mode(FEATURE_MODE_OVERRIDE);
