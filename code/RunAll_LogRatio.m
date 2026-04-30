%% RunAll_LogRatio.m
% ============================================================
% CORAL REEF BLEACHING PIPELINE - LOG-RATIO BRANCH
%
% MAIN ENTRY POINT
%   Run this file to execute the production methodology pipeline
%   using the default feature mode. The shared runner is
%   RunAll_LogRatio_Mode.m.
%
% REQUIRED LOCAL SETUP
%   1. Create config/paths_local.m from config/paths_example.m.
%   2. Point paths_local.m to local data roots outside this repo.
%   3. Ensure required local vector/raster support files exist
%      under the configured ACA/data paths.
%
% TESTED ENVIRONMENT
%   - MATLAB R2025a
%   - Windows
%
% DATA BOUNDARY
%   This repository is code only.
%   Do not add PlanetScope RAW data, processed .mat datasets,
%   generated outputs, or local archives to Git.
%
% OPTIONAL EXTERNAL DEPENDENCY
%   Coreg_AROSICS_runner.m may require a local Python/Conda
%   environment with AROSICS if that preprocessing path is used.
%
% PRIMARY METHOD SECTIONS
%   1. Preprocessing
%   2. Baseline construction
%   3. Classification
%   4. Evaluation note only (trimmed in first public release)
%
% ============================================================
clear; clc;
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(genpath(SCRIPT_DIR));
if ~exist('FEATURE_MODE_OVERRIDE', 'var') || isempty(FEATURE_MODE_OVERRIDE)
    FEATURE_MODE_OVERRIDE = 'A3';
end
RunAll_LogRatio_Mode(FEATURE_MODE_OVERRIDE);
