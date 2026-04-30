function P = paths_example()
% Example machine-local path configuration for the CoralRemote methodology pipeline.
% Copy this file to paths_local.m and edit the values for your machine.

P = struct();
P.REPO_ROOT = 'C:/path/to/CoralRemote_GitHub_ProductionCandidate';
P.DATA_ROOT = 'D:/CoralRemoteData';
P.RAW_DIR = fullfile(P.DATA_ROOT, 'raw', 'cropped_common_extent');
P.ACA_DIR = fullfile(P.DATA_ROOT, 'ancillary');
P.OUT_ROOT = fullfile(P.DATA_ROOT, 'output_logratio');
P.PREP_SOURCE_ROOT = fullfile(P.DATA_ROOT, 'prepared_source_logratio');
P.MASTERLOG_DIR = fullfile(P.OUT_ROOT, 'masterlog');
end
