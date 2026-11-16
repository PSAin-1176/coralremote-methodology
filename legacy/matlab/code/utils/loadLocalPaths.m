function P = loadLocalPaths(repoRoot)
% Load machine-local paths from config/paths_local.m.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

configDir = fullfile(repoRoot, 'config');
localPathFile = fullfile(configDir, 'paths_local.m');
if ~exist(localPathFile, 'file')
    error(['Missing config/paths_local.m. Copy config/paths_example.m to ' ...
           'config/paths_local.m and edit the paths for this machine.']);
end

addpath(configDir);
P = paths_local();

requiredFields = {'RAW_DIR','ACA_DIR','OUT_ROOT','PREP_SOURCE_ROOT','MASTERLOG_DIR'};
for i = 1:numel(requiredFields)
    fieldName = requiredFields{i};
    if ~isfield(P, fieldName) || isempty(P.(fieldName))
        error('paths_local.m must define %s.', fieldName);
    end
end
end
