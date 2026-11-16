function cfg = loadPipelineConfigPortable()
% Load the current pipeline config without machine-specific fallback paths.

global PIPELINE_CONFIG_PATH;

if ischar(PIPELINE_CONFIG_PATH) && ~isempty(PIPELINE_CONFIG_PATH) && exist(PIPELINE_CONFIG_PATH, 'file')
    cfg = load(PIPELINE_CONFIG_PATH);
    return;
end

candidates = {
    fullfile(pwd, 'pipeline_logratio_config.mat')
    fullfile(pwd, 'pipeline_config.mat')
    };

for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        PIPELINE_CONFIG_PATH = candidates{i};
        cfg = load(PIPELINE_CONFIG_PATH);
        return;
    end
end

[f, p] = uigetfile('*.mat', 'Select pipeline configuration (.mat)');
if isequal(f, 0)
    error('Pipeline configuration file was not selected.');
end

PIPELINE_CONFIG_PATH = fullfile(p, f);
cfg = load(PIPELINE_CONFIG_PATH);
end
