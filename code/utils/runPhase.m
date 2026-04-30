function runPhase(scriptDir, scriptName, logFid)
%RUNPHASE Execute a phase script in an isolated workspace.
scriptPath = fullfile(scriptDir, scriptName);
t0 = tic;
logWrite(logFid, '[START] %s  %s\n', scriptName, ...
    char(datetime('now','Format','HH:mm:ss')));
try
    runInSandbox(scriptPath);
    logWrite(logFid, '[END]   %s  elapsed=%.1fs\n', scriptName, toc(t0));
catch ME
    elapsed = toc(t0);
    lineNum = 0;
    if ~isempty(ME.stack)
        lineNum = ME.stack(1).line;
    end
    logWrite(logFid, '[ERROR] %s  elapsed=%.1fs  %s (line %d)\n  %s\n', ...
        scriptName, elapsed, ME.identifier, lineNum, ME.message);
    rethrow(ME);
end
end

function runInSandbox(scriptPath)
run(scriptPath);
end
