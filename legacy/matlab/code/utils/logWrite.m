function logWrite(fid, varargin)
%LOGWRITE Write to terminal and optional master log file.
msg = sprintf(varargin{:});
fprintf('%s', msg);
if isnumeric(fid) && fid > 0
    try
        fprintf(fid, '%s', msg);
    catch
    end
end
end
