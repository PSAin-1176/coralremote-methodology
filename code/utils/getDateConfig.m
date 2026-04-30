function dcfg = getDateConfig(DATE_CFG, dateStr)
%GETDATECONFIG Return DATE_CFG record for a specific YYYYMMDD date.
idx = find(strcmp({DATE_CFG.date}, dateStr), 1);
if isempty(idx)
    error('Date %s not found in DATE_CFG.', dateStr);
end
dcfg = DATE_CFG(idx);
end
