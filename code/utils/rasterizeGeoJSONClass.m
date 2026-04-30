function mask = rasterizeGeoJSONClass(path, H, W, R, className)
%RASTERIZEGEOJSONCLASS Rasterize only features matching a Class property.
mask = false(H, W);
if ~exist(path, 'file')
    warning('rasterizeGeoJSONClass: file not found: %s', path);
    return;
end
try
    data = jsondecode(fileread(path));
    feats = data.features;
    xmin = R.XWorldLimits(1);
    ymax = R.YWorldLimits(2);
    cw = R.CellExtentInWorldX;
    ch = R.CellExtentInWorldY;
    found = false;
    for fi = 1:numel(feats)
        prop = feats(fi).properties;
        cls = '';
        if isfield(prop, 'Class')
            cls = prop.Class;
        elseif isfield(prop, 'class')
            cls = prop.class;
        end
        if ~strcmpi(cls, className), continue; end
        found = true;
        rings = extractRings(feats(fi).geometry);
        for ri = 1:numel(rings)
            ring = rings{ri};
            [xu, yu] = latlon2utm47(ring(:, 2), ring(:, 1));
            col = round((xu - xmin) / cw) + 1;
            row = round((ymax - yu) / ch) + 1;
            valid = col >= 1 & col <= W & row >= 1 & row <= H;
            if any(valid)
                mask = mask | poly2mask(col(valid), row(valid), H, W);
            end
        end
    end
    if ~found
        warning('rasterizeGeoJSONClass: no feature with Class="%s" in %s', ...
            className, path);
    end
catch ME
    warning('rasterizeGeoJSONClass: %s', ME.message);
end
end
