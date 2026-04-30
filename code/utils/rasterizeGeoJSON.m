function mask = rasterizeGeoJSON(path, H, W, R, label)
%RASTERIZEGEOJSON Rasterize all polygons from a GeoJSON file.
if nargin < 5
    label = 'GeoJSON';
end
mask = false(H, W);
if ~exist(path, 'file')
    warning('rasterizeGeoJSON: file not found: %s', path);
    return;
end
try
    data = jsondecode(fileread(path));
    feats = data.features;
    xmin = R.XWorldLimits(1);
    ymax = R.YWorldLimits(2);
    cw = R.CellExtentInWorldX;
    ch = R.CellExtentInWorldY;
    for fi = 1:numel(feats)
        rings = extractRings(feats(fi).geometry);
        for ri = 1:numel(rings)
            ring = rings{ri};
            if size(ring, 1) < 3, continue; end
            [xu, yu] = latlon2utm47(ring(:, 2), ring(:, 1));
            col = round((xu - xmin) / cw) + 1;
            row = round((ymax - yu) / ch) + 1;
            valid = col >= 1 & col <= W & row >= 1 & row <= H;
            if any(valid)
                mask = mask | poly2mask(col(valid), row(valid), H, W);
            end
        end
    end
catch ME
    warning('rasterizeGeoJSON (%s): %s', label, ME.message);
end
end
