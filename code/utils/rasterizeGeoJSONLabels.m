function [mask, label_img] = rasterizeGeoJSONLabels(path, H, W, R, class_map, tag)
%RASTERIZEGEOJSONLABELS Rasterize labeled polygons from a GeoJSON file.
%
% class_map is an Nx2 cell array:
%   { 'HC', 1; 'SND', 2; 'RBL', 3 }

if nargin < 6
    tag = 'GeoJSONLabels';
end

mask = false(H, W);
label_img = zeros(H, W, 'uint8');

if ~exist(path, 'file')
    warning('rasterizeGeoJSONLabels: file not found: %s', path);
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
        prop = feats(fi).properties;
        if isfield(prop, 'Class')
            cstr = prop.Class;
        elseif isfield(prop, 'class')
            cstr = prop.class;
        else
            continue;
        end

        cid = [];
        for mi = 1:size(class_map, 1)
            if strcmp(class_map{mi, 1}, cstr)
                cid = class_map{mi, 2};
                break;
            end
        end
        if isempty(cid), continue; end

        rings = extractRings(feats(fi).geometry);
        for ri = 1:numel(rings)
            ring = rings{ri};
            if size(ring, 1) < 3, continue; end
            [xu, yu] = latlon2utm47(ring(:, 2), ring(:, 1));
            col = round((xu - xmin) / cw) + 1;
            row = round((ymax - yu) / ch) + 1;
            valid = col >= 1 & col <= W & row >= 1 & row <= H;
            if ~any(valid), continue; end
            pm = poly2mask(col(valid), row(valid), H, W);
            mask = mask | pm;
            label_img(pm & label_img == 0) = uint8(cid);
        end
    end

    fprintf('  %s: %d features -> %d px\n', tag, numel(feats), nnz(mask));
catch ME
    warning('rasterizeGeoJSONLabels (%s): %s', tag, ME.message);
end
end
