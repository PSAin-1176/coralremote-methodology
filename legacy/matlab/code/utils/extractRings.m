function rings = extractRings(geom)
%EXTRACTRINGS Extract polygon rings from GeoJSON geometry.
rings = {};
try
    if strcmp(geom.type, 'Polygon')
        c = geom.coordinates;
        if isnumeric(c)
            rings{end+1} = c; %#ok<AGROW>
        elseif iscell(c)
            r = c{1};
            if iscell(r), r = cell2mat(r); end
            rings{end+1} = r; %#ok<AGROW>
        end
    elseif strcmp(geom.type, 'MultiPolygon')
        c = geom.coordinates;
        if ~iscell(c)
            for pi = 1:size(c, 1)
                rings{end+1} = squeeze(c(pi, 1, :, :)); %#ok<AGROW>
            end
            return;
        end
        for pi = 1:numel(c)
            poly = c{pi};
            if isnumeric(poly)
                rings{end+1} = poly; %#ok<AGROW>
            elseif iscell(poly)
                r = poly{1};
                if iscell(r), r = cell2mat(r); end
                if isnumeric(r) && size(r, 2) >= 2
                    rings{end+1} = r; %#ok<AGROW>
                end
            end
        end
    end
catch
end
end
