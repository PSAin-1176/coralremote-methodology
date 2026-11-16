function saveUint8Tif(path, A, R)
%SAVEUINT8TIF Save uint8 raster as GeoTIFF.
if ~isa(A, 'uint8')
    A = uint8(A);
end
try
    geotiffwrite(path, A, R);
    return;
catch
end
try
    geotiffwrite(path, A, R, 'CoordRefSysCode', 32647);
    return;
catch
end
try
    [H, W, nb] = size(A);
    t = Tiff(path, 'w');
    t.setTag('Photometric', 1);
    t.setTag('PlanarConfig', 1);
    t.setTag('ImageLength', H);
    t.setTag('ImageWidth', W);
    t.setTag('BitsPerSample', 8);
    t.setTag('SampleFormat', 1);
    t.setTag('SamplesPerPixel', nb);
    t.setTag('Compression', 5);
    t.write(A);
    t.close();
catch ME
    error('saveUint8Tif: %s', ME.message);
end
end
