function saveMultiBandTif(path, A, R)
%SAVEMULTIBANDTIF Save single or multiband single-precision GeoTIFF.
if ~isa(A, 'single')
    A = single(A);
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
    t.setTag('BitsPerSample', 32);
    t.setTag('SampleFormat', 3);
    t.setTag('SamplesPerPixel', nb);
    t.setTag('Compression', 5);
    t.write(A);
    t.close();
catch ME
    error('saveMultiBandTif: %s', ME.message);
end
end
