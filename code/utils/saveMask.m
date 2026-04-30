function saveMask(path, mask, R)
%SAVEMASK Save logical mask as a uint8 GeoTIFF.
saveUint8Tif(path, uint8(mask), R);
end
