function [x, y] = latlon2utm47(lat, lon)
%LATLON2UTM47 Convert WGS84 lat/lon to UTM zone 47N.
a = 6378137.0;
f = 1 / 298.257223563;
k0 = 0.9996;
E0 = 500000;
lon0 = 99;
e2 = 2 * f - f^2;
e4 = e2^2;
e6 = e2^3;
lr = deg2rad(lat);
wr = deg2rad(lon);
dw = wr - deg2rad(lon0);
N_ = a ./ sqrt(1 - e2 .* sin(lr).^2);
T = tan(lr).^2;
C = e2 / (1 - e2) .* cos(lr).^2;
A = cos(lr) .* dw;
M = a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) .* lr ...
    - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) .* sin(2 * lr) ...
    + (15 * e4 / 256 + 45 * e6 / 1024) .* sin(4 * lr) ...
    - (35 * e6 / 3072) .* sin(6 * lr));
x = k0 .* N_ .* (A + (1 - T + C) .* A.^3 / 6 + ...
    (5 - 18 * T + T.^2 + 72 * C - 58 * e2 / (1 - e2)) .* A.^5 / 120) + E0;
y = k0 .* (M + N_ .* tan(lr) .* (A.^2 / 2 + ...
    (5 - T + 9 * C + 4 * C.^2) .* A.^4 / 24 + ...
    (61 - 58 * T + T.^2 + 600 * C - 330 * e2 / (1 - e2)) .* A.^6 / 720));
end
