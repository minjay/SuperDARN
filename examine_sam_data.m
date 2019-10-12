load('2012-02-29_00_00-2012-02-29_04_00_sam_pot_data.mat')
% it seems that lon_grd is in lon_MLT unit.

% convert latitude to colatitude.
lat_grd_mat = fliplr((90 - reshape(lat_grd, [41 180])') / 180 * pi);
lon_grd_mat = fliplr(reshape(lon_grd, [41 180])' / 180 * pi);
% lat_grd_mat is 180*41
% 0 0.0175 0.0349 ... 0.6981
% 0 0.0175 0.0349 ... 0.6981
% ...
% lon_grd_mat is 180*41
% 0 0 0 ...
% 0.0351 0.0351 0.0351 ...
% ...
% 6.2832 6.2832 6.2832 ...

% fix the MLT time label misalignment "caused by function mypolar".
lon_rot = lon_grd_mat - pi/2;
[x, y] = pol2cart(lon_rot, lat_grd_mat / pi * 180);

for i = 20
    figure
    sam_pot = fliplr(reshape(pot(i, :), [41 180])');

    vmag = linspace(min(sam_pot(:)), max(sam_pot(:)), 10);
    mypolar([0 2*pi], [0 max(lat_grd_mat(:)) / pi * 180], x, y, sam_pot, vmag);
    
    title(ts(i, :))
end

angle = 270 / 180 * pi;

[~, I] = min(abs(lon_grd_mat(:, 1) - angle));

figure
plot(lat_grd_mat(I, :), sam_pot(I, :))


angle = 90 / 180 * pi;

[~, I] = min(abs(lon_grd_mat(:, 1) - angle));

figure
plot(lat_grd_mat(I, :), sam_pot(I, :))