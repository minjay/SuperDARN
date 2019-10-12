clear

% I basically apply the same transformation to the angles here as in the
% file examine_sam_data.

% load data from different fields.
hinfo = hdf5info('OP_2012-02-20_000000.h5');
mag_lat = hdf5read('OP_2012-02-20_000000.h5', '/lat_grid');
mag_lon = hdf5read('OP_2012-02-20_000000.h5', '/lon_grid');
mag_lat = fliplr(reshape(mag_lat, 41, 180)');
mag_lon = fliplr(reshape(mag_lon, 41, 180)');
ped_cond_all = hdf5read('OP_2012-02-20_000000.h5', '/Pedersen_Conductance');
ped_cond_all = reshape(ped_cond_all, 41, 180, size(ped_cond_all, 2));

ts = hdf5read('OP_2012-02-20_000000.h5', '/ts');

% convert latitude to colatitude.
mag_colat = (90 - mag_lat) / 180 * pi;

mag_lon = mag_lon / 180 * pi;

% fix the MLT time label misalignment "caused by function mypolar".
mag_lon_rot = mag_lon - pi / 2;
[x, y] = pol2cart(mag_lon_rot, mag_colat / pi * 180);

ped_cond = fliplr(ped_cond_all(:, :, 1)');
vmag = linspace(min(ped_cond(:)), max(ped_cond(:)), 10);
mypolar([0 2*pi], [0 max(mag_colat(:)) / pi * 180], x, y, ped_cond, vmag);