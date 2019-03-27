% load data from different fields.
hinfo = hdf5info('op_cond_matuso_2012229.h5');
mag_lat = hdf5read('op_cond_matuso_2012229.h5', '/magnetic_latitude');
mag_local_time = hdf5read('op_cond_matuso_2012229.h5', '/magnetic_local_time');
% name starts from /5, /10 to /235. There are in total 47 of them.
% /5 represents 5 mins.
ped_cond = hdf5read('op_cond_matuso_2012229.h5', '/5/pedersen');

% convert latitude to colatitude.
mag_colat = (90 - mag_lat) / 180 * pi;

% mag_local_time is in the unit of MLT.
mag_local_time(mag_local_time < 0) = mag_local_time(mag_local_time < 0) + 24;
mag_long = mag_local_time / 24 * 2 * pi;

% fix the MLT time label misalignment "caused by function mypolar".
mag_long_rot = mag_long - pi / 2;
[x, y] = pol2cart(mag_long_rot, mag_colat / pi * 180);

vmag = linspace(min(ped_cond(:)), max(ped_cond(:)), 10);
mypolar([0 2*pi], [0 max(mag_colat(:)) / pi * 180], x, y, ped_cond, vmag);