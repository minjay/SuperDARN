clear
clc

% compute Joule heating rate of the entire field (sum of all components)

% load sam pot (provided by Xueling)
load('2012-02-29_00_00-2012-02-29_04_00_sam_pot_data.mat')

% convert latitude to colatitude (90 - x).
lat_grd_mat = fliplr((90 - reshape(lat_grd, [41 180])') / 180 * pi);
lon_grd_mat = fliplr(reshape(lon_grd, [41 180])' / 180 * pi);

% select first time point (do we need to draw a time series?)
time_point = 1;
% potential field from SAM
sam_pot = fliplr(reshape(pot(time_point, :), [41 180])');

% residual potential field fitted by needlets
load('vs_0_2000_tweak_252.mat')

% size = number of data points * number of replicates
post_samples_eta = post_samples.eta;
% size = number of data points * number of replicates
post_samples_c = squeeze(post_samples.c);

% stretch colatitude to the entire sphere (*4)
theta_grd_vec = lat_grd_mat(:) * 4;
lon_grd_vec = lon_grd_mat(:);

B = 2;
j_min = 2;
j_max = 3;
[~, ~, A] = get_A_ss(B, j_min, j_max, theta_grd_vec, lon_grd_vec);

knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat_grid, ~] = bspline_basismatrix(4, knots, theta_grd_vec);
b_mat_grid(:, 1) = 1;
% size = number of data points * number of replicates
std_grid_vec = exp(b_mat_grid * post_samples_eta);

% size = number of data points * number of replicates
Ac = A * post_samples_c;

% Unit is 10^-6 Volt.
% confirm the unit?
% conditional distribution E(Z(s)|Z_1,...,Z_n)
cf = reshape(mean(std_grid_vec .* Ac, 2), size(lon_grd_mat));
% convert to Volt
cf = cf / 1e6;

phi_rot = lon_grd_mat + pi/2;
[x, y] = pol2cart(phi_rot, lat_grd_mat / pi * 180);
vmag = linspace(min(cf(:)), max(cf(:)), 10);
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, cf, vmag);

HX = lat_grd_mat(1, :);
HY = lon_grd_mat(:, 1);
[FX, FY] = gradient(sam_pot + cf, HX, HY);
[FX_sam, FY_sam] = gradient(cf, HX, HY);

R = 6.5*1e6;
E_theta = -FX/R;
E_phi = -FY./(R*sin(lat_grd_mat));
E_theta_sam = -FX_sam/R;
E_phi_sam = -FY_sam./(R*sin(lat_grd_mat));

% load conductance data
ped_cond_all = hdf5read('OP_2012-02-29_000000.h5', '/Pedersen_Conductance');
ped_cond_all = reshape(ped_cond_all, 41, 180, size(ped_cond_all, 2));
ped_cond = fliplr(ped_cond_all(:, :, time_point + 1)');

% compute energy
energy = (E_theta.^2+E_phi.^2).*ped_cond;
energy_sam = (E_theta_sam.^2+E_phi_sam.^2).*ped_cond;

% drop first column (Inf)
energy = energy(:, 2:end);
energy_sam = energy_sam(:, 2:end);
lat_grd_mat = lat_grd_mat(:, 2:end);
lon_grd_mat = lon_grd_mat(:, 2:end);

phi_rot = lon_grd_mat + pi/2;
[x, y] = pol2cart(phi_rot, lat_grd_mat / 4 /pi * 180);
vmag = linspace(min(energy(:)), max(energy(:)), 10);
figure
mypolar([0 2*pi], [0 max(lat_grd_mat(:) / 4)/pi*180], x, y, energy, vmag);
vmag = linspace(min(energy_sam(:)), max(energy_sam(:)), 10);
figure
mypolar([0 2*pi], [0 max(lat_grd_mat(:) / 4)/pi*180], x, y, energy_sam, vmag);

% compute the area of each latitudinal band
theta_one = lat_grd_mat(1, :);
n_theta = length(theta_one);
area_theta = zeros(1, n_theta);
for i = 1:n_theta
    if i==1
        theta_lower = 0;
    else
        theta_lower = (theta_one(i-1)+theta_one(i))/2;
    end
    if i==n_theta
        theta_upper = pi/4;
    else
        theta_upper = (theta_one(i)+theta_one(i+1))/2;
    end
    area_theta(i) = areaquad(90-theta_lower/pi*180, -180, 90-theta_upper/pi*180, 180);
end

% compute integrated energy
tot_area = 4*pi*R^2;
energy = reshape(energy, size(lon_grd_mat));
energy_sam = reshape(energy_sam, size(lon_grd_mat));
% dividing by 1e9 is due to the unit giga
% since the areaquad function only gives a fraction of the unit sphere's area ranging from 0 to 1,
% we need to multiply back the total area
int_energy = sum(mean(energy, 1).*area_theta)*tot_area/1e9;
int_energy_sam = sum(mean(energy_sam, 1).*area_theta)*tot_area/1e9;