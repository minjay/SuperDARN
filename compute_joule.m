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
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_grd_vec, lon_grd_vec);

knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat_grid, ~] = bspline_basismatrix(4, knots, theta_grd_vec);
b_mat_grid(:, 1) = 1;
% size = number of data points * number of replicates
std_grid_vec = exp(b_mat_grid * post_samples_eta);

% size = number of data points * number of replicates
Ac = A * post_samples_c;

Ac_l2 = A(:, 1:Npix(1)) * post_samples_c(1:Npix(1), :);
Ac_l3 = A(:, (Npix(1)+1):size(A, 2)) * post_samples_c((Npix(1)+1):size(A, 2), :);

% Unit is 10^-6 Volt.
% confirm the unit?
% conditional expectation E(Z(s)|Z_1,...,Z_n)
% conditional distribution Z(s)|Z_1,...,Z_n
cf_all = reshape(std_grid_vec .* Ac, [size(lon_grd_mat) size(Ac, 2)]);
cf_l2_all = reshape(std_grid_vec .* Ac_l2, [size(lon_grd_mat) size(Ac_l2, 2)]);
cf_l3_all = reshape(std_grid_vec .* Ac_l3, [size(lon_grd_mat) size(Ac_l3, 2)]);

% convert to Volt
cf_all = cf_all / 1e6;
cf_l2_all = cf_l2_all / 1e6;
cf_l3_all = cf_l3_all / 1e6;
cf = mean(cf_all, 3);
cf_l2 = mean(cf_l2_all, 3);
cf_l3 = mean(cf_l3_all, 3);

% convert to kV
cf_kV = cf / 1e3;
cf_l2_kV = cf_l2 / 1e3;
cf_l3_kV = cf_l3 / 1e3;

phi_rot = lon_grd_mat + pi/2;
[x, y] = pol2cart(phi_rot, lat_grd_mat / pi * 180);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.075], [0.05 0.02], [0.05 0.2]);
subplot(1, 3, 1)
vmag = linspace(min(cf_l2_kV(:)), max(cf_l2_kV(:)), 10);
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, cf_l2_kV, vmag);
title('j = 2')
text(-50, -50, sprintf('Min\n%2.1f',min(cf_l2_kV(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf_l2_kV(:))),'FontName','times','Fontsize',10)
subplot(1, 3, 2)
vmag = linspace(min(cf_l3_kV(:)), max(cf_l3_kV(:)), 10);
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, cf_l3_kV, vmag);
title('j = 3')
text(-50, -50, sprintf('Min\n%2.1f',min(cf_l3_kV(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf_l3_kV(:))),'FontName','times','Fontsize',10)
subplot(1, 3, 3)
vmag = linspace(min(cf_kV(:)), max(cf_kV(:)), 10);
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, cf_kV, vmag);
title('Total')
text(-50, -50, sprintf('Min\n%2.1f',min(cf_kV(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf_kV(:))),'FontName','times','Fontsize',10)

cmax = max([max(abs(cf_l2_kV)) max(abs(cf_l3_kV)) max(abs(cf_kV))]);
caxis([-cmax cmax])
colormap(jet)
h = colorbar;
set(h, 'Position', [.85 0.1 .025 .8]);

print -painters -depsc multi_res_fitted.eps

%% Compute Joule heating rate

HX = lat_grd_mat(1, :);
HY = lon_grd_mat(:, 1);
% get numeric x-, y- gradient fields
[FX_sam, FY_sam] = gradient(sam_pot, HX, HY);
FX = zeros(size(cf_all));
FY = zeros(size(cf_all));
n_rep = size(cf_all, 3);
for t = 1:n_rep
    [FX(:, :, t), FY(:, :, t)] = gradient(sam_pot + cf_all(:, :, t), HX, HY);
end
R = 6.5*1e6;
E_theta = -FX/R;
E_phi = zeros(size(cf_all));
for t = 1:n_rep
    E_phi(:, :, t) = -FY(:, :, t)./(R*sin(lat_grd_mat));
end
E_theta_sam = -FX_sam/R;
E_phi_sam = -FY_sam./(R*sin(lat_grd_mat));

% load conductance data
ped_cond_all = hdf5read('OP_2012-02-29_000000.h5', '/Pedersen_Conductance');
ped_cond_all = reshape(ped_cond_all, 41, 180, size(ped_cond_all, 2));
% discard the first data point
ped_cond = fliplr(ped_cond_all(:, :, time_point + 1)');

% compute energy
energy = zeros(size(cf_all));
for t = 1:n_rep
    energy(:, :, t) = (E_theta(:, :, t).^2+E_phi(:, :, t).^2).*ped_cond;
end
energy_sam = (E_theta_sam.^2+E_phi_sam.^2).*ped_cond;

% drop first column (Inf)
energy = energy(:, 2:end, :);
energy_sam = energy_sam(:, 2:end);
lat_grd_mat = lat_grd_mat(:, 2:end);
lon_grd_mat = lon_grd_mat(:, 2:end);

phi_rot = lon_grd_mat + pi/2;
[x, y] = pol2cart(phi_rot, lat_grd_mat / pi * 180);
energy_to_plot = energy(:, :, 1);
vmag = linspace(min(energy_to_plot(:)), max(energy_to_plot(:)), 10);
figure
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, energy_to_plot, vmag);
vmag = linspace(min(energy_sam(:)), max(energy_sam(:)), 10);
figure
mypolar([0 2*pi], [0 max(lat_grd_mat(:))/ pi * 180], x, y, energy_sam, vmag);

% compute the area of each latitudinal band
% first column with Inf has already been dropped
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
    % 90 - x to convert colatitude to latitude
    area_theta(i) = areaquad(90-theta_lower/pi*180, -180, 90-theta_upper/pi*180, 180);
end

% compute integrated energy
tot_area = 4*pi*R^2;
% dividing by 1e9 is due to the unit giga
% since the areaquad function only gives a fraction of the unit sphere's area ranging from 0 to 1,
% we need to multiply back the total area
int_energy = zeros(n_rep, 1);
for t = 1:n_rep
    int_energy(t) = sum(mean(energy(:, :, t), 1).*area_theta)*tot_area/1e9;
end
int_energy_mean = mean(int_energy)
int_energy_sam = sum(mean(energy_sam, 1).*area_theta)*tot_area/1e9
