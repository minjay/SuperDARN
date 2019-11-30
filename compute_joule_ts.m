clear
clc

% compute Joule heating rate of the entire field (sum of all components)

% load sam pot (provided by Xueling)
load('2012-02-29_00_00-2012-02-29_04_00_sam_pot_data.mat')

% load conductance data
ped_cond_all = hdf5read('OP_2012-02-29_000000.h5', '/Pedersen_Conductance');
ped_cond_all = reshape(ped_cond_all, 41, 180, size(ped_cond_all, 2));

% convert latitude to colatitude (90 - x).
lat_grd_mat = fliplr((90 - reshape(lat_grd, [41 180])') / 180 * pi);
lon_grd_mat = fliplr(reshape(lon_grd, [41 180])' / 180 * pi);

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
% conditional expectation E(Z(s)|Z_1,...,Z_n)
% conditional distribution Z(s)|Z_1,...,Z_n
cf_all = reshape(std_grid_vec .* Ac, [size(lon_grd_mat) size(Ac, 2)]);
% convert to Volt
cf_all = cf_all / 1e6;

%% Compute Joule heating rate

n_t = size(pot, 1);
int_energy_mean = zeros(n_t, 1);
int_energy_sam = zeros(n_t, 1);
for time_point = 1:n_t
    % potential field from SAM
    sam_pot = fliplr(reshape(pot(time_point, :), [41 180])');

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
    lat_grd_mat_drop = lat_grd_mat(:, 2:end);
    lon_grd_mat_drop = lon_grd_mat(:, 2:end);

    % compute the area of each latitudinal band
    % first column with Inf has already been dropped
    theta_one = lat_grd_mat_drop(1, :);
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
    int_energy_mean(time_point) = mean(int_energy);
    int_energy_sam(time_point) = sum(mean(energy_sam, 1).*area_theta)*tot_area/1e9;
end

%% load AE index
activity_data = readtable('OMNI_Feb29_2012.csv');

% downsample data to every 2 mins
activity_data = activity_data(1:2:size(activity_data, 1), :);

% drop hr = 0, min = 0
activity_data = activity_data(2:size(activity_data, 1), :);

% select the first 4 hours
activity_data = activity_data(1:120, :);

% select column 6 (AE-index, [nT])
ae_index = activity_data.x6;

%% plot
figure
yyaxis left
int_energy_sam_ts = timeseries(int_energy_sam, ts);
plot(int_energy_sam_ts, 'k-o')
hold on
int_energy_mean_ts = timeseries(int_energy_mean, ts);
plot(int_energy_mean_ts, '-o')
legend('SAM', 'SAM + Needlet')
xlabel('Time [UT]')
ylabel('Integrated Joule heating rate [GW]')

yyaxis right
ae_index_ts = timeseries(ae_index, ts);
plot(ae_index_ts, '-o')
ylabel('AE index [nT]')
