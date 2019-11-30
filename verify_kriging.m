% load un-downsampled data
load('sd_data.mat')
% theta is the co-latitude after stretching
load('sd_theta.mat')
load('sd_phi.mat')
load('sd_k_theta.mat');
load('sd_k_phi.mat');

B = 2;
j_min = 2;
j_max = 3;

knots = [0 0 0 0 0.5 1 1 1 1]*pi;

% to avoid override, rename
all_theta = theta';
all_phi = phi';
k_theta = k_theta';
k_phi = k_phi';

load('vs_0_2000_tweak_252.mat')

Bf = compute_mag_field(glat, glon);
radius = 6356.752 + 300;
factor = Bf * radius;

%[Npix, ~, A_1, A_2] = get_A_full(B, j_min, j_max, all_theta, all_phi, k_theta, k_phi);
%save('A_mat.mat', 'A_1', 'A_2')
load('A_mat.mat')

% for R code
%save('all_theta.mat', 'all_theta')

[b_mat, ~] = bspline_basismatrix(4, knots, all_theta);
b_mat(:, 1) = 1;
eta_hat = mean(post_samples.eta, 2);
std_vec = exp(b_mat * eta_hat);

load('deriv_all_theta.mat')
b_mat_deriv = bS;
b_mat_deriv(:, 1) = 0;
std_vec_deriv = b_mat_deriv * eta_hat .* std_vec;

N = length(all_theta);
M = size(A_1, 2);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i) * A_1(i, :) - std_vec_deriv(i) * A_2(i,:);
end

c_hat = mean(squeeze(post_samples.c), 2);

Y_fitted = DA * c_hat ./ factor;

Y = Vlos_res';

load('voroni_sphere_index.mat')

figure
unsampled_index = setdiff(1:N, sampled_index);
plot(Y(unsampled_index), Y_fitted(unsampled_index), 'o')
axis('equal')
axis('tight')
hline = refline(1, 0);
hline.Color = 'r';
hline.LineWidth = 2;
hline.LineStyle = '--';
xlabel('Observed')
ylabel('Fitted')

figure
plot(Y(sampled_index), Y_fitted(sampled_index), 'o')
axis('equal')
axis('tight')
hline = refline(1, 0);
hline.Color = 'r';
hline.LineWidth = 2;
hline.LineStyle = '--';
xlabel('Observed')
ylabel('Fitted')

tau_hat = mean(1 ./ sqrt(post_samples.tau_sq_inv)) ./ factor;

errs = Y(sampled_index) - Y_fitted(sampled_index);
sqrt(mean(errs.^2))
sqrt(mean(Y(sampled_index).^2))

errs = Y(unsampled_index) - Y_fitted(unsampled_index);
sqrt(mean(errs.^2))
sqrt(mean(Y(unsampled_index).^2))
