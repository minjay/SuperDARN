clear 
clc

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
% pre-computed by the above two lines of code and saved
load('A_mat.mat')

% for R code
%save('all_theta.mat', 'all_theta')

[b_mat, ~] = bspline_basismatrix(4, knots, all_theta);
b_mat(:, 1) = 1;

post_samples_eta = post_samples.eta;
post_samples_c = squeeze(post_samples.c);

std_vec = exp(b_mat * post_samples_eta);

load('deriv_all_theta.mat')
b_mat_deriv = bS;
b_mat_deriv(:, 1) = 0;
std_vec_deriv = (b_mat_deriv * post_samples_eta) .* std_vec;

N = length(all_theta);

A1c = A_1 * post_samples_c;
A2c = A_2 * post_samples_c;

DAc = mean(std_vec .* A1c - std_vec_deriv .* A2c, 2);

Y_fitted = DAc ./ factor;

Y = Vlos_res';

load('voroni_sphere_index.mat')

% out-of-sample predictions
figure
unsampled_index = setdiff(1:N, sampled_index);
plot(Y(unsampled_index), Y_fitted(unsampled_index), 'o')
axis('equal')
axis('tight')
hline = refline(1, 0);
hline.Color = 'r';
hline.LineWidth = 2;
hline.LineStyle = '--';
xlabel('Observed','Fontsize',12)
ylabel('Fitted','Fontsize',12)

% compute Pearson correlation coefficient
corr_out_sample = corrcoef(Y(unsampled_index), Y_fitted(unsampled_index));
corr_out_sample = corr_out_sample(1, 2);
text(-600, 800, strcat('Pearson corr coef=', sprintf('%.2f', corr_out_sample)),'FontName','times','Fontsize',16);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];

% in-sample predictions
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

% compute Pearson correlation coefficient
corr_in_sample = corrcoef(Y(sampled_index), Y_fitted(sampled_index));
corr_in_sample = corr_in_sample(1, 2);
text(-800, 800, strcat('Pearson corr coef=', sprintf('%.2f', corr_in_sample)),'FontName','times','Fontsize',16);

tau_hat = mean(1 ./ sqrt(post_samples.tau_sq_inv)) ./ factor;

errs = Y(sampled_index) - Y_fitted(sampled_index);
sqrt(mean(errs.^2))
sqrt(mean(Y(sampled_index).^2))

errs = Y(unsampled_index) - Y_fitted(unsampled_index);
sqrt(mean(errs.^2))
sqrt(mean(Y(unsampled_index).^2))
