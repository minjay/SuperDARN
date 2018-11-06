% Main idea: Apply Voronoi tessellation to data points on the sphere.
% Sample each data point based on the probability proportional to the area
% of the corresponding Voronoi cell. The larger the area is, the higher
% probability this data point is sampled.
addpath('sphere_voronoi')

% load un-downsampled data
load('sd_data.mat')
% theta is the co-latitude after stretching
load('sd_theta.mat')
load('sd_phi.mat')

[x, y, z] = sph2cart(phi, pi / 2 - theta, 1);
n = length(x);
[ face_num, face ] = sphere_delaunay ( n, [x; y; z] );
v = voronoi_vertices ( n, [x; y; z], face_num, face );
area = voronoi_areas_direct ( n, [x; y; z], face_num, face, v );

rng(0)
% sample 2000 data points
sampled_index = sort(datasample(1:n, 2000, 'Weights', area, 'Replace', false));

theta_thin = theta(sampled_index);
phi_thin = phi(sampled_index);
x_thin = x(sampled_index);
y_thin = y(sampled_index);
z_thin = z(sampled_index);

fig = figure;
[sx, sy, sz] = sphere;
h = surf(sx, sy, sz);
set(h, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold on
scatter3(x_thin, y_thin, z_thin)
axis equal
saveas(fig, 'thinned_data', 'png')

fig = figure;
[sx, sy, sz] = sphere;
h = surf(sx, sy, sz);
set(h, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold on
scatter3(x, y, z)
axis equal
saveas(fig, 'raw_data', 'png')