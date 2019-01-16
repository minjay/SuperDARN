function [Npix, grid_points, A_1, A_2] = get_A_full(B, j_min, j_max, theta, phi, k_theta, k_phi)
%GET_A_full_SS   Computes the design matrices A_1 and A_2 when the grid is the symmetric
%spherical t-design.
%
%   [Npix, grid_points, A_1,A_2] = get_A_full_ss(B, j_min, j_max, theta, phi, k_theta, k_phi);
%
% Inputs:
%   B - the parameter
%   j_min - the minimal frequency
%   j_max - the maximal frequency
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%   k_theta - component of LOS vector
%   k_phi - component of LOS vector
%
% Outputs:
%   Npix - the number of grid points at each frequency,
%   (j_max-j_min+1)-by-1 vector
%   grid_points - the location of the grid points at each frequency in R^3,
%   (j_max-j_min+1)-by-1 cell
%   A_1 - the design matrix for psi_jk^(1), N-by-M matrix
%   A_2 - the design matrix for psi_jk^(2), N-by-M matrix
%
% Author: Caleb Miller, 2018 building upon Minjie Fan, 2015

%%% Constructing psi_jk^(1) -> A_1

% Get 1/sin(theta) \pdiff{(\psi_{jk})}{\phi} call it A_phi
    [~,~,A_phi] = get_A_part_ss('phi', B, j_min, j_max, theta, phi);
% Get \pdiff{(psi_{jk})}{\theta} call it A_theta
    [~,~,A_theta] = get_A_part_ss('theta', B, j_min, j_max, theta, phi);  
% Construct A_1 

    A_1 = k_theta.*sin(theta)./sin(theta/4).*A_phi-4*k_phi.*A_theta;

%%% Constructing psi_jk^(2) -> A_2

% Get \psi_{jk}, Npix, grid_points
    [Npix, grid_points, A] = get_A_ss(B, j_min, j_max, theta, phi);
% Construct A_2
    A_2 = 4*k_phi.*A;














