function Bf = compute_mag_field(glat, glon)
% Compute magnetic field at a number of locations.
% 
% Args:
%   glat: geodetic latitudes.
%   glon: geodetic longitudes.
%
% Returns:
%   Magnetic field at the provided geodetic latitudes and longitudes.

coef2015 = [-29442.0 -1501.0 4797.1 -2445.1 3012.9 -2845.6 1676.7 -641.9 1350.7 -2352.3 -115.3 1225.6 244.9 582.0 -538.4 907.6 813.7 283.3 120.4 -188.7 -334.9 180.9 70.4 -329.5 -232.6 360.1 47.3 192.4 197.0 -140.9 -119.3 -157.5 16.0 4.1 100.2 70.0 67.7 -20.8 72.7 33.2 -129.9 58.9 -28.9 -66.7 13.2 7.3 -70.9 62.6 81.6 -76.1 -54.1 -6.8 -19.5 51.8 5.7 15.0 24.4 9.4 3.4 -2.8 -27.4 6.8 -2.2 24.2 8.8 10.1 -16.9 -18.3 -3.2 13.3 -20.6 -14.6 13.4 16.2 11.7 5.7 -15.9 -9.1 -2.0 2.1 5.4 8.8 -21.6 3.1 10.8 -3.3 11.8 0.7 -6.8 -13.3 -6.9 -0.1 7.8 8.7 1.0 -9.1 -4.0 -10.5 8.4 -1.9 -6.3 3.2 0.1 -0.4 0.5 4.6 -0.5 4.4 1.8 -7.9 -0.7 -0.6 2.1 -4.2 2.4 -2.8 -1.8 -1.2 -3.6 -8.7 3.1 -1.5 -0.1 -2.3 2.0 2.0 -0.7 -0.8 -1.1 0.6 0.8 -0.7 -0.2 0.2 -2.2 1.7 -1.4 -0.2 -2.5 0.4 -2.0 3.5 -2.4 -1.9 -0.2 -1.1 0.4 0.4 1.2 1.9 -0.8 -2.2 0.9 0.3 0.1 0.7 0.5 -0.1 -0.3 0.3 -0.4 0.2 0.2 -0.9 -0.9 -0.1 0.0 0.7 0.0 -0.9 -0.9 0.4 0.4 0.5 1.6 -0.5 -0.5 1.0 -1.2 -0.2 -0.1 0.8 0.4 -0.1 -0.1 0.3 0.4 0.1 0.5 0.5 -0.3 -0.4 -0.4 -0.3 -0.8];
coef2010 = [-29496.57	-1586.42	4944.26	-2396.06	3026.34	-2708.54	1668.17	-575.73	1339.85	-2326.54	-160.4	1232.1	251.75	633.73	-537.03	912.66	808.97	286.48	166.58	-211.03	-356.83	164.46	89.4	-309.72	-230.87	357.29	44.58	200.26	189.01	-141.05	-118.06	-163.17	-0.01	-8.03	101.04	72.78	68.69	-20.9	75.92	44.18	-141.4	61.54	-22.83	-66.26	13.1	3.02	-78.09	55.4	80.44	-75	-57.8	-4.55	-21.2	45.24	6.54	14	24.96	10.46	7.03	1.64	-27.61	4.92	-3.28	24.41	8.21	10.84	-14.5	-20.03	-5.59	11.83	-19.34	-17.41	11.61	16.71	10.85	6.96	-14.05	-10.74	-3.54	1.64	5.5	9.45	-20.54	3.45	11.51	-5.27	12.75	3.13	-7.14	-12.38	-7.42	-0.76	7.97	8.43	2.14	-8.42	-6.08	-10.08	7.01	-1.94	-6.24	2.73	0.89	-0.1	-1.07	4.71	-0.16	4.44	2.45	-7.22	-0.33	-0.96	2.13	-3.95	3.09	-1.99	-1.03	-1.97	-2.8	-8.31	3.05	-1.48	0.13	-2.03	1.67	1.65	-0.66	-0.51	-1.76	0.54	0.85	-0.79	-0.39	0.37	-2.51	1.79	-1.27	0.12	-2.11	0.75	-1.94	3.75	-1.86	-2.12	-0.21	-0.87	0.3	0.27	1.04	2.13	-0.63	-2.49	0.95	0.49	-0.11	0.59	0.52	0	-0.39	0.13	-0.37	0.27	0.21	-0.86	-0.77	-0.23	0.04	0.87	-0.09	-0.89	-0.87	0.31	0.3	0.42	1.66	-0.45	-0.59	1.08	-1.14	-0.31	-0.07	0.78	0.54	-0.18	0.1	0.38	0.49	0.02	0.44	0.42	-0.25	-0.26	-0.53	-0.26	-0.79];
% coef from linear interpolation
coef2012 = coef2010 + 2 / 5 * (coef2015 - coef2010);
% altitude
alt = 300;

% type of coordinate
coord = 'geodetic';

len = length(glat);
Bf = zeros(len, 1);

for i = 1:len
  [~, ~ , Bf(i, 1)] = igrf(coef2012, glat(i), glon(i), alt, coord);
end

end