function [h] = plot_fov(x, psi, r, fov, varargin)
%plot_fov(x, psi, r, fov) Plot the sonar's field of view
%   Inputs:
%       x       Sonar position [n, e]
%       psi     Sonar heading, degrees
%       r       Sonar end range
%       theta   Included angle in the sonar's field of view, degrees
%
%   Outputs:
%       h       Handle of plotted objects

l_lim = psi - fov/2;
r_lim = psi + fov/2;
theta = (l_lim:.1:r_lim)';
arc = [x(2) + r .* sind(theta) x(1) + r .* cosd(theta)];
r_pts = (0:.1:r)';
l_line = [x(2) + r_pts.* sind(l_lim) x(1) + r_pts .* cosd(l_lim)];
r_line = [x(2) + r_pts.* sind(r_lim) x(1) + r_pts .* cosd(r_lim)];

h_a = plot(arc(:, 1), arc(:, 2), varargin{:});
hold on
h_l = plot(l_line(:, 1), l_line(:, 2), varargin{:});
h_r = plot(r_line(:, 1), r_line(:, 2), varargin{:});

h = [h_a, h_l, h_r];
end