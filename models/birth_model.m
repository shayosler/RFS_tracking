%%
clc
clear
close all

%% Import


%% birth model:
% new targets are most likely to appear at the leading edge of the sonar's
% field of view
deg_to_rad = pi / 180;
l_lim = -45 * deg_to_rad;
r_lim = 45 * deg_to_rad;
range = 40;

n_pts = 100;
r = linspace(0, range, n_pts);
theta = linspace(l_lim, r_lim, n_pts);

[R, THETA] = meshgrid(r, theta);

% probability only depends on range
mu = range;
sigma = range/4;
p = 2 * normpdf(R(:), range, range/4);
p = reshape(p, size(R));
[x, y, p_b] = pol2cart(THETA, R, p);

pcolor(y, x, p_b)
shading interp;
axis equal


%% Fit gaussian mixture
samples = 100000;
samples = 10000;
bearings = l_lim + (r_lim-l_lim).*rand(samples,1);
ranges = -abs(normrnd(0, sigma, samples, 1)) + range;
pts = [bearings ranges];
% ensure no points lay past the end range
[sn, se] = pol2cart(bearings, ranges);
figure 
plot(se, sn, '.');

% Fit gaussian mixture
num_gaussians = 100;
options = statset('MaxIter', 1000);
dist = fitgmdist([sn, se], num_gaussians, 'SharedCovariance', false, 'Options',options);
n = 0:.1:60;
e = -40:.1:40;

gamma = RFS.utils.GMRFS(dist.mu', dist.Sigma, dist.ComponentProportion);
figure
RFS.utils.plotgmphd(gamma, n, e);
axis equal
