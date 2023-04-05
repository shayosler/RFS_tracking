%% clear
clc
clear
close all

%% Set up
%seed = 52490;
seed = 1;
rng(seed);

%% Sensor
fov = 90;
range = 40;
sigma_rb = [.25 0;   % range
            0  10];   % bearing

pd0 = .6;   % probability of detecting a target
pd1 = 0.6;  % probability of detecting a "clutter" object

% clutter: poisson process
A_fov = (pi * range^2) * (fov / 360);
lambda = 2; % Expected number of clutter returns
kappa = lambda / A_fov;

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = range * sind(sigma_rb(2, 2));


%% Define targets

% Simple: line of targets heading straight at sensor
n_tgt = 1;  % Start with a single target
spacing = 4; 
v_tgt = [1, 0];
tgt_e = zeros(n_tgt, 1);
tgt_n = (35:spacing:71)';
map = [tgt_n tgt_e ones(n_tgt, 1) * pd0];

% Plot map
map_fig = figure;
plot(map(:, 2), map(:, 1), 'b*')
hold on
title 'Target Locations'
set(gca, 'Fontsize', 18)
axis equal;
axis([-30, 30, 0, 100])

%% Filtering parameters
n_states = 2;
F = eye(n_states);
N = 0;
ps = 1.0;
pd = [pd0 pd1];
bm = load('birth_model_100.mat');
birth_rate = 1.0; % Expected rate of new births
gamma = birth_rate .* bm.gamma;
H = eye(n_states);
Q = 0.1 * eye(2);
U = 4;
T = 1e-5;
Jmax = 100;
w_min = 0.5;

v_fig = figure;

%% Simulate
sim_dt = 0.1;
t_total = 400;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions
x = zeros(sim_steps, 12);
obs = cell(sim_steps, 1);
N = zeros(sim_steps, 1);
v = cell(sim_steps, 1);
v{1} = GMRFS();
j = 1;
pass = 1;
all_obs = [];
for k = 2:sim_steps
    tk = t(j);
    j = j+1;

    % Update position
    x(k, :) = x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
    n = x(k, 1);
    e = x(k, 2);
    psi = x(k, 6);

    % Get relative positions of simuated observations
    [r, b, r_true, b_true] = simulate_obs([n, e, psi], map, range, fov, lambda, sigma_rb);

    % Convert to absolute coordinates
    b = b + psi; 
    n_obs = n + r .* cosd(b);
    e_obs = e + r .* sind(b);

    % Store the observations
    obs{k} = zeros(length(r), 2);
    obs{k}(:, 1) = n_obs;
    obs{k}(:, 2) = e_obs;
    all_obs = [all_obs; obs{k}];

    [v_k_obs, N_in_k, ~] = phd_filter(v{k-1}, N_in, F, Q, ps, pd, gamma, obs{k}', H, R, kappa, U, T, Jmax, w_min);

end
