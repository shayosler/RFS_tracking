%% clear
clc
clear
close all

%% Set up
%seed = 52490;
seed = 1;
rng(seed);

%% True system
sim_dt = 1;

% Target dynamics: constant velocity model
% x = [n ndot e edot]'
F = [1  sim_dt  0  0;
     0  1       0  0;
     0  0       1  sim_dt;
     0  0       0  1]; % Constant velocity
Q = diag([.1 0 .1 0]) * 0;

% Detection and survival probabilities
pd0 = 0.9;  % Probability of detecting a "clutter" object
pd1 = 0.9;  % Probability of detecting a target
ps0 = 0.1;  % Probability of clutter survival
ps1 = 0.99; % Probability of target survival

% Sensor
sensor = RFS.sim.Sonar_RB;
sensor.fov = 90;
sensor.range = 40;
sensor.sigma_range = .25 * 1e-6;
sensor.sigma_bearing = 5 * 1e-6;
sensor.lambda = 0;
sensor.pd = pd0;

% Measurements
% True R matrix for a measurement directly ahead
% TODO: will need to be rotated for other measurements
R_0deg_true = [sensor.sigma_range 0; 0 sensor.range * sind(sensor.sigma_bearing)];

% Clutter
A_fov = (pi * sensor.range^2) * (sensor.fov / 360);

%% Define environment
% Bounds for northing and easting
min_n = 0;
max_n = 100;
min_e = -30;
max_e = 30;
bounds = [min_e max_e min_n max_n];
northings = min_n:.1:max_n;
eastings = min_e-10:.1:max_e;

% Initialize targets. 
n_tgt = 10;  
targets(n_tgt) = RFS.sim.Target_2D;
spacing = 4;
for k = 1:n_tgt
    targets(k).F = F;
    targets(k).Q = Q;
    tgt_e = min_e;
    tgt_n = 10 + 2.5 * k;
    tgt_ndot = 0;
    tgt_edot = 1;
    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
    targets(k).X = x0;
end

% Plot map
map_fig = figure;
h_map = RFS.utils.plot_targets(targets, 'b*'); %  plot(map(:, 2), map(:, 1), 'b*');
hold on
title 'Target Locations'
set(gca, 'Fontsize', 18)
axis equal;
axis(bounds)
delete(h_map)

%% l-CPHD Filtering parameters
lcphd_params = RFS.CPHD.cphd_params();
lcphd_params.Nmax = 100;      % Maximum number of tracked objects
lcphd_params.U = 4;           % Threshold for merging components during pruning
lcphd_params.T = 1e-5;        % Threshold for discarding components during pruning
lcphd_params.Jmax = 100;      % Maximum number of components to keep during pruning
lcphd_params.w_min = 0.5;     % Minimum weight to track

%% System model used for filtering
n_states = 2;
model = RFS.CPHD.cphd_model();

% Target dynamics
model.F = eye(n_states);
model.Q = .5 * eye(n_states);

% Detection and survival probabilities
model.pd0 = pd0;    % Probability of detecting a "clutter" object
model.pd1 = pd1;    % Probability of detecting a target
model.ps0 = ps0;    % Probability of clutter survival
model.ps1 = ps1;    % Probability of target survival

% Target birth model
bm = load('./models/birth_model_100.mat');
birth_rate = 0.01; % Expected rate of new births
model.gamma1 = birth_rate .* bm.gamma;

% Sensor/measurement model
sigma_rb = [.25 0;   % range
            0   1];   % bearing

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = sensor.range * sind(sigma_rb(2, 2));
H = eye(n_states);

% clutter
model.Ngamma0 = 0.5;        % Mean clutter birth rate
model.kappa = 1 / A_fov;    % Clutter is equally likely anywhere

%% Simulate
t_total = 100;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions
x = zeros(sim_steps, 12);
Xhat = cell(sim_steps, 1);
obs = cell(sim_steps, 1);
true_obs = cell(sim_steps, 1);
lambda_hat = zeros(sim_steps, 1);
lambda_hat(1) = 1;
lcphd_states(sim_steps, 1) = RFS.CPHD.cphd_state;
lcphd_states(1).v = RFS.utils.GMRFS();
lcphd_states(1).N0 = 0;
lcphd_states(1).N1 = 0;
lcphd_states(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown

% Figures for plotting
rfs_figures = struct();
rfs_figures.v_fig = figure;
rfs_figures.n_fig = figure;
rfs_figures.map_fig = map_fig;

lcpdh_figures = struct();
lcpdh_figures.rho_fig = figure;
lcpdh_figures.lambda_fig = figure;
rfs_handles = [];
for k = 2:sim_steps
    % Update vehicle position
    x(k, :) = x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
    n = x(k, 1);
    e = x(k, 2);
    psi = x(k, 6);

    % Update target positions
    existing = [targets.t_birth] <= k & [targets.t_death] > k;
    targets(existing) = targets(existing).step();

    % Get measurements of targets
    [r, b, r_true, b_true] = sensor.measure([n, e, psi], targets(existing));   

    % Convert to absolute coordinates
    b = b + psi; 
    n_obs = n + r .* cosd(b);
    e_obs = e + r .* sind(b);
    b_true = b_true + psi;
    n_obs_true = n + r_true .* cosd(b_true);
    e_obs_true = e + r_true .* sind(b_true);

    % Store the observations
    obs{k} = zeros(length(r), 2);
    obs{k}(:, 1) = n_obs;
    obs{k}(:, 2) = e_obs;
    true_obs{k} = zeros(length(r_true), 2);
    true_obs{k}(:, 1) = n_obs_true;
    true_obs{k}(:, 2) = e_obs_true;

    % Generate measurement
    measurement = RFS.CPHD.cphd_measurement();
    measurement.Z = obs{k}';
    measurement.H = H;
    measurement.R = R;

    % Update estimates
    [lcphd_states(k), Xhat{k}, lambda_hat(k)] = RFS.CPHD.lcphd_filter(lcphd_states(k-1), measurement, model, lcphd_params);
    %[v_k_obs, N_in_k, ~] = phd_filter(v{k-1}, N_in, F, Q, ps, pd, gamma, obs{k}', H, R, kappa, U, T, Jmax, w_min);

    %% Plots
    if mod(k, 1) == 0
        delete(rfs_handles);
        rfs_handles = RFS.utils.rfs_tracker_plots(rfs_figures, 'l-CPHD', k, lcphd_states(k).v, [n e psi], Xhat, sensor, obs, true_obs, targets, bounds);
    end % for k = 2:sim_steps
end

%% Final plots
% General RFS plots
delete(rfs_handles);
RFS.utils.rfs_tracker_plots(rfs_figures, 'l-CPHD', k, lcphd_states(k).v, [n e psi], Xhat, sensor, obs, true_obs, targets, bounds);

% l-CPHD specific plots
