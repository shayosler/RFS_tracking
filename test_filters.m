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
true_pd0 = 0.2; % Probability of detecting a "clutter" generator
true_pd1 = .99; % Probability of detecting a target
true_ps0 = 0.9; % Probability of clutter generator survival
true_ps1 = .9;  % Probability of target survival

% Sensor
sensor = RFS.sim.Sonar_RB;
sensor.fov = 90;
sensor.range = 40;
sensor.sigma_range = .25;
sensor.sigma_bearing = 1;
sensor.lambda = 10;
sensor.pd = true_pd1;

% Measurements
% True R matrix for a measurement directly ahead
% TODO: will need to be rotated for other measurements
R_0deg_true = [sensor.sigma_range 0; 0 sensor.range * sind(sensor.sigma_bearing)];

% Clutter
A_fov = (pi * sensor.range^2) * (sensor.fov / 360);

%% Define environment/targets
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
    tgt_edot = .25;
    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
    targets(k).X = x0;
end

% Plot map
map_fig = figure;
h_map = RFS.utils.plot_targets(targets, 'b*'); 
hold on
title 'Target Locations'
set(gca, 'Fontsize', 18)
axis equal;
axis(bounds)
delete(h_map)

%% Some models/parameters that will be used in multiple filters

% Dynamics model
n_states = 2;
model_F = eye(n_states);
model_Q = 5 * eye(n_states);

% Detection/survival probabilities
model_pd0 = true_pd0;
model_pd1 = true_pd1;
model_ps0 = true_ps0;
model_ps1 = true_ps1;

% Target birth model
%bm = load('./models/birth_model_100.mat');
bm = load('./models/birth_model_200_uniform_rb.mat');
birth_rate = 0.01; % Expected rate of new births
birth_gmrfs = birth_rate .* bm.gamma;

% Sensor/measurement model
sigma_rb = [.25 0;   % range
            0   1];   % bearing

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = sensor.range * sind(sigma_rb(2, 2));
R = R_0deg_true;
H = eye(n_states);

%% l-CPHD Filter
% Filter parameters 
lcphd_params = RFS.CPHD.cphd_params();
lcphd_params.Nmax = 100;      % Maximum number of tracked objects
lcphd_params.U = 4;           % Threshold for merging components during pruning
lcphd_params.T = 1e-5;        % Threshold for discarding components during pruning
lcphd_params.Jmax = 100;      % Maximum number of components to keep during pruning
lcphd_params.w_min = 0.5;     % Minimum weight to track

% l-CPHD system model 
lcphd_model = RFS.CPHD.cphd_model();

% Target dynamics
lcphd_model.F = model_F;
lcphd_model.Q = model_Q;

% Detection and survival probabilities
lcphd_model.pd0 = model_pd0;    % Probability of detecting a "clutter" object
lcphd_model.pd1 = model_pd1;    % Probability of detecting a target
lcphd_model.ps0 = model_ps0;    % Probability of clutter survival
lcphd_model.ps1 = model_ps1;    % Probability of target survival

% Birth model
lcphd_model.gamma1 = birth_rate .* bm.gamma;

% clutter
lcphd_model.Ngamma0 = 5;        % Mean clutter birth rate, lower -> smaller lambda_hat (or maybe just takes longer to converge?)
lcphd_model.kappa = 1 / A_fov;    % Clutter is equally likely anywhere

%% l-pD-CPHD Filter
lpdcphf_model = RFS.LPDCPHD.lpdcphd_model();

% Target dynamics
lpdcphd_model.F = model_F;
lpdcphd_model.Q = model_Q;

% Detection and survival probabilities
lpdcphd_model.pd0 = model_pd0;    % Probability of detecting a "clutter" object
lpdcphd_model.pd1 = model_pd1;    % Probability of detecting a target
lpdcphd_model.ps0 = model_ps0;    % Probability of clutter survival
lpdcphd_model.ps1 = model_ps1;    % Probability of target survival

% Target birth model
gamma_s = ones(size(birth_gmrfs.w));
gamma_t = gamma_s;
lpdcphd_model.gamma1 = RFS.utils.BGMRFS(birth_gmrfs.w, birth_gmrfs.m, birth_gmrfs.P, gamma_s, gamma_t);

% Clutter birth model
% TODO: put in something real for this
lpdcphd_model.gamma0 = RFS.utils.BMRFS(birth_gmrfs.w, gamma_s, gamma_t);

% Dilation constant for beta distributions
lpdcphd_model.kB = 1;


%% Simulate
t_total = 300;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions and storage variables
x = zeros(sim_steps, 12);
obs = cell(sim_steps, 1);
true_obs = cell(sim_steps, 1);
lambda_true = zeros(sim_steps, 1);
lambda_true(1) = sensor.lambda;

% Initial conditions and storage for l-CPHD
lcphd_Xhat = cell(sim_steps, 1);
lambda_hat = zeros(sim_steps, 1);
lambda_hat(1) = 1;
lcphd_states(sim_steps, 1) = RFS.CPHD.cphd_state();
lcphd_states(1).v = RFS.utils.GMRFS();
lcphd_states(1).N0 = 0;
lcphd_states(1).N1 = 0;
lcphd_states(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown

% Initial conditions and storage for GMPHD
gmphd_v(sim_steps) = RFS.utils.GMRFS();
gmphd_N = zeros(sim_steps, 1);
gmphd_Xhat = cell(sim_steps, 1);

% Initial conditions and storage for l-pd-CPHD
lpdcphd_Xhat = cell(sim_steps, 1);
lpdcphd_states(sim_steps, 1) = RFS.LPDCPHD.lpdcphd_state();
lpdcphd_states(1).v0 = RFS.utils.BMRFS();
lpdcphd_states(1).v1 = RFS.utils.BGMRFS();
lpdcphd_states(1).N0 = 0;
lpdcphd_states(1).N1 = 0;
lpdcphd_states(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown

% Figures for plotting
lcphd_figures = struct();
lcphd_figures.v_fig = figure;
lcphd_figures.n_fig = figure;
lcphd_figures.map_fig = figure;
lcpdh_figures.rho_fig = figure;
lcpdh_figures.lambda_fig = figure;

gmphd_figures = struct();
gmphd_figures.v_fig = figure;
gmphd_figures.n_fig = figure;
gmphd_figures.map_fig = figure;
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
    if k < 100
        sensor.lambda = 7;
    elseif k < 200
        sensor.lambda = 12;
    else
        sensor.lambda = 2;
    end
    lambda_true(k) = sensor.lambda;
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
    [lcphd_states(k), lcphd_Xhat{k}, lambda_hat(k)] = RFS.CPHD.lcphd_filter(lcphd_states(k-1), measurement, lcphd_model, lcphd_params);
    [lpdcphd_states(k), lpdcphd_Xhat{k}, ~, ~] = RFS.CPHD.lpdcphd_filter(lpdcphd_states(k-1), measurement, lpdcphd_model, lcphd_params);
    [gmphd_v(k), gmphd_N(k), gmphd_Xhat{k}] = RFS.GMPHD.phd_filter(gmphd_v(k-1), ...
        gmphd_N(k-1), ...
        lcphd_model.F, ...
        lcphd_model.Q, ...
        lcphd_model.ps1, ...
        lcphd_model.pd1, ...
        lcphd_model.gamma1, ...
        obs{k}', ...
        H, ...
        R, ...
        lcphd_model.kappa, ...
        lcphd_params.U, ...
        lcphd_params.T, ...
        lcphd_params.Jmax, ...
        lcphd_params.w_min);

    %% Plots
    if mod(k, 10) == 0 && false
        delete(rfs_handles);
        lcphd_handles = RFS.utils.rfs_tracker_plots(lcphd_figures, 'l-CPHD', k, lcphd_states(k).v, [n e psi], lcphd_Xhat, sensor, obs, true_obs, targets, bounds);
        gmphd_handles = RFS.utils.rfs_tracker_plots(gmphd_figures, 'GMPHD', k, gmphd_v(k), [n e psi], gmphd_Xhat, sensor, obs, true_obs, targets, bounds);
        rfs_handles = [lcphd_handles gmphd_handles];
    end % for k = 2:sim_steps
end

%% Final plots
% General RFS plots
delete(rfs_handles);
lcphd_handles = RFS.utils.rfs_tracker_plots(lcphd_figures, 'l-CPHD', k, lcphd_states(k).v, [n e psi], lcphd_Xhat, sensor, obs, true_obs, targets, bounds);
gmphd_handles = RFS.utils.rfs_tracker_plots(gmphd_figures, 'GMPHD', k, gmphd_v(k), [n e psi], gmphd_Xhat, sensor, obs, true_obs, targets, bounds);

% l-CPHD specific plots
RFS.CPHD.lcphd_plots(lcpdh_figures, k, targets, lcphd_Xhat, lambda_hat, lambda_true);
