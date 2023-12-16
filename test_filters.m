%% clear
clc
clear
close all

%% Set up
%seed = 52490;
seed = 1;
rng(seed);

% Select which filters to run
gmphd = true;
lcphd = true;
lpdcphd = true;
lmb = true;
almb = true;

if almb && ~lpdcphd
    error('The adaptive LMB filter requires also runnning the l-pd-CPHD filter')
end

% Simulation params
n_runs = 10;
t_total = 100;
sim_dt = 1;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% OSPA paramaters
ospa_c = 10;
ospa_p = 1;

%% True system

% Target dynamics: constant velocity model
% x = [n ndot e edot]'
F = [1  sim_dt  0  0;
     0  1       0  0;
     0  0       1  sim_dt;
     0  0       0  1]; % Constant velocity
Q = diag([.1 0 .1 0]) * .05;

% Detection and survival probabilities
true_pd0 = 0.2; % Probability of detecting a "clutter" generator
true_pd1 = .95; % Probability of detecting a target
true_ps0 = 0.9; % Probability of clutter generator survival
true_ps1 = 0.9; % Probability of target survival

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
max_n = 45;
min_e = -30;
max_e = 30;
bounds = [min_e max_e min_n max_n];
northings = min_n:.1:max_n;
eastings = min_e-10:.1:max_e;

% Number of targets. 
n_tgt = 10;  

% Plot map
map_fig = figure;

%% Some models/parameters that will be used in multiple filters
% Dynamics model
% Constant velocity
F_cv = [1  sim_dt  0  0;
        0  1       0  0;
        0  0       1  sim_dt;
        0  0       0  1];
H_cv = [1 0 0 0;
        0 0 1 0];

% Static
F_static = eye(2);
H_static = [1 0;
            0 1];

% Set model to use
model_F = F_static;
n_states = size(model_F, 2);
model_Q = 5 * eye(n_states);

% Detection/survival probabilities
model_pd0 = true_pd0;
model_pd1 = .5 * true_pd1;
model_ps0 = true_ps0;
model_ps1 = true_ps1;

% Target birth model
%bm = load('./models/birth_model_100.mat');
uniform_bm = load('./models/birth_model_200_uniform_rb.mat');
uniform_bm = uniform_bm.gamma;
centerline_bm = RFS.utils.GMRFS([23; 0], [35 0; 0 2], 1); % One component along center
left_edge_bm = RFS.utils.transform_gmrfs(centerline_bm, 0, 0, -pi/4);
center_bm = RFS.utils.GMRFS([20; 0], [40 0; 0 80], 1);

left_edge_bm_cv = RFS.utils.GMRFS([16.25; 0; -16.25; 0], [18.5 0 -16.5 0; 0 1 0 0; -16.5 0 18.5 0; 0 0 0 1], 1);
center_bm_cv = RFS.utils.GMRFS([23; 0; 0; 0], [15 0 0 0; 0 1 0 0; 0 0 15 0; 0 0 0 1], 1);

birth_rate = 0.01; % Expected rate of new births
birth_gmrfs = birth_rate .* center_bm; % uniform_bm;
birth_fig = figure;
h_birth = RFS.utils.plotgmphd(birth_gmrfs, min_n:.1:max_n, min_e:.1:max_e);
hold on
h_fov = sensor.plot_fov(0, 0, 0, 'r');
title 'Birth Model'

% Sensor/measurement model
n_meas = 2;
sigma_rb = [.25 0;   % range
            0   1];   % bearing

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = sensor.range * sind(sigma_rb(2, 2));
R = R_0deg_true;
H = H_static;

% Mean clutter birth rate
% lower -> smaller lambda_hat (or maybe just takes longer to converge?)
Ngamma0 = 5;

%% GMPHD Filter

% Initial conditions and storage for GMPHD
gmphd_results.label = 'GMPHD';
gmphd_results.est(n_runs) = struct();
gmphd_kappa = sensor.lambda / A_fov;
for i = 1:n_runs
    gmphd_results.est(i).v(sim_steps) = RFS.utils.GMRFS();
    gmphd_results.est(i).N = zeros(sim_steps, 1);
    gmphd_results.est(i).Xhat = cell(sim_steps, 1);
    gmphd_results.est(i).ospa = ones(sim_steps, 1) * ospa_c;
    gmphd_results.est(i).ospa_l = ones(sim_steps, 1) * ospa_c;
    gmphd_results.est(i).ospa_n = ones(sim_steps, 1) * ospa_c;
end
% Figures for GMPHD plots
if gmphd
    gmphd_results.figures = struct();
    gmphd_results.figures.v_fig = figure;
    gmphd_results.figures.n_fig = figure;
    gmphd_results.figures.map_fig = figure;
    gmphd_results.figures.ospa_fig = figure;
end

%% l-CPHD Filter
% Filter parameters 
lcphd_params = RFS.LCPHD.cphd_params();
lcphd_params.Nmax = 100;      % Maximum number of tracked objects
lcphd_params.U = 4;           % Threshold for merging components during pruning
lcphd_params.T = 1e-5;        % Threshold for discarding components during pruning
lcphd_params.Jmax = 100;      % Maximum number of components to keep during pruning
lcphd_params.w_min = 0.5;     % Minimum weight to track

% l-CPHD system model 
lcphd_model = RFS.LCPHD.cphd_model();

% Target dynamics
lcphd_model.F = model_F;
lcphd_model.Q = model_Q;

% Detection and survival probabilities
lcphd_model.ps0 = model_ps0;    % Probability of clutter survival
lcphd_model.ps1 = model_ps1;    % Probability of target survival
lcphd_model.pd0 = model_pd0;    % Probability of detecting a "clutter" object
lcphd_model.pd1 = model_pd1;    % Probability of detecting a target

% Birth model
lcphd_model.gamma1 = birth_gmrfs;

% clutter
lcphd_model.Ngamma0 = Ngamma0;  % Mean clutter birth rate, 
lcphd_model.kappa = 1 / A_fov;  % Clutter is equally likely anywhere

% Initial conditions
lcphd_states0(sim_steps, 1) = RFS.LCPHD.cphd_state();
lcphd_states0(1).v = RFS.utils.GMRFS();
lcphd_states0(1).N0 = 0;
lcphd_states0(1).N1 = 0;
lcphd_states0(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown

% Storage for results
lcphd_results.label = '\lambda-CPHD';
lcphd_results.est(n_runs) = struct();
lcphd_results.states = cell(n_runs, 1);
for i = 1:n_runs
    lcphd_results.states{i} = lcphd_states0;
    lcphd_results.est(i).N = zeros(sim_steps, 1);
    lcphd_results.est(i).Xhat = cell(sim_steps, 1);
    lcphd_results.est(i).ospa = ones(sim_steps, 1) * ospa_c;
    lcphd_results.est(i).ospa_l = ones(sim_steps, 1) * ospa_c;
    lcphd_results.est(i).ospa_n = ones(sim_steps, 1) * ospa_c;
    lcphd_results.est(i).lambda_hat = zeros(sim_steps, 1);
    lcphd_results.est(i).lambda_hat(1) = 1;
end

% Figures for plotting
if lcphd
    lcphd_results.figures = struct();
    lcphd_results.figures.v_fig = figure;
    lcphd_results.figures.n_fig = figure;
    lcphd_results.figures.map_fig = figure;
    lcphd_results.figures.ospa_fig = figure;
    lcphd_results.figures.rho_fig = figure;
    lcphd_results.figures.lambda_fig = figure;
end

%% l-pD-CPHD Filter
% Filter parameters 
lpdcphd_params = RFS.LPDCPHD.lpdcphd_params();
lpdcphd_params.Nmax = 100;  % Maximum number of targets to track
lpdcphd_params.U0 = .03;    % Distance threshold for merging clutter components during pruning
lpdcphd_params.T0 = 1e-5;   % Weight threshold for discarding clutter components during pruning
lpdcphd_params.Jmax0 = 100; % Maximum number of clutter components to keep
lpdcphd_params.U1 = 4;      % Distance threshold for merging target components during pruning
lpdcphd_params.T1 = 1e-5;   % Weight threshold for discarding target components during pruning
lpdcphd_params.Jmax1 = 100; % Maximum number of target components to keep
lpdcphd_params.w_min = 0.5; % Minimum weight to be included in the output state estimate

% System model
lpdcphd_model = RFS.LPDCPHD.lpdcphd_model();

% Target dynamics
lpdcphd_model.F = model_F;
lpdcphd_model.Q = model_Q;

% Detection and survival probabilities
lpdcphd_model.ps0 = model_ps0;    % Probability of clutter survival
lpdcphd_model.ps1 = model_ps1;    % Probability of target survival

% Clutter birth model: uniform detection probability, Ngamma0 expected
% clutter births
lpdcphd_model.gamma0 = RFS.utils.BMRFS(Ngamma0, 1, 1);

% Target birth model, s = 1, t = 1 -> uniform
gamma_s = ones(size(birth_gmrfs.w));
gamma_t = gamma_s;
lpdcphd_model.gamma1 = RFS.utils.BGMRFS(birth_gmrfs.w, birth_gmrfs.m, birth_gmrfs.P, gamma_s, gamma_t);

% Clutter distribution
lpdcphd_model.kappa = 1 / A_fov;  % Clutter is equally likely anywhere

% Dilation constant for beta distributions, factor by which to enlarge
% variance of beta distribution when predicting
lpdcphd_model.kB = 1.5;

% Initial conditions and storage for l-pd-CPHD
lpdcphd_states0(sim_steps, 1) = RFS.LPDCPHD.lpdcphd_state();
lpdcphd_states0(1).v0 = RFS.utils.BMRFS();
lpdcphd_states0(1).v1 = RFS.utils.BGMRFS();
lpdcphd_states0(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown
lpdcphd_states0(1).N0 = 0;
lpdcphd_states0(1).N1 = 0;
lpdcphd_states0(1).lambda = 0;
lpdcphd_states0(1).pd0 = 0;
lpdcphd_states0(1).pd1 = 0;

% Storage for results
lpdcphd_results.label = '\lambda-p_D-CPHD';
lpdcphd_results.est(n_runs) = struct();
lpdcphd_results.states = cell(1, n_runs);
for i = 1:n_runs
    lpdcphd_results.states{i} = lpdcphd_states0;
    lpdcphd_results.est(i).N = zeros(sim_steps, 1);
    lpdcphd_results.est(i).Xhat = cell(sim_steps, 1);
    lpdcphd_results.est(i).ospa = ones(sim_steps, 1) * ospa_c;
    lpdcphd_results.est(i).ospa_l = ones(sim_steps, 1) * ospa_c;
    lpdcphd_results.est(i).ospa_n = ones(sim_steps, 1) * ospa_c;
    lpdcphd_results.est(i).lambda_hat = zeros(sim_steps, 1);
    lpdcphd_results.est(i).lambda_hat(1) = 1;
    lpdcphd_results.est(i).pd0_hat = zeros(sim_steps, 1);
    lpdcphd_results.est(i).pd1_hat = zeros(sim_steps, 1);
end

% l-pd-CPHD figures
if lpdcphd
    lpdcphd_results.figures = struct();
    lpdcphd_results.figures.v_fig = figure;
    lpdcphd_results.figures.n_fig = figure;
    lpdcphd_results.figures.map_fig = figure;
    lpdcphd_results.figures.ospa_fig = figure;
    lpdcphd_results.figures.rho_fig = figure;
    lpdcphd_results.figures.lambda_fig = figure;
    lpdcphd_results.figures.pd0_fig = figure;
    lpdcphd_results.figures.pd1_fig = figure;
end

%% Vo LMB Filter

%filter parameters
lmb_params.T_max= 100;                  %maximum number of tracks
lmb_params.track_threshold= 1e-3;       %threshold to prune tracks

lmb_params.H_bth= 5;                    %requested number of birth components/hypotheses (for LMB to GLMB casting before update)
lmb_params.H_sur= 1000;                 %requested number of surviving components/hypotheses (for LMB to GLMB casting before update)
lmb_params.H_upd= 1000;                 %requested number of updated components/hypotheses (for GLMB update)
lmb_params.H_max= 1000;                 %cap on number of posterior components/hypotheses (not used yet)
lmb_params.hyp_threshold= 1e-15;        %pruning threshold for components/hypotheses (not used yet)

lmb_params.L_max= 10;                   %limit on number of Gaussians in each track 
lmb_params.elim_threshold= 1e-5;        %pruning threshold for Gaussians in each track 
lmb_params.merge_threshold= 4;          %merging threshold for Gaussians in each track

lmb_params.P_G= 0.9999999;                           %gate size in percentage
lmb_params.gamma= chi2inv(lmb_params.P_G, n_states);   %inv chi^2 dn gamma value
lmb_params.gate_flag= 1;                             %gating on or off 1/0

% model parameters
lmb_model.x_dim= n_states;   %dimension of state vector
lmb_model.z_dim= 2;   %dimension of observation vector

% dynamical model parameters (CV model)
lmb_model.T= sim_dt;                                     %sampling period
%lmb_model.A0= [ 1 lmb_model.T; 0 1 ];                         %transition matrix                     
lmb_model.F= model_F;
%lmb_model.B0= [ (lmb_model.T^2)/2; lmb_model.T ];
%lmb_model.B= [ lmb_model.B0 zeros(2,1); zeros(2,1) lmb_model.B0 ];
%lmb_model.sigma_v = 5;
lmb_model.Q = model_Q;   %process noise covariance

% survival/death parameters
lmb_model.P_S= model_ps1;
lmb_model.Q_S= 1-lmb_model.P_S;

% birth parameters (LMB birth model, single component only)
lmb_birth_rate = .75;
lmb_birth_gmrfs = lmb_birth_rate .* center_bm;
lmb_birth_fig = figure;
%h_lmb_birth = RFS.utils.plotgmphd(lmb_birth_gmrfs, min_n:.1:max_n, min_e:.1:max_e);
lmb_model.T_birth= lmb_birth_gmrfs.J;         %no. of LMB birth terms
lmb_model.L_birth= ones(lmb_model.T_birth,1);                                       %no of Gaussians in each LMB birth term
lmb_model.r_birth= lmb_birth_gmrfs.w;                                                   %prob of birth for each LMB birth term
lmb_model.w_birth= num2cell(ones(lmb_birth_gmrfs.J, 1));                                %weights of GM for each LMB birth term
lmb_model.m_birth= mat2cell(lmb_birth_gmrfs.m, n_states, ones(1, lmb_birth_gmrfs.J))';       %means of GM for each LMB birth term
%lmb_model.B_birth= cell(lmb_model.T_birth,1);                                      %std of GM for each LMB birth term
lmb_model.P_birth= cell(lmb_model.T_birth,1);                                       %cov of GM for each LMB birth term
for j = 1:lmb_birth_gmrfs.J
    lmb_model.P_birth{j} = lmb_birth_gmrfs.P(:, :, j);
end

% observation model parameters (noisy x/y only)
lmb_model.H= H;    %observation matrix
%lmb_model.D= diag([ 10; 10 ]); 
lmb_model.R= R;              %observation noise covariance

% detection parameters
lmb_model.P_D= .98;   %probability of detection in measurements
lmb_model.Q_D= 1-lmb_model.P_D; %probability of missed detection in measurements

% clutter parameters
lmb_model.lambda_c= 10;                             %poisson average rate of uniform clutter (per scan)
%lmb_model.range_c= [ -1000 1000; -1000 1000 ];      %uniform clutter region
lmb_model.pdf_c= 1 / A_fov; %uniform clutter density

%initial prior
%tt_lmb_update= cell(0,1);

% Storage for results
lmb_results.label = 'LMB';
lmb_results.est(n_runs) = struct();
lmb_results.tt_lmb_update = cell(1, n_runs);
for i = 1:n_runs
    lmb_results.tt_lmb_update{i} = cell(0, 1);
    lmb_results.est(i).N = zeros(sim_steps, 1);
    lmb_results.est(i).Xhat = cell(sim_steps, 1);
    lmb_results.est(i).L = cell(sim_steps, 1);
    lmb_results.est(i).ospa = ones(sim_steps, 1) * ospa_c;
    lmb_results.est(i).ospa_l = ones(sim_steps, 1) * ospa_c;
    lmb_results.est(i).ospa_n = ones(sim_steps, 1) * ospa_c;
end

% Figures for lmb
if lmb
    lmb_results.figures = struct();
    lmb_results.figures.v_fig = figure;
    lmb_results.figures.n_fig = figure;
    lmb_results.figures.map_fig = figure;
    lmb_results.figures.ospa_fig = figure;
end

%% Adaptive Vo LMB

% Storage for results
almb_results.label = 'Adaptive LMB';
almb_results.est(n_runs) = struct();
almb_results.tt_lmb_update = cell(1, n_runs);
for i = 1:n_runs
    almb_results.tt_lmb_update{i} = cell(0, 1);
    almb_results.est(i).N = zeros(sim_steps, 1);
    almb_results.est(i).Xhat = cell(sim_steps, 1);
    almb_results.est(i).L = cell(sim_steps, 1);
    almb_results.est(i).ospa = ones(sim_steps, 1) * ospa_c;
    almb_results.est(i).ospa_l = ones(sim_steps, 1) * ospa_c;
    almb_results.est(i).ospa_n = ones(sim_steps, 1) * ospa_c;
end

% Figures for adaptive lmb
if almb
    almb_results.figures = struct();
    almb_results.figures.v_fig = figure;
    almb_results.figures.n_fig = figure;
    almb_results.figures.map_fig = figure;
    almb_results.figures.ospa_fig = figure;
end

%% Simulate

% Initial conditions and storage variables
obs = cell(1, n_runs);
truth(n_runs).targets = struct();
t_gmphd = 0;
t_lcphd = 0;
t_lpdcphd = 0;
t_lmb = 0;
t_almb = 0;
for i = 1:n_runs
    truth(i).x = zeros(sim_steps, 12);
    truth(i).vis_tgts = cell(sim_steps, 1); % True locations of visible targets
    truth(i).lambda = zeros(sim_steps, 1);  % True clutter rate
    truth(i).pd0 = zeros(sim_steps, 1);     % True clutter detection probability
    truth(i).pd1 = zeros(sim_steps, 1);     % True target detection probability
    truth(i).N = zeros(sim_steps, 1);       % True target cardinality    
    truth(i).lambda(1) = sensor.lambda;         
    truth(i).pd0(1) = true_pd0;
    truth(i).pd1(1) = true_pd1;

    obs{i} = cell(sim_steps, 1);
end

rfs_handles = [];
initial_targets = generate_targets(F, Q, bounds, n_tgt, sim_steps, sensor);

% Run multiple simulation runs
for r = 1:n_runs
    % Generate new targets
    %targets = generate_targets(F, Q, bounds, n_tgt, sim_steps, sensor);

    % Reset targets
    targets = initial_targets;

    %figure(map_fig);
    %h_map = RFS.utils.plot_targets(targets, 'b*');
    %hold on
    %title 'Target Locations'
    %set(gca, 'Fontsize', 18)
    %axis equal;
    %axis(bounds)
    %delete(h_map)

    % Simulate
    for k = 2:sim_steps
        % Update vehicle position
        truth(r).x(k, :) = truth(r).x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
        n = truth(r).x(k, 1);
        e = truth(r).x(k, 2);
        psi = truth(r).x(k, 6);

        % Update target positions
        existing = [targets.t_birth] <= k & [targets.t_death] > k;
        targets(existing) = targets(existing).step();

        % Update system parameters
        %if k < 100
        %    sensor.lambda = 7;
        %elseif k < 200
        %    sensor.lambda = 12;
        %else
        %    sensor.lambda = 2;
        %end
        truth(r).lambda(k) = sensor.lambda;
        truth(r).pd0(k) = true_pd0;
        truth(r).pd1(k) = true_pd1;

        % Get measurements of targets
        [rng, brng, r_true, b_true] = sensor.measure([n, e, psi], targets(existing));

        % Convert to absolute coordinates
        brng = brng + psi;
        n_obs = n + rng .* cosd(brng);
        e_obs = e + rng .* sind(brng);
        b_true = b_true + psi;
        n_obs_true = n + r_true .* cosd(b_true);
        e_obs_true = e + r_true .* sind(b_true);

        % Store the observations
        obs{r}{k} = zeros(length(rng), 2);
        obs{r}{k}(:, 1) = n_obs;
        obs{r}{k}(:, 2) = e_obs;
        truth(r).vis_tgts{k} = zeros(2, length(r_true));
        truth(r).vis_tgts{k}(1, :) = n_obs_true;
        truth(r).vis_tgts{k}(2, :) = e_obs_true;
        truth(r).N(k) = length(r_true);

        % Generate measurement
        measurement = RFS.LCPHD.cphd_measurement();
        measurement.Z = obs{r}{k}';
        measurement.H = H;
        measurement.R = R;

        % Run filters
        if gmphd
            tic
            [gmphd_results.est(r).v(k), gmphd_results.est(r).N(k), gmphd_results.est(r).Xhat{k}] = RFS.GMPHD.phd_filter(gmphd_results.est(r).v(k-1), ...
                gmphd_results.est(r).N(k-1), ...
                model_F, ...
                model_Q, ...
                model_ps1, ...
                model_pd1, ...
                birth_gmrfs, ...
                obs{r}{k}', ...
                H, ...
                R, ...
                gmphd_kappa, ...
                lcphd_params.U, ...
                lcphd_params.T, ...
                lcphd_params.Jmax, ...
                lcphd_params.w_min);
            t_gmphd = t_gmphd + toc;
            [gmphd_results.est(r).ospa(k), gmphd_results.est(r).ospa_l(k), gmphd_results.est(r).ospa_n(k)] = ospa_dist(gmphd_results.est(r).Xhat{k}, ...
                truth(r).vis_tgts{k}, ...
                ospa_c, ...
                ospa_p);
        end

        if lcphd
            tic
            [lcphd_results.states{r}(k), lcphd_results.est(r).Xhat{k}, lcphd_results.est(r).lambda_hat(k)] = RFS.LCPHD.lcphd_filter(lcphd_results.states{r}(k-1), measurement, lcphd_model, lcphd_params);
            t_lcphd = t_lcphd + toc;
            lcphd_results.est(r).N(k) = size(lcphd_results.est(r).Xhat{k}, 2); 
            [lcphd_results.est(r).ospa(k), lcphd_results.est(r).ospa_l(k), lcphd_results.est(r).ospa_n(k)] = ospa_dist(lcphd_results.est(r).Xhat{k}, ...
                truth(r).vis_tgts{k}, ...
                ospa_c, ...
                ospa_p);
        end

        if lpdcphd
            tic
            [lpdcphd_results.states{r}(k), lpdcphd_results.est(r).Xhat{k}] = RFS.LPDCPHD.lpdcphd_filter(lpdcphd_results.states{r}(k-1), measurement, lpdcphd_model, lpdcphd_params);
            t_lpdcphd = t_lpdcphd + toc;
            lpdcphd_results.est(r).lambda_hat(k) = lpdcphd_results.states{r}(k).lambda;
            lpdcphd_results.est(r).N(k) = size(lpdcphd_results.est(r).Xhat{k}, 2);
            lpdcphd_results.est(r).pd0_hat(k) = lpdcphd_results.states{r}(k).v0.w' * lpdcphd_results.states{r}(k).pd0 / sum(lpdcphd_results.states{r}(k).v0.w);
            lpdcphd_results.est(r).pd1_hat(k) = lpdcphd_results.states{r}(k).v1.w' * lpdcphd_results.states{r}(k).pd1 / sum(lpdcphd_results.states{r}(k).v1.w);
            [lpdcphd_results.est(r).ospa(k), lpdcphd_results.est(r).ospa_l(k), lpdcphd_results.est(r).ospa_n(k)] = ospa_dist(lpdcphd_results.est(r).Xhat{k}, ...
                truth(r).vis_tgts{k}, ...
                ospa_c, ...
                ospa_p);
        end

        if lmb
            tic
            [lmb_results.tt_lmb_update{r}, lmb_results.est(r).Xhat{k}, lmb_results.est(r).N(k), lmb_results.est(r).L{k}] = vo.lmb.jointlmb_filter(lmb_results.tt_lmb_update{r}, lmb_model, lmb_params, measurement, k);
            t_lmb = t_lmb + toc;
            [lmb_results.est(r).ospa(k), lmb_results.est(r).ospa_l(k), lmb_results.est(r).ospa_n(k)] = ospa_dist(lmb_results.est(r).Xhat{k}, ...
                truth(r).vis_tgts{k}, ...
                ospa_c, ...
                ospa_p);
        end

        if almb            
            almb_model = lmb_model;
            almb_model.lambda_c = lpdcphd_results.states{r}(k).lambda;
            almb_model.P_D = lpdcphd_results.est(r).pd1_hat(k);
            tic
            [almb_results.tt_lmb_update{r}, almb_results.est(r).Xhat{k}, almb_results.est(r).N(k), almb_results.est(r).L{k}] = vo.lmb.jointlmb_filter(almb_results.tt_lmb_update{r}, almb_model, lmb_params, measurement, k);
            t_almb = t_almb + toc;
            [almb_results.est(r).ospa(k), almb_results.est(r).ospa_l(k), almb_results.est(r).ospa_n(k)] = ospa_dist(almb_results.est(r).Xhat{k}, ...
                truth(r).vis_tgts{k}, ...
                ospa_c, ...
                ospa_p);
        end        

        fprintf("Run %d/%d, t = %d/%d\n", r, n_runs, k, sim_steps);
    end % Simulation run

    % Save true target tracks
    truth(r).targets = targets;

end % All simulation runs
run_complete = true;
%% Save data
data_dir = "./data/";
now = datetime;
timestampstr = string(yyyymmdd(now)) + "T" ...
    + string(hour(now)) ...
    + string(minute(now)) ...
    + string(second(now));
filename = data_dir + timestampstr;
if gmphd
    filename = filename + "_gmphd";
end
if lcphd
    filename = filename + "_lcphd";
end
if lpdcphd
    filename = filename + "_lpdcphd";
end
if lmb
    filename = filename + "_lmb";
end
if almb
    filename = filename + "_almb";
end
filename = filename + ".mat";
save(filename);
fprintf('Saved data to %s\n', filename);

%% Load data
if exist('run_complete', 'var') == 0
    clc
    clear
    close all
    data_dir = "./data/";
    to_load = "20231216T6958.0626_gmphd_lcphd_lpdcphd_lmb_almb.mat";
    load_path = data_dir + to_load;
    fprintf('Loading data from %s\n', load_path)
    load(data_dir + to_load);
    clear 'run_complete'
    dbstop 'error'
end

%% Plots from individual runs
% General RFS plots
runs_to_plot = []; % Runs to plot
for r = runs_to_plot
    delete(rfs_handles);
    rfs_handles = [];
    if gmphd
        gmphd_handles = RFS.utils.rfs_tracker_plots(gmphd_results.figures, ...
            gmphd_results.label, ...
            k, ...
            [n e psi], ...
            sensor, ...
            obs{r}, ...
            truth(r).vis_tgts, ...
            truth(r).targets, ...
            bounds, ...
            gmphd_results.est(r).v(end), ...
            gmphd_results.est(r).Xhat, ...
            gmphd_results.est(r).ospa, ...
            gmphd_results.est(r).ospa_l, ...
            gmphd_results.est(r).ospa_n);
        rfs_handles = [rfs_handles gmphd_handles];
    end

    if lcphd
        lcphd_handles = RFS.utils.rfs_tracker_plots(lcphd_results.figures, ...
            lcphd_results.label, ...
            k, ...
            [n e psi], ...
            sensor, ...
            obs{r}, ...
            truth(r).vis_tgts, ...
            truth(r).targets, ...
            bounds, ...
            lcphd_results.states{r}(end).v, ...
            lcphd_results.est(r).Xhat, ...
            lcphd_results.est(r).ospa, ...
            lcphd_results.est(r).ospa_l, ...
            lcphd_results.est(r).ospa_n);    
        % l-CPHD specific plots
        lcphd_h = RFS.LCPHD.lcphd_plots(lcphd_results.figures, lcphd_results.label, k, truth(r).targets, truth(r), lcphd_results.est(r).Xhat, lcphd_results.est(r).lambda_hat);
        rfs_handles = [rfs_handles lcphd_h lcphd_handles];
    end

    if lpdcphd
        lpdcphd_handles = RFS.utils.rfs_tracker_plots(lpdcphd_results.figures, ...
            lpdcphd_results.label, ...
            k, ...
            [n e psi], ...
            sensor, ...
            obs{r}, ...
            truth(r).vis_tgts, ...
            truth(r).targets, ...
            bounds, ...
            lpdcphd_results.states{r}(end).v1, ...
            lpdcphd_results.est(r).Xhat, ...
            lpdcphd_results.est(r).ospa, ...
            lpdcphd_results.est(r).ospa_l, ...
            lpdcphd_results.est(r).ospa_n);

        lpdcphd_h1 = RFS.LCPHD.lcphd_plots(lpdcphd_results.figures, lpdcphd_results.label, k, truth(r).targets, truth(r), lpdcphd_results.est(r).Xhat, lpdcphd_results.est(r).lambda_hat);
        lpdcphd_h2 = RFS.LPDCPHD.lpdcphd_plots(lpdcphd_results.figures, lpdcphd_results.label, k, truth(r).targets, truth(r), lpdcphd_results.est(r).Xhat, lpdcphd_results.states{r});
        rfs_handles = [rfs_handles lpdcphd_handles lpdcphd_h1 lpdcphd_h2];
    end

    if lmb
        lmb_handles = RFS.utils.rfs_tracker_plots(lmb_results.figures, ...
            lmb_results.label, ...
            k, ...
            [n e psi], ...
            sensor, ...
            obs{r}, ...
            truth(r).vis_tgts, ...
            truth(r).targets, ...
            bounds, ...
            RFS.utils.GMRFS(), ... % No easy intensity plot
            lmb_results.est(r).Xhat, ...
            lmb_results.est(r).ospa, ...
            lmb_results.est(r).ospa_l, ...
            lmb_results.est(r).ospa_n);
        rfs_handles = [rfs_handles lmb_handles];
    end

    if almb
        almb_handles = RFS.utils.rfs_tracker_plots(almb_results.figures, ...
            almb_results.label, ...
            k, ...
            [n e psi], ...
            sensor, ...
            obs{r}, ...
            truth(r).vis_tgts, ...
            truth(r).targets, ...
            bounds, ...
            RFS.utils.GMRFS(), ... % No easy intensity plot
            almb_results.est(r).Xhat, ...
            almb_results.est(r).ospa, ...
            almb_results.est(r).ospa_l, ...
            almb_results.est(r).ospa_n);
        rfs_handles = [rfs_handles almb_handles];
    end
end

%% Runtime stats
total_steps = n_runs * sim_steps;
fprintf("GMPHD: Total time = %0.2f, mean time = %0.2f\n", t_gmphd, t_gmphd / total_steps);
fprintf("l-CPHD: Total time = %0.2f, mean time = %0.2f\n", t_lcphd, t_lcphd / total_steps);
fprintf("l-pd-CPHD: Total time = %0.2f, mean time = %0.2f\n", t_lpdcphd, t_lpdcphd / total_steps);
fprintf("LMB: Total time = %0.2f, mean time = %0.2f\n", t_lmb, t_lmb/ total_steps);
fprintf("ALMB: Total time = %0.2f, mean time = %0.2f\n", t_almb, t_almb / total_steps);

%% Aggregate plots and data
RFS.utils.aggregate_plots(sim_steps, truth, gmphd_results, lcphd_results, lpdcphd_results, lmb_results, almb_results);

%% Target generation function
function targets = generate_targets(F, Q, bounds, n_tgt, sim_steps, sensor)
%targets = generate_targets(F, Q, bounds, n_tgts, sim_steps)
% Generate n targets
min_e = bounds(1);
max_e = bounds(2);
min_n = bounds(3);
max_n = bounds(4);

targets(n_tgt) = RFS.sim.Target_2D;
% Line of targets moving from left to right
spacing = 2.5;
for k = 1:n_tgt
    targets(k).F = F;
    targets(k).Q = Q;
    tgt_e = min_e;
    tgt_n = 10 + spacing * k;
    tgt_ndot = 0;
    tgt_edot = .25;
    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
    targets(k).X = x0;
end

% Targets initialized from a random location, moving in a random direction
vel = 0.2;
for k = 1:n_tgt
    % Target dynamics
    targets(k).F = F;
    targets(k).Q = Q;

    % Choose a random target starting location
    tgt_brng = (rand() * sensor.fov / 2) - sensor.fov/4;
    tgt_rng = (rand() * sensor.range / 2) + sensor.range/4;
    tgt_psi = rand() * 2 * pi;

    tgt_e = tgt_rng * sind(tgt_brng);
    tgt_n = tgt_rng * cosd(tgt_brng);
    tgt_ndot = vel * cos(tgt_psi);
    tgt_edot = vel * sin(tgt_psi);
    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
    targets(k).X = x0;

    % Target birth and death times
    targets(k).t_birth = randi([0, round(sim_steps / 5)]);
end

% Line of targets stationary in center of field of view
%for k = 1:n_tgt
%    targets(k).F = F;
%    targets(k).Q = Q;
%    tgt_e = 0;
%    tgt_n = 10 + spacing * k;
%    tgt_ndot = 0;
%    tgt_edot = 0;
%    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
%    targets(k).X = x0;
%end
end
