%% clear
clc
clear
close all

%% Set up
%seeds = 52490;
%seeds = [1 2 3 4 5 6 7 8 9 10];
seeds = 1;
rng(seeds(1));

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
t_total = 200;
sim_dt = 1;
t = 0:sim_dt:t_total;
sim_steps = numel(t);
n_runs = length(seeds);

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
Q = diag([.1 0 .1 0]) * .1;

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

% Initialize targets. 
n_tgt = 10;  
targets(n_tgt) = RFS.sim.Target_2D;
% Line of targets moving from left to right
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

% Line of targets stationary in center of field of view
%for k = 1:n_tgt
%    targets(k).F = F;
%    targets(k).Q = Q;
%    tgt_e = 0;
%    tgt_n = 10 + 2.5 * k;
%    tgt_ndot = 0;
%    tgt_edot = 0;
%    x0 = [tgt_n tgt_ndot tgt_e tgt_edot];
%    targets(k).X = x0;
%end

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
F_cv = [1  sim_dt  0  0;
        0  1       0  0;
        0  0       1  sim_dt;
        0  0       0  1]; % Constant velocity
F_static = eye(n_states);
model_F = eye(n_states);
model_Q = 5 * eye(n_states);

% Detection/survival probabilities
model_pd0 = true_pd0;
model_pd1 = true_pd1;
model_ps0 = true_ps0;
model_ps1 = true_ps1;

% Target birth model
%bm = load('./models/birth_model_100.mat');
uniform_bm = load('./models/birth_model_200_uniform_rb.mat');
uniform_bm = uniform_bm.gamma;
centerline_bm = RFS.utils.GMRFS([23; 0], [35 0; 0 2], 1); % One component along center
left_edge_bm = RFS.utils.transform_gmrfs(centerline_bm, 0, 0, -pi/4);

birth_rate = 0.01; % Expected rate of new births
birth_gmrfs = birth_rate .* centerline_bm; % uniform_bm;
birth_fig = figure;
h_birth = RFS.utils.plotgmphd(birth_gmrfs, [min_n:.1:max_n], [min_e:.1:max_e]);

% Sensor/measurement model
sigma_rb = [.25 0;   % range
            0   1];   % bearing

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = sensor.range * sind(sigma_rb(2, 2));
R = R_0deg_true;
H = eye(n_states);

% Mean clutter birth rate
% lower -> smaller lambda_hat (or maybe just takes longer to converge?)
Ngamma0 = 5;

%% GMPHD Filter

% Initial conditions and storage for GMPHD
gmphd_results.label = 'GMPHD';
gmphd_results.est(r_runs)
gmphd_v(sim_steps) = RFS.utils.GMRFS();
gmphd_N = zeros(sim_steps, 1);
gmphd_Xhat = cell(sim_steps, 1);
gmphd_ospa = zeros(sim_steps, 1);

% Figures for GMPHD plots
if gmphd
    gmphd_figures = struct();
    gmphd_figures.v_fig = figure;
    gmphd_figures.n_fig = figure;
    gmphd_figures.map_fig = figure;
    gmphd_figures.ospa_fig = figure;
end

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
lcphd_model.ps0 = model_ps0;    % Probability of clutter survival
lcphd_model.ps1 = model_ps1;    % Probability of target survival
lcphd_model.pd0 = model_pd0;    % Probability of detecting a "clutter" object
lcphd_model.pd1 = model_pd1;    % Probability of detecting a target

% Birth model
lcphd_model.gamma1 = birth_gmrfs;

% clutter
lcphd_model.Ngamma0 = Ngamma0;  % Mean clutter birth rate, 
lcphd_model.kappa = 1 / A_fov;  % Clutter is equally likely anywhere

% Initial conditions and storage for l-CPHD
lcphd_Xhat = cell(sim_steps, 1);
lcphd_lambda_hat = zeros(sim_steps, 1);
lcphd_lambda_hat(1) = 1;
lcphd_ospa = zeros(sim_steps, 1);
lcphd_states(sim_steps, 1) = RFS.CPHD.cphd_state();
lcphd_states(1).v = RFS.utils.GMRFS();
lcphd_states(1).N0 = 0;
lcphd_states(1).N1 = 0;
lcphd_states(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown

% Figures for plotting
if lcphd
    lcphd_figures = struct();
    lcphd_figures.v_fig = figure;
    lcphd_figures.n_fig = figure;
    lcphd_figures.map_fig = figure;
    lcphd_figures.ospa_fig = figure;
    lcphd_figures.rho_fig = figure;
    lcphd_figures.lambda_fig = figure;
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
lpdcphd_Xhat = cell(sim_steps, 1);
lpdcphd_ospa = zeros(sim_steps, 1);
lpdcphd_states(sim_steps, 1) = RFS.LPDCPHD.lpdcphd_state();
lpdcphd_states(1).v0 = RFS.utils.BMRFS();
lpdcphd_states(1).v1 = RFS.utils.BGMRFS();
lpdcphd_states(1).rho = ones(lcphd_params.Nmax + 1, 1) ./ (lcphd_params.Nmax + 1); % initial cardinality distribution is unknown
lpdcphd_states(1).N0 = 0;
lpdcphd_states(1).N1 = 0;
lpdcphd_states(1).lambda = 0;
lpdcphd_states(1).pd0 = 0;
lpdcphd_states(1).pd1 = 0;

% l-pd-CPHD figures
if lpdcphd
    lpdcphd_figures = struct();
    lpdcphd_figures.v_fig = figure;
    lpdcphd_figures.n_fig = figure;
    lpdcphd_figures.map_fig = figure;
    lpdcphd_figures.ospa_fig = figure;
    lpdcphd_figures.rho_fig = figure;
    lpdcphd_figures.lambda_fig = figure;
    lpdcphd_figures.pd0_fig = figure;
    lpdcphd_figures.pd1_fig = figure;
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
lmb_birth_rate = .5;
lmb_birth_gmrfs = lmb_birth_rate .* centerline_bm;
lmg_birth_fig = figure;
h_lmb_birth = RFS.utils.plotgmphd(lmb_birth_gmrfs, min_n:.1:max_n, min_e:.1:max_e);
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
tt_lmb_update= cell(0,1);

% Data storage
lmb_est.X = cell(sim_steps, 1);
lmb_est.N = zeros(sim_steps, 1);
lmb_est.L = cell(sim_steps, 1);
lmb_ospa = zeros(sim_steps, 1);

% Figures for lmb
if lmb
    lmb_figures = struct();
    lmb_figures.v_fig = figure;
    lmb_figures.n_fig = figure;
    lmb_figures.map_fig = figure;
    lmb_figures.ospa_fig = figure;
end

%% Adaptive Vo LMB

%initial prior
tt_almb_update= cell(0,1);

% Data storage
almb_est.X = cell(sim_steps, 1);
almb_est.N = zeros(sim_steps, 1);
almb_est.L = cell(sim_steps, 1);
almb_ospa = zeros(sim_steps, 1);

% Figures for lmb
if almb
    almb_figures = struct();
    almb_figures.v_fig = figure;
    almb_figures.n_fig = figure;
    almb_figures.map_fig = figure;
    almb_figures.ospa_fig = figure;
end

%% Simulate

% Initial conditions and storage variables
x = zeros(sim_steps, 12);
obs = cell(sim_steps, 1);
truth.vis_tgts = cell(sim_steps, 1);    % True locations of visible targets
truth.lambda = zeros(sim_steps, 1);     % True clutter rate
truth.pd0 = zeros(sim_steps, 1);        % True clutter detection probability
truth.pd1 = zeros(sim_steps, 1);        % True target detection probability
truth.lambda(1) = sensor.lambda;
truth.pd0(1) = true_pd0;
truth.pd1(1) = true_pd1;

rfs_handles = [];

% Run one simulator run for each seed
for s = 1:length(seeds)
    seed = seeds(s);
    rng(seed);

    % Simulate
    for k = 2:sim_steps
        % Update vehicle position
        x(k, :) = x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
        n = x(k, 1);
        e = x(k, 2);
        psi = x(k, 6);

        % Update target positions
        existing = [targets.t_birth] <= k & [targets.t_death] > k;
        targets(existing) = targets(existing).step();

        % Update system parameters
        if k < 100
            sensor.lambda = 7;
        elseif k < 200
            sensor.lambda = 12;
        else
            sensor.lambda = 2;
        end
        truth.lambda(k) = sensor.lambda;
        truth.pd0(k) = true_pd0;
        truth.pd1(k) = true_pd1;

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
        truth.vis_tgts{k} = zeros(2, length(r_true));
        truth.vis_tgts{k}(1, :) = n_obs_true;
        truth.vis_tgts{k}(2, :) = e_obs_true;

        % Generate measurement
        measurement = RFS.CPHD.cphd_measurement();
        measurement.Z = obs{k}';
        measurement.H = H;
        measurement.R = R;

        % Update estimates
        if gmphd
            [gmphd_v(k), gmphd_N(k), gmphd_Xhat{k}] = RFS.GMPHD.phd_filter(gmphd_v(k-1), ...
                gmphd_N(k-1), ...
                model_F, ...
                model_Q, ...
                model_ps1, ...
                model_pd1, ...
                birth_gmrfs, ...
                obs{k}', ...
                H, ...
                R, ...
                lcphd_model.kappa, ...
                lcphd_params.U, ...
                lcphd_params.T, ...
                lcphd_params.Jmax, ...
                lcphd_params.w_min);
            gmphd_ospa(k) = ospa_dist(gmphd_Xhat{k}, truth.vis_tgts{k}, ospa_c, ospa_p);
        end

        if lcphd
            [lcphd_states(k), lcphd_Xhat{k}, lcphd_lambda_hat(k)] = RFS.CPHD.lcphd_filter(lcphd_states(k-1), measurement, lcphd_model, lcphd_params);
            lcphd_ospa(k) = ospa_dist(lcphd_Xhat{k}, truth.vis_tgts{k}, ospa_c, ospa_p);
        end

        if lpdcphd
            [lpdcphd_states(k), lpdcphd_Xhat{k}] = RFS.LPDCPHD.lpdcphd_filter(lpdcphd_states(k-1), measurement, lpdcphd_model, lpdcphd_params);
            lpdcphd_ospa(k) = ospa_dist(lpdcphd_Xhat{k}, truth.vis_tgts{k}, ospa_c, ospa_p);
        end

        if lmb
            [tt_lmb_update, lmb_est.X{k}, lmb_est.N(k), lmb_est.L{k}] = vo.lmb.jointlmb_filter(tt_lmb_update, lmb_model, lmb_params, measurement, k);
            lmb_ospa(k) = ospa_dist(lmb_est.X{k}, truth.vis_tgts{k}, ospa_c, ospa_p);
        end

        if almb
            almb_model = lmb_model;
            almb_model.lambda_c = lpdcphd_states(k).lambda;
            almb_model.P_D = lpdcphd_states(k).v1.w' * lpdcphd_states(k).pd1 / sum(lpdcphd_states(k).v1.w);
            [tt_almb_update, almb_est.X{k}, almb_est.N(k), almb_est.L{k}] = vo.lmb.jointlmb_filter(tt_almb_update, almb_model, lmb_params, measurement, k);
            almb_ospa(k) = ospa_dist(almb_est.X{k}, truth.vis_tgts{k}, ospa_c, ospa_p);
        end        

        fprintf("t = %0.2f\n", t(k));

        %% Plots
        if mod(k, 10) == 0 && false
            delete(rfs_handles);
            if gmphd
                gmphd_handles = RFS.utils.rfs_tracker_plots(gmphd_figures, ...
                    'GMPHD', ...
                    k, ...
                    [n e psi], ...
                    sensor, ...
                    obs, ...
                    truth.vis_tgts, ...
                    targets, ...
                    bounds, ...
                    gmphd_v(k), ...
                    gmphd_Xhat, ...
                    gmphd_ospa);
                rfs_handles = [rfs_handles gmphd_handles];
            end
            if lcphd
                lcphd_handles = RFS.utils.rfs_tracker_plots(lcphd_figures, ...
                    'l-CPHD', ...
                    k, ...
                    [n e psi], ...
                    sensor, ...
                    obs, ...
                    truth.vis_tgts, ...
                    targets, ...
                    bounds, ...
                    lcphd_states(k).v, ...
                    lcphd_Xhat, ...
                    lcphd_ospa);
                rfs_handles = [rfs_handles lcphd_handles];
            end
            if lpdcphd
                lpdcphd_handles = RFS.utils.rfs_tracker_plots(lpdcphd_figures, ...
                    'l-pd-CPHD', ...
                    k, ...
                    [n e psi], ...
                    sensor, ...
                    obs, ...
                    truth.vis_tgts, ...
                    targets, ...
                    bounds, ...
                    lpdcphd_states(k).v1, ...
                    lpdcphd_Xhat, ...
                    lpdcphd_ospa);
                rfs_handles = [rfs_handles lpdcphd_handles];
            end
        end % for k = 2:sim_steps
    end % Simulation run

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
filename = filename + ".mat";
save(filename);
fprintf('Saved data to %s\n', filename);

%% Load data
if exist('run_complete', 'var') == 0
    clc
    clear
    close all
    data_dir = "./data/";
    to_load = "20231214T22142.7108_gmphd_lcphd_lpdcphd.mat";
    load_path = data_dir + to_load;
    fprintf('Loading data from %s\n', load_path)
    load(data_dir + to_load);
    clear 'run_complete'
    dbstop 'error'
end

%% Final plots
% General RFS plots
delete(rfs_handles);
if gmphd
    gmphd_handles = RFS.utils.rfs_tracker_plots(gmphd_figures, ...
        'GMPHD', ...
        k, ...
        [n e psi], ...
        sensor, ...
        obs, ...
        truth.vis_tgts, ...
        targets, ...
        bounds, ...
        gmphd_v(k), ...
        gmphd_Xhat, ...
        gmphd_ospa);
end

if lcphd
    lcphd_handles = RFS.utils.rfs_tracker_plots(lcphd_figures, ...
        'l-CPHD', ...
        k, ...
        [n e psi], ...
        sensor, ...
        obs, ...
        truth.vis_tgts, ...
        targets, ...
        bounds, ...
        lcphd_states(k).v, ...
        lcphd_Xhat, ...
        lcphd_ospa);    % l-CPHD specific plots
    lcphd_h = RFS.CPHD.lcphd_plots(lcphd_figures, 'l-CPHD', k, targets, truth, lcphd_Xhat, lcphd_lambda_hat);
end

if lpdcphd
    lpdcphd_handles = RFS.utils.rfs_tracker_plots(lpdcphd_figures, ...
        'l-pd-CPHD', ...
        k, ...
        [n e psi], ...
        sensor, ...
        obs, ...
        truth.vis_tgts, ...
        targets, ...
        bounds, ...
        lpdcphd_states(k).v1, ...
        lpdcphd_Xhat, ...
        lpdcphd_ospa);

    lpdcphd_h1 = RFS.CPHD.lcphd_plots(lpdcphd_figures, 'l-pd-CPHD', k, targets, truth, lpdcphd_Xhat, [lpdcphd_states.lambda]);
    lpdcphd_h2 = RFS.LPDCPHD.lpdcphd_plots(lpdcphd_figures, 'l-pd-CPHD', k, targets, truth, lpdcphd_Xhat, lpdcphd_states);
end

if lmb
    lmb_handles = RFS.utils.rfs_tracker_plots(lmb_figures, ...
        'LMB', ...
        k, ...
        [n e psi], ...
        sensor, ...
        obs, ...
        truth.vis_tgts, ...
        targets, ...
        bounds, ...
        RFS.utils.GMRFS(), ... % No easy intensity plot
        lmb_est.X, ...
        lmb_ospa);
end

if almb
    almb_handles = RFS.utils.rfs_tracker_plots(almb_figures, ...
        'Adaptive-LMB', ...
        k, ...
        [n e psi], ...
        sensor, ...
        obs, ...
        truth.vis_tgts, ...
        targets, ...
        bounds, ...
        RFS.utils.GMRFS(), ... % No easy intensity plot
        almb_est.X, ...
        almb_ospa);
end
