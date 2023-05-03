%% clear
clc
clear
close all

%% Set up
%seed = 52490;
seed = 1;
rng(seed);

%% Filtering parameters
params = cphd_params();
params.Nmax = 32;
params.U = 4;
params.T = 1e-5;
params.Jmax = 100;
params.w_min = 0.5;

%% System model
n_states = 2;
model = cphd_model();

% Target dynamics
model.F = eye(n_states);
model.Q = .5 * eye(n_states);

% Detection and survival probabilities
model.pd0 = 0.6;  % Probability of detecting a "clutter" object
model.pd1 = 1.0;  % Probability of detecting a target
model.ps0 = 0.1;  % Probability of clutter survival
model.ps1 = 0.99;  % Probability of target survival

% Target birth model
bm = load('../models/birth_model_100.mat');
birth_rate = 0.01; % Expected rate of new births
model.gamma1 = birth_rate .* bm.gamma;
model.rho_gamma1 = poisspdf(0:1:params.Nmax, birth_rate);

% Sensor/measurement model
fov = 90;
range = 40;
sigma_rb = [.25 0;   % range
            0  10];   % bearing
sigma_rb = [0.001 0;
            0 0.001];
% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = range * sind(sigma_rb(2, 2));

H = eye(n_states);

% clutter: poisson 
A_fov = (pi * range^2) * (fov / 360);
lambda_true = 0; %2; % Expected number of clutter returns
model.Ngamma0 = 0.5; % Mean clutter birth rate
model.kappa = 1 / A_fov; % Clutter is equally likely anywhere

%% Define environment
% Bounds for northing and easting
min_n = 0;
max_n = 100;
min_e = -20;
max_e = 100;
northings = min_n:.1:max_n;
eastings = min_e-10:.1:max_e;

% Simple targets: line of targets heading straight at sensor
n_tgt = 10;  % Start with a single target
spacing = 4; 
v_n = 0;
v_e = 0;
v_tgt = repmat([v_n, v_e], n_tgt, 1);
tgt_e = zeros(n_tgt, 1); % + min_e;
tgt_n = linspace(35, 10, n_tgt)';
map = [tgt_n tgt_e ones(n_tgt, 1) * model.pd1]; %FIXME: testing perfect detection

% Plot map
map_fig = figure;
h_map = plot(map(:, 2), map(:, 1), 'b*');
hold on
title 'Target Locations'
set(gca, 'Fontsize', 18)
axis equal;
axis([-30, 30, 0, 100])
delete(h_map)

%% Simulate
sim_dt = 1;
t_total = 100;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions
x = zeros(sim_steps, 12);
Xhat = cell(sim_steps, 1);
obs = cell(sim_steps, 1);
true_obs = cell(sim_steps, 1);
lambda = zeros(sim_steps, 1);
lambda(1) = 1;
states(sim_steps, 1) = cphd_state;
states(1).v = GMRFS();
states(1).N0 = 0;
states(1).N1 = 0;
states(1).rho = ones(params.Nmax + 1, 1) ./ (params.Nmax + 1); % initial cardinality distribution is unknown
j = 1;
pass = 1;

v_fig = figure;
n_fig = figure;
rho_fig = figure;
for k = 2:sim_steps
    tk = t(j);
    j = j+1;

    % Update vehicle position
    x(k, :) = x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
    n = x(k, 1);
    e = x(k, 2);
    psi = x(k, 6);

    % Get relative positions of simuated observations
    [r, b, r_true, b_true] = simulate_sonar_obs([n, e, psi], map, range, fov, lambda_true, sigma_rb);
    map(:, 1:2) = map(:, 1:2) + v_tgt .* sim_dt;

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
    %all_obs = [all_obs; obs{k}];   

    measurement = cphd_measurement();
    measurement.Z = obs{k}';
    measurement.H = H;
    measurement.R = R;

    % Update estimates
    [states(k), Xhat{k}, lambda(k)] = cphd_filter(states(k-1), measurement, model, params);
    %[v_k_obs, N_in_k, ~] = phd_filter(v{k-1}, N_in, F, Q, ps, pd, gamma, obs{k}', H, R, kappa, U, T, Jmax, w_min);

    %% Plots
    if mod(k, 100) == 0
        handles = [];

        % Current target locations
        figure(map_fig);
        h_map = plot(map(:, 2), map(:, 1), 'b*');
        handles = [handles h_map];

        % observations
        figure(map_fig);
        h_obs = [];
        n_trackers = zeros(k, 1);
        for kk = 1:k
            n_trackers(kk) = size(Xhat{kk}, 2);
            obs_k = obs{kk};
            if ~isempty(obs_k)
                h_obs = [h_obs plot(obs_k(:, 2), obs_k(:, 1), '.')];
            end
        end
        h_obs = [h_obs plot(obs{k}(:, 2), obs{k}(:, 1), 'o')];
        handles = [handles h_obs];

        % Plot estimates
        figure(map_fig);
        trackers = Xhat{k};
        if ~isempty(trackers)
            h_tracker = plot(trackers(2, :), trackers(1, :), 'go', 'LineWidth', 2, 'MarkerSize', 6);
            title(['Map and Observations After t = ' num2str(k)])
            legend('Targets', 'Measurements', 'Tracked Objects')
            set(gca, 'Fontsize', 18)
            axis equal;
            axis([-30, 30, 0, 100])
            handles = [handles h_tracker];
        end
        
        % Plot current fov
        figure(map_fig);
        h_fov = plot_sonar_fov([n, e], psi, range, fov, 'b');
        handles = [handles h_fov];

        % Plot estimated number of targets/clutter
        figure(n_fig)
        N0 = [states(1:k).N0];
        N1 = [states(1:k).N1];
        plot(N0, 'LineWidth', 2)
        hold on
        plot(N1, 'LineWidth', 2)
        plot(n_trackers, 'LineWidth', 2);
        title 'Cardinality'
        xlabel 'Time step'
        legend('Clutter Generators', 'Targets', 'Tracked Objects')
        set(gca, 'Fontsize', 18)

        % Cardinality distribution
        figure(rho_fig);
        bar(0:length(states(k).rho) - 1, states(k).rho);
        set(gca, 'Fontsize', 18)
        title 'Hybrid Cardinality Distribution'
        xlabel 'N'
        ylabel '\rho'
        
        % Plot current intensity
        figure(v_fig);
        h_v = plotgmphd(states(k).v, northings, eastings);
        handles = [handles h_v];
        title(['PHD Intensity After t = ' num2str(k)])
        set(gca, 'Fontsize', 18)
        axis equal;
        axis([-30, 30, 0, 100])
        colorbar


        % Clear plotted stuff
        delete(handles);
        clf(n_fig);
        break;
    end

end
