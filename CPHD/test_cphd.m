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
params.Nmax = 1000;
params.U = 4;
params.T = 1e-5;
params.Jmax = 100;
params.w_min = 0.5;

%% System model
n_states = 2;
model = cphd_model();

% Target dynamics
model.F = eye(n_states);
model.Q = 0.1 * eye(n_states);

% Detection and survival probabilities
model.pd0 = 0.6;  % Probability of detecting a "clutter" object
model.pd1 = 0.6;  % Probability of detecting a target
model.ps0 = 0.1;  % Probability of clutter survival
model.ps1 = 0.9;  % Probability of target survival

% Target birth model
bm = load('birth_model_100.mat');
birth_rate = 1.0; % Expected rate of new births
model.gamma1 = birth_rate .* bm.gamma;
model.rho_gamma = poisspdf(0:1:params.Nmax, birth_rate);

% clutter: poisson 
A_fov = (pi * range^2) * (fov / 360);
lambda = 2; % Expected number of clutter returns
model.Ngamma = lambda;
kappa = 1 / A_fov; % Clutter is equally likely anywhere

% Sensor/measurement model
fov = 90;
range = 40;
sigma_rb = [.25 0;   % range
            0  10];   % bearing
% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = range * sind(sigma_rb(2, 2));
H = eye(n_states);

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

%% Simulate
sim_dt = 0.1;
t_total = 400;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions
x = zeros(sim_steps, 12);
obs = cell(sim_steps, 1);
true_obs = cell(sim_steps, 1);
lambda = zeros(sim_steps, 1);
states = cell(sim_steps, 1);
states{1} = cphd_state();
states{1}.v = GMRFS();
states{1}.N0 = 0;
states{1}.N1 = 0;
states{1}.rho = ones(params.Nmax + 1) ./ params.Nmax + 1; % initial cardinality distribution is unknown
j = 1;
pass = 1;

v_fig = figure;
for k = 2:sim_steps
    tk = t(j);
    j = j+1;

    % Update vehicle position
    x(k, :) = x(1, :); %grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);
    n = x(k, 1);
    e = x(k, 2);
    psi = x(k, 6);

    % Get relative positions of simuated observations
    [r, b, r_true, b_true] = simulate_sonar_obs([n, e, psi], map, range, fov, lambda, sigma_rb);

    % Convert to absolute coordinates
    b = b + psi; 
    n_obs = n + r .* cosd(b);
    e_obs = e + r .* sind(b);
    b_true = b_true + psi;
    n_obs_true = n + r_true .* cosd(b_true);
    e_obs_true = e + r .* sind(b_true);

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

    [states{k}, Xhat, lambda(k)] = phd_filter(states{k-1}, measurement, model, params);
    %[v_k_obs, N_in_k, ~] = phd_filter(v{k-1}, N_in, F, Q, ps, pd, gamma, obs{k}', H, R, kappa, U, T, Jmax, w_min);

    %% Plots
    if true
        % observations
        %figure(map_fig);
        %h_obs = [];
        %h_obs = plot(all_obs(:, 2), all_obs(:, 1), 'r.');
        for kk = 1:k
            obs_k = obs{kk};
            if ~isempty(obs_k)
                h_obs = [h_obs plot(obs_k(:, 2), obs_k(:, 1), '.')];
            end
        end

        % Plot estimates
        figure(map_fig);
        h_tracker = plot(Xhat(2, :), Xhat(1, :), 'go','LineWidth', 2, 'MarkerSize', 6);
        title(['Map and Observations After Pass ' num2str(pass)])
        legend('Targets', 'Measurements', 'Tracked Objects')
        set(gca, 'Fontsize', 18)
        axis equal;
        axis([-30, 30, 0, 100])
        

        % Plot current fov
        figure(map_fig);
        h_fov = plot_fov([n, e], psi, range, fov, 'b');
        
        % Plot current intensity
        figure(v_fig);
        h_v = plotgmphd(v{k}, northings, eastings);
        title(['PHD Intensity After Pass ' num2str(pass)])
        set(gca, 'Fontsize', 18)
        axis equal;
        axis([-30, 30, 0, 100])
        
        colorbar

        % Clear plotted stuff
        delete(h_v);
        delete(h_fov);
        delete(h_tracker);
        delete(h_obs);
    end

end
