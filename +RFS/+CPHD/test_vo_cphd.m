%% clear
%clc
%clear
close all

%% Set up
%seed = 52490;
seed = 1;
rng(seed);

%% Set up model from paper
vo_model= robust.lcphd.gms.gen_model;
vo_truth= robust.lcphd.gms.gen_truth(vo_model);
vo_meas=  robust.lcphd.gms.gen_meas(vo_model,vo_truth);
vo_est=   robust.lcphd.gms.run_filter(vo_model, vo_meas);
handles=  robust.lcphd.gms.plot_results(vo_model, vo_truth, vo_meas, vo_est);

%% Filtering parameters (from paper)
params = RFS.CPHD.cphd_params();
params.Nmax = 300;
params.U = 4;
params.T = 1e-5;
params.Jmax = 100;
params.w_min = 0.5;

%% System model
n_states = 4;
model = RFS.CPHD.cphd_model();

% Target dynamics
model.F = vo_model.F;
model.Q = vo_model.Q;

% Detection and survival probabilities
model.pd0 = vo_model.clutter_P_D;  % Probability of detecting a "clutter" object
model.pd1 = vo_model.P_D;  % Probability of detecting a target
model.ps0 = vo_model.clutter_P_S;  % Probability of clutter survival
model.ps1 = vo_model.P_S;  % Probability of target survival

% Target birth model
model.gamma1 = RFS.utils.GMRFS(vo_model.m_birth, vo_model.P_birth, vo_model.w_birth);

% Sensor/measurement model
R = vo_model.R;
H = vo_model.H;

% clutter: poisson 
A_fov = (vo_model.range_c(1, 2) - vo_model.range_c(1, 1)) * (vo_model.range_c(2, 2) - vo_model.range_c(2, 1));
lambda_true = 50; %2; % Expected number of clutter returns
model.Ngamma0 = vo_model.lambda_cb;
model.kappa = 1 / A_fov; % Clutter is equally likely anywhere

%% Define environment
% Bounds for northing and easting
min_n = vo_model.range_c(1, 1);
max_n = vo_model.range_c(1, 2);
min_e = vo_model.range_c(2, 1);
max_e = vo_model.range_c(2, 2);
northings = min_n:1:max_n;
eastings = min_e-10:1:max_e;

%% Simulate
sim_dt = 1;
t_total = 100;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions
x = zeros(sim_steps, 12);
Xhat = cell(sim_steps, 1);
lambda = zeros(sim_steps, 1);
lambda(1) = 1;
states(sim_steps, 1) = RFS.CPHD.cphd_state;
states(1).v = RFS.utils.GMRFS([0.1;0;0.1;0], diag([1 1 1 1]).^2, eps);
states(1).N0 = round((size(vo_meas.Z{1},2)-vo_model.P_D*sum(vo_model.w_birth))/vo_model.clutter_P_D);
states(1).N1 = 0;
states(1).rho = [zeros(states(1).N0-1,1);1;zeros(params.Nmax-states(1).N0+1,1)];
j = 1;
pass = 1;

map_fig = figure;
v_fig = figure;
n_fig = figure;
rho_fig = figure;
for k = 2:sim_steps
    tk = t(j);
    j = j+1;

    % Simulated observations
    measurement = RFS.CPHD.cphd_measurement();
    measurement.Z = vo_meas.Z{k-1};
    measurement.H = H;
    measurement.R = R;

    % Update estimates
    [states(k), Xhat{k}, lambda(k)] = RFS.CPHD.lcphd_filter(states(k-1), measurement, model, params);

    %% Plots
    if mod(k, 100) == 0
        handles = [];

        % Current target locations
%         figure(map_fig);
%         h_map = plot(map(:, 2), map(:, 1), 'b*');
%         handles = [handles h_map];
% 
%         % observations
%         figure(map_fig);
%         h_obs = [];

        n_trackers = zeros(k, 1);
        for kk = 1:k
            n_trackers(kk) = size(Xhat{kk}, 2);
            %obs_k = obs{kk};
            %if ~isempty(obs_k)
            %    h_obs = [h_obs plot(obs_k(:, 2), obs_k(:, 1), '.')];
            %end
        end
        %h_obs = [h_obs plot(obs{k}(:, 2), obs{k}(:, 1), 'o')];
        %handles = [handles h_obs];

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
        %figure(map_fig);
        %h_fov = plot_sonar_fov([n, e], psi, range, fov, 'b');
        %handles = [handles h_fov];

        % Plot estimated number of targets/clutter
        figure(n_fig)
        N0 = [states(1:k).N0];
        N1 = [states(1:k).N1];
        plot(N0, 'LineWidth', 2)
        hold on
        plot(N1, 'LineWidth', 2)
        plot(n_trackers, 'LineWidth', 2);
        plot(vo_truth.N(1:k), '--', 'LineWidth', 2)
        title 'Cardinality'
        xlabel 'Time step'
        legend('Clutter Generators', 'Targets', 'Tracked Objects', 'True Number of Targets')
        set(gca, 'Fontsize', 18)

        % Cardinality distribution
        figure(rho_fig);
        bar(0:length(states(k).rho) - 1, states(k).rho);
        set(gca, 'Fontsize', 18)
        title 'Hybrid Cardinality Distribution'
        xlabel 'N'
        ylabel '\rho'
        
        % Plot current intensity
        %figure(v_fig);
        %h_v = plotgmphd(states(k).v, northings, eastings);
        %handles = [handles h_v];
        %title(['PHD Intensity After t = ' num2str(k)])
        %set(gca, 'Fontsize', 18)
        %axis equal;
        %axis([-30, 30, 0, 100])
        %colorbar


        % Clear plotted stuff
        delete(handles);
        clf(n_fig);
    end

end
