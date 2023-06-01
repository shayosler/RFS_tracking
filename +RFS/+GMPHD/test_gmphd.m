%% Generate a dataset of simulated observations
clc
clear
close all

%% General setup parameters
%seed = 52490;
seed = 1;
rng(seed);

%% Sensor
fov = 90;
range = 40;
sigma_rb = [.25 0;   % range
            0  10];   % bearing
p_mine = .6; %0.3;   % probability of detecting a mine
p_other = p_mine; %.6; %0.1;  % probability of detecting a non-mine target as a mine

% clutter: poisson process
A_fov = (pi * range^2) * (fov / 360);
lambda = 2; % Expected number of clutter returns
kappa = lambda / A_fov;

% TODO: R should be calculated for each measurement
R = sigma_rb;
R(2, 2) = range * sind(sigma_rb(2, 2));


%% Define the map
% Bounds for northing and easting
min_n = 0;
max_n = 100;
min_e = -20;
max_e = 100;

northings = min_n:.1:max_n;
eastings = min_e-10:.1:max_e;

% Number of mines and other sonar targets
n_mines = 10;
n_other = 10;

% Randomly place targets
mines = [min_n + (max_n - min_n).*rand(n_mines,1), ...
    min_n + (max_e - min_e).*rand(n_mines,1), ...
    ones(n_mines, 1) * p_mine];
mines = [mines; 39, 0, p_mine];
other = [min_n + (max_n - min_n).*rand(n_other,1), ...
    min_n + (max_e - min_e).*rand(n_other,1), ...
    ones(n_other, 1) * p_other];
map = [mines; other];

% Plot map
map_fig = figure;
plot(map(:, 2), map(:, 1), 'b*')
%plot(mines(:, 2), mines(:, 1), 'b*')
hold on
%plot(other(:, 2), other(:, 1), 'r*')
title 'Target Locations'
legend('Targets', 'Measurements', 'Tracked Objects')
set(gca, 'Fontsize', 18)
axis equal;
axis([-30, 30, 0, 100])

%% Path
%   u       Forward velocity
%   l       Line length (m)
%   d       line separation (m)
%   psi0    Initial line direction
%   n0      Northing of start point    
%   e0      Easting of start point
%   z0;     Depth to run at
u = 0.5;
l = 100;
d = 10;
psi0 = 0;
n0 = 0;
e0 = 0;
z0 = 0;
lines = 10;
t_total = (l / u * lines) + (lines - 1) * (pi * d / 2) / u;

%% PHD Filtering parameters
n_states = 2;
F = eye(n_states);
N = 0;
ps = 1.0;
pd = p_mine;
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
sim_dt = 2;
t = 0:sim_dt:t_total;
sim_steps = floor(t_total / sim_dt);
sim_steps = 1000;

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

    % just run the same line over and over
    if x(k-1, 1) > 90
        j = 1;
    end
    tk = t(j);
    j = j+1;

    % Update position
    x(k, :) = grid_trajectory(tk, u, l, d, psi0, n0, e0, z0);

    n = x(k, 1);
    e = x(k, 2);
    psi = x(k, 6);

    % Get relative positions of simuated observations
    [r, b] = simulate_obs([n, e, psi], map, range, fov, lambda, sigma_rb);

    % Approximate measurement error.
    % TODO: use a gaussian mixture for this
    %R = sigma_rb;
    %R(2, 2) = r * sind(sigma_rb(2, 2));
    % TODO: rotate

    % Convert to absolute coordinates
    b = b + psi; 
    n_obs = n + r .* cosd(b);
    e_obs = e + r .* sind(b);

    % Store the observations
    obs{k} = zeros(length(r), 2);
    obs{k}(:, 1) = n_obs;
    obs{k}(:, 2) = e_obs;
    all_obs = [all_obs; obs{k}];
    %waitforbuttonpress;

    % Update birth model
    gamma = transform_gmrfs(bm.gamma, n, e, psi * pi / 180);
    
    % Filter through GMPHD
    % Only filter on the subset of the RFS that is within the FOV
    visible = in_fov([n; e; psi], v{k-1}.m, range, -fov/2, fov/2);
    if v{k-1}.J > 0
        v_in_fov = GMRFS(v{k-1}.m(:, visible), v{k-1}.P(:, :, visible), v{k-1}.w(visible));
        v_out_fov = GMRFS(v{k-1}.m(:, ~visible), v{k-1}.P(:, :, ~visible), v{k-1}.w(~visible));
    else
        v_in_fov = v{k-1};
        v_out_fov = GMRFS();
    end
    %v_in_fov = v{k-1};
    N_in = sum(v_in_fov.w);
    N_out = sum(v_out_fov.w);
    [v_k_obs, N_in_k, ~] = phd_filter(v_in_fov, N_in, F, Q, ps, pd, gamma, obs{k}', H, R, kappa, U, T, Jmax, w_min);
    v_k_unpruned = v_k_obs + v_out_fov;
    N(k) = sum(v_k_unpruned.w); % TODO: this might not be right

    %% Prune
    v_k = prune_gmphd(v_k_unpruned, T, U, Jmax);

    % Renormalize the weights in v_k so that N_k stays the same
    % TODO: this wasn't done in the paper, but it seems like a step that should
    % happen
    %sum_pruned = sum(v_k.w);
    %v_k = N(k) / sum_pruned .* v_k;

    %% Extract a state estimate from the RFS
    Xhat = [];
    for i = 1:v_k.J
        if v_k.w(i) > w_min
            Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
        end
    end
    v{k} = v_k;

    %% Plots
    if(x(k, 1) > 90)
        % observations
        figure(map_fig);
        h_obs = [];
        h_obs = plot(all_obs(:, 2), all_obs(:, 1), 'r.');
        %for kk = 1:k
        %    obs_k = obs{kk};
        %    if ~isempty(obs_k)
        %        h_obs = [h_obs plot(obs_k(:, 2), obs_k(:, 1), '.')];
        %    end
        %end

        % Plot estimates
        figure(map_fig);
        h_tracker = plot(Xhat(2, :), Xhat(1, :), 'go','LineWidth', 2, 'MarkerSize', 6);
        title(['Map and Observations After Pass ' num2str(pass)])
        legend('Targets', 'Measurements', 'Tracked Objects')
        set(gca, 'Fontsize', 18)
        axis equal;
        axis([-30, 30, 0, 100])
        

        % Plot current fov
        %figure(map_fig);
        %h_fov = plot_fov([n, e], psi, range, fov, 'b');
        
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
        %delete(h_fov);
        delete(h_tracker);
        delete(h_obs);
        pass = pass + 1;
    end
end

%% Plots 

% observations
figure(map_fig);
h_obs = plot(all_obs(:, 2), all_obs(:, 1), 'r.');
%for k = 1:sim_steps
%    obs_k = obs{k};
%    if ~isempty(obs_k)
%        h_obs = [h_obs plot(obs_k(:, 2), obs_k(:, 1), '.')];
%    end
%end

% Plot estimates
figure(map_fig);
h_tracker = plot(Xhat(2, :), Xhat(1, :), 'go','LineWidth', 2, 'MarkerSize', 6);
title(['Map and Observations After Pass ' num2str(pass)])
legend('Targets', 'Measurements', 'Tracked Objects')
set(gca, 'Fontsize', 18)
axis equal;
axis([-30, 30, 0, 100])

% Plot current fov
%figure(map_fig);
%h_fov = plot_fov([n, e], psi, range, fov, 'b');
%axis equal;

figure(v_fig);
h_v = plotgmphd(v{k}, northings, eastings);
title(['PHD Intensity After Pass ' num2str(pass)])
set(gca, 'Fontsize', 18)
axis equal;
axis([-30, 30, 0, 100])



