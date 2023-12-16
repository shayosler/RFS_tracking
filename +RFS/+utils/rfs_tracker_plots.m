function handles = rfs_tracker_plots(figures, prefix, k, x, sensor, z, z_true, targets, bounds, v, Xhat, ospa, ospa_l, ospa_n)
%rfs_tracker_plots Generate plots common to all trackers
% Inputs:
%   prefix
%   t
%   v
%   x           Vehicle position [n e psi]
%   Xhat
%   z
%   ztrue
%   targets
%   bounds      Plot boundaries: [min_e max_e min_n max_n]

% Validate
if ~isvector(targets)
    error('targets must be 1xN or Nx1');
end

% Plot target trajectories
handles = [];
if ~isfield(figures, 'map_fig') || isempty(figures.map_fig)
    figures.map_fig = figure;
end
figure(figures.map_fig);
for n = 1:length(targets)
    traj = targets(n).get_trajectory();
    h = plot(traj(3, :), traj(1, :), 'b');
    hold on
    handles = [handles h];
end

% Plot observations and tracker positions at each time step
%k = length(Xhat);
h_obs = [];
h_trackers = [];
n_trackers = zeros(k, 1);
n_true = zeros(k, 1);
for n = 1:k
    n_trackers(n) = size(Xhat{n}, 2);
    n_true(n) = size(z_true{n}, 2);
    z_n = z{n};
    if ~isempty(z_n)
        h_obs = [h_obs plot(z_n(:, 2), z_n(:, 1), '.')];
        hold on
    end
    if ~isempty(Xhat{n})
        h_trackers = [h_trackers plot(Xhat{n}(2, :), Xhat{n}(1,:), 'go')];
        hold on
    end
end

% Plot final target locations?
tgt_pos = targets.get_position();
h_map = plot(tgt_pos(2, :), tgt_pos(1, :), 'b*');
hold on
handles = [handles h_map];

% Plot final observations
h_obs = [h_obs plot(z{k}(:, 2), z{k}(:, 1), 'bo')];
handles = [handles h_obs h_trackers];

% Plot final position estimates
%trackers = Xhat{k};
%if ~isempty(trackers)
%    h_tracker = plot(trackers(2, :), trackers(1, :), 'go', 'LineWidth', 2, 'MarkerSize', 6);
%    handles = [handles h_tracker];
%end

% Plot current fov
n = x(1);
e = x(2);
psi = x(3);
h_fov = sensor.plot_fov(n, e, psi, 'b');
handles = [handles h_fov];

title([prefix ' Map and Observations After t = ' num2str(k)])
%legend('Targets', 'Measurements', 'Tracked Objects')
set(gca, 'Fontsize', 18)
axis equal;
axis(bounds)

% Plot current intensity
if ~isfield(figures, 'v_fig') || isempty(figures.v_fig)
    figures.v_fig = figure;
end
figure(figures.v_fig);
northings = bounds(3):.1:bounds(4);
eastings = bounds(1):.1:bounds(2);
h_v = RFS.utils.plotgmphd(v, northings, eastings);
handles = [handles h_v];
title([prefix ' PHD Intensity After t = ' num2str(k)])
set(gca, 'Fontsize', 18)
axis equal;
axis(bounds)
colorbar


% Plot true vs estimated number of targets
if ~isfield(figures, 'n_fig') || isempty(figures.n_fig)
    figures.n_fig = figure;
end
figure(figures.n_fig);
h_nx = plot(n_trackers, 'g', 'LineWidth', 2);
hold on
h_ntrue = plot(n_true, 'b--', 'LineWidth', 2);
handles = [handles h_nx h_ntrue];
title([prefix ' Cardinalities after t = ' num2str(k)])
xlabel 'Time step'
ylabel 'Number'
legend('Tracked Objects', 'True Targets')
set(gca, 'Fontsize', 18)

% Plot OSPA
if ~isfield(figures, 'ospa_fig') || isempty(figures.ospa_fig)
    figures.ospa_fig = figure;
end
figure(figures.ospa_fig);
subplot(3, 1, 1)
h_ospa = plot(ospa(2:end), 'LineWidth', 2);
xlabel 'Time step'
ylabel 'OSPA Distance'
title([prefix ' OSPA Metrics'])
set(gca, 'Fontsize', 18)

subplot(3, 1, 2)
h_ospa_l = plot(ospa_l(2:end), 'LineWidth', 2);
xlabel 'Time step'
ylabel 'OSPA Location'
set(gca, 'Fontsize', 18)

subplot(3, 1, 3)
h_ospa_n = plot(ospa_n(2:end), 'LineWidth', 2);
xlabel 'Time step'
ylabel 'OSPA Cardinality'
set(gca, 'Fontsize', 18)
handles = [handles h_ospa h_ospa_l h_ospa_n];

%% print out some statistices
n_err = n_true - n_trackers;
rms_n_err = rms(n_err);
mean_n_err = mean(n_err);
std_n_err = std(n_err);
mean_ospa = mean(ospa(2:end));
fprintf('%s RMS Cardinality Error: %0.3f\n', prefix, rms_n_err);
fprintf('%s Mean Cardinality Error: %0.3f\n', prefix, mean_n_err);
fprintf('%s Std Cardinality Error: %0.3f\n', prefix, std_n_err);
fprintf('%s Mean OSPA: %0.3f\n', prefix, mean_ospa);

end