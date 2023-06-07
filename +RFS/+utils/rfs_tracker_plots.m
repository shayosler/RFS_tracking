function handles = rfs_tracker_plots(figures, prefix, k, v, x, Xhat, sensor, z, z_true, targets, bounds)
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
    n_true(n) = size(z_true{n}, 1);
    z_n = z{n};
    if ~isempty(z_n)
        h_obs = [h_obs plot(z_n(:, 2), z_n(:, 1), '.')];
        hold on
    end
    if ~isempty(Xhat{n})
        h_trackers = [h_trackers plot(Xhat{n}(2, :), Xhat{n}(1,:), 'o')];
        hold on
    end
end

% Plot final target locations?
tgt_pos = targets.get_position();
h_map = plot(tgt_pos(2, :), tgt_pos(1, :), 'b*');
hold on
handles = [handles h_map];

% Plot final observations
h_obs = [h_obs plot(z{k}(:, 2), z{k}(:, 1), 'o')];
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
figure(figures.n_fig);
h_nx = plot(n_trackers, 'LineWidth', 2);
hold on
h_ntrue = plot(n_true, '--', 'LineWidth', 2);
handles = [handles h_nx h_ntrue];
title([prefix ' Cardinalities after t = ' num2str(k)])
xlabel 'Time step'
ylabel 'Number'
legend('Tracked Objects', 'True Targets')
set(gca, 'Fontsize', 18)



end