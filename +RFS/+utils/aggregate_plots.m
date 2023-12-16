function [h] = aggregate_plots(sim_steps, truth, gmphd, lcphd, lpdcphd, lmb, almb)
%[h] = aggregate_plots(gmphd, lcphd, lpdcphd, lmb, almb)
% Generate plots from aggregate data sets

h = [];

% Plot styling
width = {'LineWidth', 2};
styles = {'g', 'r', 'm', 'c', 'y'};

% Figures and truth data
n_fig = figure;
h_ntrue = plot(mean([truth.N], 2), 'b--', width{:});
h = [h h_ntrue];
hold on

n_err_fig = figure;
o_fig = figure;

% Extract data to a format we can do math on it
estimates = {gmphd, lcphd, lpdcphd, lmb, almb};
N = cell(1, length(estimates));
N_err = cell(1, length(estimates));
ospa = cell(1, length(estimates));
ospa_l = cell(1, length(estimates));
ospa_n = cell(1, length(estimates));
labels = cell(1, length(estimates));
for k = 1:length(estimates)
    labels{k} = estimates{k}.label;
    est = estimates{k}.est;
    N{k} = [est.N];
    N_err{k} = [truth.N] - [est.N];
    ospa{k} = [est.ospa];
    ospa_l{k} = [est.ospa_l];
    ospa_n{k} = [est.ospa_n];

    % Cardinality
    figure(n_fig);
    h_nx = plot(mean([est.N], 2), styles{k}, width{:});
    h = [h h_nx];

    % Cardinality error
    figure(n_err_fig)
    plot(mean(N_err{k}, 2), styles{k}, width{:});
    hold on

    % OSPA
    figure(o_fig);
    subplot(3, 1, 1)
    h_ospa = plot(mean(ospa{k}, 2), styles{k}, width{:});
    hold on

    subplot(3, 1, 2)
    h_ospa_l = plot(mean(ospa_l{k}, 2), styles{k}, width{:});
    hold on

    subplot(3, 1, 3)
    h_ospa_n = plot(mean(ospa_n{k}, 2), styles{k}, width{:});
    hold on
    h = [h h_ospa h_ospa_l h_ospa_n];
end % for each estimate

% Legends, formatting

% Cardinality
figure(n_fig);
legend([{'Truth'} labels]);
title('Cardinality');
xlabel('Time');

% Cardinality error
figure(n_err_fig);
legend(labels);
title('Cardinality Error');
xlabel('Time');

% Ospa
figure(o_fig);
subplot(3, 1, 1)
legend(labels);
title('OSPA Metrics')
xlabel 'Time step'
ylabel 'OSPA Distance'
set(gca, 'Fontsize', 18)
subplot(3, 1, 2)
xlabel 'Time step'
ylabel 'OSPA Location'
set(gca, 'Fontsize', 18)
subplot(3, 1, 3)
xlabel 'Time step'
ylabel 'OSPA Cardinality'
set(gca, 'Fontsize', 18)

% Print stats table

% Lambda estimate from lcphd

% Lambda and pD estimates from lpdcphd


end