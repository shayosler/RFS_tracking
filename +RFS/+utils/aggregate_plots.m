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
n_estimates = length(estimates);
N = cell(1, length(estimates));
N_err = cell(1, length(estimates));
ospa = cell(1, length(estimates));
ospa_l = cell(1, length(estimates));
ospa_n = cell(1, length(estimates));
labels = cell(1, length(estimates));
stats_cols = {'RMS $N_{err}$', 'Mean $N_{err}$', 'Std $N_{err}$', 'Mean OSPA', 'Mean OSPA Loc.', 'Mean OSPA Card.'};
n_cols = length(stats_cols);
stats = zeros(n_estimates, n_cols);
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
    plot(abs(mean(N_err{k}, 2)), styles{k}, width{:});
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

    rms_n_err = rms(N_err{k}, "all");
    mean_n_err = mean(mean(N_err{k}));
    std_n_err = std(N_err{k}, 0, 'all');
    mean_ospa = mean(mean(ospa{k}));
    mean_ospa_l = mean(mean(ospa_l{k}));
    mean_ospa_n = mean(mean(ospa_n{k}));
    fprintf('%s RMS Cardinality Error: %0.3f\n', estimates{k}.label, rms_n_err);
    fprintf('%s Mean Cardinality Error: %0.3f\n', estimates{k}.label, mean_n_err);
    fprintf('%s Std Cardinality Error: %0.3f\n', estimates{k}.label, std_n_err);
    fprintf('%s Mean OSPA: %0.3f\n', estimates{k}.label, mean_ospa);
    fprintf('%s Mean OSPA Loc.: %0.3f\n', estimates{k}.label, mean_ospa_l);
    fprintf('%s Mean OSPA Card.: %0.3f\n', estimates{k}.label, mean_ospa_n);
    stats(k, :) = [rms_n_err mean_n_err std_n_err mean_ospa mean_ospa_l mean_ospa_n ];
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
ylabel 'OSPA Loc.'
set(gca, 'Fontsize', 18)
subplot(3, 1, 3)
xlabel 'Time step'
ylabel 'OSPA Card.'
set(gca, 'Fontsize', 18)

% Print stats table
% first row, black col then stats labels
indent = "  ";

row1 = indent;
col_align = " c ";
for c = 1:n_cols
    row1 = row1 + " & " + stats_cols{c};
    col_align = col_align + "c ";
end
fprintf("\\begin{tabular}{%s}\n", col_align);
row1 = row1 + "\\";
fprintf("%s\n", row1);
for r = 1:n_estimates
    row = indent + "$" + labels{r} + "$ ";
    for c = 1:n_cols
        row = row + "& " + num2str(stats(r, c), 3) + " ";
    end
    row = row + "\\";
    fprintf("%s\n", row);
end
fprintf("\\end{tabular}\n");

% Lambda estimate from lcphd
figure
h_lambda_hat = plot(mean([lcphd.est.lambda_hat], 2), 'g', 'LineWidth', 2);
hold on
h_lambda_true = plot(truth(1).lambda,  'b--', 'LineWidth', 2);
title(lcphd.label + " \lambda Estimate")
xlabel 'Time'
ylabel '\lambda'
legend('Estimate', 'Truth')
set(gca, 'Fontsize', 18)
h = [h h_lambda_hat h_lambda_true];

% Lambda estimate from lpdcphd
figure
h_lambda_hat = plot(mean([lpdcphd.est.lambda_hat], 2), 'g', 'LineWidth', 2);
hold on
h_lambda_true = plot(truth(1).lambda,  'b--', 'LineWidth', 2);
title(lpdcphd.label + " \lambda Estimate")
xlabel 'Time'
ylabel '\lambda'
legend('Estimate', 'Truth')
set(gca, 'Fontsize', 18)
h = [h h_lambda_hat h_lambda_true];

% Lambda and pD estimates from lpdcphd
% Plot true vs estimated detection probability
figure
h_pd0_mean = plot(mean([lpdcphd.est.pd0_hat],2), 'g', 'LineWidth', 2);
hold on
h_pd0_true = plot(truth(1).pd0,  'b--', 'LineWidth', 2);
title('p_D^{(0)} Estimate')
xlabel 'Time'
ylabel 'p_D^{(0)}'
legend('Estimate', 'Truth')
set(gca, 'Fontsize', 18)
h = [h h_pd0_mean h_pd0_true];

figure
h_pd1_mean = plot(mean([lpdcphd.est.pd1_hat],2), 'g', 'LineWidth', 2);
hold on
h_pd1_true = plot(truth(1).pd1,  'b--', 'LineWidth', 2);
title(' p_D^{(1)} Estimate')
xlabel 'Time'
ylabel 'p_D^{(1)}'
legend('Estimate', 'Truth')
set(gca, 'Fontsize', 18)
h = [h h_pd1_mean h_pd1_true];


end