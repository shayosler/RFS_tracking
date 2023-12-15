function [h] = lpdcphd_plots(figures, prefix, k, targets, truth, trackers, states)
%lpdcphd_plots Plots specific to the lambda-pd-CPHD filter

% plot true number of targets, estimated number of targets, and estimated
num = length(states);
pd0_mean = zeros(num, 1);
pd1_mean = zeros(num, 1);
pd0_wmean = zeros(num, 1);
pd1_wmean = zeros(num, 1);

for j = 1:num
    pd0_mean(j) = mean(states(j).pd0);
    pd1_mean(j) = mean(states(j).pd1);
    
    if states(j).v0.J == 0
        pd0_wmean(j) = NaN;
    else
        pd0_wmean(j) = states(j).v0.w' * states(j).pd0 / sum(states(j).v0.w);
    end
    if states(j).v1.J == 0
        pd1_wmean(j) = NaN;
    else
        pd1_wmean(j) = states(j).v1.w' * states(j).pd1 / sum(states(j).v1.w);
    end
end

% Plot true vs estimated detection probability
figure(figures.pd0_fig)
h_pd0_mean = plot(pd0_mean, 'g', 'LineWidth', 2);
hold on
h_pd0_wmean = plot(pd0_wmean, 'r', 'LineWidth', 2);
h_pd0_true = plot(truth.pd0,  'b--', 'LineWidth', 2);
title([prefix ' p_D^{(0)} after t = ' num2str(k)])
xlabel 'Time step'
ylabel 'p_D^{(0)}'
legend('Estimate', 'Estimate (weighted)', 'Truth')
set(gca, 'Fontsize', 18)

figure(figures.pd1_fig)
h_pd1_mean = plot(pd1_mean, 'g', 'LineWidth', 2);
hold on
h_pd1_wmean = plot(pd1_wmean, 'r', 'LineWidth', 2);
h_pd1_true = plot(truth.pd1,  'b--', 'LineWidth', 2);
title([prefix ' p_D^{(1)} after t = ' num2str(k)])
xlabel 'Time step'
ylabel 'p_D^{(1)}'
legend('Estimate', 'Estimate (weighted)', 'Truth')
set(gca, 'Fontsize', 18)


% Plot cardinality distribution
h = [h_pd0_mean h_pd0_wmean h_pd0_true h_pd1_mean h_pd1_wmean h_pd1_true];
end