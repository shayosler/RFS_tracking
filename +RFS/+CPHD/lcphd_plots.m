function [h] = lcphd_plots(figures, prefix, k, targets, truth, trackers, lambda)
%lcphd_plots Plots specific to the lambda-CPHD filter

% plot true number of targets, estimated number of targets, and estimated
% number of clutter generators

% Plot true vs estimated clutter rate (lambda)
figure(figures.lambda_fig)
h_lambda_hat = plot(lambda, 'g', 'LineWidth', 2);
hold on
h_lambda_true = plot(truth.lambda,  'b--', 'LineWidth', 2);
title([prefix ' \lambda after t = ' num2str(k)])
xlabel 'Time step'
ylabel '\lambda'
legend('Estimate', 'True')
set(gca, 'Fontsize', 18)

% Plot cardinality distribution
h = [h_lambda_hat h_lambda_true];
end