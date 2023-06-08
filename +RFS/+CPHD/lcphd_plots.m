function [h] = lcphd_plots(figures, k, targets, trackers, lambda, lambda_true)
%lcphd_plots Plots specific to the lambda-CPHD filter

% plot true number of targets, estimated number of targets, and estimated
% number of clutter generators

% Plot true vs estimated clutter rate (lambda)
figure(figures.lambda_fig)
plot(lambda, 'g', 'LineWidth', 2);
hold on
plot(lambda_true,  'b--', 'LineWidth', 2);
title(['l-CPHD \lambda after t = ' num2str(k)])
xlabel 'Time step'
ylabel '\lambda'
legend('Estimate', 'True')
set(gca, 'Fontsize', 18)

% Plot cardinality distribution
end