function [v_k, N_k, Xhat] = cphd_filter(v, N, F, Q, ps, pd, gamma, Z, H, R, kappa, U, T, Jmax, w_min)
% Cardinalized PHD filter for tracking in 2d
% Inputs:
%   v       Current "intensity" of the RFS being estimated
%   N       N estimated number of targets
%   F       State transition matrix x[k+1] = Fx[k]
%   Q       Process noise, 2x2
%   ps      Survival probability, state independent
%   pd      Detection probability
%   gamma   birth model intensity
%   Z       Measurement set 2xJz
%   H       Measurement model z = Hx
%   R       Measurement noise covariance, 2x2
%   kappa   Clutter intensity,
%   U       Threshold for merging components during pruning
%   T       Threshold for discarding components during pruning
%   Jmax    Maximum number of components to keep
%   w_min   Minimum weight to be included in the output state estimate
%   

%if isscalar(pd)
%    pd = @(m) pd;
%end

%% Prediction

%% Update


% TODO: pruning and estimate extraction has been moved out for now
% I think that once there is functionality to only filter over the
% current FOV implemented in this function I can move it back in
v_k = v_k_unpruned;
Xhat = [];
% %% Prune
% v_k = prune_gmphd(v_k_unpruned, T, U, Jmax);
% 
% % Renormalize the weights in v_k so that N_k stays the same
% % TODO: this wasn't done in the paper, but it seems like a step that should
% % happen
% sum_pruned = sum(v_k.w);
% v_k = N_k / sum_pruned .* v_k;
% 
% %% Extract a state estimate from the RFS
% Xhat = [];
% for i = 1:v_k.J
%     if v_k.w(i) > w_min
%         Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
%     end
% end

end