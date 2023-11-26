function [mu, P] = merge_gm(w, mu_in, P_in)
%merge_gm Merge components of a gaussian mixture into a single gaussian
%
% Inputs 
%   w       Weights of components to merge, Jx1
%   mu_in   Means of components to merge, NxJ
%   P_in    Covariances of components to merge, NxJ
%
% Outputs
%   mu      Mean of gaussian representing merged components
%   P       Covariance of gaussian representing merged components

% Check sizes
if size(w, 2) ~= 1
    error('w must be a column vector');
end
J = size(w, 1);

if size(mu_in, 2) ~= J
    error('mu must represent the same number of components as w')
end

w_merged = sum(w);
dim = size(mu_in, 1);
mu = (1/w_merged) .* sum(repmat(w', dim, 1) .* mu_in, 2);

P = zeros(size(P_in(:, :, 1)));
for i = 1:length(w)
    P = P + w(i)*(P_in(:, :, i) + (mu - mu_in(:, i))*(mu - mu_in(:, i))');
end
P = P ./ w_merged;
end