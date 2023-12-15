function [w_merged, mu_merged, P_merged] = merge_gm(w, mu_in, P_in)
%merge_gm Merge components of a gaussian mixture into a single gaussian
%
% Inputs 
%   w       Weights of components to merge, Jx1
%   mu_in   Means of components to merge, NxJ
%   P_in    Covariances of components to merge, NxNxJ
%
% Outputs
%   w_merged    New weight of merged gaussian
%   mu_merged   Mean of gaussian representing merged components, Nx1
%   P_merged    Covariance of gaussian representing merged components, NxN

% Check sizes
if size(w, 2) ~= 1
    error('w must be a column vector');
end
J = size(w, 1);

if size(mu_in, 2) ~= J
    error('mu must represent the same number of components as w')
end
dim = size(mu_in, 1);

if size(P_in, 1) ~= dim || size(P_in, 2) ~= dim
    error('P must have the same dimensionality as mu')
end
if size(P_in, 3) ~= J
    error('mu must represent the same number of components as w')
end

w_merged = sum(w);
mu_merged = (1/w_merged) .* sum(repmat(w', dim, 1) .* mu_in, 2);

P_merged = zeros(size(P_in(:, :, 1)));
for i = 1:length(w)
    P_merged = P_merged + w(i)*(P_in(:, :, i) + (mu_merged - mu_in(:, i))*(mu_merged - mu_in(:, i))');
end
P_merged = P_merged ./ w_merged;
end