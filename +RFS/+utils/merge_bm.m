function [w_merged, s_merged, t_merged] = merge_bm(w, s, t)
%merge_bm Merge components of a beta mixture into a single beta
% distribution
%
% Inputs
%   w
%   s
%   t
%
% Outputs
%   w_merged
%   s_merged
%   t_merged

% Check sizes
if size(w, 2) ~= 1
    error('w must be a column vector');
end
J = size(w, 1);
if ~isvector(s) || size(s, 1) ~= J
    error('s must be a column vector with the same number of elements as w')
end
if ~isvector(t) || size(t, 1) ~= J
    error('t must be a column vector with the same number of elements as w')
end

w_merged = sum(w);
sigsq_beta = s .* t ./( (s + t).^2 .* (s + t + 1) );
mu_beta = s ./ (s + t);
sigsq_beta_merged = 1 / w_merged * (w' * sigsq_beta);
mu_beta_merged = 1 / w_merged * (w' * mu_beta);

com = (mu_beta_merged * (1 - mu_beta_merged) / sigsq_beta_merged) - 1;
s_merged = com * mu_beta_merged;
t_merged = com * (1 - mu_beta_merged);

end