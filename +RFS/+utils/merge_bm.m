function [s_merged, t_merged] = merge_bm(w, s, t)
%merge_bm Merge components of a beta mixture into a single beta
% distribution
%
% Inputs
%   w
%   s
%   t
%
% Outputs
%   s_merged
%   t_merged

sigsq_beta = s .* t ./( (s + t).^2 .* (s + t + 1) );
mu_beta = s ./ (s + t);
sigsq_beta_merged = 1 / wl * (w' * sigsq_beta);
mu_beta_merged = 1 / wl * (w' * mu_beta);

com = (mu_beta_merged * (1 - mu_beta_merged) / sigsq_beta_merged) - 1;
s_merged = com * mu_beta_merged;
t_merged = com * (1 - mu_beta_merged);

end