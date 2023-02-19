function [v_out] = transform_gmrfs(dist, n, e, psi)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

% Rotate
rot = erm_b2e([0, 0, psi]);
rot = rot(1:2, 1:2);
mu = rot * dist.m;
P = zeros(size(dist.P));
for j = 1:dist.J
    P(:, :, j) = rot * dist.P(:, :, j) * rot';
end

% Translate
mu = mu + repmat([n; e], 1, dist.J);
v_out = GMRFS(mu, P, dist.w);

end