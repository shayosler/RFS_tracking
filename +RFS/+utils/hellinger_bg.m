function d = hellinger_bg(m1, P1, s1, t1, m2, P2, s2, t2)
%HELLINGER Hellinger distance between beta gaussian distributions
% d = HELLINGER(m1, P1, s1, t1, m2, P2, s2, t2) calculates the Hellinger 
% distance between two beta-gaussian distributions defined by parameters
% m, P, s, and t
% 
% Inputs:
%
% Outputs:

% TODO: make this accept a vector for each of m, P, s, t
% Beta component of distance
b = beta(mean([s1, s2]), mean([t1, t2])) / sqrt(beta(s1, t1) * beta(s2, t2));

% gaussian component of distance
z = zeros(size(m1));
g = sqrt(mvnpdf(z, m1 - m2, P1 + P2)) / det(8 * pi * (inv(P1) + inv(P2)))^.25;

% Compute distance
d = sqrt(1 - b * g);

end