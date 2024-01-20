function d = hellinger_beta(s1, t1, s2, t2)
%HELLINGER_BETA Square Hellinger distance between beta distributions
% d = HELLINGER_BETA(s1, t1, s2, t2) calculates the Squared Hellinger 
% distance between two beta distributions defined by parameters
% s, and t
% 
% Inputs:
%
% Outputs:

d = 1 - (beta(mean([s1, s2]), mean([t1, t2])) / sqrt(beta(s1, t1) * beta(s2, t2)));

end
