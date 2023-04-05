function [r, b] = simulate_obs(x, map, end_range, fov, lambda, sigma)
% [r, b] = simulate_obs(x, map, end_range, fov, fp, fn)
% Simulate sonar detections of targets. Sonar is assumed to be oriented
% such that it looks directly forward, and has a symmetric field of view
% Inputs:
%   x           Current vehicle position [x, y, psi]
%   map         Map of object locations and the probability of producing a 
%               detection from the object when its inside the FOV
%               Nx3 [x1, y1, p1; x2, y2, p2; ...; xN, yN, pN]
%   end_range   Current sonar end range
%   fov         Total sonar field of view, degrees
%   lambda      Expected number of false positives. False positives are
%               Treated as a Poisson process.
%               False positives are sampled uniform randomly over the 
%               range/bearing space in the field of view
%   sigma       Covariance matrix for range and bearing noise
% Outputs:
%   r         Ranges to observed targets, Mx1
%   b         Bearings to observed targets, degrees relative to current
%             vehicle heading, Mx1

% Validate args
if numel(x) ~= 3
    error('x must be of length 3')
end

% Calculate range/bearing to each target
n_obj = size(map, 1);
offset = map(:, 1:2) - repmat([x(1) x(2)], [n_obj, 1]);
r = sqrt(sum(offset.^2, 2));

% Determine which targets are within the field of view
% and which ones are actually detected
b = atan2d(offset(:, 2), offset(:, 1)) - x(3);
in_fov = r < end_range & abs(b) < fov/2;
detected = map(:, 3) > rand(n_obj, 1);
r = r(in_fov & detected);
b = b(in_fov & detected);

% Corrupt observations
noise = mvnrnd([0, 0], sigma, length(r));
r = r + noise(:, 1);
b = b + noise(:, 2);

% False positives
% normalize fp
if lambda > 0
    n_fp = poissrnd(lambda);
    b_fp = (rand(n_fp, 1)) * fov - (fov / 2);
    r_fp = rand(n_fp, 1) * end_range;
    % Output is observed objects and random false positives
    r = [r; r_fp];
    b = [b; b_fp];
end