classdef Sonar_RB
    % A forward looking sonar the returns range/bearing measurements of
    % targets corrupted by 0 mean gaussian noise. The sensor will also
    % produce false positive "clutter" returns that are modelled as a 
    % poisson process evenly distributed over the field of view
    properties
        range (1, 1) double {mustBeReal} = 30           % End range of the sonar
        fov (1, 1) double {mustBeReal} = 90             % Total sonar field of view, degrees
        sigma_range (1, 1) double {mustBeReal} = 0      % Range noise standard deviation
        sigma_bearing (1, 1) double {mustBeReal} = 0    % Bearing noise standard deviation
        lambda (1, 1) double {mustBeReal} = 0           % Mean number of clutter returns
        pd (1, 1) double {mustBeReal} = 1               % Probability that this sonar detects a visible target
    end % properties

    methods
        function [r, b, r_true, b_true] = measure(this, x, targets)
            % [r, b, r_true, b_true] = measure(x, targets)
            % Simulate sonar detections of targets. Sonar is assumed to be oriented
            % such that it looks directly forward, and has a symmetric field of view
            % Inputs:
            %   x           Current vehicle position [x, y, psi]
            %   targets     List of all targets Kx1
            % Outputs:
            %   r         "Measured" ranges to observed targets and clutter, Mx1
            %   b         "Measured" bearings to observed targets and clutter, degrees relative to current
            %             vehicle heading, Mx1
            %   r_true    True ranges to all visible targets, Kx1
            %   b_true    True bearings to all visible targets, Kx1

            % Validate args
            if numel(x) ~= 3
                error('x must be of length 3')
            end
            if ~isvector(targets)
                error('targets must be 1xN or Nx1')
            end

            % Generate map of all existing targets
            map = zeros(length(targets), 1); % [n e pd; ...]
            for k = 1:length(targets)
                tgt = targets(k);
                map(k, :) = [tgt.state(1) tgt.state(3) tgt.pd * this.pd];
            end

            % Calculate range/bearing to each target
            n_obj = size(map, 1);
            offset = map(:, 1:2) - repmat([x(1) x(2)], [n_obj, 1]);
            r_true = sqrt(sum(offset.^2, 2));
            b_true = atan2d(offset(:, 2), offset(:, 1)) - x(3);

            % Determine which targets are within the field of view
            % and which ones are actually detected
            in_fov = r_true < this.range & abs(b_true) < this.fov/2;
            detected = map(:, 3) > rand(n_obj, 1);
            r = r_true(in_fov & detected);
            b = b_true(in_fov & detected);
            r_true = r_true(in_fov);
            b_true = b_true(in_fov);

            % Corrupt observations
            noise = mvnrnd([0, 0], diag([this.sigma_range.^2 this.sigma_bearing.^2]), length(r));
            r = r + noise(:, 1);
            b = b + noise(:, 2);

            % False positives
            % normalize fp
            if this.lambda > 0
                n_fp = poissrnd(this.lambda);
                b_fp = (rand(n_fp, 1)) * this.fov - (this.fov / 2);
                r_fp = rand(n_fp, 1) * end_range;
                % Output is observed objects and random false positives
                r = [r; r_fp];
                b = [b; b_fp];
            end
        end % measure()

    end % methods

end