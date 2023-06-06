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
        function in = in_fov(this, x, y)
            %in_fov Test if a point is in the field of view
            %
            % Inputs:
            %   x       X coordinate of the point, sensor frame (fwd)
            %   y       Y coordinate of the point, sensor frame (right)
            %
            % Outputs:
            %   in      Whether or not the point (x, y) is visible
            r = sqrt(x.^2 + y.^2);
            b = atan2d(y, x);
            in = r < this.range & abs(b) < this.fov / 2;
        end

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
            %   r_true    True ranges to all visible targets, Nx1
            %   b_true    True bearings to all visible targets, Nx1

            r = [];
            b = [];
            r_true = [];
            b_true = [];
            % Measure targets
            if ~isempty(targets)
                % Validate args
                if numel(x) ~= 3
                    error('x must be of length 3')
                end
                if ~isvector(targets)
                    error('targets must be 1xN or Nx1')
                end

                % Generate map containing data we care about for of all targets
                state = [targets.X]';
                map = [state(:, 1) state(:, 2) [targets.pd]'.*this.pd];

                % Calculate range/bearing to each target
                n_obj = size(map, 1);
                offset = map(:, 1:2) - repmat([x(1) x(2)], [n_obj, 1]);
                r_true = sqrt(sum(offset.^2, 2));
                b_true = atan2d(offset(:, 2), offset(:, 1)) - x(3);

                % Determine which targets are within the field of view
                % and which ones are actually detected
                in_fov = r_true < this.range & abs(b_true) < this.fov/2;
                detected = map(:, 3) >= rand(n_obj, 1);
                r = r_true(in_fov & detected);
                b = b_true(in_fov & detected);
                r_true = r_true(in_fov);
                b_true = b_true(in_fov);

                % Corrupt observations
                noise = mvnrnd([0, 0], diag([this.sigma_range.^2 this.sigma_bearing.^2]), length(r));
                r = r + noise(:, 1);
                b = b + noise(:, 2);
            end

            % False positives (clutter)
            if this.lambda > 0
                % Generate clutter evenly over field of view
                % To do this generate clutter as if it were over a
                % rectangle enclosing the field of view, then only take the
                % points that fall in the actual field of view.

                % Calculate expected mean clutter rate over the enclosing
                % rectangle, scale actual FOV lambda by ratio of areas
                fov_min_x = 0;
                fov_max_x = this.range;
                fov_max_y = this.range * cosd(this.fov / 2);
                fov_min_y = -fov_max_y;
                A_enclosing = (fov_max_y - fov_min_y) * (fov_max_x - fov_min_x);
                A_fov = pi * this.range.^2 * this.fov / 360;

                enclosing_lambda = this.lambda * A_enclosing / A_fov;
                n_fp = poissrnd(enclosing_lambda);

                % Generate random clutter over enclosing rectangle,
                % then select only the points inside the FOV
                x_fp = fov_min_x + rand(n_fp, 1) * (fov_max_x - fov_min_x);
                y_fp = fov_min_y + rand(n_fp, 1) * (fov_max_y - fov_min_y);
                r_clut = sqrt(x_fp.^2 + y_fp.^2);
                b_clut = atan2d(y_fp, x_fp);
                visible_clutter = this.in_fov(x_fp, y_fp);
                r_clut = r_clut(visible_clutter);
                b_clut = b_clut(visible_clutter);
                r = [r; r_clut];
                b = [b; b_clut];
            end

            if this.lambda > 0 && false
                n_fp = poissrnd(this.lambda);
                b_fp = (rand(n_fp, 1)) * this.fov - (this.fov / 2);
                r_fp = rand(n_fp, 1) * this.range;
                % Output is observed objects and random false positives
                r = [r; r_fp];
                b = [b; b_fp];
            end
        end % measure()

    end % methods

end