classdef Sonar_XY
    % A forward looking sonar the returns X/Y measurements of
    % targets positions relative to the vehicle position, corrupted by 0 
    % mean gaussian noise. The sensor will also produce false positive 
    % "clutter" returns that are modelled as a poisson process evenly 
    % distributed over the field of view
    properties
        range (1, 1) double {mustBeReal} = 30   % End range of the sonar
        fov (1, 1) double {mustBeReal} = 90     % Total sonar field of view, degrees
        sigma_x (1, 1) double {mustBeReal} = 0  % X noise standard deviation
        sigma_y (1, 1) double {mustBeReal} = 0  % Y noise standard deviation
        lambda (1, 1) double {mustBeReal} = 0   % Mean number of clutter returns
        pd (1, 1) double {mustBeReal} = 1       % Probability that this sonar detects a visible target
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

        function [zx, zy, zx_true, zy_true] = measure(this, x, targets)
            % [zx, zy, zx_true, zy_true] = measure(this, x, targets)
            % Simulate sonar detections of targets. Sonar is assumed to be oriented
            % such that it looks directly forward, and has a symmetric field of view
            % Inputs:
            %   x           Current vehicle position [x, y, psi]
            %   targets     List of all targets Kx1
            % Outputs:
            %   zx      "Measured" X distance to observed targets and clutter, Mx1
            %   zy      "Measured" Y distance to observed targets and clutter, Mx1
            %   r_true  True X distance to to all visible targets, Nx1
            %   b_true  True Y distance to to all visible targets, Nx1
            
            zx_true = [];
            zy_true = [];
            zx = [];
            zy = [];   
            % "Measure" targets
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
                zx_true = map(in_fov, 1);
                zy_true = map(in_fov, 2);
                detected = map(:, 3) > rand(n_obj, 1);
                zx = map(detected & in_fov, 1);
                zy = map(detected & in_fov, 2);

                % Corrupt observations
                noise = normrnd([0, 0], [this.sigma_x this.sigma_y], length(zx));
                zx = zx + noise(:, 1);
                zx = zx + noise(:, 2);
            end            

            % False positives
            % normalize fp
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
                visible_clutter = this.in_fov(x_fp, y_fp);
                zx = [zx; x_fp(visible_clutter)];
                zy = [zy; y_fp(visible_clutter)];
            end
        end % measure()

    end % methods

end