classdef Target_2D
    % A target occupying 2 spacial dimensions

    properties
        X   (4, 1) double {mustBeReal} = zeros(4, 1) % Current target state [n ndot e edot]'
        pd  (1, 1) double {mustBeReal} = 1  % Target detection probability, scaling factor applied to 
                                            % detection probability from sensor
        F   (4, 4) double {mustBeReal} = eye(4)     % Target dynamics                                 
        Q   (4, 4) double {mustBeReal} = zeros(4)   % Target process noise
        t_birth (1, 1) double {mustBeReal} = 0      % Target birth time
        t_death (1, 1) double {mustBeReal} = Inf    % Targer death time

    end % properties

    methods
        function [n, e] = get_position(this)
            %get_position Get the target position
            % [n, e] = get_position
            % Outputs:
            %   n       Target northing
            %   e       Target easting
            %
            % P = get_position
            % Outputs:
            %   P       Target position [northing, easting]
            northing = this.X(1);
            easting = this.X(3);
            if nargout == 1
                n = [northing, easting];
            else
                n = northing;
                e = easting;
            end
            return 
        end

        function obj = step(this)
            %step(this) Propagate the target one time step forward
            obj = this;
            obj.X = this.F * this.X + mvnrnd([0, 0, 0, 0], sqrtm(this.Q));
        end

    end % methods

end