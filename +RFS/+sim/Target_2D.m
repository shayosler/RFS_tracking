classdef Target_2D
    % A target occupying 2 spatial dimensions

    properties
        X   (4, 1) double {mustBeReal} = zeros(4, 1) % Current target state [n ndot e edot]'
        pd  (1, 1) double {mustBeReal} = 1  % Target detection probability, scaling factor applied to 
                                            % detection probability from sensor
        F   (4, 4) double {mustBeReal} = eye(4)     % Target dynamics                                 
        Q   (4, 4) double {mustBeReal} = zeros(4)   % Target process noise
        t_birth (1, 1) double {mustBeReal} = 0      % Target birth time
        t_death (1, 1) double {mustBeReal} = Inf    % Targer death time

    end % properties

    properties (Access = private)
        trajectory
    end

    methods
        function [n, e] = get_position(this)
            %get_position Get the target position(s)
            % [n, e] = get_position
            % Outputs:
            %   n       Target northing
            %   e       Target easting
            %
            % P = get_position
            % Outputs:
            %   P       Target position [northing; easting]
            state = [this.X];

            northing = state(1, :);
            easting = state(3, :);
            if nargout == 1
                n = [northing; easting];
            else
                n = northing;
                e = easting;
            end
            return 
        end

        function objs = step(objs)
            %step(this) Propagate the target one time step forward
            for i = 1:numel(objs)
                obj = objs(i);
                obj.trajectory = [obj.trajectory obj.X];
                obj.X = objs(i).F * objs(i).X + mvnrnd([0, 0, 0, 0], objs(i).Q)';
                objs(i) = obj;
            end            
        end

        function traj = get_trajectory(this)
            %get_trajectory Get the target's trajectory
            traj = [this.trajectory this.X];
        end

    end % methods

end