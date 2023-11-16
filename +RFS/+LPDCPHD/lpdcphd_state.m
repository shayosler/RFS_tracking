classdef lpdcphd_state
    %LPDCPHD_STATE State of a lambda-pD-CPHD filter at some time

    properties
        v0  RFS.utils.BMRFS     % "Intensity" of the clutter RFS
        v1  RFS.utils.BGMRFS    % Intensity of the target RFS
        N0                      % Estimated number of clutter generators TODO: do we need this and v0?
        N1          % Estimated number of targets
        rho         % Hybrid cardinality distribution, 1xn
                    % Assumes distribution is (potentially) nonzero on 
                    % some range [0,n], but definitely 0 elsewhere
        lambda      % Estimated clutter rate
        pd1         % Estimated target detection probability
    end

    methods
        function obj = cphd_state(v0, v1, N0, N1, rho, lambda, pd1)
            if nargin == 5                
                obj.v0 = v0;
                obj.v1 = v1;
                obj.N0 = N0;
                obj.N1 = N1;
                obj.rho = rho;
                obj.lambda = lambda;
                obj.pd1 = pd1;
            elseif nargin == 0
                obj.v0 = RFS.utils.BMRFS();
                obj.v1 = RFS.utils.BGMRFS();
                obj.N0 = 0;
                obj.N1 = 0;
                obj.rho = [];
                obj.lambda = 0;
                obj.pd1 = 0;
            else
                error('Invalid number of arguments')
            end
        end
    end
end