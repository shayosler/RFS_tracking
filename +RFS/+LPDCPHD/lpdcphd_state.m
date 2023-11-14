classdef lpdcphd_state
    %LPDCPHD_STATE State of a lambda-pD-CPHD filter at some time

    properties
        v0  RFS.utils.BMRFS     % "Intensity" of the clutter RFS
        v1  RFS.utils.BGMRFS    % Intensity of the target RFS
        N0                      % Estimated number of clutter generators
        N1          % Estimated number of targets
        rho         % Hybrid cardinality distribution, 1xn
                    % Assumes distribution is (potentially) nonzero on 
                    % some range [0,n], but definitely 0 elsewhere
    end

    methods
        function obj = cphd_state(v0, v1, N0, N1, rho)
            if nargin == 5                
                obj.v0 = v0;
                obj.v1 = v1;
                obj.N0 = N0;
                obj.N1 = N1;
                obj.rho = rho;
            elseif nargin == 0
                obj.v0 = RFS.utils.BMRFS();
                obj.v1 = RFS.utils.BGMRFS();
                obj.N0 = 0;
                obj.N1 = 0;
                obj.rho = [];
            else
                error('Invalid number of arguments')
            end
        end
    end
end