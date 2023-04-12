classdef cphd_state
    %CPHD_STATE State of a Cardinalized Probability Hypothesis Density
    % filter at some time

    properties
        v           % Current "intensity" of the RFS being estimated
        N0          % Estimated number of clutter generators
        N1          % Estimated number of targets
        rho         % Hybrid cardinality distribution, 1xn
                    % Assumes distribution is (potentially) nonzero on 
                    % some range [0,n], but definitely 0 elsewhere
    end

    methods
        function obj = cphd_state(v, N0, N1, rho)
            if nargin == 4                
                obj.v = v;
                obj.N0 = N0;
                obj.N1 = N1;
                obj.rho = rho;
            elseif nargin == 0
                obj.v = GMRFS();
                obj.N0 = 0;
                obj.N1 = 0;
                obj.rho = [];
            else
                error('Invalid number of arguments')
            end
        end
    end
end