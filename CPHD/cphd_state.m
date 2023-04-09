classdef cphd_state
    %CPHD_STATE State of a Cardinalized Probability Hypothesis Density
    % filter at some time

    properties
        v           % Current "intensity" of the RFS being estimated
        N0          % Estimated number of clutter generators
        rho         % Hybrid cardinality distribution, 1xn
                    % Assumes distribution is (potentially) nonzero on 
                    % some range [0,n], but definitely 0 elsewhere
        
        Z           % Measurement set UxJz
        H           % Measurement model, z = Hx, UxN
        R           % Measurement noise covariance, UxU
    end

    methods

    end
end