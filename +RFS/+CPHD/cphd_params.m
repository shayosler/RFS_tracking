classdef cphd_params
    %CPHD_PARAMS Parameters controlling the behavior of a Cardinalized 
    % Probability Hypothesis Density filter

    properties
        Nmax        % Maximum number of targets to track
        U           % Threshold for merging components during pruning        
        T           % Threshold for discarding components during pruning
        Jmax        % Maximum number of components to keep
        w_min       % Minimum weight to be included in the output state estimate
    end

    methods

    end
end