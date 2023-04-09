classdef cphd_params
    %CPHD_PARAMS Parameters controlling the behavior of a Cardinalized 
    % Probability Hypothesis Density filter

    properties
        T_merge     % Threshold for merging components during pruning        
        T_discard   % Threshold for discarding components during pruning
        Jmax        % Maximum number of components to keep
        w_min       % Minimum weight to be included in the output state estimate
    end

    methods

    end
end