classdef lpdcphd_params
    %LPDCPHD_PARAMS Parameters controlling the behavior of a Cardinalized 
    % Probability Hypothesis Density filter

    properties
        Nmax        % Maximum number of targets to track
        U0          % Distance threshold for merging clutter components during pruning        
        T0          % Weight threshold for discarding clutter components during pruning
        Jmax0       % Maximum number of clutter components to keep        
        U1          % Distance threshold for merging target components during pruning        
        T1          % Weight threshold for discarding target components during pruning
        Jmax1       % Maximum number of target components to keep
        w_min       % Minimum weight to be included in the output state estimate
    end

    methods

    end
end