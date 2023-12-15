classdef cphd_model
    %CPHD_MODEL System model for a Cardinalized Probability Hypothesis 
    % Density filter

    properties
        F           % State transition matrix x[k+1] = Fx[k], NxN
        Q           % Process noise, NxN
        ps0         % Clutter generator survival probability
        ps1         % Target survival probability
        pd0         % Clutter generator detection probability
        pd1         % Target detection probability
        Ngamma0     % Mean number of clutter generator births
        gamma1      % Target birth model, filter assumes birth process is Poisson
        kappa       % Spatial likelihood of clutter
    end

    methods

    end
end