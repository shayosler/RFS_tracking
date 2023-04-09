classdef cphd_model
    %CPHD_MODEL System model for a Cardinalized Probability Hypothesis 
    % Density filter

    properties
        F           % State transition matrix x[k+1] = Fx[k], NxN
        Q           % Process noise, NxN
        ps0         % Clutter survival probability
        ps1         % Target survival probability
        pd0         % Clutter detection probability
        pd1         % Target detection probability
        Ngamma      % Mean number of clutter births (clutter birth assumed to be poisson distributed)
        gamma1      % Target birth model
        rho_gamma   % Target birth model cardinality distribution
        kappa       % Clutter intensity??
    end

    methods

    end
end