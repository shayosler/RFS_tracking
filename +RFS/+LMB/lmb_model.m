classdef lmb_model
    %LMB_MODEL System model for a Labeled Multi-Bernoulli filter

    properties
        F           % State transition matrix x[k+1] = Fx[k], NxN
        Q           % Process noise, NxN
        kB          % "Dilation" constant used to increase variance of beta
                    % distributions during prediction. Must be > 0 
                    % (usually > 1)
        ps0         % Clutter generator survival probability
        ps1         % Target survival probability
        gamma0  RFS.utils.BMRFS     % Clutter birth model
        gamma1  RFS.utils.BGMRFS    % Target birth model, filter assumes birth process is Poisson
        kappa       % Spatial likelihood of clutter
    end

    methods

    end
end