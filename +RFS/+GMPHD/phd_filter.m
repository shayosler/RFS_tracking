function [v_k, N_k, Xhat] = phd_filter(v, N, F, Q, ps, pd, gamma, Z, H, R, kappa, U, T, Jmax, w_min)
% PHD filter for tracking in 2d, assumes spawn probability = 0;
% Inputs:
%   v       Current "intensity" of the RFS being estimated
%   N       N estimated number of targets
%   F       State transition matrix x[k+1] = Fx[k]
%   Q       Process noise, 2x2
%   ps      Survival probability, state independent
%   pd      Detection probability
%   gamma   birth model intensity
%   Z       Measurement set 2xJz
%   H       Measurement model z = Hx
%   R       Measurement noise covariance, 2x2
%   kappa   Clutter intensity,
%   U       Threshold for merging components during pruning
%   T       Threshold for discarding components during pruning
%   Jmax    Maximum number of components to keep
%   w_min   Minimum weight to be included in the output state estimate
%   

%if isscalar(pd)
%    pd = @(m) pd;
%end

%% Prediction
% predicted intensity =
%   survival intensity
%   + spawn intensity (= 0)
%   + birth intensity (= gamma)

% survival intensity
% propagate each element in the RFS forward
P_k = v.P;
m_k = v.m;
m_kk1 = zeros(size(v.m));
P_kk1 = zeros(size(v.P));
Jk = v.J;
for j = 1:Jk
    m_kk1(:, j) = F*m_k(:, j);
    P_kk1(:, :, j) = Q + F*P_k(:, :, j)*F';
end

% Full prediction is sum of births and propagation
v_kk1 = RFS.utils.GMRFS(m_kk1, P_kk1, ps.*v.w) + gamma;
P_kk1 = v_kk1.P;
m_kk1 = v_kk1.m;
N_kk1 = ps * N + sum(gamma.w);

%% Update
% new estimated states is false detections + weighted sum over detection updates
%  = (1-pd)v_kk1 + sum over detections

% Compute updated covariance matrix and kalman gains, which don't depend
% on individual measurements 
% TODO: if R and/or H get updated so they can vary with each measurement 
% in Z then this will need to get moved into the other loop
Jz = size(Z, 2);
Jkk1 = v_kk1.J;
n_obs = size(Z, 1); % dimensionality of each observation
n_states = size(m_kk1, 1); % dimensionality of state vector
P_kk = zeros(size(P_k));
K = zeros(n_states, n_obs, Jkk1);
for j = 1:Jkk1
    P_kk1j = P_kk1(:, :, j);
    % Kalman gain:
    Kj = P_kk1j * (H' / (H * P_kk1j * H' + R));
    K(:, :, j) = Kj;

    % Updated estimate covariances
    %P_kk(:, :, j) = (eye(size(P_kk1j)) - K(:, :, j)*H) * P_kk1j;
    I = eye(size(P_kk1j));
    P_kk(:, :, j) = (I - Kj * H) * P_kk1j * (I - Kj * H)' + Kj * R * Kj';
end

sum_vD = RFS.utils.GMRFS();
for jz = 1:Jz
    % Apply each measurement to each element of the RFS, basically
    % generating a new gaussian mixture for each measurement
    z = Z(:, jz);

    % Sum over all of the likelihoods of z given predicted m_kk1
    % weighted by wk
    sum_likelihood = 0;
    for l = 1:Jkk1
        m_kk1l = m_kk1(:, l);
        P_kk1l = P_kk1(:, :, l);
        wkl = v_kk1.w(l);
        sum_likelihood = sum_likelihood + wkl * mvnpdf(z, H*m_kk1l, R + H*P_kk1l*H');
    end

    m_kkz = zeros(size(m_kk1));
    wkz = zeros(Jkk1, 1);
    for j = 1:Jkk1
        Kj = K(:, :, j);
        m_kk1j = m_kk1(:, j);

        % Updated mean based on this measurement
        m_kkz(:, j) = m_kk1j + Kj * (z - H*m_kk1j);
       
        % likelihood of z
        % TODO: should the q(z) calculation be using the predicted m_ and 
        % P_ or the updated ones?
        P_kk1j = P_kk1(:, :, j);
        qz = mvnpdf(z, H*m_kk1j, R + H*P_kk1j*H');

        % Clutter likelihood
        % TODO: per measurement clutter probabilities?
        kappaz = kappa;

        % Updated weights
        wkz(j) = (pd * v_kk1.w(j) * qz) / (kappaz + pd * sum_likelihood);
        if wkz(j) >= 1
            warning('Suspicious weight')
        end
    end
    sum_vD = sum_vD + RFS.utils.GMRFS(m_kkz, P_kk, wkz);
end

v_k_unpruned = (1 - pd) .* v_kk1 + sum_vD;
N_k = N_kk1 * (1-pd) + sum(sum_vD.w);

% TODO: pruning and estimating can be moved out to allow filtering only ove
% the current FOV
%v_k = v_k_unpruned;
%Xhat = [];

%% Prune
v_k = RFS.utils.prune_gmphd(v_k_unpruned, T, U, Jmax);

% Renormalize the weights in v_k so that N_k stays the same
% TODO: this wasn't done in the paper, but it seems like a step that should
% happen
sum_pruned = sum(v_k.w);
v_k = N_k / sum_pruned .* v_k;

%% Extract a state estimate from the RFS
Xhat = [];
for i = 1:v_k.J
    if v_k.w(i) > w_min
        Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
    end
end

end