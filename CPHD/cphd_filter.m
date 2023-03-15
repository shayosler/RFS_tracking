function [v_k, N_k, Xhat] = cphd_filter(v, N, rho, F, Q, ps, pd, gamma, rho_gamma, Z, H, R, kappa, U, T, Jmax, w_min)
% Cardinalized PHD filter for tracking in 2d, with unknown clutter model
% Inputs:
%   v       Current "intensity" of the RFS being estimated
%   N       N estimated number of clutter generators
%   rho     Hybrid cardinality distribution, 1xn
%           Assumes distribution is (potentially) nonzero on some range [0,
%           n], but definitely 0 elsewhere
%   F       State transition matrix x[k+1] = Fx[k]
%   Q       Process noise, 2x2
%   ps      Survival probability, state independent [clutter survival; target survival]
%   pd      Detection probability [clutter detection prob; target detection prob]
%   gamma   birth model intensity {mean number of clutter births (scalar); target birth model}
%   rhogamma birth model cardinality distribution
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
%if isscalar(ps)
%    ps = @(m) ps;
%end

% separate out augmented states
Ngamma0 = gamma(1);
gamma1 = gamma(2);
ps0 = ps(1);
ps1 = ps(2);
pd0 = pd(1);
pd1 = pd(2);


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
v_kk1 = GMRFS(m_kk1, P_kk1, ps1.*v.w) + gamma;
P_kk1 = v_kk1.P;
m_kk1 = v_kk1.m;

% Predict number of clutter generators: clutter birth + clutter survival
N_kk1 = Ngamma0 + ps0 * N;

% Predict hybrid cardinality distribution
phi = (ps1 * sum(v.w + ps0 * N)) / (sum(v.w + N));
rho_kk1 = zeros(size(rho));

% TODO: I think this should work, double check for off by one errors
for ndd = 1:size(rho, 2)
    rho_kk1_n = 0;
    for j = 1:ndd
        if j - ndd > size(rho_gamma, 2)
            continue;
        end

        binom_sum = 0;
        for l = j:size(rho, 2)
            binom_sum = binom_sum + nchoosek(l, j) * rho(l) * (1-phi)^(l-j) * phi^j;
        end
        rho_kk1_n = rho_kk1_n + rho_gamma(ndd - j) * binom_sum;
    end
    rho_kk1(ndd) = rho_kk1_n;
end


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

sum_vD = GMRFS();
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
        wkz(j) = (pd1 * v_kk1.w(j) * qz) / (pd0 * kappaz * N_kk1 + pd1 * sum_likelihood);
        if wkz(j) >= 1 || wkz(j) <= 0
            warning('Suspicious weight')
        end
    end
    sum_vD = sum_vD + GMRFS(m_kkz, P_kk, wkz);
end

v_k_unpruned = (1 - pd) .* v_kk1 + sum_vD;
N_k = N_kk1 * (1-pd) + sum(sum_vD.w);

% TODO: pruning and estimate extraction has been moved out for now
% I think that once there is functionality to only filter over the
% current FOV implemented in this function I can move it back in
v_k = v_k_unpruned;
Xhat = [];
% %% Prune
% v_k = prune_gmphd(v_k_unpruned, T, U, Jmax);
% 
% % Renormalize the weights in v_k so that N_k stays the same
% % TODO: this wasn't done in the paper, but it seems like a step that should
% % happen
% sum_pruned = sum(v_k.w);
% v_k = N_k / sum_pruned .* v_k;
% 
% %% Extract a state estimate from the RFS
% Xhat = [];
% for i = 1:v_k.J
%     if v_k.w(i) > w_min
%         Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
%     end
% end

end

function P = permnj(n, j)
% P = permnj(n, j) Compute a permutation coefficient
% P = n! / (n -j)!
P = nchoosek(n, j) * factorial(j);
end
