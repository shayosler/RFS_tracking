function [state_k, Xhat, lambda] = cphd_filter(state, measurement, model, params)
% Cardinalized PHD filter for tracking with unknown clutter model
% Inputs:
%   state       cphd_state object containing the previous state of the system
%   measurement cphd_measurement object containing the measurement(s) for 
%               the current time step
%   model       cphd_model object defining the model for the system
%   params      cphd_params object defining the system parameters
%   
% Outputs:
%   state_k     cphd_state object containing the estimated current state of
%               the system
%   lambda      Estimated clutter rate

%if isscalar(pd)
%    pd = @(m) pd;
%end
%if isscalar(ps)
%    ps = @(m) ps;
%end

% Extract state
v = state.v;
N0 = state.N0;
rho = state.rho;

% Extract model data
Ngamma0 = model.Ngamma0;
gamma1 = model.gamma1;
rho_gamma1 = model.rho_gamma1;
ps0 = model.ps0;
ps1 = model.ps1;
pd0 = model.pd0;
pd1 = model.pd0;
Q = model.Q;
F = model.F;

% Extract params
Jmax = params.Jmax;
T = params.T;
U = params.U;
w_min = params.w_min;

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
v_kk1 = GMRFS(m_kk1, P_kk1, ps1.*v.w) + gamma1;
P_kk1 = v_kk1.P;
m_kk1 = v_kk1.m;

% Predict number of clutter generators: clutter birth + clutter survival
N0_kk1 = Ngamma0 + ps0 * N0;

% Predict hybrid cardinality distribution
phi = (ps1 * sum(v.w + ps0 * N0)) / (sum(v.w + N0));
rho_kk1 = zeros(size(rho));

% TODO: I think this should work, double check for off by one errors
% rho(j) = P(N==j-1
for ndd = 1:size(rho, 2)
    rho_kk1_n = 0;
    for j = 1:ndd
        if j - ndd > size(rho_gamma1, 2)
            continue;
        end

        binom_sum = 0;
        for l = j:size(rho, 2)
            binom_sum = binom_sum + nchoosek(l, j) * rho(l) * (1-phi)^(l-j) * phi^j;
        end
        rho_kk1_n = rho_kk1_n + rho_gamma1(ndd - j) * binom_sum;
    end
    rho_kk1(ndd) = rho_kk1_n;
end

%% Update

% Extract measurement
Z = measurement.Z;
R = measurement.R;
H = measurement.H;

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

% Update cardinality distribution
sum_w_kk1 = 0;
for l = 1:Jkk1
    sum_w_kk1 = sum_w_kk1 + v_kk1.w(l);
end
phi_kk1 = 1 - (pd1 * sum_w_kk1 + pd0 * N0_kk1) / (sum_w_kk1 + N0_kk1);
% I'm pretty sure the denominator of this equation in the paper is
% just a normalizing term
rho_k = zeros(size(rho_kk1));
psi0_k = zeros(size(rho_kk1));
psi1_k = zeros(size(rho_kk1));
for ndd = 1:size(rho_kk1, 2)
    psi0_k = psi_ddot(0, phi_kk1, Jz, ndd);
    psi1_k = psi_ddot(1, phi_kk1, Jz, ndd);
    if ndd < Jz
        rho_k(ndd) = 0;
    else
        rho_k(ndd) = rho_kk1(ndd) * psi_ddot(0, phi_kk1, Jz, ndd);
    end
end
rho_k = rho_k ./ sum(rho_k);

% Update intensity RFS
sum_vD = GMRFS();
sum_Nk = 0;
for jz = 1:Jz
    % Apply each measurement to each element of the RFS, basically
    % generating a new gaussian mixture for each measurement
    z = Z(:, jz);

    % Clutter likelihood
    % TODO: per measurement clutter probabilities?
    % TODO: non-uniform clutter probabilities
    kappaz = kappa;

    % Sum over all of the likelihoods of z given predicted m_kk1
    % weighted by wk
    sum_likelihood = 0;
    for l = 1:Jkk1
        m_kk1l = m_kk1(:, l);
        P_kk1l = P_kk1(:, :, l);
        wkl = v_kk1.w(l);
        sum_likelihood = sum_likelihood + wkl * mvnpdf(z, H*m_kk1l, R + H*P_kk1l*H');
    end

    % Denominator used for normalizing weights
    weight_denom = (pd0 * kappaz * N0_kk1 + pd1 * sum_likelihood);

    sum_Nk = sum_Nk + (pd0 * kappaz) / weight_denom;

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

        % Updated weights
        wkz(j) = (pd1 * v_kk1.w(j) * qz) / weight_denom;
        if wkz(j) >= 1 || wkz(j) <= 0
            warning('Suspicious weight')
        end

    end
    sum_vD = sum_vD + GMRFS(m_kkz, P_kk, wkz);
end

% This ratio of inner products gets reused:
psi_rho_ratio = sum(psi1_k .* rho_kk1) / sum(psi0_k .* rho_kk1);
psi_rho_weights = psi_rho_ratio / (sum_w_kk1 + N0_kk1);

% qd = 1-pd
v_k_unpruned = ((1 - pd1) * psi_rho_weights) .*  v_kk1 + sum_vD;
N0_k = N0_kk1 * ((1 - pd0) * psi_rho_weights + sum_Nk);

%% Prune
v_k = prune_gmphd(v_k_unpruned, T, U, Jmax); 
 
%% Extract a state estimate from the RFS
Xhat = [];
for i = 1:v_k.J
    if v_k.w(i) > w_min
        Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
    end
end

%% Outputs
state_k = cphd_state();
state_k.N0 = N0_k;
state_k.v = v_k;
state_k.rho = rho_k;

lambda = N0_k * pd0;

end


%% Helper Functions

function P = permnj(n, j)
% P = permnj(n, j) Compute a permutation coefficient
% P = n! / (n -j)!
P = nchoosek(n, j) * factorial(j);
end

function out = psi_ddot(u, phi, Jz, n)
% PSI double dot function. I think that |Zk| in the paper
% is referring the the cardinality of the observation RFS
Zk_u = Jz + u;

if n < Jz + u
    out = 0;
    return;
end

out = permnj(n, Zk_u) * (phi ^ (n - Zku));
end
