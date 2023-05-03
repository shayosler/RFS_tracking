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
if ~isvector(rho)
    error 'state.rho must by a 1xN or Nx1 vector'
end

% Extract model data
Ngamma0 = model.Ngamma0;
gamma1 = model.gamma1;
rho_gamma1 = model.rho_gamma1;
if ~isvector(rho_gamma1)
    error 'model.rho_gamma1 must by a 1xN or Nx1 vector'
end
ps0 = model.ps0;
ps1 = model.ps1;
pd0 = model.pd0;
pd1 = model.pd1;
Q = model.Q;
F = model.F;
kappa = model.kappa;

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
v_kk1 = gamma1 + GMRFS(m_kk1, P_kk1, ps1.*v.w); %GMRFS(m_kk1, P_kk1, ps1.*v.w) + gamma1; % Order swapped to make comparison to vo easier
P_kk1 = v_kk1.P;
m_kk1 = v_kk1.m;

% Predict number of clutter generators: clutter birth + clutter survival
N0_kk1 = Ngamma0 + ps0 * N0;

% Predict hybrid cardinality distribution
if sum(v.w) + N0 == 0
    phi = 0;
else
    phi = (ps1 * sum(v.w) + ps0 * N0) / (sum(v.w) + N0);
end

% Variables as named by vo:
w_update = v.w;
Nc_update = N0;
vo_model.P_S = ps1;
vo_model.clutter_P_S = ps0;
filter.N_max = length(rho) - 1;
vo_model.w_birth = model.gamma1.w;
vo_model.lambda_cb = model.Ngamma0;
cdn_update = rho;

%---cardinality prediction
%surviving cardinality distribution
survive_cdn_predict = zeros(filter.N_max+1,1);
survival_factor= sum(w_update)/(sum(w_update)+Nc_update)*vo_model.P_S + Nc_update/(sum(w_update)+Nc_update)*vo_model.clutter_P_S;
if isnan(survival_factor), survival_factor=0; end %catch the degernerate zero case

% TODO: this if block added in to help prevent NaNs - so
if survival_factor == 0
    log_survival_factor = log(eps(0));
else
    log_survival_factor = log(survival_factor);
end

for j=0:filter.N_max
    idxj=j+1;
    terms= zeros(filter.N_max+1,1);
    for ell=j:filter.N_max
        idxl= ell+1;
        terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log_survival_factor+(ell-j)*log(1-survival_factor))*cdn_update(idxl);
    end
    survive_cdn_predict(idxj) = sum(terms);
end

%predicted cardinality= convolution of birth and surviving cardinality distribution
cdn_predict = zeros(filter.N_max+1,1);
lambda_b = sum(vo_model.w_birth)+vo_model.lambda_cb;
for n=0:filter.N_max
    idxn=n+1;
    terms= zeros(filter.N_max+1,1);
    for j=0:n
        idxj= j+1;
        terms(idxj)= exp(-sum(lambda_b)+(n-j)*log(lambda_b)-sum(log(1:n-j)))*survive_cdn_predict(idxj);
    end
    cdn_predict(idxn) = sum(terms);
end
%normalize predicted cardinality distribution
cdn_predict = cdn_predict/sum(cdn_predict);
rho_kk1 = cdn_predict;

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
    I = eye(size(P_kk1j));
    P_kk(:, :, j) = (I - Kj * H) * P_kk1j * (I - Kj * H)' + Kj * R * Kj';
end

% Phi_kk1 is 1 - ratio of expected number of detected objects (targets + clutter)
% to estimated  total number of objects. ie its the expected fraction of 
% total objects that are "missed" by the observation
sum_w_kk1 = 0;
for l = 1:Jkk1
    sum_w_kk1 = sum_w_kk1 + v_kk1.w(l);
end
phi_kk1 = 1 - (pd1 * sum_w_kk1 + pd0 * N0_kk1) / (sum_w_kk1 + N0_kk1);

rho_k = zeros(size(rho_kk1));
psi0_k = zeros(size(rho_kk1));
psi1_k = zeros(size(rho_kk1));
for n = 1:length(rho_kk1)
    ndd = n - 1;
    psi0_k(n) = psi_ddot(0, phi_kk1, Jz, ndd);
    psi1_k(n) = psi_ddot(1, phi_kk1, Jz, ndd);
end

% Update cardinality distribution
% I'm pretty sure the denominator of this equation in the paper is
% just a normalizing term
denom = sum(rho_kk1 .* psi0_k);
for n = 1:length(rho_kk1)
    ndd = n - 1;
    if ndd < Jz
        rho_k(n) = 0;
    else
        rho_k(n) = (rho_kk1(n) * psi0_k(n)) / denom;
    end
end
if sum(rho_k) ~= 0
    rho_k = rho_k ./ sum(rho_k);
end

% Update intensity RFS
sum_vD = GMRFS();
sum_Nk = 0;
for jz = 1:Jz
    % Apply each measurement to each element of the RFS, basically
    % generating a new gaussian mixture for each measurement
    z = Z(:, jz);

    % Clutter likelihood
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

    % Update each component of the RFS based on the current measurement
    m_kkz = zeros(size(m_kk1));
    wkz = zeros(Jkk1, 1);
    for j = 1:Jkk1
        Kj = K(:, :, j);
        m_kk1j = m_kk1(:, j);

        % Kalman update mean based on this measurement
        m_kkz(:, j) = m_kk1j + Kj * (z - H*m_kk1j);
       
        % likelihood of z given current RFS component
        P_kk1j = P_kk1(:, :, j);
        qz = mvnpdf(z, H*m_kk1j, R + H*P_kk1j*H');

        % Updated weights
        wkz(j) = (pd1 * v_kk1.w(j) * qz) / weight_denom;
        if wkz(j) >= 1 || wkz(j) < 0
            warning('Suspicious weight')
        end

    end
    sum_vD = sum_vD + GMRFS(m_kkz, P_kk, wkz);
end

% This ratio of inner products gets reused:
if sum(psi0_k .* rho_kk1) == 0
    psi_rho_ratio = 0;
else
    psi_rho_ratio = sum(psi1_k .* rho_kk1) / sum(psi0_k .* rho_kk1);
end
psi_rho_weights = psi_rho_ratio / (sum_w_kk1 + N0_kk1);

% qd = 1-pd
v_k_unpruned = (1 - pd1) * psi_rho_weights .*  v_kk1 + sum_vD;
N0_k = N0_kk1 * ((1 - pd0) * psi_rho_weights + sum_Nk); 

%% Prune
v_k = prune_gmphd(v_k_unpruned, T, U, Jmax); 
 
%% Outputs
state_k = cphd_state();
state_k.N0 = N0_k;
state_k.N1 = sum(v_k.w);
state_k.v = v_k;
state_k.rho = rho_k;

% Estimated clutter rate
lambda = N0_k * pd0;

%% Extract a state estimate from the RFS
% TODO: Take the round(N1) highest components? How then to handle
% components with a rounded weight > 1, e.g. if w = 2 then that
% theoretically means that component represents 2 tracked objects (ie two
% objects at the same location)
Xhat = [];
for i = 1:v_k.J
    if v_k.w(i) > w_min
        Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
    end
end

end


%% Helper Functions

function P = permnj(n, j)
% P = permnj(n, j) Compute a permutation coefficient
% P = n! / (n -j)!

% "naive" implementation
%P = nchoosek(n, j) * factorial(j);

% Implementation from vo
P = exp(sum(log(1:n))-sum(log(1:n-j)));
end

function out = psi_ddot(u, phi, Jz, ndd)
% PSI double dot function. I think that |Zk| in the paper
% is referring the the cardinality of the observation RFS
if ~isscalar(u) || ~isscalar(phi) || ~isscalar(Jz) || ~isscalar(ndd)
    error('All inputs must be scalars')
end

Zk_u = Jz + u;

if ndd < Jz + u
    out = 0;
    return;
end

% Calculate using approximation of permutation coefficient and
% log likelihoods to avoid numerical issues
out = exp(sum(log(1:ndd))-sum(log(1:ndd - Zk_u))+(ndd - Zk_u)*log(phi));

end
