function [state_k, Xhat, lambda_hat, pd_hat] = lpdcphd_filter(state, measurement, model, params)
% Cardinalized PHD filter for tracking with unknown clutter model and 
% unknown detection profile
% (lambda-pD-CPHD filter)
% Based on: https://ieeexplore.ieee.org/document/5730505
% Inputs:
%   state       lpdcphd_state object containing the previous state of the system
%   measurement lpdcphd_measurement object containing the measurement(s) for 
%               the current time step
%   model       lpdcphd_model object defining the model for the system
%   params      lpdcphd_params object defining the system parameters
%   
% Outputs:
%   state_k     lpdcphd_state object containing the estimated current state of
%               the system
%   Xhat        Estimated target locations
%   lambda_hat  Estimated clutter rate
%   pd_hat      Estimated target detection rate

%if isscalar(pd)
%    pd = @(m) pd;
%end
%if isscalar(ps)
%    ps = @(m) ps;
%end

% Extract state
v0 = state.v0;
v1 = state.v1;
rho = state.rho;
if ~isvector(rho)
    error 'state.rho must by a 1xN or Nx1 vector'
end

% Extract model data:
% 0 suffix -> clutter
% 1 suffix -> targets
% delt suffix -> part of augmented state space associated with detection
% prob
gamma0 = model.gamma0;  % Clutter generator birth RFS
gamma1 = model.gamma1;  % Target birth model, filter assumes birth process is Poisson
ps0 = model.ps0;        % Clutter generator survival probability
ps1 = model.ps1;        % Target survival probability
%pd0 = model.pd0;        % Clutter generator detection probability
%pd1 = model.pd1;        % Target detection probability
Q = model.Q;            % Process noise, NxN
F = model.F;            % State transition matrix x[k+1] = Fx[k], NxN
kappa = model.kappa;    % Spatial likelihood of clutter
k_beta = model.kB;      % "Dilation" constant for beta distribution

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

%% Predict target intensity
% Propagate each gaussian mixture element in the  target RFS forward (survival intensity)
P_k = v1.P;
m_k = v1.m;
m_kk1 = zeros(size(v1.m));
P_kk1 = zeros(size(v1.P));
Jk = v1.J;
for j = 1:Jk
    m_kk1(:, j) = F*m_k(:, j);
    P_kk1(:, :, j) = Q + F*P_k(:, :, j)*F';
end

% Predict beta part of distribution
mu_beta = v1.s ./ (v1.s + v1.t);
sigsq_beta = k_beta * (v1.s .* v1.t) ./ ( (v1.s + v1.t).^2 .* (v1.s + v1.t + 1));
mm = ((mu_beta .* (1 - mu_beta)) / sigsq_beta) - 1;
s_kk1 = mm .* mu_beta;
t_kk1 = mm .* (1 - mu_beta);

% Predicted target RFS:
v1_kk1 = gamma1 + RFS.utils.BGMRFS(ps1.*v1.w, m_kk1, P_kk1, s_kk1, t_kk1);
P_kk1 = v1_kk1.P;
m_kk1 = v1_kk1.m;

%% Predict clutter intensity
v0_kk1 = ps0 .* v0 + gamma0;

%% Predict hybrid cardinality distribution
% Survival factor
if sum(v1.w) + sum(v0.w) == 0
    phi = 0;
else
    phi = (ps1 * sum(v1.w) + ps0 * sum(v0.w)) / (sum(v1.w) + sum(v0.w));
end

% Variables as named by vo:
%w_update = v.w;
%Nc_update = N0;
%vo_model.P_S = ps1;
%vo_model.clutter_P_S = ps0;
%vo_model.w_birth = model.gamma1.w;
%vo_model.lambda_cb = model.Ngamma0;
%cdn_update = rho;

%---cardinality prediction
%surviving cardinality distribution
N_max = length(rho) - 1;
survive_cdn_predict = zeros(N_max+1, 1);

% Prevent NaNs 
if phi == 0
    log_phi = log(eps(0));
else
    log_phi = log(phi);
end

for j=0:N_max
    idxj=j+1;
    terms= zeros(N_max+1,1);
    for ell=j:N_max
        idxl= ell+1;
        terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log_phi+(ell-j)*log(1-phi))*rho(idxl);
    end
    survive_cdn_predict(idxj) = sum(terms);
end

% predicted cardinality = convolution of birth and surviving cardinality distribution
cdn_predict = zeros(N_max+1,1);
lambda_b = sum(model.gamma1.w) + sum(model.gamma0.w); % TODO: check this
for n=0:N_max
    idxn=n+1;
    terms= zeros(N_max+1,1);
    for j=0:n
        idxj= j+1;
        terms(idxj)= exp(-sum(lambda_b)+(n-j)*log(lambda_b)-sum(log(1:n-j)))*survive_cdn_predict(idxj);
    end
    cdn_predict(idxn) = sum(terms);
end
% Normalize predicted cardinality distribution
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
Jkk1 = v1_kk1.J;
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
    sum_w_kk1 = sum_w_kk1 + v1_kk1.w(l);
end
sum_w0_kk1 = sum(v0_kk1.w);
sum_w1_kk1 = sum(v1_kk1.w);

d0 = v0_kk1.s ./ (v0_kk1.s + v0_kk1.t);
d1 = v1_kk1.s ./ (v1_kk1.s + v1_kk1.t);

% TODO: check dimensions of w and d for this. Want dot product
phi_kk1 = 1 - (v1_kk1.w * d1' + v0_kk1.w * d0') / (sum_w1_kk1 + sum_w0_kk1);

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
sum_vD = RFS.utils.GMRFS();
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
        wkl = v1_kk1.w(l);
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
        wkz(j) = (pd1 * v1_kk1.w(j) * qz) / weight_denom;
        if wkz(j) >= 1 || wkz(j) < 0
            warning('Suspicious weight')
        end

    end
    sum_vD = sum_vD + RFS.utils.GMRFS(m_kkz, P_kk, wkz);
end

% This ratio of inner products gets reused:
if sum(psi0_k .* rho_kk1) == 0
    psi_rho_ratio = 0;
else
    psi_rho_ratio = sum(psi1_k .* rho_kk1) / sum(psi0_k .* rho_kk1);
end
psi_rho_weights = psi_rho_ratio / (sum_w_kk1 + N0_kk1);

% qd = 1-pd
v_k_unpruned = (1 - pd1) * psi_rho_weights .*  v1_kk1 + sum_vD;
N0_k = N0_kk1 * ((1 - pd0) * psi_rho_weights + sum_Nk); 

%% Prune
v_k = RFS.utils.prune_gmphd(v_k_unpruned, T, U, Jmax); 
 
%% Outputs
state_k = RFS.CPHD.cphd_state();
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
% for i = 1:v_k.J
%     if v_k.w(i) > w_min
%         Xhat = [Xhat repmat(v_k.m(:, i), 1, round(v_k.w(i)))];
%     end
% end

% Duplicate any components that represent more than 1 target
% Remove any components that don't represent any targets
for i = 1:v_k.J
    n_tgts = round(v_k.w(i));

    % Don't actually remove any?
    if n_tgts < 1
        n_tgts = 1;
    end
    Xhat = [Xhat repmat(v_k.m(:, i), 1, n_tgts)];
end

% Take up to the rounded cardinality estimate targets
if size(Xhat, 2) > round(state_k.N1)
    Xhat = Xhat(:, 1:round(state_k.N1));
end

end


%% Helper Functions

function P = permnj(n, j)
% P = permnj(n, j) Compute a permutation coefficient
% P = n! / (n -j)!

% "naive" implementation, doesn't work with large values of n, j
%P = nchoosek(n, j) * factorial(j);

% Implementation from vo
P = exp(sum(log(1:n))-sum(log(1:n-j)));
end % P()

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

end % psi_ddot()
