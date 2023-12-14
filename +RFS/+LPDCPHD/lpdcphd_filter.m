function [state_k, Xhat] = lpdcphd_filter(state, measurement, model, params)
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
    P_pred = Q + F*P_k(:, :, j)*F';
    P_kk1(:, :, j) = RFS.utils.make_symmetric(P_pred);
end

% Predict beta part of distribution
mu_beta = v1.s ./ (v1.s + v1.t);
sigsq_beta = k_beta * (v1.s .* v1.t) ./ ( (v1.s + v1.t).^2 .* (v1.s + v1.t + 1));
mm = ((mu_beta .* (1 - mu_beta)) ./ sigsq_beta) - 1;
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

% surviving cardinality distribution
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
    P_kkj = (I - Kj * H) * P_kk1j * (I - Kj * H)' + Kj * R * Kj';
    P_kk(:, :, j) = RFS.utils.make_symmetric(P_kkj);
end

% Detection probabilities for prediction RFSs
d0_kk1 = v0_kk1.s ./ (v0_kk1.s + v0_kk1.t);
d1_kk1 = v1_kk1.s ./ (v1_kk1.s + v1_kk1.t);

% Phi_kk1 is 1 - ratio of expected number of detected objects (targets + clutter)
% to estimated  total number of objects. ie its the expected fraction of 
% total objects that are "missed" by the observation
sum_w0_kk1 = sum(v0_kk1.w);
sum_w1_kk1 = sum(v1_kk1.w);
phi_kk1 = 1 - (v1_kk1.w' * d1_kk1 + v0_kk1.w' * d0_kk1) / (sum_w1_kk1 + sum_w0_kk1);

% Calculate psi0 and psi1. Not sure exactly what they represent
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
rho_k = zeros(size(rho_kk1));
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

% Calculate updated weights for the prediction component
wM_k_denom = sum_w0_kk1 + sum_w1_kk1;
psi_ratio = (psi1_k' * rho_kk1) / (psi0_k' * rho_kk1);
wM0_k = v0_kk1.w .* (beta(v0_kk1.s, v0_kk1.t + 1)./beta(v0_kk1.s, v0_kk1.t)) * psi_ratio / wM_k_denom;
wM1_k = v1_kk1.w .* (beta(v1_kk1.s, v1_kk1.t + 1)./beta(v1_kk1.s, v1_kk1.t)) * psi_ratio / wM_k_denom;

% Update intensity RFS. Update is weighted sum of prediction RFS, and an
% RFS derived from each measurement
v0_z = RFS.utils.BMRFS();
v1_z = RFS.utils.BGMRFS();

for jz = 1:Jz
    % Apply each measurement to each element of the RFS, basically
    % generating a new gaussian mixture for each measurement
    z = Z(:, jz);

    % Update each component of the RFS based on the current measurement
    m_kz = zeros(size(m_kk1));
    qz = zeros(Jkk1, 1);
    for j = 1:Jkk1
        Kj = K(:, :, j);
        m_kk1j = m_kk1(:, j);
        P_kk1j = P_kk1(:, :, j);

        % Kalman update mean based on this measurement
        m_kz(:, j) = m_kk1j + Kj * (z - H*m_kk1j);
       
        % likelihood of z given current RFS component
        qz(j) = mvnpdf(z, H*m_kk1j, R + H*P_kk1j*H');
    end

    % Spatial clutter likelihood
    kappaz = kappa; % TODO: support non-uniform clutter likelihood

    % Denominator used for normalizing weights
    weight_denom = v0_kk1.w' * d0_kk1 * kappaz + sum(v1_kk1.w .* d1_kk1 .* qz);

    % Updated weights
    wD0_k = v0_kk1.w .* (beta(v0_kk1.s + 1, v0_kk1.t) ./ beta(v0_kk1.s, v0_kk1.t)) .* kappaz ./ weight_denom;
    wD1_k = v1_kk1.w .* (beta(v1_kk1.s + 1, v1_kk1.t) ./ beta(v1_kk1.s, v1_kk1.t)) .* qz ./ weight_denom;
    if any(wD0_k >= 1) || any(wD0_k < 0) || any(wD1_k >= 1) || any(wD1_k < 0)
        warning('Suspicious weight')
    end
    v0_z = v0_z + RFS.utils.BMRFS(wD0_k, v0_kk1.s + 1, v0_kk1.t);   % TODO: avoid extra allocations by pre-allocation vectors for all components then filling them during the loop, and then only creating one RFS at the end
    v1_z = v1_z + RFS.utils.BGMRFS(wD1_k, m_kz, P_kk, v1_kk1.s + 1, v1_kk1.t);
end

v0_k_unpruned = RFS.utils.BMRFS(wM0_k, v0_kk1.s, v0_kk1.t + 1) + v0_z;
v1_k_unpruned = RFS.utils.BGMRFS(wM1_k, v1_kk1.m, v1_kk1.P, v1_kk1.s, v1_kk1.t + 1) + v1_z;

%% Prune
% TODO: different pruning parameters for clutter only beta mixture vs 
% beta-gaussian mixture representing targets? There definitely needs to be 
% because prune_bmrfs uses Hellinger distance which is bounded on [0, 1]
% but prune_bgmrfs uses a likelihood-like distance which is unbounded
v0_k = RFS.utils.prune_bmrfs(v0_k_unpruned, T, U, Jmax); 
v1_k = RFS.utils.prune_bgmrfs(v1_k_unpruned, T, U, Jmax);
 
%% Outputs

% Calculate detection probabilities for pruned update RFSs
pd0 = v0_k.s ./ (v0_k.s + v0_k.t);
pd1 = v1_k.s ./ (v1_k.s + v1_k.t);

% Estimated clutter rate
lambda = v0_k.w * d0_kk1';

% Updated state
state_k = RFS.LPDCPHD.lpdcphd_state();
state_k.v0 = v0_k;
state_k.v1 = v1_k;
state_k.rho = rho_k;
state_k.N0 = sum(v0_k.w);
state_k.N1 = sum(v1_k.w);
state_k.lambda = lambda;
state_k.pd0 = pd0;
state_k.pd1 = pd1;


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
for i = 1:v1_k.J
    n_tgts = round(v1_k.w(i));

    % Don't actually remove any?
    if n_tgts < 1
        n_tgts = 1;
    end
    Xhat = [Xhat repmat(v1_k.m(:, i), 1, n_tgts)];
end

% Take up to the rounded cardinality estimate targets
% TODO: are the components ordered by weight at this point?
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
