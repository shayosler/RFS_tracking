%% clear
clc
clear
close all

%% Set up
seed = 1;
rng(seed);

%% Run VO model
vo_model= robust.jointcphd.gms.gen_model;
vo_truth= robust.jointcphd.gms.gen_truth(vo_model);
vo_meas=  robust.jointcphd.gms.gen_meas(vo_model,vo_truth);
%vo_est=   robust.jointcphd.gms.run_filter(vo_model,vo_meas);
%vo_handles= robust.jointcphd.gms.plot_results(vo_model,vo_truth,vo_meas,vo_est);

%% Filtering parameters (from paper)
params = RFS.LPDCPHD.lpdcphd_params();
params.Nmax = 100;
params.U = 4;
params.T = 1e-5;
params.Jmax = 500;
params.w_min = 0.5;

%% Define environment
% Bounds for northing and easting
min_n = vo_model.range_c(1, 1);
max_n = vo_model.range_c(1, 2);
min_e = vo_model.range_c(2, 1);
max_e = vo_model.range_c(2, 2);
northings = min_n:1:max_n;
eastings = min_e-10:1:max_e;

%% Set up model for my implementation
n_states = 4;
A_fov = (vo_model.range_c(1, 2) - vo_model.range_c(1, 1)) * (vo_model.range_c(2, 2) - vo_model.range_c(2, 1));
model = RFS.LPDCPHD.lpdcphd_model();

% Target dynamics
model.F = vo_model.F;
model.Q = vo_model.Q;

model.kB = 1.1; % TODO: find actual value

% Detection and survival probabilities
model.ps0 = 0.9; 
model.ps1 = .99;

% Birth models
model.gamma0 = RFS.utils.BMRFS(1, vo_model.u_cb, vo_model.v_cb); % TODO: confirm weight
model.gamma1 = RFS.utils.BGMRFS(vo_model.w_birth, vo_model.m_birth, vo_model.P_birth, vo_model.u_b, vo_model.v_b);

model.kappa = 1 / A_fov; % Clutter is equally likely anywhere

% Sensor/measurement model
R = vo_model.R;
H = vo_model.H;

% clutter: poisson 
Ngamma0 = vo_model.clutter_Nt;
lambda_true = vo_model.clutter_Nt * vomodel.clutter_P_D;

%% Simulate
sim_dt = vo_model.T;
t_total = 100;
t = 0:sim_dt:t_total;
sim_steps = numel(t);

% Initial conditions and storage for l-pd-CPHD
lpdcphd_Xhat = cell(sim_steps, 1);
lpdcphd_states(sim_steps, 1) = RFS.LPDCPHD.lpdcphd_state();

PD_init = 0.5;
Nc_init = round((size(meas.Z{1},2)-PD_init*sum(vo_model.w_birth))/PD_init);
lpdcphd_states(1).v0 = RFS.utils.BMRFS(Nc_init, 1, 1);
lpdcphd_states(1).v1 = RFS.utils.BGMRFS(eps, [0.1;0;0.1;0], diag([1 1 1 1]).^2, 1, 1);
lpdcphd_states(1).rho = [zeros(Nc_init-1,1);1;zeros(filter.N_max-Nc_init+1,1)];

lpdcphd_states(1).N0 = 1;
lpdcphd_states(1).N1 = 1;
lpdcphd_states(1).lambda = 0;
lpdcphd_states(1).pd0 = 0;
lpdcphd_states(1).pd1 = 0;

for k = 2:sim_steps
    tk = t(j);
    j = j+1;

    % Simulated observations
    measurement = RFS.LPDCPHD.lpdcphd_measurement();
    measurement.Z = vo_meas.Z{k-1};
    measurement.H = H;
    measurement.R = R;

    [lpdcphd_states(k), lpdcphd_Xhat{k}] = RFS.LPDCPHD.lpdcphd_filter(lpdcphd_states(k-1), measurement, model, params);

end






%% Run my implementation