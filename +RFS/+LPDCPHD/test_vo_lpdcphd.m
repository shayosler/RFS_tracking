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
vo_est=   robust.jointcphd.gms.run_filter(vo_model,vo_meas);
vo_handles= robust.jointcphd.gms.plot_results(vo_model,vo_truth,vo_meas,vo_est);
