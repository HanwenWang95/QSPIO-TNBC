%% In silico Virtual Clinical Trial (Paramster Sensitivity Analysis PSA)
% Script for setting up and running in silico clinical trial
clear
close all
sbioreset

%% Create the model
immune_oncology_model_TNBC

%% Define dosing

% Default dose regimen for atezolizumab is 1200 mg Q3W
% Default dose regimen for nab-paclitaxel is 100 mg/m2 Q3/4W

% dose_schedule = [];
% dose_schedule = schedule_dosing({'atezolizumab'});
% dose_schedule = schedule_dosing({'nabPaclitaxel'});
dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_schedule',[0,28,15]);

% Sequential Therapy
% dose_schedule = schedule_dosing({'atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30]);
% dose_schedule = schedule_dosing({'atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30]);
% dose_schedule = schedule_dosing({'atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30]);

% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_schedule',[0,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_schedule',[0,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_schedule',[14,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_schedule',[14,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_schedule',[14,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_schedule',[14,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_schedule',[28,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_schedule',[28,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_schedule',[28,28,15]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_schedule',[28,28,15]);

% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',125, 'nabPaclitaxel_schedule',[28,21,20]);

% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[0,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[14,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel'}, 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[0,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[14,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[28,21,20]);
% dose_schedule = schedule_dosing({'nabPaclitaxel','atezolizumab'}, 'atezolizumab_dose',840/80, 'atezolizumab_schedule',[28,14,30], 'nabPaclitaxel_dose',260, 'nabPaclitaxel_schedule',[28,21,20]);


%% Define and Prepare the Input and Output Parameters
% Load Previously Generated Virtual Patients if applicable
% load('VP.mat')

% Generate New Parameter Sets (New Virtual Patient Cohort)
n_PSA = 1300;

% Set distributions of the selected parameters for random sampling
params_in  = PSA_param_in_TNBC;
% Set boundaries for model species (only include those added in the present model)
params_out = PSA_param_out(model);
% Add selected model outputs to sensitivity analysis
params_in  = PSA_param_obs(params_in);
% Randomly generate parameter sets using Latin-Hypercube Sampling
params_in  = PSA_setup(model,params_in,n_PSA);

% Save Virtual Patient Cohort
save('VP.mat', 'params_in', 'params_out')

%% Run Batch Simulations
warning('off','all')
% dbstop if warning
% dbclear all

sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

%% Postprocess
% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

% Prepare the data for the rest of the analysis
params_out = PSA_prep(simDataPSA,simDataPSApost,params_out);

% Save and print data of interest (by assigning a unique code name for the trial)
sprint_data(simDataPSA, simDataPSApost, params_in, params_out, 'combo')

%% Perform and plot different types of analysis

% Partial Rank Correlation Coefficients
% PSA_PRCC(params_in,params_out,'plausible')
% PSA_PRCC(params_in,params_out)

% t-SNE Analysis
% PSA_tSNE(params_in,params_out,'plausible')
% PSA_tSNE(params_in,params_out,'patient')

% eFAST

% Principle Component Analysis

%% Plot Results
% close all
% Plot percent change in size and RECIST
% PSA_plot_RECIST(simDataPSA,simDataPSApost,params_out)

% Tumor Size
% PSA_plot_TumSize(simDataPSA,simDataPSApost,params_out)

% Kaplan-Meier Progression-free survival (PFS)
% PSA_plot_KaplanMeier(simDataPSA,simDataPSApost,params_out)

%% Waterfall plots
% Plot percent change in size using waterfall plots for a parameter.
% Waterfall plots can be plotted using either end tumor sizes or best overall responses.
%
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'n_T1_clones')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_P1_d1')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_C_resist')

% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.2 .4 .7]) % blue
% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.2 .6 .4]) % green
% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.8 .6 .2]) % yellow/orange
