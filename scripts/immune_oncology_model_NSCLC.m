%% Immune Oncology Model
% Script for setting up and running the immune oncology model in simbiology
clear
close all
sbioreset

%% Add use-defined units into both Simbiology and Symbolic libraries
% Add 'cell' unit to SimBiology and Symbolic Toolboxes
if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

% Add 'mU' unit to SimBiology and Symbolic Toolboxes
mU_unit = sbiounit('mU','mole/liter');
sbioaddtolibrary(mU_unit);
% Symbolic Unit
u = symunit;
try u.mU;
catch
    newUnit('mU',u.molarity);
end

%% Setup Parameters
% Setup Parameters
params_in     = parameters_NSCLC;
params_out    = load_parameters(params_in);

%% Create the SimBiology Model
% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % [days]
time_step = 1; % [days] 0.01 days ~ 15 mins
end_time = 400; % [days]
absolute_tolerance = 1e-9;
relative_tolerance = 1e-6;
% solver = 'ode15s';
solver = 'sundials';
% Model Object
time = start_time:time_step:end_time;
model = simbio_init(model_name,time,solver,absolute_tolerance,relative_tolerance,params_out);
% Maximal simulation time
config = getconfigset(model);
set(config, 'MaximumWallClock', 120)
% set(config, 'MaximumNumberOfLogs', 1)
set(config.SolverOptions, 'AbsoluteToleranceScaling', false)

%% Add Modules to the Model
% Cancer Modules
model = cancer_module(model,'C1',params_out,1); % 4.7e6
% model = cancer_module(model,'C2',params_out,0); % chemotherapy-resistant clone
% T cell Modules
model = Treg_module(model,params_out);
model = Teff_module(model,'1',params_out,{'C1'});
% APC Module
model = APC_module(model,params_out);
% Antigen Modules
antigenCP = create_antigen({'C1'},5.4e-13,'antigenID',0);
model = antigen_module(model,'0',params_out,antigenCP);
antigen   = create_antigen({'C1'},5.4e-13,'antigenID',1);
model = antigen_module(model,'1',params_out,antigen);
% PK Modules
params_aPD1    = pk_parameters('pembrolizumab');
params_aPDL1   = pk_parameters('durvalumab');
params_aCTLA4  = pk_parameters('tremelimumab');
model = pk_module(model,'aPD1',params_aPD1);
model = pk_module(model,'aPDL1',params_aPDL1,'n');
model = pk_module(model,'aCTLA4' ,params_aCTLA4);
% Checkpoint Modules
model = checkpoint_module(model,params_out,'T','C1');
model = checkpoint_module(model,params_out,'T','APC');
% ADCC Module (use in ipilimumab therapy; in development)
% model = Treg_ADCC_module(model,params_out);

% QSPIO-TNBC Modules
model = Th_module(model,params_out);
model = MDSC_module(model,params_out,{'C1'},'inostat',0,'drugName','entinostat');
% model = nabpaclitaxel_module(model,params_out);
model = macrophage_module(model,params_out,{'C1'},'aCD47',0); % PK module in development for aCD47

%% Setup Dosing
dose_schedule = [];
% dose_schedule = schedule_dosing({'atezolizumab'});
% dose_schedule = schedule_dosing({'nabPaclitaxel'});

% dbstop if warning
% simData = sbiosimulate(model,[],[],dose_schedule);
% figure
% % semilogy(0:4000,simData.Data(:,66))
% hold on
% plot(0:4000,(simData.Data(:,204).*6./pi).^(1/3))
% plot(0:4000,simData.Data(:,22))
% plot([find(simData.Data(:,198)>=77,1) find(simData.Data(:,198)>=77,1)], [1e-10 1e10])
% ylim([1 1e5])
% legend('PDL1','Vt')

%% Initialize and Run the Model (should run with realistic baseline parameters)
% (should be commented out when conducting in silico virtual clinical trial)
% tic
% [model,success,simDataInit] = initial_conditions(model);
% % [model,success] = initial_conditions(model);
% toc
%
% % Generate a list of parameters and species for debug
% modelComp = listModelComp(model);
%
% % Run Simulation
% if (success)
%     tic
%     simData = sbiosimulate(model,[],[],dose_schedule);
%     toc
% else
%     simData = simDataInit;
%     disp('Tumour did not reach specified initial tumour diameter with current parameters');
% end
%
% VDT = 74.*log(2)./(log(simData.Data(75,204))-log(simData.Data(1,204)))

%% Plots
% Plot diagnostics
% if (success)
%    diagnostic_plot(simData,model);
%    diagnostic_plot_H(simData);
%    diagnostic_plot_KPR(simData,model);
% end
