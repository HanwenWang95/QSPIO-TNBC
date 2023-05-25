% Model Settings
model_name = 'data_model';
start_time = 0.0; % [days]
time_step = .1; % [days] 0.01 days ~ 15 mins
end_time = 365; % [days]
tol_abs = 1e-12;
tol_rel = 1e-6;
solver = 'ode15s';
% Model Object
time = start_time:time_step:end_time;

% Model Object
model = sbiomodel(model_name);
config = getconfigset(model);
options = get(config,'CompileOptions');
set(options,'UnitConversion',true);
set(config,'TimeUnits','day');
set(config,'SolverType',solver);
set(config.SolverOptions,'OutputTimes',time);
% set(config.SolverOptions,'MaxStep',0.01);
set(config,'StopTime',time(end));
set(config.SolverOptions,'AbsoluteTolerance',tol_abs);
set(config.SolverOptions,'RelativeTolerance',tol_rel);

% Setup Compartments
comp_1 = addcompartment(model,'V_C',3.51,'CapacityUnits','liter');
comp_2 = addcompartment(model,'V_P',3.45,'CapacityUnits','liter');

% Add Species
% Central
A_C = addspecies(comp_1,'A',0.0,'InitialAmountUnits','microgram/milliliter');
% Peripheral
A_P = addspecies(comp_2,'A',0.0,'InitialAmountUnits','microgram/milliliter');

% Add Parameters
MW = addparameter(model,'MW',1.463e8,'ValueUnits','milligram/mole'); % durva
    set(MW,'Notes',['Drug molecular weight']);

q_P = addparameter(model,'q_P',0.476,'ValueUnits','liter/day');
    set(q_P,'Notes',['Volumetric flow rate of the drug between central and peripheral compartment']);

k_cl = addparameter(model,'k_cl',0.232,'ValueUnits','liter/day');
    set(k_cl,'Notes',['Clearance rate of the drug from central compartment']);
k_cln = addparameter(model,'k_cln',0.824,'ValueUnits','milligram/day');
    set(k_cln,'Notes',['Non-linear clearance rate of the drug from central compartment']);
Kcl = addparameter(model,'Kcl',0.344,'ValueUnits','milligram/liter');
    set(Kcl,'Notes',['M-M constant for non-linear clearance']);

% Add Reactions
% Diffusive Transport: Central to Peripheral
reaction = addreaction(model,'V_C.A <-> V_P.A'); % q_P = perm_CP_Ab*P_ratio_Nivo*S_CP*Peri, A = Nivo/f_vol
    set(reaction,'ReactionRate','q_P*(V_C.A/V_C - V_P.A/V_P)');
    set(reaction,'Notes',['Drug diffusive transport to peripheral compartment']);

% Clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cl*V_C.A/V_C');
    set(reaction,'Notes','Drug clearance from central compartment');

% Non-linear clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cln*V_C.A/(V_C.A+Kcl)');
    set(reaction,'Notes','Drug clearance from central compartment');

data_model = model;
