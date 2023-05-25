% Model Settings
model_name = 'model';
start_time = 0; % [days]
time_step = 1; % [days] 0.01 days ~ 15 mins
end_time = 21; % [days]
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
comp_C = addcompartment(model,'V_C',5,'CapacityUnits','liter');
comp_P = addcompartment(model,'V_P',60,'CapacityUnits','liter');
comp_T = addcompartment(model,'V_T',20,'CapacityUnits','milliliter');
comp_LN = addcompartment(model,'V_LN',1,'CapacityUnits','milliliter');

% Add Species
% Central
A_C = addspecies(comp_C,'A',0.0,'InitialAmountUnits','molarity');
% Peripheral
A_P = addspecies(comp_P,'A',0.0,'InitialAmountUnits','molarity');
% Tumour
A_T = addspecies(comp_T,'A',0.0,'InitialAmountUnits','molarity');
% Lymph Node
A_LN = addspecies(comp_LN,'A',0.0,'InitialAmountUnits','molarity');

% Add Parameters
q_P = addparameter(model,'q_P',0,'ValueUnits','liter/second');
    set(q_P,'Notes',['Volumetric flow rate of the drug between central and peripheral compartment']);
q_T = addparameter(model,'q_T',8.52e-5,'ValueUnits','milliliter/second');
    set(q_T,'Notes',['Volumetric flow rate of the drug between central and tumor compartment']);
q_LN = addparameter(model,'q_LN',3.25e-6,'ValueUnits','milliliter/second');
    set(q_LN,'Notes',['Volumetric flow rate of the drug between central and TDLN compartment']);
q_LD = addparameter(model,'q_LD',0.0015,'ValueUnits','1/minute');
    set(q_LD,'Notes',['Rate of lymphatic drainage of the drug from TDLN to central compartment']);
k_cl = addparameter(model,'k_cl',0,'ValueUnits','liter/day');
    set(k_cl,'Notes',['Clearance rate of the drug from central compartment']);
k_cln = addparameter(model,'k_cln',0,'ValueUnits','nanomole/day');
    set(k_cln,'Notes',['Non-linear clearance rate of the drug from central compartment']);
Kcl = addparameter(model,'Kcl',0,'ValueUnits','nanomole/liter');
    set(Kcl,'Notes',['M-M constant for non-linear clearance']);

gamma_C = addparameter(model,'gamma_C',.55,'ValueUnits','dimensionless');
    set(gamma_C,'Notes',['Volume fraction of interstitial space available to the drug in central compartment']);
gamma_P = addparameter(model,'gamma_P',.06,'ValueUnits','dimensionless');
    set(gamma_P,'Notes',['Volume fraction of interstitial space available to the drug in peripheral compartment']);
gamma_T = addparameter(model,'gamma_T',0.522,'ValueUnits','dimensionless');
    set(gamma_T,'Notes',['Volume fraction of interstitial space available to the drug in tumor compartment']);
gamma_LN = addparameter(model,'gamma_LN',0.2,'ValueUnits','dimensionless');
    set(gamma_LN,'Notes',['Volume fraction of interstitial space available to the drug in TDLN compartment']);

% Add Reactions
% Diffusive Transport: Central to Peripheral
reaction = addreaction(model,'V_C.A <-> V_P.A'); % q_P = perm_CP_Ab*P_ratio_Nivo*S_CP*Peri, A = Nivo/f_vol
    set(reaction,'ReactionRate','q_P*(V_C.A/gamma_C - V_P.A/gamma_P)');
    set(reaction,'Notes',['Drug diffusive transport to peripheral compartment']);
% Diffusive Transport: Central to Tumour
reaction = addreaction(model,'V_C.A <-> V_T.A'); % q_T = perm_CT_Ab*S_CT*Tum
    set(reaction,'ReactionRate','q_T*(V_C.A/gamma_C - V_T.A/gamma_T)');
    set(reaction,'Notes',['Drug diffusive transport to tumor compartment']);
% Diffusive Transport: Central to Lymph Node
reaction = addreaction(model,'V_C.A <-> V_LN.A');
    set(reaction,'ReactionRate','q_LN*(V_C.A/gamma_C - V_LN.A/gamma_LN)');
    set(reaction,'Notes',['Drug diffusive transport to lymph node compartment']);
% Convective Transport: Tumour to Lymph Node
reaction = addreaction(model,'V_T.A -> V_LN.A');
    set(reaction,'ReactionRate','q_LD*V_T*V_T.A/gamma_T');
    set(reaction,'Notes',['Drug convective transport from tumor to lymph node']);
% Convective Transport: Lymph Node to Central
reaction = addreaction(model,'V_LN.A -> V_C.A');
    set(reaction,'ReactionRate','q_LD*V_T*V_LN.A/gamma_LN');
    set(reaction,'Notes',['Drug convective transport from lymph node to central']);
% Clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cl*V_C.A');
    set(reaction,'Notes','Drug clearance from central compartment');

reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cln*V_C.A/(V_C.A+Kcl)');
    set(reaction,'Notes','Drug clearance from central compartment');

fit_model = model;
