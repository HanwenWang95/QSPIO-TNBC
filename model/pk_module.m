% Pharmacokinetic Module
%
% Models pharmacokinetics in four compartments
%
% Inputs: model        -- SimBiology model object with four compartments: C,P,T,LN
%         species_name -- name of drug [must be unique]
%         params       -- object containing model parameter
%                         - q_P--rate of diffusive transport C<->P
%                         - q_T--rate of diffusive transport C<->T
%                         - q_LN--rate of diffusive transport C<->LN
%                         - q_LD--rate of convective transport T->LN
%                         - k_cl--clearence rate from central
%                         - gamma_C--volume fraction in C
%                         - gamma_P--volume fraction in P
%                         - gamma_T--volume fraction in T
%                         - gamma_LN--volume fraction in LN
%                         - k_a1--buccal absorption rate
%                         - k_a2--GI absorption rate
%                         - k_cln--Vmax of non-linear clearance from central compartment
%                         - Kc--drug concentration when non-linear clearance is .5*Vmax
%
% Outputs: model -- SimBiology model object with new pk module


function model = pk_module(model,species_name,params,dose_type)

% Get Compartments
comp_C = model.Compartment(1);
comp_P = model.Compartment(2);
comp_T = model.Compartment(3);
comp_LN = model.Compartment(4);

% Add Species
% Central
A_C = addspecies(comp_C,'A',0.0,'InitialAmountUnits','nanomolarity');
    set(A_C,'Notes',['Concentration of ' species_name ' in central compartment']);
% Peripheral
A_P = addspecies(comp_P,'A',0.0,'InitialAmountUnits','nanomolarity');
    set(A_P,'Notes',['Concentration of ' species_name ' in peripheral compartment']);
% Tumour
A_T = addspecies(comp_T,'A',0.0,'InitialAmountUnits','nanomolarity');
    set(A_T,'Notes',['Concentration of ' species_name ' in tumour compartment']);
% Lymph Node
A_LN = addspecies(comp_LN,'A',0.0,'InitialAmountUnits','nanomolarity');
    set(A_LN,'Notes',['Concentration of ' species_name ' in lymph node compartment']);
% Oral Administration
if (nargin==4 && dose_type == 'o')

    A_Buccal = addspecies(comp_C,'A_Buccal',0.0,'InitialAmountUnits','nanomolarity');
        set(A_Buccal,'Notes',['Concentration of ' species_name ' from buccal absorption']);
    A_GI = addspecies(comp_C,'A_GI',0.0,'InitialAmountUnits','nanomolarity');
        set(A_GI,'Notes',['Concentration of ' species_name ' from GI absorption']);
    k_a1 = addparameter(model,'k_a1',params.k_a1.Value,'ValueUnits',params.k_a1.Units);
        set(k_a1,'Notes',['Absorption rate of ' species_name ' from buccal to central compartment' params.k_a1.Notes]);
    k_a2 = addparameter(model,'k_a2',params.k_a2.Value,'ValueUnits',params.k_a2.Units);
        set(k_a2,'Notes',['Absorption rate of ' species_name ' from GI to central compartment' params.k_a2.Notes]);
    k_cln = addparameter(model,'k_cln',params.k_cln.Value,'ValueUnits',params.k_cln.Units);
        set(k_cln,'Notes',['Non-linear clearance rate of ' species_name ' from central compartment' params.k_cln.Notes]);
    Kc = addparameter(model,'Kc',params.Kc.Value,'ValueUnits',params.Kc.Units);
        set(Kc,'Notes',['Half-maximal concentration of ' species_name ' in central compartment for nonlinear clearance' params.Kc.Notes]);
    % add parameters for dose schedule
    lagP = addparameter(model,'lagP','Value',params.lagP.Value,'ValueUnits',params.lagP.Units);
        set(lagP,'Notes',['Lag time of ' species_name ' absorption into central compartment' params.lagP.Notes]);
    durP = addparameter(model,'durP','Value',params.durP.Value,'ValueUnits',params.durP.Units);
        set(durP,'Notes',['Duration of zero-order absorption of ' species_name ' into central compartment' params.durP.Notes]);

    Dose_GI = addspecies(model.Compartment(1),'Dose_GI',0,'InitialAmountUnits','nanomolarity');
        set(Dose_GI,'Notes',['Fraction of ' species_name ' dose for GI absorption']);
    k_GI = addparameter(model,'k_GI',params.k_GI.Value,'ValueUnits',params.k_GI.Units);
        set(k_GI,'Notes',['Absorption rate of ' species_name ' to gastro-intestinal compartment' params.k_GI.Notes]);

    reaction = addreaction(model,'V_C.Dose_GI -> V_C.A_GI');
        set(reaction,'ReactionRate','k_GI*V_C.Dose_GI');
        set(reaction,'Notes',['First-order absorption of ' species_name ' to gastro-intestinal compartment']);
    reaction = addreaction(model,'V_C.A_Buccal -> V_C.A');
        set(reaction,'ReactionRate','k_a1*V_C.A_Buccal');
        set(reaction,'Notes',['Absorption of ' species_name ' from buccal to central comparment']);
    reaction = addreaction(model,'V_C.A_GI -> V_C.A');
        set(reaction,'ReactionRate','k_a2*V_C.A_GI');
        set(reaction,'Notes',['Absorption of ' species_name ' from gastro-intestinal to central comparment']);
    reaction = addreaction(model,'V_C.A -> null');
        set(reaction,'ReactionRate','k_cln*V_C.A/(V_C.A + Kc)');
        set(reaction,'Notes',['Nonlinear clearance of ' species_name ' from central compartment']);

    rename(A_Buccal,[species_name '_Buccal']);
    rename(A_GI,[species_name '_GI']);
    rename(k_a1,['k_a1_' species_name]);
    rename(k_a2,['k_a2_' species_name]);
    rename(k_cln,['k_cln_' species_name]);
    rename(Kc,['Kc_' species_name]);
    rename(lagP,['lagP_' species_name]);
    rename(durP,['durP_' species_name]);
    rename(k_GI,['k_GI_' species_name]);
end

% Add Parameters
q_P = addparameter(model,'q_P',params.q_P.Value,'ValueUnits',params.q_P.Units);
    set(q_P,'Notes',['Capillary filtration rate of ' species_name ' between central and peripheral compartment' params.q_P.Notes]);
q_T = addparameter(model,'q_T',params.q_T.Value,'ValueUnits',params.q_T.Units);
    set(q_T,'Notes',['Capillary filtration rate of ' species_name ' between central and tumor compartment' params.q_T.Notes]);
q_LN = addparameter(model,'q_LN',params.q_LN.Value,'ValueUnits',params.q_LN.Units);
    set(q_LN,'Notes',['Capillary filtration rate of ' species_name ' between central and TDLN compartment' params.q_LN.Notes]);
q_LD = addparameter(model,'q_LD',params.q_LD.Value,'ValueUnits',params.q_LD.Units);
    set(q_LD,'Notes',['Rate of lymphatic drainage of ' species_name ' from TDLN to central compartment' params.q_LD.Notes]);
k_cl = addparameter(model,'k_cl',params.k_cl.Value,'ValueUnits',params.k_cl.Units);
    set(k_cl,'Notes',['Clearance rate of ' species_name ' from central compartment' params.k_cl.Notes]);
gamma_C = addparameter(model,'gamma_C',params.gamma_C.Value,'ValueUnits',params.gamma_C.Units);
    set(gamma_C,'Notes',['Volume fraction of interstitial space available to ' species_name ' in central compartment' params.gamma_C.Notes]);
gamma_P = addparameter(model,'gamma_P',params.gamma_P.Value,'ValueUnits',params.gamma_P.Units);
    set(gamma_P,'Notes',['Volume fraction of interstitial space available to ' species_name ' in peripheral compartment' params.gamma_P.Notes]);
gamma_T = addparameter(model,'gamma_T',params.gamma_T.Value,'ValueUnits',params.gamma_T.Units);
    set(gamma_T,'Notes',['Volume fraction of interstitial space available to ' species_name ' in tumor compartment' params.gamma_T.Notes]);
gamma_LN = addparameter(model,'gamma_LN',params.gamma_LN.Value,'ValueUnits',params.gamma_LN.Units);
    set(gamma_LN,'Notes',['Volume fraction of interstitial space available to ' species_name ' in TDLN compartment' params.gamma_LN.Notes]);

% Add Reactions
% Diffusive Transport: Central to Peripheral
reaction = addreaction(model,'V_C.A <-> V_P.A'); % q_P = perm_CP_Ab*P_ratio_Nivo*S_CP*Peri, A = Nivo/f_vol
    set(reaction,'ReactionRate','q_P*(V_C.A/gamma_C - V_P.A/gamma_P)');
    set(reaction,'Notes',[species_name ' diffusive transport to peripheral compartment']);
% Diffusive Transport: Central to Tumour
reaction = addreaction(model,'V_C.A <-> V_T.A'); % q_T = perm_CT_Ab*S_CT*Tum
    set(reaction,'ReactionRate','q_T*(V_C.A/gamma_C - V_T.A/gamma_T)');
    set(reaction,'Notes',[species_name ' diffusive transport to tumor compartment']);
% Diffusive Transport: Central to Lymph Node
reaction = addreaction(model,'V_C.A <-> V_LN.A');
    set(reaction,'ReactionRate','q_LN*(V_C.A/gamma_C - V_LN.A/gamma_LN)');
    set(reaction,'Notes',[species_name ' diffusive transport to lymph node compartment']);
% Convective Transport: Tumour to Lymph Node
reaction = addreaction(model,'V_T.A -> V_LN.A');
    set(reaction,'ReactionRate','q_LD*V_T*V_T.A/gamma_T');
    set(reaction,'Notes',[species_name ' convective transport from tumor to lymph node']);
% Convective Transport: Lymph Node to Central
reaction = addreaction(model,'V_LN.A -> V_C.A');
    set(reaction,'ReactionRate','q_LD*V_T*V_LN.A/gamma_LN');
    set(reaction,'Notes',[species_name ' convective transport from lymph node to central']);
% Clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cl*V_C.A');
    set(reaction,'Notes','Drug clearance from central compartment');

if (nargin==4 && dose_type == 'n')

k_cln = addparameter(model,'k_cln',params.k_cln.Value,'ValueUnits',params.k_cln.Units);
    set(k_cln,'Notes',['Non-linear clearance rate of ' species_name ' from central compartment' params.k_cln.Notes]);
Kc = addparameter(model,'Kc',params.Kc.Value,'ValueUnits',params.Kc.Units);
    set(Kc,'Notes',['Half-maximal concentration of ' species_name ' in central compartment for nonlinear clearance' params.Kc.Notes]);

% Non-Linear Clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cln*V_C.A/(V_C.A + Kc)');
    set(reaction,'Notes',['Nonlinear clearance of ' species_name ' from central compartment']);

rename(k_cln,['k_cln_' species_name]);
rename(Kc,['Kc_' species_name]);

end

% Rename Objects with 'species_name'
rename(A_C,species_name);
rename(A_P,species_name);
rename(A_T,species_name);
rename(A_LN,species_name);

rename(q_P,['q_P_' species_name]);
rename(q_T,['q_T_' species_name]);
rename(q_LN,['q_LN_' species_name]);
rename(q_LD,['q_LD_' species_name]);
rename(k_cl,['k_cl_' species_name]);
rename(gamma_C,['gamma_C_' species_name]);
rename(gamma_P,['gamma_P_' species_name]);
rename(gamma_T,['gamma_T_' species_name]);
rename(gamma_LN,['gamma_LN_' species_name]);
