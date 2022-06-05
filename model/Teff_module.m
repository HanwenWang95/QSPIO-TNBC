% Teff Module
%
% Models Teff transport and activation by APCs [Use before antigen corresponding module]
%
% Inputs: model        -- SimBiology model object with four compartments
%         ID           -- T cell-antigen ID number [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%         cancer_types -- cell array of strings containing names of cancer types Teff kill
%
% Outputs: model -- SimBiology model object with new Tcell module


function model = Teff_module(model,ID,params,cancer_types)

% Species Names
species_name = ['T' ID];
antigen = ['P' ID];

% Add Species
% Naive T cells
nCD8_C  = addspecies(model.Compartment(1),'nT',params.nCD8_C.Value/params.nCD8_div.Value,'InitialAmountUnits','cell');
    set(nCD8_C,'Notes',['Number of naive ' species_name ' cells in the central compartment']);
nCD8_P  = addspecies(model.Compartment(2),'nT',params.nCD8_P.Value/params.nCD8_div.Value,'InitialAmountUnits','cell');
    set(nCD8_P,'Notes',['Number of naive ' species_name ' cells in the peripheral compartment']);
% nCD8_T  = addspecies(model.Compartment(3),'nT',0,'InitialAmountUnits','cell');
%     set(nCD8_T,'Notes',['Number of naive ' species_name ' cells in the tumor compartment']);
nCD8_LN = addspecies(model.Compartment(4),'nT',params.nCD8_LN.Value/params.nCD8_div.Value,'InitialAmountUnits','cell');
    set(nCD8_LN,'Notes',['Number of naive ' species_name ' cells in the lymph node']);
% Activated T cells
aT = addspecies(model.Compartment(4),'aT',0,'InitialAmountUnits','cell');
    set(aT,'Notes',['Number of activated ' species_name ' cells in the lymph node']);
% aT_T = addspecies(model.Compartment(3),'aT',0,'InitialAmountUnits','cell');
%     set(aT_T,'Notes',['Number of activated ' species_name ' cells in the tumor comparment']);
% Mature T cells
T_C = addspecies(model.Compartment(1),'T',0,'InitialAmountUnits','cell');
    set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
T_P = addspecies(model.Compartment(2),'T',0,'InitialAmountUnits','cell');
    set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);
T_T = addspecies(model.Compartment(3),'T',0,'InitialAmountUnits','cell');
    set(T_T,'Notes',['Number of ' species_name ' cells in the tumor compartment']);
T_LN = addspecies(model.Compartment(4),'T',0,'InitialAmountUnits','cell');
    set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);

% Determine if first call
first_call = true;
try % add IL2 if it does not exist yet
IL2 = addspecies(model.Compartment(4),'IL2',1.9e-4,'InitialAmountUnits','nanomolarity'); % PMID: 21774806
    set(IL2,'Notes','Concentration of IL2 in the lymph node compartment');
% IL2_T = addspecies(model.Compartment(3),'IL2',1e-18,'InitialAmountUnits','molarity');
%     set(IL2_T,'Notes','Concentration of IL2 in the tumor compartment');
catch
    first_call = false;
end
% Add Hill Functions for APC/mAPC
H_APC = 'H_mAPC';

% Add Parameters
nCD8_div = addparameter(model,'nCD8_div',params.nCD8_div.Value,'ValueUnits',params.nCD8_div.Units);
    set(nCD8_div,'Notes',['T Cell Diversity ',params.nCD8_div.Notes]);
n_clones_tum = addparameter(model,'n_clones_tum',params.n_clones_tum.Value,'ValueUnits',params.n_clones_tum.Units);
    set(n_clones_tum,'Notes',['Number of T cell clones ' params.n_clones_tum.Notes]);
q_nCD8_LN_in = addparameter(model,'q_nCD8_LN_in',params.q_nCD8_LN_in.Value,'ValueUnits',params.q_nCD8_LN_in.Units);
    set(q_nCD8_LN_in,'Notes',['Rate of naive T cell transport into the LN ' params.q_nCD8_LN_in.Notes]);
q_CD8_LN_out = addparameter(model,'q_CD8_LN_out',params.q_CD8_LN_out.Value,'ValueUnits',params.q_CD8_LN_out.Units);
    set(q_CD8_LN_out,'Notes',['Rate of activated T cell transport out of the LN ' params.q_CD8_LN_out.Notes]);
q_nCD8_LN_out = addparameter(model,'q_nCD8_LN_out',params.q_nCD8_LN_out.Value,'ValueUnits',params.q_nCD8_LN_out.Units);
    set(q_nCD8_LN_out,'Notes',['Rate of naive T cell transport out of the LN ' params.q_nCD8_LN_out.Notes]);
k_nCD8_act = addparameter(model,'k_nCD8_act',params.k_nCD8_act.Value,'ValueUnits',params.k_nCD8_act.Units);
    set(k_nCD8_act,'Notes',[species_name ' activation rate ' params.k_nCD8_act.Notes]);
k_CD8_pro = addparameter(model,'k_CD8_pro',params.k_CD8_pro.Value,'ValueUnits',params.k_CD8_pro.Units);
    set(k_CD8_pro,'Notes',[species_name ' proliferation rate ' params.k_CD8_pro.Notes]);
k_CD8_death = addparameter(model,'k_CD8_death',params.k_CD8_death.Value,'ValueUnits',params.k_CD8_death.Units);
    set(k_CD8_death,'Notes',[species_name ' death rate ' params.k_CD8_death.Notes]);
q_CD8_P_in = addparameter(model,'q_CD8_P_in',params.q_CD8_P_in.Value,'ValueUnits',params.q_CD8_P_in.Units);
    set(q_CD8_P_in,'Notes',['rate of ' species_name ' transport into the peripheral compartment ' params.q_CD8_P_in.Notes]);
q_CD8_P_out = addparameter(model,'q_CD8_P_out',params.q_CD8_P_out.Value,'ValueUnits',params.q_CD8_P_out.Units);
    set(q_CD8_P_out,'Notes',['rate of ' species_name ' transport out of the peripheral compartment ' params.q_CD8_P_out.Notes]);
q_CD8_T_in = addparameter(model,'q_CD8_T_in',params.q_CD8_T_in.Value,'ValueUnits',params.q_CD8_T_in.Units);
    set(q_CD8_T_in,'Notes',['rate of ' species_name ' transport into the tumor compartment ' params.q_CD8_T_in.Notes]);
q_nCD8_P_in = addparameter(model,'q_nCD8_P_in',params.q_nCD8_P_in.Value,'ValueUnits',params.q_nCD8_P_in.Units);
    set(q_nCD8_P_in,'Notes',['rate of n' species_name ' transport into the peripheral compartment ' params.q_nCD8_P_in.Notes]);
% q_nCD8_T_in = addparameter(model,'q_nCD8_T_in',params.q_nCD8_T_in.Value,'ValueUnits',params.q_nCD8_T_in.Units);
%     set(q_nCD8_T_in,'Notes',['rate of n' species_name ' transport into the tumor compartment ' params.q_nCD8_T_in.Notes]);
q_nCD8_P_out = addparameter(model,'q_nCD8_P_out',params.q_nCD8_P_out.Value,'ValueUnits',params.q_nCD8_P_out.Units);
    set(q_nCD8_P_out,'Notes',['rate of n' species_name ' transport out of the peripheral compartment ' params.q_nCD8_P_out.Notes]);
Q_nCD8_thym = addparameter(model,'Q_nCD8_thym',params.Q_nCD8_thym.Value,'ValueUnits',params.Q_nCD8_thym.Units);
    set(Q_nCD8_thym,'Notes',['Thymic output of n',species_name,' ', params.Q_nCD8_thym.Notes]);
k_nCD8_pro = addparameter(model,'k_nCD8_pro',params.k_nCD8_pro.Value,'ValueUnits',params.k_nCD8_pro.Units);
    set(k_nCD8_pro,'Notes',['rate of n',species_name,' proliferation ', params.k_nCD8_pro.Notes]);
K_nT_pro = addparameter(model,'K_nT_pro',params.K_nT_pro.Value,'ValueUnits',params.K_nT_pro.Units);
    set(K_nT_pro,'Notes',['half-maximal peripheral proliferation of of n',species_name,' ', params.K_nT_pro.Notes]);
k_nT_death = addparameter(model,'k_nT_death',params.k_nT_death.Value,'ValueUnits',params.k_nT_death.Units);
    set(k_nT_death,'Notes',['rate of n',species_name,' death ', params.k_nT_death.Notes]);

% IL2 Parameters
if (first_call)
    p = addparameter(model,'k_IL2_deg',params.k_IL2_deg.Value,'ValueUnits',params.k_IL2_deg.Units);
        set(p,'Notes',['rate of IL2 degradation ' params.k_IL2_deg.Notes]);
    p = addparameter(model,'k_IL2_cons',params.k_IL2_cons.Value,'ValueUnits',params.k_IL2_cons.Units);
        set(p,'Notes',['rate of IL2 consumption by T cells ' params.k_IL2_cons.Notes]);
    p = addparameter(model,'k_IL2_sec',params.k_IL2_sec.Value,'ValueUnits',params.k_IL2_sec.Units);
        set(p,'Notes',['rate of IL2 secretion from T cells ' params.k_IL2_sec.Notes]);
    p = addparameter(model,'IL2_50',params.IL2_50.Value,'ValueUnits',params.IL2_50.Units);
        set(p,'Notes',['T cell activation half-maximal IL2 concentration ' params.IL2_50.Notes]);
    p = addparameter(model,'IL2_50_Treg',params.IL2_50_Treg.Value,'ValueUnits',params.IL2_50_Treg.Units);
        set(p,'Notes',['Treg activation half-maximal IL2 concentration ' params.IL2_50_Treg.Notes]);
    p = addparameter(model,'N0',params.N0.Value,'ValueUnits',params.N0.Units);
        set(p,'Notes',['numer of activated T cell generation by TCR signaling only' params.N0.Notes]);
    p = addparameter(model,'N_costim',params.N_costim.Value,'ValueUnits',params.N_costim.Units);
        set(p,'Notes',['numer of activated T cell generation by co-stimulatory signaling only ' params.N_costim.Notes]);
    p = addparameter(model,'N_IL2_CD8',params.N_IL2_CD8.Value,'ValueUnits',params.N_IL2_CD8.Units);
        set(p,'Notes',['maximum number of activated T cell generations due to IL2 ' params.N_IL2_CD8.Notes]);
    p = addparameter(model,'N_IL2_CD4',params.N_IL2_CD4.Value,'ValueUnits',params.N_IL2_CD4.Units);
        set(p,'Notes',['maximum number of activated T cell generations due to IL2 ' params.N_IL2_CD4.Notes]);
end

k_Tcell = addparameter(model,'k_Tcell',params.k_Tcell.Value,'ValueUnits',params.k_Tcell.Units);
    set(k_Tcell,'Notes',['Rate of T cell exhaustion by cancer cells ' params.k_Tcell.Notes]);
k_C_Tcell = addparameter(model,'k_C_Tcell',params.k_C_Tcell.Value,'ValueUnits',params.k_C_Tcell.Units);
    set(k_C_Tcell,'Notes',['Rate of cancer cell death by T cells ' params.k_C_Tcell.Notes]);

try % only add once
    k_Treg = addparameter(model,'k_Treg',params.k_Treg.Value,'ValueUnits',params.k_Treg.Units);
        set(k_Treg,'Notes',['Rate of T cell death by Tregs ' params.k_Treg.Notes]);
catch
end
% Antigen Default Hill Function
p = addparameter(model,['H_' antigen],1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes',['Hill function of tumor antigen ' antigen]);

% Add Reactions
% Thymic output of naive T cell
reaction = addreaction(model,'null -> V_C.nT');
    set(reaction,'ReactionRate','Q_nCD8_thym/nCD8_div');
    set(reaction,'Notes','Thymic output of naive T cell to blood');
% Naive T cell proliferation in the peripheral compartment and the LN compartment
reaction = addreaction(model,'null -> V_P.nT');
    set(reaction,'ReactionRate','k_nCD8_pro/nCD8_div*V_P.nT/(K_nT_pro/nCD8_div+V_P.nT)');
    set(reaction,'Notes','Naive T cell proliferation in the peripheral compartment');
% reaction = addreaction(model,'null -> V_C.nT');
%     set(reaction,'ReactionRate','k_nCD8_pro/nCD8_div*V_C.nT/(K_nT_pro/nCD8_div+V_C.nT)');
%     set(reaction,'Notes','Naive T cell proliferation in the central compartment');
% reaction = addreaction(model,'null -> V_T.nT');
%     set(reaction,'ReactionRate','k_nCD8_pro/nCD8_div*V_T.nT/(K_nT_pro/nCD8_div+V_T.nT)');
%     set(reaction,'Notes','Naive T cell proliferation in the Tumor compartment');
reaction = addreaction(model,'null -> V_LN.nT');
    set(reaction,'ReactionRate','k_nCD8_pro/nCD8_div*V_LN.nT/(K_nT_pro/nCD8_div+V_LN.nT)');
    set(reaction,'Notes','Naive T cell proliferation in the TDLN compartment');
% Naive T cell death in the peripheral compartment and the central compartment
reaction = addreaction(model,'V_P.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_P.nT');
    set(reaction,'Notes','Naive T cell death in the peripheral compartment');
reaction = addreaction(model,'V_C.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_C.nT');
    set(reaction,'Notes','Naive T cell death in the central compartment');
reaction = addreaction(model,'V_LN.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_LN.nT');
    set(reaction,'Notes','Naive T cell death in the TDLN compartment');
% reaction = addreaction(model,'V_T.nT -> null');
%     set(reaction,'ReactionRate','k_nT_death*V_T.nT');
%     set(reaction,'Notes','Naive T cell death in the tumor compartment');
% Naive T cell transport into and out of the peripheral compartment
reaction = addreaction(model,'V_C.nT -> V_P.nT');
    set(reaction,'ReactionRate','q_nCD8_P_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the peripheral compartment');
reaction = addreaction(model,'V_P.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_nCD8_P_out*V_P.nT');
    set(reaction,'Notes','Naive T cell exit from the peripheral compartment');
% Naive T cell transport into and out of the lymph node
reaction = addreaction(model,'V_C.nT -> V_LN.nT');
    set(reaction,'ReactionRate','q_nCD8_LN_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the lymph node');
reaction = addreaction(model,'V_LN.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_nCD8_LN_out*V_LN.nT');
    set(reaction,'Notes','Naive T cell exit from the lymph node');
% Naive T cell transport into the tumor compartment
% reaction = addreaction(model,'V_C.nT -> V_T.nT');
%     set(reaction,'ReactionRate','q_nCD8_T_in*V_T*V_C.nT');
%     set(reaction,'Notes','Naive T cell entry into the tumor compartment');

% Naive T cell activation
reaction = addreaction(model,'V_LN.nT -> null');
    set(reaction,'ReactionRate',['k_nCD8_act*' H_APC '*H_' antigen '*V_LN.nT']);
    set(reaction,'Notes','Naive T cell activation');
reaction = addreaction(model,'null -> V_LN.aT');
    set(reaction,'ReactionRate',['k_nCD8_act*' H_APC '*H_' antigen '*V_LN.nT*n_clones_tum']);
    set(reaction,'Notes','Naive T cell activation');
% Activated T cell proliferation
reaction = addreaction(model,'V_LN.aT -> null');
    set(reaction,'ReactionRate','k_CD8_pro/N_aT*V_LN.aT');
    set(reaction,'Notes',['a' species_name ' cell proliferation']);
reaction = addreaction(model,'null -> V_LN.T');
    set(reaction,'ReactionRate','k_CD8_pro/N_aT*2^(N_aT)*V_LN.aT'); % *(1-H_PD1_APC)
    set(reaction,'Notes',['a' species_name ' cell proliferation']);

% T cell Death
reaction = addreaction(model,'V_C.T -> null');
    set(reaction,'ReactionRate','k_CD8_death*V_C.T');
    set(reaction,'Notes','T cell death in the central compartment');
reaction = addreaction(model,'V_P.T -> null');
    set(reaction,'ReactionRate','k_CD8_death*V_P.T');
    set(reaction,'Notes','T cell death in the peripheral compartment');
reaction = addreaction(model,'V_LN.T -> null');
    set(reaction,'ReactionRate','k_CD8_death*V_LN.T');
    set(reaction,'Notes','T cell death in the lymph node compartment');
reaction = addreaction(model,'V_T.T -> V_T.T1_exh');
    set(reaction,'ReactionRate','k_CD8_death*V_T.T');
    set(reaction,'Notes','T cell death in the tumor compartment');

% T cell clearance upon Ag clearance
reaction = addreaction(model,'V_T.T -> V_T.T1_exh');
    set(reaction,'ReactionRate','k_cell_clear*V_T.T*(Kc_rec/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','T cell clearance upon antigen clearance');

% T cell death from Treg
reaction = addreaction(model,'V_T.T -> V_T.T1_exh');
    set(reaction,'ReactionRate','k_Treg*V_T.T*Tregs_/(V_T.T+Tregs_+cell)');
    set(reaction,'Notes','T cell death from Tregs');
% T cell exhaustion from cancer
reaction = addreaction(model,'V_T.T -> V_T.T1_exh');
    set(reaction,'ReactionRate','k_Tcell*V_T.T*C_total/(C_total+V_T.T+cell)*H_PD1_C1');
    set(reaction,'Notes','T cell death from cancer');
% T cell exhaustion from APCs
% reaction = addreaction(model,'V_T.T -> V_T.T_exh');
%     set(reaction,'ReactionRate','k_Tcell*V_T.T*V_T.mAPC/(V_T.mAPC+V_T.APC+V_T.T1+V_T.T0+cell)*H_PD1_APC');
%     set(reaction,'Notes','T cell death from mAPCs');


% T cell transport
% Central & Peripheral
reaction = addreaction(model,'V_C.T -> V_P.T');
    set(reaction,'ReactionRate','q_CD8_P_in*V_C.T');
    set(reaction,'Notes','T cell transport into the peripheral compartment');
reaction = addreaction(model,'V_P.T -> V_C.T');
    set(reaction,'ReactionRate','q_CD8_P_out*V_P.T');
    set(reaction,'Notes','T cell transport out of the peripheral compartment');
% Central & tumor
reaction = addreaction(model,'V_C.T -> V_T.T');
    set(reaction,'ReactionRate','q_CD8_T_in*V_T*V_C.T*(C_total^2/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','T cell transport into the tumor compartment');
% Central & LN
reaction = addreaction(model,'V_LN.T -> V_C.T');
    set(reaction,'ReactionRate','q_CD8_LN_out*V_LN.T');
    set(reaction,'Notes','T cell transport out of the lymph node compartment');
% IL2 Reactions
if (first_call)
    % IL2 Degradation
    reaction = addreaction(model,'V_LN.IL2 -> null');
        set(reaction,'ReactionRate','k_IL2_deg*V_LN.IL2*V_LN');
        set(reaction,'Notes','IL2 degradation');
    % IL2 Consumption
    reaction = addreaction(model,'V_LN.IL2 -> null');
        set(reaction,'ReactionRate','k_IL2_cons*V_LN.T1*V_LN.IL2/(IL2_50+V_LN.IL2)');
        set(reaction,'Notes','IL2 consumption by T cells');
end
% IL2 Secretion by Activated T Cells
reaction = addreaction(model,'null -> V_LN.IL2');
    set(reaction,'ReactionRate','k_IL2_sec*V_LN.aT');
    set(reaction,'Notes','IL2 secretion from activated T cells');

% Add Rules
if (first_call)
    % Set Number of Activated T Cell Generations
    p = addparameter(model,'N_aT',1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes',['Number of Activated CD8+ T Cell Generations (see Rules)']);

    addrule(model,'N_aT = N0 + N_costim*H_CD28_APC + N_IL2_CD8*V_LN.IL2/(IL2_50+V_LN.IL2)','repeatedAssignment');

    p = addparameter(model,'N_aT0',1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes',['Number of Activated Treg Generations (see Rules)']);

    addrule(model,'N_aT0 = N0 + N_costim*H_CD28_APC + N_IL2_CD4*V_LN.IL2/(IL2_50+V_LN.IL2)','repeatedAssignment');
end

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update Total T Cells in tumor (Rule 3)
Tcell_rule = model_rules(3);
rule = get(Tcell_rule,'Rule');
set(Tcell_rule,'Rule',[rule '+V_T.' species_name]);

% Update Total T Cells in LN (Rule 4)
Tcell_rule = model_rules(4);
rule = get(Tcell_rule,'Rule');
set(Tcell_rule,'Rule',[rule '+V_LN.' species_name]);

K_T_C = addparameter(model,'K_T_C',params.K_T_C.Value,'ValueUnits',params.K_T_C.Units);
    set(K_T_C,'Notes',['Dependence of Teff killing rate on Teff/C ratio ',params.K_T_C.Notes]);
K_T_Treg = addparameter(model,'K_T_Treg',params.K_T_Treg.Value,'ValueUnits',params.K_T_Treg.Units);
    set(K_T_Treg,'Notes',['Dependence of Teff killing rate on Teff/Treg ratio ',params.K_T_Treg.Notes]);
% Update Cancer Killing by T Cells (Rule 5)
if (exist('cancer_types','var'))
    for i = 1:length(cancer_types)
        reaction = addreaction(model,['V_T.' cancer_types{i} ' -> V_T.C_x']);
            set(reaction,'ReactionRate',['k_C_Tcell*V_T.' cancer_types{i} '*V_T.T/(K_T_C*V_T.' cancer_types{i} '+V_T.T+cell)*V_T.T/(V_T.T+K_T_Treg*Tregs_+cell)*(1-H_TGFb_Teff)*(1-H_PD1_C1)']);
            set(reaction,'Notes','Cancer cell killing by T cells');
        rule = get(model_rules(5),'Rule');
            set(model_rules(5),'Rule',[rule,'+(k_' cancer_types{i} '_death+k_' cancer_types{i} '_therapy)*V_T.' cancer_types{i} '+k_C_Tcell*V_T.' cancer_types{i} '*V_T.T/(K_T_C*V_T.' cancer_types{i} '+V_T.T+cell)*V_T.T/(V_T.T+K_T_Treg*Tregs_+cell)*(1-H_TGFb_Teff)*(1-H_PD1_C1)']);
    end
end

% Rename Objects with 'species_name'
rename(nCD8_div,['div_' species_name]);
rename(n_clones_tum,['n_' species_name '_clones']);
rename(nCD8_C,['n' species_name]);
rename(nCD8_P,['n' species_name]);
% rename(nCD8_T,['n' species_name]);
rename(nCD8_LN,['n' species_name]);
rename(aT,['a' species_name]);
% rename(aT_T,['a' species_name]);
rename(T_C,species_name);
rename(T_P,species_name);
rename(T_T,species_name);
rename(T_LN,species_name);
rename(q_CD8_LN_out,['q_' species_name '_LN_out']);
rename(k_nCD8_act,['k_' species_name '_act']);
rename(k_CD8_pro,['k_' species_name '_pro']);
rename(k_CD8_death,['k_' species_name '_death']);
rename(q_CD8_P_in,['q_' species_name '_P_in']);
rename(q_CD8_P_out,['q_' species_name '_P_out']);
rename(q_CD8_T_in,['q_' species_name '_T_in']);
rename(q_nCD8_P_in,['q_n' species_name '_P_in']);
% rename(q_nCD8_T_in,['q_n' species_name '_T_in']);
rename(q_nCD8_LN_in,['q_n' species_name '_LN_in']);
rename(q_nCD8_P_out,['q_n' species_name '_P_out']);
rename(q_nCD8_LN_out,['q_n' species_name '_LN_out']);
rename(Q_nCD8_thym,['Q_n' species_name '_thym']);
rename(k_nCD8_pro,['k_n' species_name '_pro']);
rename(K_nT_pro,['K_n' species_name '_pro']);
rename(k_nT_death,['k_n' species_name '_death']);
rename(k_Tcell,['k_' species_name]);
rename(k_C_Tcell,['k_C_' species_name]);

warning('off','SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');
