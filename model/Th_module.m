% Helper T Cell Module
%
% Models Th activation and transport
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing the default parameters
%
% Outputs: model -- SimBiology model object with new Th module


function model = Th_module(model,params)

% nT0 is naive CD4+, T0 is CD4+/Treg
species_name = ['Th'];
n_clones = 'n_T1_clones'; % number of tumor neoantigen clones

% Add Species
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

s = addspecies(model.Compartment(3),'TGFb',0,'InitialAmountUnits','nanomolarity');
    set(s,'Notes',['TGFb in the tumor comparment']);
s = addspecies(model.Compartment(3),'IFNg',0,'InitialAmountUnits','nanomolarity');
    set(s,'Notes',['IFNg in the tumor comparment']);

% Add Hill Functions for APC/mAPC
H_APC = 'H_APCh';

% Add Parameters
k_Th_act = addparameter(model,'k_Th_act',params.k_Th_act.Value,'ValueUnits',params.k_Th_act.Units);
    set(k_Th_act,'Notes',['Rate of T helper cell activation ' params.k_Th_act.Notes]);
k_Th_Treg = addparameter(model,'k_Th_Treg',params.k_Th_Treg.Value,'ValueUnits',params.k_Th_Treg.Units);
    set(k_Th_Treg,'Notes',[species_name ' differentiation rate to Treg ' params.k_Th_Treg.Notes]);
k_TGFb_Tsec = addparameter(model,'k_TGFb_Tsec',params.k_TGFb_Tsec.Value,'ValueUnits',params.k_TGFb_Tsec.Units); % 1.2e-10 12193750; 23393376 27589056 9697990
    set(k_TGFb_Tsec,'Notes',['TGF secretion rate by Treg ' params.k_TGFb_Tsec.Notes]);
k_TGFb_deg = addparameter(model,'k_TGFb_deg',params.k_TGFb_deg.Value,'ValueUnits',params.k_TGFb_deg.Units);
    set(k_TGFb_deg,'Notes',['TGF degradtion rate ' params.k_TGFb_deg.Notes]);

TGFb_50 = addparameter(model,'TGFb_50',params.TGFb_50.Value,'ValueUnits',params.TGFb_50.Units);
    set(TGFb_50,'Notes',['Half-Maximal TGFb level for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization ' params.TGFb_50.Notes]);
TGFb_50_Teff = addparameter(model,'TGFb_50_Teff',params.TGFb_50_Teff.Value,'ValueUnits',params.TGFb_50_Teff.Units);
    set(TGFb_50_Teff,'Notes',['Half-Maximal TGFb level for CD8 T cell inhibition ' params.TGFb_50_Teff.Notes]);
Kc_rec = addparameter(model,'Kc_rec',params.Kc_rec.Value,'ValueUnits',params.Kc_rec.Units);
    set(Kc_rec,'Notes',['Half-Maximal cancer cell number for T cell recruitment ' params.Kc_rec.Notes]);
TGFbase = addparameter(model,'TGFbase',params.TGFbase.Value,'ValueUnits',params.TGFbase.Units);
    set(TGFbase,'Notes',['Baseline TGFb level in breast tumor ' params.TGFbase.Notes]);
k_IFNg_sec = addparameter(model,'k_IFNg_sec',params.k_IFNg_sec.Value,'ValueUnits',params.k_IFNg_sec.Units);
    set(k_IFNg_sec,'Notes',['IFNg secretion rate by T helper cell ' params.k_IFNg_sec.Notes]);
k_IFNg_deg = addparameter(model,'k_IFNg_deg',params.k_IFNg_deg.Value,'ValueUnits',params.k_IFNg_deg.Units);
    set(k_IFNg_deg,'Notes',['IFNg degradation rate ' params.k_IFNg_deg.Notes]);
IFNg_50_ind = addparameter(model,'IFNg_50_ind',params.IFNg_50_ind.Value,'ValueUnits',params.IFNg_50_ind.Units);
    set(IFNg_50_ind,'Notes',['Half-Maximal IFNg level for PD-L1 induction ' params.IFNg_50_ind.Notes]);

p = addparameter(model,'H_TGFb',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of TGFb for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization');
addrule(model,'H_TGFb = V_T.TGFb/(V_T.TGFb+TGFb_50)','repeatedAssignment');
p = addparameter(model,'H_TGFb_Teff',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of TGFb for Teff inhibition and chemoresistance development');
addrule(model,'H_TGFb_Teff = V_T.TGFb/(V_T.TGFb+TGFb_50_Teff)','repeatedAssignment');

% Add Reactions
% Naive T cell activation
reaction = addreaction(model,'V_LN.nT0 -> null');
    set(reaction,'ReactionRate',['k_Th_act*' H_APC '*H_P1*V_LN.nT0']);
    set(reaction,'Notes','Naive T cell activation');
reaction = addreaction(model,'null -> V_LN.aT');
    set(reaction,'ReactionRate',['k_Th_act*' H_APC '*H_P1*V_LN.nT0*' n_clones]);
    set(reaction,'Notes','Naive T cell activation');
% Activated T cell proliferation
reaction = addreaction(model,'V_LN.aT -> null');
    set(reaction,'ReactionRate','(k_T0_pro/N_aTh)*V_LN.aT');
    set(reaction,'Notes',['a' species_name ' cell proliferation']);
reaction = addreaction(model,'null -> V_LN.T');
    set(reaction,'ReactionRate','(k_T0_pro/N_aTh)*2^N_aTh*V_LN.aT'); % *(1-H_PD1_APC)
    set(reaction,'Notes',['a' species_name ' cell proliferation']);
% Differentiation between Treg and Th Cells
reaction = addreaction(model,'V_T.T -> V_T.T0');
    set(reaction,'ReactionRate','k_Th_Treg*V_T.T*H_TGFb');
    set(reaction,'Notes',['Differentiation between Treg and Th Cells']);

% T cell Death
reaction = addreaction(model,'V_C.T -> null');
    set(reaction,'ReactionRate','k_T0_death*V_C.T');
    set(reaction,'Notes','T cell death in the central compartment');
reaction = addreaction(model,'V_P.T -> null');
    set(reaction,'ReactionRate','k_T0_death*V_P.T');
    set(reaction,'Notes','T cell death in the peripheral compartment');
reaction = addreaction(model,'V_LN.T -> null');
    set(reaction,'ReactionRate','k_T0_death*V_LN.T');
    set(reaction,'Notes','T cell death in the lymph node compartment');
reaction = addreaction(model,'V_T.T -> V_T.Th_exh');
    set(reaction,'ReactionRate','k_T0_death*V_T.T');
    set(reaction,'Notes','T cell death in the tumor compartment');

% T cell clearance upon Ag clearance
reaction = addreaction(model,'V_T.T -> V_T.Th_exh');
    set(reaction,'ReactionRate','k_cell_clear*V_T.T*(Kc_rec/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','T cell clearance upon antigen clearance');

% T cell transport
% Central & Peripheral
reaction = addreaction(model,'V_C.T -> V_P.T');
    set(reaction,'ReactionRate','q_T0_P_in*V_C.T');
    set(reaction,'Notes','T cell transport into the peripheral compartment');
reaction = addreaction(model,'V_P.T -> V_C.T');
    set(reaction,'ReactionRate','q_T0_P_out*V_P.T');
    set(reaction,'Notes','T cell transport out of the peripheral compartment');

% Central & tumor
reaction = addreaction(model,'V_C.T -> V_T.T');
    set(reaction,'ReactionRate','q_T0_T_in*V_T*V_C.T*(C_total^2/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','T cell transport into the tumor compartment');
% Central & LN
reaction = addreaction(model,'V_LN.T -> V_C.T');
    set(reaction,'ReactionRate','q_T0_LN_out*V_LN.T');
    set(reaction,'Notes','T cell transport out of the lymph node compartment');

% IL2 Reactions
% IL2 Secretion by Activated T Cells
reaction = addreaction(model,'null -> V_LN.IL2');
   set(reaction,'ReactionRate','k_IL2_sec*V_LN.aT');
   set(reaction,'Notes','IL2 secretion from activated T cells');

% TGFb Reactions
% TGFb secretion and consumption by cancer cells
reaction = addreaction(model,'null -> V_T.TGFb');
    set(reaction,'ReactionRate','k_TGFb_deg*(TGFbase - V_T.TGFb)*V_T');
    set(reaction,'Notes','TGFb secretion by triple-negative breast cancer cells');
% TGFb Secretion by Activated Treg Cells
reaction = addreaction(model,'null -> V_T.TGFb');
   set(reaction,'ReactionRate','k_TGFb_Tsec*V_T.T0');
   set(reaction,'Notes','TGFb secretion from activated T cells in tumor');


% IFNg Secretion by CD4 T helper Cells
reaction = addreaction(model,'null -> V_T.IFNg');
    set(reaction,'ReactionRate','k_IFNg_sec*V_T.T');
    set(reaction,'Notes','IFNg secretion from T helper cells in tumor');
% IFNg Degradation
reaction = addreaction(model,'V_T.IFNg -> null');
    set(reaction,'ReactionRate','k_IFNg_deg*V_T.IFNg');
    set(reaction,'Notes','IFNg degradtion in tumor');

% Add Rules
% Set Number of Activated T Cell Generations
p = addparameter(model,'N_aTh',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes',['Number of Activated T Helper Cell Generations (see Rules)']);
addrule(model,'N_aTh = N0 + N_costim*H_CD28_APC + N_IL2_CD4*V_LN.IL2/(IL2_50+V_LN.IL2)','repeatedAssignment');

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

% Rename Objects with 'species_name'
rename(aT,['a' species_name]);
% rename(aT_T,['a' species_name]);
rename(T_C,species_name);
rename(T_P,species_name);
rename(T_T,species_name);
rename(T_LN,species_name);

warning('off','SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');
