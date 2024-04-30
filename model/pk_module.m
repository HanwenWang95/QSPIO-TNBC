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

if species_name == "aCD47"

    % Add Endosomal Compartment
    comp_endo = addcompartment(model,'V_endo',params.V_endo_tot.Value,'CapacityUnits',params.V_endo_tot.Units);
        set(comp_endo,'Notes',['Endosomal compartment (endo) for aCD47 PK' params.V_endo_tot.Notes]);

    % Add Species
    A_endo = addspecies(comp_endo,'aCD47',0.0,'InitialAmountUnits','nanomolarity');
        set(A_endo,'Notes',['Concentration of ' species_name ' in endosomal compartment ']);
    CD47 = addspecies(comp_C,'CD47',0.0,'InitialAmountUnits','nanomolarity'); 
        set(CD47,'Notes',['Concentration of CD47 on RBC in central compartment ']);
    CD47_aCD47 = addspecies(comp_C,'CD47_aCD47',0.0,'InitialAmountUnits','nanomolarity');
        set(CD47_aCD47,'Notes','Concentration of CD47-aCD47 complex in central compartment ');
    CD47_aCD47_CD47 = addspecies(comp_C,'CD47_aCD47_CD47',0.0,'InitialAmountUnits','nanomolarity');
        set(CD47_aCD47_CD47,'Notes','Concentration of CD47-aCD47-CD47 complex in central compartment ');
    FcRn = addspecies(comp_endo,'FcRn',params.FcRn_tot.Value,'InitialAmountUnits',params.FcRn_tot.Units);
        set(FcRn,'Notes',['Concentration of FcRn in endosomal compartment ' params.FcRn_tot.Notes]);
    FcRn_aCD47 = addspecies(comp_endo,'FcRn_aCD47',0.0,'InitialAmountUnits','nanomolarity');
        set(FcRn_aCD47,'Notes','Concentration of FcRn-aCD47 complex in endosomal compartment ');

    % Add Parameters
    RBC_CD47   = addparameter(model,'RBC_CD47',params.RBC_CD47.Value,'ValueUnits',params.RBC_CD47.Units);
        set(RBC_CD47,'Notes',['CD47 expression on RBCs ' params.RBC_CD47.Notes]);
    kint_CD47   = addparameter(model,'kint_CD47',params.kint_CD47.Value,'ValueUnits',params.kint_CD47.Units);
        set(kint_CD47,'Notes',['Internalization rate of bound CD47 on RBCs ' params.kint_CD47.Notes]);

    kint_aCD47 = addparameter(model,'kint_aCD47',params.kint_aCD47.Value,'ValueUnits',params.kint_aCD47.Units);
        set(kint_aCD47,'Notes',['Internalization rate of anti-CD47 antibody by endothelial cells ' params.kint_aCD47.Notes]);
    krec_aCD47 = addparameter(model,'krec_aCD47',params.krec_aCD47.Value,'ValueUnits',params.krec_aCD47.Units);
        set(krec_aCD47,'Notes',['Recycling rate of anti-CD47 antibody from endothelial cell back into plasma ' params.krec_aCD47.Notes]);
    koff_FcRn_aCD47 = addparameter(model,'koff_FcRn_aCD47',params.koff_FcRn_aCD47.Value,'ValueUnits',params.koff_FcRn_aCD47.Units);
        set(koff_FcRn_aCD47,'Notes',['Dissociation rate of aCD47-FcRn complex' params.koff_FcRn_aCD47.Notes]);
    kd_FcRn_aCD47   = addparameter(model,'kd_FcRn_aCD47',params.kd_FcRn_aCD47.Value,'ValueUnits',params.kd_FcRn_aCD47.Units);
        set(kd_FcRn_aCD47,'Notes',['Binding affinity of anti-CD47 antibody to FcRn ' params.kd_FcRn_aCD47.Notes]);
    kdeg_aCD47 = addparameter(model,'kdeg_aCD47',params.kdeg_aCD47.Value,'ValueUnits',params.kdeg_aCD47.Units);
        set(kdeg_aCD47,'Notes',['Degradation rate of anti-CD47 antibody in the endosomal compartment ' params.kdeg_aCD47.Notes]);
    kd_CD47_aCD47   = addparameter(model,'kd_CD47_aCD47',params.kd_CD47_aCD47.Value,'ValueUnits',params.kd_CD47_aCD47.Units);
        set(kd_CD47_aCD47,'Notes',['Binding affinity of anti-CD47 antibody and CD47 ' params.kd_CD47_aCD47.Notes]);
    kon_CD47_aCD47   = addparameter(model,'kon_CD47_aCD47',params.koff_CD47_aCD47.Value/(2*params.kd_CD47_aCD47.Value),'ValueUnits','1/minute/nanomolarity'); % KD = koff/(2*kon)
        set(kon_CD47_aCD47,'Notes',['Association rate of anti-CD47 antibody and CD47 ']);
    koff_CD47_aCD47   = addparameter(model,'koff_CD47_aCD47',params.koff_CD47_aCD47.Value,'ValueUnits',params.koff_CD47_aCD47.Units);
        set(koff_CD47_aCD47,'Notes',['Dissociation rate of anti-CD47 antibody and CD47 ' params.koff_CD47_aCD47.Notes]);
    Chi_CD47_aCD47_3D   = addparameter(model,'Chi_CD47_aCD47_3D',params.Chi_CD47_aCD47_3D.Value,'ValueUnits',params.Chi_CD47_aCD47_3D.Units);
        set(Chi_CD47_aCD47_3D,'Notes',['aCD47 antibody cross-arm binding strength ' params.Chi_CD47_aCD47_3D.Notes]);

    addrule(model,'V_C.CD47 = RBC_CD47/V_C','initialAssignment');

    % Add Reactions
    reaction = addreaction(model,'null -> V_C.CD47');
        set(reaction,'ReactionRate','kdeg_CD47*RBC_CD47/gamma_C');
        set(reaction,'Notes','Production of CD47 on red blood cell ');
    reaction = addreaction(model,'V_C.CD47 -> null');
        set(reaction,'ReactionRate','kdeg_CD47*V_C.CD47/gamma_C');
        set(reaction,'Notes','Internalization of CD47 on red blood cell ');
    reaction = addreaction(model,'V_C.CD47_aCD47 -> null');
        set(reaction,'ReactionRate','kint_CD47*V_C.CD47_aCD47/gamma_C');
        set(reaction,'Notes','Internalization of CD47-aCD47 on red blood cell ');
    reaction = addreaction(model,'V_C.CD47_aCD47_CD47 -> null');
        set(reaction,'ReactionRate','kint_CD47*V_C.CD47_aCD47_CD47/gamma_C');
        set(reaction,'Notes','Internalization of CD47-aCD47-CD47 on red blood cell ');

    reaction = addreaction(model,'V_C.A + V_C.CD47 <-> V_C.CD47_aCD47');
        set(reaction,'ReactionRate','2*kon_CD47_aCD47*V_C.A/gamma_C*V_C.CD47/gamma_C - koff_CD47_aCD47*V_C.CD47_aCD47/gamma_C');
        set(reaction,'Notes','Binding and unbinding of anti-CD47 antibody to CD47 on red blood cell ');
    reaction = addreaction(model,'V_C.CD47_aCD47 + V_C.CD47 <-> V_C.CD47_aCD47_CD47');
        set(reaction,'ReactionRate','Chi_CD47_aCD47_3D*kon_CD47_aCD47*V_C.CD47_aCD47/gamma_C*V_C.CD47/gamma_C - 2*koff_CD47_aCD47*V_C.CD47_aCD47_CD47/gamma_C');
        set(reaction,'Notes','Binding and unbinding of CD47-aCD47 complex to CD47 on red blood cell ');
    reaction = addreaction(model,'V_C.A -> V_endo.aCD47');
        set(reaction,'ReactionRate','kint_aCD47*V_C.A/gamma_C');
        set(reaction,'Notes','Internalization of anti-CD47 antibody by endothelial cells ');
    reaction = addreaction(model,'V_endo.aCD47 + V_endo.FcRn <-> V_endo.FcRn_aCD47');
        set(reaction,'ReactionRate','koff_FcRn_aCD47/kd_FcRn_aCD47*V_endo.aCD47*V_endo.FcRn - koff_FcRn_aCD47*V_endo.FcRn_aCD47');
        set(reaction,'Notes','Binding and unbinding of anti-CD47 antibody to neonatal Fc receptor in the endosomal compartment ');
    reaction = addreaction(model,'V_endo.FcRn_aCD47 -> V_C.A + V_endo.FcRn');
        set(reaction,'ReactionRate','krec_aCD47*V_endo.FcRn_aCD47');
        set(reaction,'Notes','Recycling of anti-CD47 antibody from endothelial cell back into plasma ');
    reaction = addreaction(model,'V_endo.aCD47 -> null');
        set(reaction,'ReactionRate','kdeg_aCD47*V_endo.aCD47');
        set(reaction,'Notes','Degradation of anti-CD47 antibody in the endosomal compartment ');
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
