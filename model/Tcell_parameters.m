% Function to generate cancer module parameters from default physical parameters
%
% Inputs: params_in  -- object containing the default parameters
%                       - V_C--volume of the central compartment
%                       - V_P--volume of the peripheral compartment
%                       - gamma_P--volume fraction of blood vessels in the peripheral compartment
%                       - gamma_T--volume fraction of blood vessels in the tumour compartment
%                       - D_LN--lymph node diameter
%                       - nLNs--number of lymph nodes
%                       - n_clones_tum|n_clones_slf--number of T cell clones
%                       - nT_div_CD4|nT_div_CD8--naive T cell diversity
%                       - rho_nTCD4|rho_nTCD8--density of naive T cells
%                       - k_nTCD4|k_nTCD8--rate of naive T cell activation
%                       - k_pro_CD4|k_pro_CD8--rate of activated T cell proliferation
%                       - k_CD4_death|k_CD8_death--rate of mature T cell death
%                       - q_CD4_P_out|q_CD8_P_out--rate of T cell transport P->C
%                       - k_CD4_mig|k_CD8_mig--rate of T cell migration
%                       - rho_adh--adhesion site density
%                       - q_CD4_LN_out|q_CD8_LN_out--rate of T cell transport LN->C
%                       - k_Tcell--rate of cancer-T cell killing
%                       - k_Treg--reate of T cell death by Tregs
%                       - k_IL2_deg--rate of IL2 degradtion
%                       - k_IL2_cons--rate of IL2 consumption
%                       - k_IL2_sec--rate of IL2 secretion
%                       - IL2_50--IL2 concentration for half-maximal T cell activation
%                       - IL2_50_Treg--IL2 concentration for half-maximal Treg activation
%                       - N0--baseline number of activated T cell generations (no IL2)
%                       - N_IL2--maximum number of generations added due to IL2
%
% Output: params_out -- object containing parameters
%                       - n_clones--number of T cell clones
%                       - q_LN_in--rate of naive T cells into the LN
%                       - q_LN_nT_out--rate of naive T cells out of the LN
%                       - k_act--rate of naive T cell activation
%                       - k_pro--rate of activated T cell proliferation
%                       - k_death--rate of mature T cell death
%                       - q_P_in--rate of T cell transport C->P
%                       - q_P_out--rate of T cell transport P->C
%                       - q_T_in--rate of T cell transport C->T
%                       - q_LN_out--rate of T cell transport LN->C
%                       - k_Tcell--rate of cancer-T cell killing
%                       - k_Treg--reate of T cell death by Tregs
%                       - k_IL2_deg--rate of IL2 degradtion
%                       - k_IL2_cons--rate of IL2 consumption
%                       - k_IL2_sec--rate of IL2 secretion
%                       - IL2_50--IL2 concentration for half-maximal T cell activation
%                       - IL2_50_Treg--IL2 concentration for half-maximal Treg activation
%                       - N0--baseline number of activated T cell generations (no IL2)
%                       - N_IL2--maximum number of generations added due to IL2
%         isTreg      -- boolean indicating if model parameters should be calculated for Tregs


function params_out = Tcell_parameters(params_in,isTreg)

% Get T cell Type Specific Parameters
if (isTreg)
    % Number of T Cell Clones
    params_out.n_clones = params_in.n_clones_slf;
    % T Cell Activation Rate
    params_out.k_act = params_in.k_nTCD4;
    % T Cell Proliferation Rate
    params_out.k_pro = params_in.k_pro_CD4;
    % T Cell Decay Rate
    params_out.k_death = params_in.k_CD4_death;
    % Transport P->C
    params_out.q_P_out = params_in.q_CD4_P_out;
    % Transport LN->C
    params_out.q_LN_out = params_in.q_CD4_LN_out;
    % Transmigration Rate
    k_mig = params_in.k_CD4_mig;
else % Teff Th
    % Number of T Cell Clones
    params_out.n_clones = params_in.n_clones_tum;
    % T Cell Activation Rate
    params_out.k_act = params_in.k_nTCD8;
    % T Cell Proliferation Rate
    params_out.k_pro = params_in.k_pro_CD8;
    % T cell Decay Rate
    params_out.k_death = params_in.k_CD8_death;
    % Transport P->C
    params_out.q_P_out = params_in.q_CD8_P_out;
    % Transport LN->C
    params_out.q_LN_out = params_in.q_CD8_LN_out;
    % Transmigration Rate
    k_mig = params_in.k_CD8_mig;
end

% Calculate Lymph Node Volume
V_LN = params_in.nLNs*4/3*pi*(params_in.D_LN/2)^3;

% T cell Transport C->P
params_out.q_P_in = k_mig*params_in.rho_adh*params_in.gamma_P*params_in.V_P;
params_out.q_P_in.Notes = ['calculated based on T cell transmigration rate ' k_mig.Notes ' and T cell adhesion density ' params_in.rho_adh.Notes];

% T cell Transport C->T
params_out.q_T_in = k_mig*params_in.rho_adh*params_in.gamma_T;
params_out.q_T_in.Notes = ['calculated based on T cell transmigration rate ' k_mig.Notes ' and T cell adhesion density ' params_in.rho_adh.Notes];

% Cancer cell death by T cells
params_out.k_C_Tcell = params_in.k_C_Tcell;
% T cell exhaustion rate by cancer cells
params_out.k_Tcell = params_in.k_Tcell;
% T cell death rate by Tregs
params_out.k_Treg = params_in.k_Treg;

% IL2 Parameters
% IL2 Degradation Rate
params_out.k_IL2_deg = params_in.k_IL2_deg;
% IL2 Consumption Rate
params_out.k_IL2_cons = params_in.k_IL2_cons;
% IL2 Secretion Rate
params_out.k_IL2_sec = params_in.k_IL2_sec;
% IL2 Concentration for Half-Maximal T Cell Proliferation
params_out.IL2_50 = params_in.IL2_50;
% IL2 Concentration for Half-Maximal Treg Proliferation
params_out.IL2_50_Treg = params_in.IL2_50_Treg;
% Baseline Number of Activated T Cell Generations by TCR Signaling
params_out.N0 = params_in.N0;
% Baseline Number of Activated T Cell Generations by Co-stimulatory Signaling
params_out.N_costim = params_in.N_costim;
% Additional Number of Activated T Cell Generations Due to IL2
params_out.N_IL2_CD8 = params_in.N_IL2_CD8;
% Additional Number of Activated T Cell Generations Due to IL2
params_out.N_IL2_CD4 = params_in.N_IL2_CD4;

%% Naive T cell dynamics
if (isTreg)
  % Thymic output of naive CD4 T Cells into the blood
  params_out.Q_nT_thym = params_in.Q_CD4_thym;
  % Rate of naive T cell proliferation
  params_out.k_nT_pro = params_in.k_nT_pro_CD4;
  % Naive CD4 Transport P->C
  params_out.q_nT_P_out = params_in.q_nCD4_P_out;
  % Naive T Cell Lymph Node Exit Rate
  params_out.q_LN_nT_out = params_in.q_nCD4_LN_out;
  % Naive T Cell Lymph Node Entry Rate
  params_out.q_LN_in = params_in.q_CD4_LN_in;
  % Naive CD4+ T Cell Density in Blood
  params_out.nT_C = params_in.rho_nTCD4*params_in.V_C;
  % Naive CD8+ T Cell Density in Periphery estimated based on Braber 2012, PMID: 22365666
  params_out.nT_P = params_out.nT_C*50;
  % Naive CD4+ T Cell Density in TDLN calculated based on the steady-state naive T cell level in healthy individuals
  params_out.nT_LN = params_out.nT_C*0.042;
  % Naive CD4+ T Cell Diversity
  params_out.div = params_in.nT_div_CD4;
else
  % Thymic output of naive CD8 T Cells into the blood
  params_out.Q_nT_thym = params_in.Q_CD8_thym;
  % Rate of naive T cell proliferation
  params_out.k_nT_pro = params_in.k_nT_pro_CD8;
  % Naive CD8 Transport P->C
  params_out.q_nT_P_out = params_in.q_nCD8_P_out;
  % Naive T Cell Lymph Node Exit Rate
  params_out.q_LN_nT_out = params_in.q_nCD8_LN_out;
  % Naive T Cell Lymph Node Entry Rate
  params_out.q_LN_in = params_in.q_CD8_LN_in;
  % Naive CD8+ T Cell Density in Blood
  params_out.nT_C = params_in.rho_nTCD8*params_in.V_C;
  % Naive CD8+ T Cell Density in Periphery estimated based on Braber 2012, PMID: 22365666
  params_out.nT_P = params_out.nT_C*50;
  % Naive CD8+ T Cell Density in TDLN calculated based on the steady-state naive T cell level in healthy individuals
  params_out.nT_LN = params_out.nT_C*0.050;
  % Naive CD8+ T Cell Diversity
  params_out.div = params_in.nT_div_CD8;
end

% Naive T Cell Transport C->P
params_out.q_P_nT_in = params_in.k_nT_mig*params_in.rho_adh*params_in.gamma_P*params_in.V_P;
% Naive T Cell Transport C->T
params_out.q_T_nT_in = params_in.k_nT_mig*params_in.rho_adh*params_in.gamma_T;
% Naive T cell density for half-maximal peripheral proliferation
params_out.K_nT_pro = params_in.K_nT_pro;
% Rate of naive T cell death
params_out.k_nT_death = params_in.k_nT_death;
