%% Function to generate model parameters in breast cancer
%
% Output: params -- object containing parameters

function params = parameters_TNBC

% Define Cell Dimension
params.cell.Value = 1;
params.cell.Units = 'cell';
params.cell.Notes = 'A unit cell for calculation';
% Avogadro's Number
params.N_avg.Value = 6.0221409e23;
params.N_avg.Units = 'molecule/mole';
params.N_avg.Notes = '(Avogadro 1811)';

%% Compartment Volume Parameters
% Central Compartment Volume
params.V_C.Value = 5;
params.V_C.Units = 'liter';
params.V_C.Notes = '(estimated)';
% Peripheral Compartment Volume
params.V_P.Value = 60;
params.V_P.Units = 'liter';
params.V_P.Notes = '(estimated)';
% Cancer-Free Tumour Volume
params.V_Tmin.Value = 1e-6;
params.V_Tmin.Units = 'milliliter';
params.V_Tmin.Notes = '(estimated)';
% Number of Lymph Nodes
params.nLNs.Value = 17;
params.nLNs.Units = 'dimensionless';
params.nLNs.Notes = '(estimated)';
% Lymph Node Diameter
params.D_LN.Value = 5;
params.D_LN.Units = 'millimeter';
params.D_LN.Notes = '(Schmidt 2007, PMID: 17724531)';
% Lymph Nodes Compartment Volume
params.V_LN.Value = [];
params.V_LN.Units = 'milliliter';
params.V_LN.Notes = ['Volume of lumped tumor-draining lymph nodes calculated based on number of LNs ' params.nLNs.Notes ' and average LN diameter ' params.D_LN.Notes];
params.V_LN.Factors = ["nLNs", "D_LN"];
params.V_LN.Equation = 'p(1)*4/3*pi*(p(2)/2)^3';
% Cancer Cell Diameter
params.D_cell.Value = 17;
params.D_cell.Units = 'micrometer';
params.D_cell.Notes = '(Abramczyk 2015, PMID: 25730442)';
% T Cell Diameter
params.D_Tcell.Value = 6.94;
params.D_Tcell.Units = 'micrometer';
params.D_Tcell.Notes = '(Chapman 1981, PMID: 6975780)';
% Volume of a Cancer Cell
params.vol_cell.Value = [];
params.vol_cell.Units = 'micrometer^3/cell';
params.vol_cell.Notes = ['Volume of a cancer cell calculated based on cancer cell diameter ' params.D_cell.Notes];
params.vol_cell.Factors = ["D_cell", "cell"];
params.vol_cell.Equation = '4/3*pi*(p(1)/2)^3/p(2)';
% Volume of a T cell
params.vol_Tcell.Value = [];
params.vol_Tcell.Units = 'micrometer^3/cell';
params.vol_Tcell.Notes = ['Volume of a T cell calculated based on the average T cell diameter ' params.D_Tcell.Notes];
params.vol_Tcell.Factors = ["D_Tcell", "cell"];
params.vol_Tcell.Equation = '4/3*pi*(p(1)/2)^3/p(2)';
% Dead Cell Clearance Rate
params.k_cell_clear.Value = 0.1;
params.k_cell_clear.Units = '1/day';
params.k_cell_clear.Notes = '(estimated)';
% Tumor Cell Volume Fraction
params.Ve_T.Value = 0.37;
params.Ve_T.Units = 'dimensionless';
params.Ve_T.Notes = '(Finley 2012, PMID: 22547351)';

%% Cancer Parameters
% Growth Rate
params.k_C_growth.Value = 0.0072;
params.k_C_growth.Units = '1/day';
params.k_C_growth.Notes = '(Ryu 2014, PMID: 24895040; Desai 2006, PMID: 16489089)';
% Death Rate
params.k_C_death.Value = 0.0001;
params.k_C_death.Units = '1/day';
params.k_C_death.Notes = '(Palsson 2013, PMID: 24074340)';
% Maximal Tumor Capacity
params.C_max.Value = 2.7e4;
params.C_max.Units = 'cell';
params.C_max.Notes = '(Desai 2006, PMID: 16489089)';
% Initial Tumour Diameter
params.initial_tumour_diameter.Value = 2.5; % 1 - 5 cm
params.initial_tumour_diameter.Units = 'centimeter';
params.initial_tumour_diameter.Notes = '(varied)';
% Cancer Cell Density
params.rho_cell.Value = 2.06e8; % 1.67e8 - 2.06e8
params.rho_cell.Units = 'cell/milliliter';
params.rho_cell.Notes = '(Del Monte 2009, PMID: 19176997; Barnes 2016 PMID: 26332194)';

%% T Cell Parameters
% Tumor-specific T cell clone number (TMB)
params.n_clones_tum.Value = 63;
params.n_clones_tum.Units = 'dimensionless';
params.n_clones_tum.Notes = '(Morisaki 2021 PMID: 34193879; Narang 2019, PMID: 30832597)';
% Self-antigen-specific T cell clone number
params.n_clones_slf.Value = 63;
params.n_clones_slf.Units = 'dimensionless';
params.n_clones_slf.Notes = '(estimated)';
% Maximum Rate of CD8+ T Cell Activation by mAPCs
params.k_nCD8_act.Value = 23;
params.k_nCD8_act.Units = '1/day';
params.k_nCD8_act.Notes = '(De Boer and Perelson 1995, PMID: 7475092; Robertson-Tessi 2012, PMID: 22051568)';
% Maximum Rate of CD4+ Treg Activation by APCs
params.k_nCD4_act.Value = 5;
params.k_nCD4_act.Units = '1/day';
params.k_nCD4_act.Notes = '(De Boer and Perelson 1995, PMID: 7475092; Robertson-Tessi 2012, PMID: 22051568)';
% Rate of CD8+ T Cell Proliferation
params.k_CD8_pro.Value = 1.0;
params.k_CD8_pro.Units = '1/day';
params.k_CD8_pro.Notes = '(Marchingo 2014, PMID: 25430770)';
% Rate of CD4+ T Cell Proliferation
params.k_CD4_pro.Value = 1.0;
params.k_CD4_pro.Units = '1/day';
params.k_CD4_pro.Notes = '(Marchingo 2014, PMID: 25430770)';
% Rate of CD8 T Cell Decay
params.k_CD8_death.Value = 0.01;
params.k_CD8_death.Units = '1/day';
params.k_CD8_death.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Rate of CD4 T Cell Decay
params.k_CD4_death.Value = 0.01;
params.k_CD4_death.Units = '1/day';
params.k_CD4_death.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Activated CD8 Rate of Transmigration
params.k_CD8_mig.Value = 5.8e-12;
params.k_CD8_mig.Units = '1/minute/cell';
params.k_CD8_mig.Notes = '(Zhu 1996, PMID: 8706023)';
% Activated CD4 Rate of Transmigration
params.k_CD4_mig.Value = 5.8e-12;
params.k_CD4_mig.Units = '1/minute/cell';
params.k_CD4_mig.Notes = '(Zhu 1996, PMID: 8706023)';
% T cell Adhesion Density
params.rho_adh.Value = 5e8; % 5e8-1e9
params.rho_adh.Units = 'cell/centimeter^3';
params.rho_adh.Notes = '(Zhu 1996, PMID: 8706023)';
% Peripheral Vascular Volume Fractions
params.gamma_P.Value = 0.014;
params.gamma_P.Units = 'dimensionless';
params.gamma_P.Notes = '(Finley 2012, PMID: 22547351)';
% Tumour Vascular Volume Fractions
params.gamma_T.Value = 0.02;
params.gamma_T.Units = 'dimensionless';
params.gamma_T.Notes = '(Finley 2012, PMID: 22547351; Stamatelos 2014, PMID: 24342178)';
% Activated CD8+ T Cell Transport C->P
params.q_CD8_P_in.Value = [];
params.q_CD8_P_in.Units = '1/minute';
params.q_CD8_P_in.Notes = ['calculated based on T cell transmigration rate ' params.k_CD8_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_CD8_P_in.Factors = ["k_CD8_mig", "rho_adh", "gamma_P", "V_P"];
params.q_CD8_P_in.Equation = 'p(1)*p(2)*p(3)*p(4)';
% Activated CD4+ T Cell Transport C->P
params.q_CD4_P_in.Value = [];
params.q_CD4_P_in.Units = '1/minute';
params.q_CD4_P_in.Notes = ['calculated based on T cell transmigration rate ' params.k_CD4_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_CD4_P_in.Factors = ["k_CD4_mig", "rho_adh", "gamma_P", "V_P"];
params.q_CD4_P_in.Equation = 'p(1)*p(2)*p(3)*p(4)';
% Activated CD8+ T Cell Transport C->T
params.q_CD8_T_in.Value = [];
params.q_CD8_T_in.Units = '1/(centimeter^3*minute)';
params.q_CD8_T_in.Notes = ['calculated based on T cell transmigration rate ' params.k_CD8_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_CD8_T_in.Factors = ["k_CD8_mig", "rho_adh", "gamma_T"];
params.q_CD8_T_in.Equation = 'p(1)*p(2)*p(3)';
% Activated CD4+ T Cell Transport C->T
params.q_CD4_T_in.Value = [];
params.q_CD4_T_in.Units = '1/(centimeter^3*minute)';
params.q_CD4_T_in.Notes = ['calculated based on T cell transmigration rate ' params.k_CD4_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_CD4_T_in.Factors = ["k_CD4_mig", "rho_adh", "gamma_T"];
params.q_CD4_T_in.Equation = 'p(1)*p(2)*p(3)';
% Activated CD8+ T Cell Transport P->C
params.q_CD8_P_out.Value = 24;
params.q_CD8_P_out.Units = '1/day';
params.q_CD8_P_out.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Activated CD4+ T Cell Transport P->C
params.q_CD4_P_out.Value = 24;
params.q_CD4_P_out.Units = '1/day';
params.q_CD4_P_out.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Activated CD8+ T Cell Transport LN->C
params.q_CD8_LN_out.Value = 24.0;
params.q_CD8_LN_out.Units = '1/day';
params.q_CD8_LN_out.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Activated CD4+ T Cell Transport LN->C
params.q_CD4_LN_out.Value = 24.0;
params.q_CD4_LN_out.Units = '1/day';
params.q_CD4_LN_out.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';
% Death Rate Due to T Cells
params.k_C_Tcell.Value = 0.9;
params.k_C_Tcell.Units = '1/day';
params.k_C_Tcell.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Rate of T Cell Exhaustion by Cancer Cells
params.k_Tcell.Value = 0.1;
params.k_Tcell.Units = '1/day';
params.k_Tcell.Notes = '(estimated)';
% Rate of T Cell Death by Tregs
params.k_Treg.Value = 0.1;
params.k_Treg.Units = '1/day';
params.k_Treg.Notes = '(estimated)';

%% IL2 Parameters
% Degradation Rate
params.k_IL2_deg.Value = 0.2;
params.k_IL2_deg.Units = '1/minute';
params.k_IL2_deg.Notes = '(Lotze 1985, PMID: 3871099)';
% Maximum Consumption Rate by T Cells
params.k_IL2_cons.Value = 6.0e-6;
params.k_IL2_cons.Units = 'nanomole/cell/hour';
params.k_IL2_cons.Notes = '(Lotze 1985, PMID: 3871099)';
% Secretion Rate from Activated T Cells
params.k_IL2_sec.Value = 3.0e-5;
params.k_IL2_sec.Units = 'nanomole/cell/hour';
params.k_IL2_sec.Notes = '(Han 2012, PMID: 22160692; Thurley 2015, PMID: 25923703; Liu 2001)';
% IL2 Concentration for Half-Maximal T Cell Proliferation
params.IL2_50.Value = 0.32;
params.IL2_50.Units = 'nanomolarity';
params.IL2_50.Notes = '(Marchingo 2014, PMID: 25430770)';
% IL2 Concentration for Half-Maximal Treg Proliferation
params.IL2_50_Treg.Value = 0.32; % 0.0032
params.IL2_50_Treg.Units = 'nanomolarity';
params.IL2_50_Treg.Notes = '(Wang and Smith 1987, PMID: 3116143)';
% Baseline Number of Activated T Cell Generations
params.N0.Value = 2;
params.N0.Units = 'dimensionless';
params.N0.Notes = '(Marchingo 2014, PMID: 25430770)';
% Baseline Number of Activated T Cell Generations for co-stimulation
params.N_costim.Value = 3;
params.N_costim.Units = 'dimensionless';
params.N_costim.Notes = '(Marchingo 2014, PMID: 25430770)';
% Additional Number of Activated CD8+ T Cell Generations Due to IL2
params.N_IL2_CD8.Value = 11;
params.N_IL2_CD8.Units = 'dimensionless';
params.N_IL2_CD8.Notes = '(Marchingo 2014, PMID: 25430770)';
% Additional Number of Activated CD4+ T Cell Generations Due to IL2
params.N_IL2_CD4.Value = 8.5;
params.N_IL2_CD4.Units = 'dimensionless';
params.N_IL2_CD4.Notes = '(Marchingo 2014, PMID: 25430770)';

%% Naive T Cell Parameters
% Naive CD8+ T Cell Density in Blood
params.rho_nCD8.Value = 5.05e5; % 505.0, 410.5–825.4 cells/uL
params.rho_nCD8.Units = 'cell/milliliter';
params.rho_nCD8.Notes = '(Autissier 2010, PMID: 20099249)';
% Naive CD4+ T Cell Density in Blood
params.rho_nCD4.Value = 8.6e5; % 860.7, 520.8–1173.2 cells/uL
params.rho_nCD4.Units = 'cell/milliliter';
params.rho_nCD4.Notes = '(Autissier 2010, PMID: 20099249)';
% Number of Naive CD8+ T Cell in Blood
params.nCD8_C.Value = [];
params.nCD8_C.Units = 'cell';
params.nCD8_C.Notes = '';
params.nCD8_C.Factors = ["rho_nCD8", "V_C"];
params.nCD8_C.Equation = 'p(1)*p(2)';
% Number of Naive CD4+ T Cell in Blood
params.nCD4_C.Value = [];
params.nCD4_C.Units = 'cell';
params.nCD4_C.Notes = '';
params.nCD4_C.Factors = ["rho_nCD4", "V_C"];
params.nCD4_C.Equation = 'p(1)*p(2)';
% Number of Naive CD8+ T Cell in Peripheral
params.nCD8_P.Value = [];
params.nCD8_P.Units = 'cell';
params.nCD8_P.Notes = ' estimated assuming a 50x T cell number in peripheral compartment (Braber 2012, PMID: 22365666)';
params.nCD8_P.Factors = ["nCD8_C"];
params.nCD8_P.Equation = 'p(1)*50';
% Number of Naive CD4+ T Cell in Peripheral
params.nCD4_P.Value = [];
params.nCD4_P.Units = 'cell';
params.nCD4_P.Notes = ' estimated assuming a 50x T cell number in peripheral compartment (Braber 2012, PMID: 22365666)';
params.nCD4_P.Factors = ["nCD4_C"];
params.nCD4_P.Equation = 'p(1)*50';
% Number of Naive CD8+ T Cell in Lymph Nodes
params.nCD8_LN.Value = [];
params.nCD8_LN.Units = 'cell';
params.nCD8_LN.Notes = ' calculated based on the steady-state naive T cell level in healthy individuals (Autissier 2010, PMID: 20099249)';
params.nCD8_LN.Factors = ["nCD8_C"];
params.nCD8_LN.Equation = 'p(1)*0.05';
% Number of Naive CD4+ T Cell in Lymph Nodes
params.nCD4_LN.Value = [];
params.nCD4_LN.Units = 'cell';
params.nCD4_LN.Notes = ' calculated based on the steady-state naive T cell level in healthy individuals (Autissier 2010, PMID: 20099249)';
params.nCD4_LN.Factors = ["nCD4_C"];
params.nCD4_LN.Equation = 'p(1)*0.042';

% Naive CD8+ T Cell Diversity
params.nCD8_div.Value = 1.11e6;
params.nCD8_div.Units = 'dimensionless';
params.nCD8_div.Notes = '(Robins 2009, PMID: 19706884)';
% Naive CD4+ T Cell Diversity
params.nCD4_div.Value = 1.16e6;
params.nCD4_div.Units = 'dimensionless';
params.nCD4_div.Notes = '(Robins 2009, PMID: 19706884)';
% Naive T Rate of Transmigration
params.k_nT_mig.Value = 4.2e-13;
params.k_nT_mig.Units = '1/minute/cell';
params.k_nT_mig.Notes = '(Zhu 1996, PMID: 8706023)';
% Naive T cell density for half-maximal peripheral proliferation
params.K_nT_pro.Value = 1e9;
params.K_nT_pro.Units = 'cell';
params.K_nT_pro.Notes = '(Braber 2012, PMID: 22365666)';
% Rate of naive T cell death
params.k_nT_death.Value = 0.002;
params.k_nT_death.Units = '1/day';
params.k_nT_death.Notes = '(Braber 2012, PMID: 22365666)';
% Thymic output of naive CD4+ T Cells into the blood
params.Q_nCD4_thym.Value = 7e7; % median age at diagnosis of TNBC: 59
params.Q_nCD4_thym.Units = 'cell/day';
params.Q_nCD4_thym.Notes = '(Bains 2009, PMID: 19179300; Yeh 2017, PMID: 28912973)';
% Thymic output of naive CD8+ T Cells into the blood
params.Q_nCD8_thym.Value = 3.5e7;
params.Q_nCD8_thym.Units = 'cell/day';
params.Q_nCD8_thym.Notes = '(Bains 2009, PMID: 19179300; Ye 2002, PMID: 11994448)';
% Rate of naive CD4+ T Cell proliferation
params.k_nCD4_pro.Value = 3.2e8;
params.k_nCD4_pro.Units = 'cell/day';
params.k_nCD4_pro.Notes = '(Braber 2012, PMID: 22365666)';
% Rate of naive CD8+ T Cell proliferation
params.k_nCD8_pro.Value = 3.2e8;
params.k_nCD8_pro.Units = 'cell/day';
params.k_nCD8_pro.Notes = '(Braber 2012, PMID: 22365666)';
% Naive CD8+ T Cell Transport C->P
params.q_nCD8_P_in.Value = [];
params.q_nCD8_P_in.Units = '1/minute';
params.q_nCD8_P_in.Notes = ['calculated based on T cell transmigration rate ' params.k_nT_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_nCD8_P_in.Factors = ["k_nT_mig", "rho_adh", "gamma_P", "V_P"];
params.q_nCD8_P_in.Equation = 'p(1)*p(2)*p(3)*p(4)';
% Naive CD4+ T Cell Transport C->P
params.q_nCD4_P_in.Value = [];
params.q_nCD4_P_in.Units = '1/minute';
params.q_nCD4_P_in.Notes = ['calculated based on T cell transmigration rate ' params.k_nT_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_nCD4_P_in.Factors = ["k_nT_mig", "rho_adh", "gamma_P", "V_P"];
params.q_nCD4_P_in.Equation = 'p(1)*p(2)*p(3)*p(4)';
% Naive CD8+ T Cell Transport C->T
params.q_nCD8_T_in.Value = [];
params.q_nCD8_T_in.Units = '1/(centimeter^3*minute)';
params.q_nCD8_T_in.Notes = ['calculated based on T cell transmigration rate ' params.k_nT_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_nCD8_T_in.Factors = ["k_nT_mig", "rho_adh", "gamma_T"];
params.q_nCD8_T_in.Equation = 'p(1)*p(2)*p(3)';
% Naive CD4+ T Cell Transport C->T
params.q_nCD4_T_in.Value = [];
params.q_nCD4_T_in.Units = '1/(centimeter^3*minute)';
params.q_nCD4_T_in.Notes = ['calculated based on T cell transmigration rate ' params.k_nT_mig.Notes ' and T cell adhesion density ' params.rho_adh.Notes];
params.q_nCD4_T_in.Factors = ["k_nT_mig", "rho_adh", "gamma_T"];
params.q_nCD4_T_in.Equation = 'p(1)*p(2)*p(3)';
% Naive CD8+ T Cell Transport P->C
params.q_nCD8_P_out.Value = 5.1;
params.q_nCD8_P_out.Units = '1/day';
params.q_nCD8_P_out.Notes = ' calculated based on steady-state naive T cell density in healthy individuals (Autissier 2010, PMID: 20099249)';
% Naive CD4+ T Cell Transport P->C
params.q_nCD4_P_out.Value = 5.1;
params.q_nCD4_P_out.Units = '1/day';
params.q_nCD4_P_out.Notes = ' calculated based on steady-state naive T cell density in healthy individuals (Autissier 2010, PMID: 20099249)';
% Naive CD8+ T Cell Lymph Node Exit Rate
params.q_nCD8_LN_out.Value = 1.8;
params.q_nCD8_LN_out.Units = '1/day';
params.q_nCD8_LN_out.Notes = '(Mandl 2012, PMID: 23071319)';
% Naive CD4+ T Cell Lymph Node Exit Rate
params.q_nCD4_LN_out.Value = 2.88;
params.q_nCD4_LN_out.Units = '1/day';
params.q_nCD4_LN_out.Notes = '(Mandl 2012, PMID: 23071319)';
% Naive CD8+ T Cell Lymph Node Entry Rate
params.q_nCD8_LN_in.Value = 0.076;
params.q_nCD8_LN_in.Units = '1/day';
params.q_nCD8_LN_in.Notes = ' calculated based on the steady-state naive T cell level in healthy individuals (Autissier 2010, PMID: 20099249)';
% Naive CD4+ T Cell Lymph Node Entry Rate
params.q_nCD4_LN_in.Value = 0.1;
params.q_nCD4_LN_in.Units = '1/day';
params.q_nCD4_LN_in.Notes = ' calculated based on the steady-state naive T cell level in healthy individuals (Autissier 2010, PMID: 20099249)';


%% APC Module Parameters
% Rate of APC Maturation
params.k_APC_mat.Value = 1.5;
params.k_APC_mat.Units = '1/day';
params.k_APC_mat.Notes = '(Chen 2014, PMID: 25184733)';
% Rate of APC Migration
params.k_APC_mig.Value = 4.0;
params.k_APC_mig.Units = '1/day';
params.k_APC_mig.Notes = '(Russo Halin 2016, PMID: 26876174)';
% Rate of APC Death
params.k_APC_death.Value = 0.01;
params.k_APC_death.Units = '1/day';
params.k_APC_death.Notes = '(Marino and Kirschner 2004, PMID: 15038983)';
% Rate of mAPC Death
params.k_mAPC_death.Value = 0.02;
params.k_mAPC_death.Units = '1/day';
params.k_mAPC_death.Notes = '(Marino and Kirschner 2004, PMID: 15038983)';
% APC Density in Tumour
params.APC0_T.Value = 4.0e5;
params.APC0_T.Units = 'cell/milliliter';
params.APC0_T.Notes = '(Lavin 2017, PMID: 28475900)';
% APC Density in LN
params.APC0_LN.Value = 1.2e6;
params.APC0_LN.Units = 'cell/milliliter';
params.APC0_LN.Notes = '(Catron 2004, PMID: 15357945)';
% Cytokine Time Constant
params.k_c.Value = 2.0;
params.k_c.Units = '1/day';
params.k_c.Notes = '(Chen 2014, PMID: 25184733)';
% Baseline Cytokine Concentration
params.c0.Value = 1.0e-9;
params.c0.Units = 'molarity';
params.c0.Notes = '(Chen 2014, PMID: 25184733)';
% Cytokine Concentration for Half-Maximal APC Maturation
params.c50.Value = 1.0e-9;
params.c50.Units = 'molarity';
params.c50.Notes = '(Chen 2014, PMID: 25184733)';
% Concentration of Cytokines Released by Cancer Cell Death
params.DAMPs.Value = 1.34e-14;
params.DAMPs.Units = 'mole/cell';
params.DAMPs.Notes = '(Milo 2013, PMID: 24114984, Ponomarenko 2016, PMID: 27298622)';
% Maximum Number of T Cells an APC can Interact with
params.n_sites_APC.Value = 10;
params.n_sites_APC.Units = 'dimensionless';
params.n_sites_APC.Notes = '(De Boer and Perelson 1995, PMID: 7475092)';

% Antigen Module Parameters
% Number of MHC Molecule Types
params.N_MHC.Value = 1;
params.N_MHC.Units = 'dimensionless';
params.N_MHC.Notes = '';
% Total Amount of MHC
params.n_MHC_T.Value = 2e6;
params.n_MHC_T.Units = 'molecule';
params.n_MHC_T.Notes = '(Chen 2014, PMID: 25184733)';
% Total Amount of MHC per Area
params.MHC_T.Value = [];
params.MHC_T.Units = 'molecule';
params.MHC_T.Notes = ['calculated based on the total amount of MHC ' params.n_MHC_T.Notes ...
                    ' and the total surface area'];
params.MHC_T.Factors = ["n_MHC_T", "A_endo", "N_endo", "A_s"];
params.MHC_T.Equation = 'p(1)/(p(2)*p(3) + p(4))';
% Rate of MHC Internalization
params.kin.Value = 14.4;
params.kin.Units = '1/day';
params.kin.Notes = '(Chen 2014, PMID: 25184733)';
% Rate of MHC Externalization
params.kout.Value = 28.8;
params.kout.Units = '1/day';
params.kout.Notes = '(Chen 2014, PMID: 25184733)';
% Number of Endosomal Vesicles per Cell
params.N_endo.Value = 10;
params.N_endo.Units = 'dimensionless';
params.N_endo.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Endosomal Volume
params.V_endo.Value = 4.0e-17;
params.V_endo.Units = 'liter';
params.V_endo.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Endosomal Surface Area
params.A_endo.Value = 1.5;
params.A_endo.Units = 'micrometer^2';
params.A_endo.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Endosomal Surface Area
params.A_s.Value = 900.0;
params.A_s.Units = 'micrometer^2';
params.A_s.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Endosomal Volume
params.V_e.Value = [];
params.V_e.Units = 'liter';
params.V_e.Notes = ['calculated based on the volume of a single endosome ' params.V_endo.Notes ...
                        ' and the number of endosomes per cell ' params.N_endo.Notes];
params.V_e.Factors = ["V_endo", "N_endo"];
params.V_e.Equation = 'p(1)*p(2)';
% Endosomal Surface Area
params.A_e.Value = [];
params.A_e.Units = 'micrometer^2';
params.A_e.Notes = ['calculated based on the surface area of a single endosome ' params.A_endo.Notes ...
                        ' and the number of endosomes per cell ' params.N_endo.Notes];
params.A_e.Factors = ["A_endo", "N_endo"];
params.A_e.Equation = 'p(1)*p(2)';
% Surface Area of T Cells
params.A_Tcell.Value = [];
params.A_Tcell.Units = 'micrometer^2';
params.A_Tcell.Notes = ['calculated based on the average T cell diameter ' params.D_Tcell.Notes];
params.A_Tcell.Factors = ["D_Tcell"];
params.A_Tcell.Equation = '4*pi*(p(1)/2)^2';
% Surface Area of Cancer Cells
params.A_cell.Value = [];
params.A_cell.Units = 'micrometer^2';
params.A_cell.Notes = ['calculated based on the average cancer cell diameter ' params.D_cell.Notes];
params.A_cell.Factors = ["D_cell"];
params.A_cell.Equation = '4*pi*(p(1)/2)^2';
% Surface Area of APC Cells
params.A_APC.Value = 900.0;
params.A_APC.Units = 'micrometer^2';
params.A_APC.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Rate of Antigen Uptake
params.k_up.Value = 14.4;
params.k_up.Units = '1/day/cell';
params.k_up.Notes = '(Chen 2014, PMID: 25184733)';
% Rate of Extracellular Antigen Degradation
params.k_xP_deg.Value = 2.0;
params.k_xP_deg.Units = '1/day';
params.k_xP_deg.Notes = '(Palsson 2013, PMID: 24074340)';
% Rate of Endosomal Antigen Degradation
params.k_P_deg.Value = 17.28;
params.k_P_deg.Units = '1/day';
params.k_P_deg.Notes = '(Chen 2014, PMID: 25184733)';
% Rate of Endosomal Epitope Degradation
params.k_p_deg.Value = 144.0;
params.k_p_deg.Units = '1/day';
params.k_p_deg.Notes = '(Chen 2014, PMID: 25184733)';
% Rate of Antigen Binding
params.k_on.Value = 1.44e5;
params.k_on.Units = '1/day/molarity';
params.k_on.Notes = '(Agrawal and Linderman 1996, PMID: 8944895)';
% Number of Epitope Molecules for Half-Maximal T Cell Activation
params.N_p_50.Value = 1e-3;
params.N_p_50.Units = 'molecule';
params.N_p_50.Notes = '(Kimachi 1997, PMID: 9464819)';
% TCR-pMHC Concentration for Half-Maximal T Cell Activation
params.p_50.Value = [];
params.p_50.Units = 'molecule/micrometer^2';
params.p_50.Notes = ['calculated based on the number of molecules for half-maximal T cell activation '...
                         params.N_p_50.Notes ' and synapse surface area'];
params.p_50.Factors = ["N_p_50","A_syn"];
params.p_50.Equation = 'p(1)/p(2)';

% TCR signal - kinetic proofreading with limited signaling
% Number of TCR molecules on naive T cells
params.TCR_tot_abs.Value = 15708;
params.TCR_tot_abs.Units = 'molecule';
params.TCR_tot_abs.Notes = '(Lever 2014, PMID: 25145757)';
% TCR molecules density on naive T cells
params.TCR_tot.Value = [];
params.TCR_tot.Units = 'molecule/micrometer^2';
params.TCR_tot.Notes = ['calculated based on the number of molecules for TCR '...
                         params.TCR_tot_abs.Notes ' and T cell surface area'];
params.TCR_tot.Factors = ["TCR_tot_abs","A_Tcell"];
params.TCR_tot.Equation = 'p(1)/p(2)';

% Rate of modification of TCRs
params.k_TCR_p.Value = 1;
params.k_TCR_p.Units = '1/second';
params.k_TCR_p.Notes = '(Lever 2014, PMID: 25145757)';
% Unbinding rate of ag/MHC to TCR
params.k_TCR_off.Value = 1;
params.k_TCR_off.Units = '1/second';
params.k_TCR_off.Notes = '(Lever 2014, PMID: 25145757)';
% binding rate of ag/MHC to TCR
params.k_TCR_on.Value = 1e-0; %3.1831e-5;
params.k_TCR_on.Units = '1/(second*molecule/micrometer^2)';
params.k_TCR_on.Notes = '(Lever 2014, PMID: 25145757)';
% Rate of modification of TCR that leads to non-signaling
params.phi_TCR.Value = 0.09;
params.phi_TCR.Units = '1/second';
params.phi_TCR.Notes = '(Lever 2014, PMID: 25145757)';
% Number of intermediate steps
params.N_TCR.Value = 10;
params.N_TCR.Units = 'dimensionless';
params.N_TCR.Notes = '(Lever 2014, PMID: 25145757)';

%% Checkpoint Module Parameters
% Surface area of the synapse
params.A_syn.Value = 37.8;
params.A_syn.Units = 'micrometer^2';
params.A_syn.Notes = '(Jansson 2005, PMID: 16034096)';
% The synapse gap distance (kd2D = kd3D*d_syn)
params.d_syn.Value = 3.0;
params.d_syn.Units = 'nanometer';
params.d_syn.Notes = '(Jansson 2005, PMID: 16034096)';
% Expression rate of PDL1 on tumor cells
params.k_out_PDL1.Value = 5e4;
params.k_out_PDL1.Units = 'molecule/day';
params.k_out_PDL1.Notes = '(Mimura 2018, PMID: 29034543)';
% Degradation rate of PDL1 on tumor cells
params.k_in_PDL1.Value = 1;
params.k_in_PDL1.Units = '1/day';
params.k_in_PDL1.Notes = '(Hsu 2018, PMID: 30442814)';
% PD1/PDL1 Parameters
% PD1-PDL1 kd
params.kd_PD1_PDL1.Value = 8.2;
params.kd_PD1_PDL1.Units = 'micromolarity';
params.kd_PD1_PDL1.Notes = '(Cheng 2013, PMID: 23417675)';
% PD1-PDL1 kon
params.kon_PD1_PDL1_3D.Value = 0.18;
params.kon_PD1_PDL1_3D.Units = '1/(micromolarity*second)';
params.kon_PD1_PDL1_3D.Notes = '(Cheng 2013, PMID: 23417675)';
% PD1-PDL1 kon for 2D
params.kon_PD1_PDL1.Value = [];
params.kon_PD1_PDL1.Units = '';
params.kon_PD1_PDL1.Notes = ['kon for PD1-PDL1 ' params.kon_PD1_PDL1_3D.Notes];
params.kon_PD1_PDL1.Factors = ["kon_PD1_PDL1_3D","d_syn"];
params.kon_PD1_PDL1.Equation = 'p(1)/p(2)';
% PD1-PDL2 kd
params.kd_PD1_PDL2.Value = 2.3;
params.kd_PD1_PDL2.Units = 'micromolarity';
params.kd_PD1_PDL2.Notes = '(Cheng 2013, PMID: 23417675)';
% PD1-PDL2 kon
params.kon_PD1_PDL2_3D.Value = 0.24;
params.kon_PD1_PDL2_3D.Units = '1/(micromolarity*second)';
params.kon_PD1_PDL2_3D.Notes = '(Cheng 2013, PMID: 23417675)';
% PD1-PDL2 kon for 2D
params.kon_PD1_PDL2.Value = [];
params.kon_PD1_PDL2.Units = '';
params.kon_PD1_PDL2.Notes = '(Cheng 2013, PMID: 23417675)';
params.kon_PD1_PDL2.Factors = ["kon_PD1_PDL2_3D","d_syn"];
params.kon_PD1_PDL2.Equation = 'p(1)/p(2)';
% PD1-aPD1 kd
params.kd_PD1_aPD1.Value = 4;
params.kd_PD1_aPD1.Units = 'nanomolarity';
params.kd_PD1_aPD1.Notes = '(Wang 2014, PMID: 24872026; Brown 2020, PMID: 32134960)';
% PD1-aPD1 kon
params.kon_PD1_aPD1.Value = 6.7e5;
params.kon_PD1_aPD1.Units = '1/(molarity*second)';
params.kon_PD1_aPD1.Notes = '(Wang 2014, PMID: 24872026; Brown 2020, PMID: 32134960)';
% PD1-aPD1 Chi (antibody cross-arm binding strength)
params.Chi_PD1_aPD1_3D.Value = 100;
params.Chi_PD1_aPD1_3D.Units = 'dimensionless';
params.Chi_PD1_aPD1_3D.Notes = '(estimated based on Wang 2014, PMID: 24872026)';
% PD1-aPD1 Chi for 2D
params.Chi_PD1_aPD1.Value = [];
params.Chi_PD1_aPD1.Units = '';
params.Chi_PD1_aPD1.Notes = '(estimated based on Wang 2014, PMID: 24872026)';
params.Chi_PD1_aPD1.Factors = ["Chi_PD1_aPD1_3D","d_syn","N_avg"];
params.Chi_PD1_aPD1.Equation = 'p(1)/(p(2)*p(3))';
% PDL1-aPDL1 kd
params.kd_PDL1_aPDL1.Value = 0.4;
params.kd_PDL1_aPDL1.Units = 'nanomolarity';
params.kd_PDL1_aPDL1.Notes = '(Wang 2014, PMID: 24872026)';
% PDL1-aPDL1 kon
params.kon_PDL1_aPDL1.Value = 4.3e5;
params.kon_PDL1_aPDL1.Units = '1/(molarity*second)';
params.kon_PDL1_aPDL1.Notes = '(Wang 2014, PMID: 24872026)';
% PDL1-aPDL1 Chi (antibody cross-arm binding strength)
params.Chi_PDL1_aPDL1_3D.Value = 100;
params.Chi_PDL1_aPDL1_3D.Units = 'dimensionless';
params.Chi_PDL1_aPDL1_3D.Notes = '(estimated)';
% PDL1-aPDL1 Chi for 2D
params.Chi_PDL1_aPDL1.Value = [];
params.Chi_PDL1_aPDL1.Units = '';
params.Chi_PDL1_aPDL1.Notes = '(Wang 2014, PMID: 24872026)';
params.Chi_PDL1_aPDL1.Factors = ["Chi_PDL1_aPDL1_3D","d_syn","N_avg"];
params.Chi_PDL1_aPDL1.Equation = 'p(1)/(p(2)*p(3))';

% PD1-PDL1 koff
params.koff_PD1_PDL1.Value = [];
params.koff_PD1_PDL1.Units = '1/second';
params.koff_PD1_PDL1.Notes = ['calculated based on the measured kd and kon ' params.kd_PD1_PDL1.Notes];
params.koff_PD1_PDL1.Factors = ["kon_PD1_PDL1_3D","kd_PD1_PDL1"];
params.koff_PD1_PDL1.Equation = 'p(1)*p(2)';
% PD1-PDL2 koff
params.koff_PD1_PDL2.Value = [];
params.koff_PD1_PDL2.Units = '1/second';
params.koff_PD1_PDL2.Notes = ['calculated based on the measured kd and kon ' params.kd_PD1_PDL2.Notes];
params.koff_PD1_PDL2.Factors = ["kon_PD1_PDL2_3D","kd_PD1_PDL2"];
params.koff_PD1_PDL2.Equation = 'p(1)*p(2)';
% PD1-aPD1 koff
params.koff_PD1_aPD1.Value = [];
params.koff_PD1_aPD1.Units = '1/second';
params.koff_PD1_aPD1.Notes = ['calculated based on the measured kd and kon ' params.kd_PD1_aPD1.Notes];
params.koff_PD1_aPD1.Factors = ["kon_PD1_aPD1","kd_PD1_aPD1"];
params.koff_PD1_aPD1.Equation = 'p(1)*p(2)';
% PDL1-aPDL1 koff
params.koff_PDL1_aPDL1.Value = [];
params.koff_PDL1_aPDL1.Units = '1/second';
params.koff_PDL1_aPDL1.Notes = ['calculated based on the measured kd and kon ' params.kd_PDL1_aPDL1.Notes];
params.koff_PDL1_aPDL1.Factors = ["kon_PDL1_aPDL1","kd_PDL1_aPDL1"];
params.koff_PDL1_aPDL1.Equation = 'p(1)*p(2)';

% PD1 Expression on T Cells
params.T8_PD1.Value = 3.1e3*20*.45;
params.T8_PD1.Units = 'molecule';
params.T8_PD1.Notes = '(Cheng 2014, PMID: 23417675; Mkrtichyan 2012, PMID: 22837483; Gatalica 2014, PMID: 25392179)';
% PDL1 Expression on T Cells
params.T8_PDL1.Value = 9.3e3*20*.45;
params.T8_PDL1.Units = 'molecule';
params.T8_PDL1.Notes = '(Cheng 2014, PMID: 23417675; Mkrtichyan 2012, PMID: 22837483; Gatalica 2014, PMID: 25392179)';
% Average Baseline PDL1 Expression on Tumor/Immune cells in tumor
params.C_PDL1.Value = 8e4*20/6*.45;
params.C_PDL1.Units = 'molecule';
params.C_PDL1.Notes = '(Cheng 2014, PMID: 23417675; Mkrtichyan 2012, PMID: 22837483; Shin 2017, PMID: 27903500; Gatalica 2014, PMID: 25392179)';
% Average Baseline PDL1 Expression on mAPCs in TDLN
params.APC_PDL1.Value = 8e4*20/6*.45;
params.APC_PDL1.Units = 'molecule';
params.APC_PDL1.Notes = '(Cheng 2014, PMID: 23417675; Mkrtichyan 2012, PMID: 22837483; Shin 2017, PMID: 27903500; Gatalica 2014, PMID: 25392179)';
% PD1/PDL1 Concentration for Half-Maximal T Cell Killing
params.PD1_50.Value = 6;
params.PD1_50.Units = 'molecule/micrometer^2';
params.PD1_50.Notes = '(Jafarnejad 2019, PMID: 31236847)';
% Hill Coefficient for PD1/PDL1
params.n_PD1.Value = 2;
params.n_PD1.Units = 'dimensionless';
params.n_PD1.Notes = '(estimated)';
% PDL2/PDL1 ratio in tumor
params.r_PDL2C.Value = 0.07;
params.r_PDL2C.Units = 'dimensionless';
params.r_PDL2C.Notes = '(Cheng 2014, PMID: 23417675)';
% PDL2/PDL1 ratio on mAPCs
params.r_PDL2APC.Value = 0.07;
params.r_PDL2APC.Units = 'dimensionless';
params.r_PDL2APC.Notes = '(Cheng 2014, PMID: 23417675)';
% Number of folds increase of PDL1 expression by IFNg
params.r_PDL1_IFNg.Value = 6;
params.r_PDL1_IFNg.Units = 'dimensionless';
params.r_PDL1_IFNg.Notes = '(Shin 2017, PMID: 27903500)';
% PD1/PDL2 Concentration for Half-Maximal T Cell Killing
% params.PD1_PDL2_50.Value = 170/6;
% params.PD1_PDL2_50.Units = 'dimensionless';
% params.PD1_PDL2_50.Notes = '(Latchman 2001, PMID: 11224527)';

%% CD28/CTLA4/CD80/CD86 Parameters
% CD28-CD80 kd
params.kd_CD28_CD80.Value = 4.0;
params.kd_CD28_CD80.Units = 'micromolarity';
params.kd_CD28_CD80.Notes = '(van der Merwe 1997, PMID: 9053440)';
% CD28-CD80 kon
params.kon_CD28_CD80_3D.Value = 0.4;
params.kon_CD28_CD80_3D.Units = '1/(micromolarity*second)';
params.kon_CD28_CD80_3D.Notes = '(van der Merwe 1997, PMID: 9053440)';
% CD28-CD80 kon for 2D
params.kon_CD28_CD80.Value = [];
params.kon_CD28_CD80.Units = '';
params.kon_CD28_CD80.Notes = '(van der Merwe 1997, PMID: 9053440)';
params.kon_CD28_CD80.Factors = ["kon_CD28_CD80_3D","d_syn"];
params.kon_CD28_CD80.Equation = 'p(1)/p(2)';
% CD28-CD86 kd
params.kd_CD28_CD86.Value = 20.0;
params.kd_CD28_CD86.Units = 'micromolarity';
params.kd_CD28_CD86.Notes = '(Collins 2002, PMID: 12196291)';
% CD28-CD86 kon
params.kon_CD28_CD86_3D.Value = 1.4;
params.kon_CD28_CD86_3D.Units = '1/(micromolarity*second)';
params.kon_CD28_CD86_3D.Notes = '(Collins 2002, PMID: 12196291)';
% CD28-CD86 kon for 2D
params.kon_CD28_CD86.Value = [];
params.kon_CD28_CD86.Units = '';
params.kon_CD28_CD86.Notes = '(Collins 2002, PMID: 12196291)';
params.kon_CD28_CD86.Factors = ["kon_CD28_CD86_3D","d_syn"];
params.kon_CD28_CD86.Equation = 'p(1)/p(2)';
% CTLA4-CD80 kd
params.kd_CTLA4_CD80.Value = 0.42;
params.kd_CTLA4_CD80.Units = 'micromolarity';
params.kd_CTLA4_CD80.Notes = '(van der Merwe 1997, PMID: 9053440)';
% CTLA4-CD80 kon
params.kon_CTLA4_CD80_3D.Value = 1.0;
params.kon_CTLA4_CD80_3D.Units = '1/(micromolarity*second)';
params.kon_CTLA4_CD80_3D.Notes = '(van der Merwe 1997, PMID: 9053440)';
% CTLA4-CD80 kon for 2D
params.kon_CTLA4_CD80.Value = [];
params.kon_CTLA4_CD80.Units = '';
params.kon_CTLA4_CD80.Notes = '(van der Merwe 1997, PMID: 9053440)';
params.kon_CTLA4_CD80.Factors = ["kon_CTLA4_CD80_3D","d_syn"];
params.kon_CTLA4_CD80.Equation = 'p(1)/p(2)';
% CTLA4-CD86 kd
params.kd_CTLA4_CD86.Value = 2.6;
params.kd_CTLA4_CD86.Units = 'micromolarity';
params.kd_CTLA4_CD86.Notes = '(Jansson 2005, PMID: 16034096; Collins 2002, PMID: 12196291)';
% CTLA4-CD86 kon
params.kon_CTLA4_CD86_3D.Value = 2.0;
params.kon_CTLA4_CD86_3D.Units = '1/(micromolarity*second)';
params.kon_CTLA4_CD86_3D.Notes = '(Jansson 2005, PMID: 16034096; Collins 2002, PMID: 12196291)';
% CTLA4-CD86 kon for 2D
params.kon_CTLA4_CD86.Value = [];
params.kon_CTLA4_CD86.Units = '';
params.kon_CTLA4_CD86.Notes = '(Jansson 2005, PMID: 16034096; Collins 2002, PMID: 12196291)';
params.kon_CTLA4_CD86.Factors = ["kon_CTLA4_CD86_3D","d_syn"];
params.kon_CTLA4_CD86.Equation = 'p(1)/p(2)';
% CD80-PDL1 kd
params.kd_CD80_PDL1.Value = 18.8;
params.kd_CD80_PDL1.Units = 'micromolarity';
params.kd_CD80_PDL1.Notes = '(Cheng 2014, PMID: 23417675; Butte 2008, PMID: 18585785)';
% CD80-PDL1 kon
params.kon_CD80_PDL1_3D.Value = 0.32;
params.kon_CD80_PDL1_3D.Units = '1/(micromolarity*second)';
params.kon_CD80_PDL1_3D.Notes = '(Cheng 2014, PMID: 23417675)';
% CD80-PDL1 kon for 2D
params.kon_CD80_PDL1.Value = [];
params.kon_CD80_PDL1.Units = '';
params.kon_CD80_PDL1.Notes = '(Cheng 2014, PMID: 23417675)';
params.kon_CD80_PDL1.Factors = ["kon_CD80_PDL1_3D","d_syn"];
params.kon_CD80_PDL1.Equation = 'p(1)/p(2)';
% CD80-CD80 kd
params.kd_CD80_CD80.Value = 17;
params.kd_CD80_CD80.Units = 'micromolarity';
params.kd_CD80_CD80.Notes = '(Ikemizu 2000, PMID: 10661405)';
% CD80-CD80 kon
params.kon_CD80_CD80_3D.Value = 5.9e-4;
params.kon_CD80_CD80_3D.Units = '1/(micromolarity*second)';
params.kon_CD80_CD80_3D.Notes = '(Ikemizu 2000, PMID: 10661405)';
% CD80-CD80 kon for 2D
params.kon_CD80_CD80.Value = [];
params.kon_CD80_CD80.Units = '';
params.kon_CD80_CD80.Notes = '(Ikemizu 2000, PMID: 10661405)';
params.kon_CD80_CD80.Factors = ["kon_CD80_CD80_3D","d_syn"];
params.kon_CD80_CD80.Equation = 'p(1)/p(2)';
% CTLA4-aCTLA4 kd
params.kd_CTLA4_aCTLA4.Value = 18.2;
params.kd_CTLA4_aCTLA4.Units = 'nanomolarity';
params.kd_CTLA4_aCTLA4.Notes = '(Wang 2014, PMID: 24872026)';
% CTLA4-aCTLA4 kon
params.kon_CTLA4_aCTLA4.Value = 3.83e5;
params.kon_CTLA4_aCTLA4.Units = '1/(molarity*second)';
params.kon_CTLA4_aCTLA4.Notes = '(Wang 2014, PMID: 24872026)';
% CTLA4-aCTLA4 Chi (antibody cross-arm binding strength)
params.Chi_CTLA4_aCTLA4_3D.Value = 100;
params.Chi_CTLA4_aCTLA4_3D.Units = 'dimensionless';
params.Chi_CTLA4_aCTLA4_3D.Notes = '(estimated)';
% CTLA4-aCTLA4 Chi for 2D
params.Chi_CTLA4_aCTLA4.Value = [];
params.Chi_CTLA4_aCTLA4.Units = '';
params.Chi_CTLA4_aCTLA4.Notes = '(estimated)';
params.Chi_CTLA4_aCTLA4.Factors = ["Chi_CTLA4_aCTLA4_3D","d_syn","N_avg"];
params.Chi_CTLA4_aCTLA4.Equation = 'p(1)/(p(2)*p(3))';

% CD28-CD80 koff
params.koff_CD28_CD80.Value = [];
params.koff_CD28_CD80.Units = '1/second';
params.koff_CD28_CD80.Notes = ['calculated based on the measured kd and kon ' params.kd_CD28_CD80.Notes];
params.koff_CD28_CD80.Factors = ["kon_CD28_CD80_3D","kd_CD28_CD80"];
params.koff_CD28_CD80.Equation = 'p(1)*p(2)';
% CD28-CD86 koff
params.koff_CD28_CD86.Value = [];
params.koff_CD28_CD86.Units = '1/second';
params.koff_CD28_CD86.Notes = ['calculated based on the measured kd and kon ' params.kd_CD28_CD86.Notes];
params.koff_CD28_CD86.Factors = ["kon_CD28_CD86_3D","kd_CD28_CD86"];
params.koff_CD28_CD86.Equation = 'p(1)*p(2)';
% CTLA4-CD80 koff
params.koff_CTLA4_CD80.Value = [];
params.koff_CTLA4_CD80.Units = '1/second';
params.koff_CTLA4_CD80.Notes = ['calculated based on the measured kd and kon ' params.kd_CTLA4_CD80.Notes];
params.koff_CTLA4_CD80.Factors = ["kon_CTLA4_CD80_3D","kd_CTLA4_CD80"];
params.koff_CTLA4_CD80.Equation = 'p(1)*p(2)';
% CTLA4-CD86 koff
params.koff_CTLA4_CD86.Value = [];
params.koff_CTLA4_CD86.Units = '1/second';
params.koff_CTLA4_CD86.Notes = ['calculated based on the measured kd and kon ' params.kd_CTLA4_CD86.Notes];
params.koff_CTLA4_CD86.Factors = ["kon_CTLA4_CD86_3D","kd_CTLA4_CD86"];
params.koff_CTLA4_CD86.Equation = 'p(1)*p(2)';
% CD80-PDL1 koff
params.koff_CD80_PDL1.Value = [];
params.koff_CD80_PDL1.Units = '1/second';
params.koff_CD80_PDL1.Notes = ['calculated based on the measured kd and kon ' params.kd_CD80_PDL1.Notes];
params.koff_CD80_PDL1.Factors = ["kon_CD80_PDL1_3D","kd_CD80_PDL1"];
params.koff_CD80_PDL1.Equation = 'p(1)*p(2)';
% CD80-CD80 koff
params.koff_CD80_CD80.Value = [];
params.koff_CD80_CD80.Units = '1/second';
params.koff_CD80_CD80.Notes = ['calculated based on the measured kd and kon ' params.kd_CD80_CD80.Notes];
params.koff_CD80_CD80.Factors = ["kon_CD80_CD80_3D","kd_CD80_CD80"];
params.koff_CD80_CD80.Equation = 'p(1)*p(2)';
% CTLA4-aCTLA4 koff
params.koff_CTLA4_aCTLA4.Value = [];
params.koff_CTLA4_aCTLA4.Units = '1/second';
params.koff_CTLA4_aCTLA4.Notes = ['calculated based on the measured kd and kon ' params.kd_CTLA4_aCTLA4.Notes];
params.koff_CTLA4_aCTLA4.Factors = ["kon_CTLA4_aCTLA4","kd_CTLA4_aCTLA4"];
params.koff_CTLA4_aCTLA4.Equation = 'p(1)*p(2)';

% CD28 Expression on T Cells
params.T8_CD28.Value = 9200*20;
params.T8_CD28.Units = 'molecule';
params.T8_CD28.Notes = '(Jansson 2005, PMID: 16034096)';
% CTLA4 Expression in T Cell synapse
params.T8_CTLA4.Value = 400*20;
params.T8_CTLA4.Units = 'molecule';
params.T8_CTLA4.Notes = '(Jansson 2005, PMID: 16034096)';
% CD80 Expression on Cancer Cells
params.C_CD80.Value = 4000*20;
params.C_CD80.Units = 'molecule';
params.C_CD80.Notes = '(Jansson 2005, PMID: 16034096)';
% CD86 Expression on Cancer Cells
params.C_CD86.Value = 43000*20;
params.C_CD86.Units = 'molecule';
params.C_CD86.Notes = '(Jansson 2005, PMID: 16034096)';
% CD80 Expression on mAPCs
params.APC_CD80.Value = 4000*20;
params.APC_CD80.Units = 'molecule';
params.APC_CD80.Notes = '(Jansson 2005, PMID: 16034096)';
% CD86 Expression on mAPCs
params.APC_CD86.Value = 43000*20;
params.APC_CD86.Units = 'molecule';
params.APC_CD86.Notes = '(Jansson 2005, PMID: 16034096)';
% CD28/CD80/CD86 Concentration for Half-Maximal T Cell Killing
params.CD28_CD8X_50.Value = 200;
params.CD28_CD8X_50.Units = 'molecule/micrometer^2';
params.CD28_CD8X_50.Notes = 'estimated assuming T cell activation is inhibited by 70% at baseline CTLA-4 expression (Parry 2005, PMID: 16227604)';
% Hill Coefficient for CD28/CD80/CD86
params.n_CD28_CD8X.Value = 2;
params.n_CD28_CD8X.Units = 'dimensionless';
params.n_CD28_CD8X.Notes = '(estimated)';

%% Antibody Transport Properties (parameters are fitted and stored in pk_parameters.m)
% Permeability
% Central <-> Peripheral
params.Dk_P.Value = 3.0e-8;
params.Dk_P.Units = 'centimeter/second';
params.Dk_P.Notes = '(Finley 2012)';
% Central <-> Tumour
params.Dk_T.Value = 3.0e-7;
params.Dk_T.Units = 'centimeter/second';
params.Dk_T.Notes = '(Finley 2012)';
% Central <-> LN
params.Dk_LN.Value = 3.0e-8;
params.Dk_LN.Units = 'centimeter/second';
params.Dk_LN.Notes = '(Padera 2017)';
% Surface Area Density
% Peripheral
params.SA_P.Value = 105000.0;
params.SA_P.Units = 'centimeter^2/liter';
params.SA_P.Notes = '(Finley 2012)';
% Tumour
params.SA_T.Value = 105000.0;
params.SA_T.Units = 'centimeter^2/liter';
params.SA_T.Notes = '(Finley 2012)';
% LN
params.SA_LN.Value = 105000.0;
params.SA_LN.Units = 'centimeter^2/liter';
params.SA_LN.Notes = '(Finley 2012)';
% Lymphatic Drainage
params.Q_LD.Value = 0.0015;
params.Q_LD.Units = '1/minute';
params.Q_LD.Notes = '(Zhu 1996, PMID: 8706023)';

%% MDSC Module Parameters
% MDSC recruitment rate
params.k_MDSC_mig.Value = 1.1e4;
params.k_MDSC_mig.Units = 'cell/(milliliter*day)';
params.k_MDSC_mig.Notes = '(Diaz-Montero 2009, PMID: 18446337; Huang 2006, PMID: 17257744)';
% rate of MDSC death
params.k_MDSC_death.Value = 0.015;
params.k_MDSC_death.Units = '1/day';
params.k_MDSC_death.Notes = '(Lai 2018, PMID: 29735668)';
% half-maximal entinostat concentration for anti-proliferative effect on tumor cells
params.IC50_inostat_C.Value = 3.74e-7;
params.IC50_inostat_C.Units = 'molarity';
params.IC50_inostat_C.Notes = '(Bouchain 2002, PMID: 12593661; Lee 2001, PMID: 11221885)';

% rate of CCL2 degradtion
params.k_CCL2_deg.Value = 0.06;
params.k_CCL2_deg.Units = '1/hour';
params.k_CCL2_deg.Notes = '(Tanimoto 2007, PMID: 18089573)';
% rate of NO degradtion
params.k_NO_deg.Value = 135;
params.k_NO_deg.Units = '1/day';
params.k_NO_deg.Notes = '(Hakim 1996, PMID: 8953625)';
% rate of ArgI degradtion
params.k_ArgI_deg.Value = 0.173;
params.k_ArgI_deg.Units = '1/day';
params.k_ArgI_deg.Notes = '(Schimke 1964, PMID: 14257612)';
% rate of CCL2 secretion
params.k_CCL2_sec.Value = 1.7e-12; % normal: 1.065e-11; TNBC: (14.2 +/- 6)*10^-11
params.k_CCL2_sec.Units = 'nanomole/cell/day';
params.k_CCL2_sec.Notes = '(Huang 2007, PMID: 17257744; Dutta 2018, PMID: 29594759)';
% rate of NO secretion
params.k_NO_sec.Value = 4.8e-7;
params.k_NO_sec.Units = 'nanomole/cell/day';
params.k_NO_sec.Notes = '(Serafini 2008, PMID: 18593947)';
% rate of ArgI secretion
params.k_ArgI_sec.Value = 1.4e-2;
params.k_ArgI_sec.Units = '(mU*microliter)/cell/day';
params.k_ArgI_sec.Notes = '(Serafini 2008, PMID: 18593947)';
% half-maximal entinostat concentration for NO inhibition
params.IC50_inostat_NO.Value = 0.56; % <1 nM
params.IC50_inostat_NO.Units = 'nanomolarity';
params.IC50_inostat_NO.Notes = '(Choo 2010, PMID: 20421217)';
% half-maximal entinostat concentration for CCL2 inhibition
params.IC50_inostat_CCL2.Value = 1.2; % <2 nM
params.IC50_inostat_CCL2.Units = 'nanomolarity';
params.IC50_inostat_CCL2.Notes = '(Choo 2013, PMID: 24241152)';
% half-maximal entinostat concentration for Arg I inhibition
params.IC50_inostat_ArgI.Value = 500;
params.IC50_inostat_ArgI.Units = 'nanomolarity';
params.IC50_inostat_ArgI.Notes = '(Orillian 2017, PMID: 28698201)';

% rate of ArgI/NO-induced T cell death Park 2018, PMID: 29491381
% half-maximal ArgI concentration for CTL inhibition
params.ArgI_50_Teff.Value = 61.7;
params.ArgI_50_Teff.Units = 'mU';
params.ArgI_50_Teff.Notes = '(Serafini 2008, PMID: 18593947)';
% half-maximal NO concentration for CTL inhibition
params.NO_50_Teff.Value = 0.75;
params.NO_50_Teff.Units = 'nanomolarity';
params.NO_50_Teff.Notes = '(Serafini 2008, PMID: 18593947)';
% Half-maximal CCL2 concentration for macrophage/MDSC recruitment
params.CCL2_50.Value = 0.23;
params.CCL2_50.Units = 'nanomolarity';
params.CCL2_50.Notes = '(Opalek 2007, PMID: 17888174; Ernst 1994, PMID: 8144933)';
% half-maximal ArgI concentration for Treg expension
params.ArgI_50_Treg.Value = 22.1;
params.ArgI_50_Treg.Units = 'mU';
params.ArgI_50_Treg.Notes = '(estimated based on Serafini 2008, PMID: 18593947)';

%% T Helper Cell Module parameters
% T helper cell activation rate
params.k_Th_act.Value = 10;
params.k_Th_act.Units = '1/day';
params.k_Th_act.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Th differentiation rate to Treg
params.k_Th_Treg.Value = 0.022;
params.k_Th_Treg.Units = '1/day';
params.k_Th_Treg.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Secretion rate of TGFb by Treg
params.k_TGFb_Tsec.Value = 1.2e-10;
params.k_TGFb_Tsec.Units = 'nanomole/cell/day';
params.k_TGFb_Tsec.Notes = '(Liyanage 2002, PMID: 12193750)';
% Degradation rate of TGFb
params.k_TGFb_deg.Value = 14.3;
params.k_TGFb_deg.Units = '1/day';
params.k_TGFb_deg.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Half-Maximal TGFb level for differentiation of Th to Treg / chemoresistance development / M1-to-M2 polarization
params.TGFb_50.Value = 0.07;
params.TGFb_50.Units = 'nanomolarity';
params.TGFb_50.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Half-Maximal TGFb level for CD8+ T cell inhibition
params.TGFb_50_Teff.Value = 0.14;
params.TGFb_50_Teff.Units = 'nanomolarity';
params.TGFb_50_Teff.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Half-Maximal cancer cell number for T cell recruitment
params.Kc_rec.Value = 2.02e7;
params.Kc_rec.Units = 'cell^2';
params.Kc_rec.Notes = '(Phillis 2006, PMID: 16153659)';
% Baseline TGFb level in breast tumor
params.TGFbase.Value = 0.016;
params.TGFbase.Units = 'nanomolarity';
params.TGFbase.Notes = '(Ivanović 2003, PMID: 12751375; Manuyakorn 2008, PMID: 19317336; Panis 2013, PMID: 23393376)';
% IFNg secretion rate by T helper cell
params.k_IFNg_sec.Value = 5e-13;
params.k_IFNg_sec.Units = 'nanomole/cell/day';
params.k_IFNg_sec.Notes = '(Autenshlyus 2021, PMID: 34059727; Borj 2017, doi:10.15171/ijbms.2017.05)';
% IFNg degradation rate
params.k_IFNg_deg.Value = 7.68;
params.k_IFNg_deg.Units = '1/day';
params.k_IFNg_deg.Notes = '(Hofstra 1998, PMID: 9806748)';
% Half-Maximal IFNg level for PD-L1 induction
params.IFNg_50_ind.Value = 2.96;
params.IFNg_50_ind.Units = 'picomolarity';
params.IFNg_50_ind.Notes = '(Shin 2017, PMID: 27903500)';

%% Paclitaxel Module Parameters
% Max clearance rate from central compartment
params.Vmcl.Value = 8070;
params.Vmcl.Units = 'microgram/hour';
params.Vmcl.Notes = '(Chen 2014, PMID: 24719309)';
% Half-max conc. of clearance rate from central compartment
params.Kcl.Value = 40.2;
params.Kcl.Units = 'microgram/liter';
params.Kcl.Notes = '(Chen 2014, PMID: 24719309)';
% Transport rate b/t central and P compartment
params.Q2.Value = 41.6;
params.Q2.Units = 'liter/hour';
params.Q2.Notes = '(Chen 2014, PMID: 24719309)';
% Max clearance rate from central compartment
params.Vmt.Value = 325000;
params.Vmt.Units = 'microgram/hour';
params.Vmt.Notes = '(Chen 2014, PMID: 24719309)';
% Half-max conc. of clearance rate from central compartment
params.Kt.Value = 4260;
params.Kt.Units = 'microgram/liter';
params.Kt.Notes = '(Chen 2014, PMID: 24719309)';
% Body surface area
params.BSA.Value = 1.9;
params.BSA.Units = 'meter^2';
params.BSA.Notes = '(Chen 2014, PMID: 24719309)';

%% Paclitaxel PD Parameters
% Tumour to plasma concentration ratio of nab-paclitaxel
params.r_nabp.Value = 1.5;
params.r_nabp.Units = 'dimensionless';
params.r_nabp.Notes = '(Yang 2007, PMID: 17828616; Rajeshkumar 2016, PMID: 27441498)';
% Cancer cell killing rate by nab-paclitaxel
params.k_C_nabp.Value = 0.06;
params.k_C_nabp.Units = '1/hour';
params.k_C_nabp.Notes = '(Choi 2012, PMID: 23023313; Volk 2008, PMID: 18516298; Kozłowska 2018, PMID: 29769198)';
% Half-maximal nab-paclitaxel concentration for cancer cell killing
params.IC50_nabp.Value = 92;
params.IC50_nabp.Units = 'nanomolarity';
params.IC50_nabp.Notes = '(Yang 2013, PMID: 23180760)';
% Molecular weight of nab-paclitaxel
params.MW_nabp.Value = 853.9;
params.MW_nabp.Units = 'gram/mole';
params.MW_nabp.Notes = '(calculated)';
% Limited tumor penetration of nab-paclitaxel
params.Kc_nabp.Value = 8e7;
params.Kc_nabp.Units = 'cell';
params.Kc_nabp.Notes = '(estimated based on Kyle 2007, PMID: 17473214)';

% Rate of chemo-resistance development to nab-paclitaxel
params.k_C_resist.Value = 1e-4;
params.k_C_resist.Units = '1/day';
params.k_C_resist.Notes = '(estimated)';
% Number of folds increase of nab-paclitaxel IC50 in resistant cancer clones
params.r_resist.Value = 100;
params.r_resist.Units = 'dimensionless';
params.r_resist.Notes = '(Nemcova-Furstova 2016, PMID: 27664577)';

% Initial Tumor Capacity
params.K0.Value = 1.03e4;
params.K0.Units = 'cell';
params.K0.Notes = '(Desai 2006, PMID: 16489089)';
% Tumour vasculature growth rate
params.k_K_g.Value = 4.12;
params.k_K_g.Units = '1/day';
params.k_K_g.Notes = '(Desai 2006, PMID: 16489089)';
% Tumour vasculature inhibition rate
params.k_K_d.Value = 3.4e-3;
params.k_K_d.Units = '1/day';
params.k_K_d.Notes = '(Desai 2006, PMID: 16489089; Hahnfeldt 1999, PMID: 10519381)';
% Secretion rate of angiogenic factors induced by nab-paclitaxel
params.k_vas_nabp.Value = 2.8e-4;
params.k_vas_nabp.Units = 'picogram/cell/day';
params.k_vas_nabp.Notes = '(Volk 2008, PMID: 18516298)';
% Half-maximal conc. of nab-paclitaxel on angiogenic factor induction
params.IC50_nabp_vas.Value = 5.2;
params.IC50_nabp_vas.Units = 'nanomolarity';
params.IC50_nabp_vas.Notes = '(Volk 2008, PMID: 18516298)';
% Secretion rate of angiogenic factors by cancer cells
params.k_vas_Csec.Value = 1.1e-4;
params.k_vas_Csec.Units = 'picogram/cell/day';
params.k_vas_Csec.Notes = '(Volk 2008, PMID: 18516298)';
% Degradation rate of angiogenic factors
params.k_vas_deg.Value = 16.6;
params.k_vas_deg.Units = '1/day';
params.k_vas_deg.Notes = '(Finley 2011, PMID: 22104283)';
% Half-maximal conc. of angiogenic factor on tumor capacity growth
params.c_vas_50.Value = 1.07e3; 
params.c_vas_50.Units = 'picogram/milliliter';
params.c_vas_50.Notes = '(Desai 2006, PMID: 16489089)';
% Inhibition rate of maximal tumor capacity by nab-paclitaxel
params.k_K_nabp.Value = 17.8;
params.k_K_nabp.Units = '1/hour/molarity';
params.k_K_nabp.Notes = '(Mollard 2017, PMID: 28416742)';


%% Properties for CTLA4-mediated Treg depletion (ADCC)
% CTLA4-Ipi Concentration for Half-Maximal Treg death
params.Treg_CTLA4_50.Value = 1000;
params.Treg_CTLA4_50.Units = 'molecule';
params.Treg_CTLA4_50.Notes = '(estimated)';
% Hill Coefficient for CTLA4-Ipi
params.n_Treg_CTLA4.Value = 2;
params.n_Treg_CTLA4.Units = 'dimensionless';
params.n_Treg_CTLA4.Notes = '(estimated)';
% CTLA4 Expression on Tregs
params.Treg_CTLA4_tot.Value = 5000;
params.Treg_CTLA4_tot.Units = 'molecule';
params.Treg_CTLA4_tot.Notes = '(Jansson 2005, PMID: 16034096)';
% Ipi ADCC (antibody-dependent cellular cytotoxicity) rate of Treg
params.k_CTLA4_ADCC.Value = 0.1;
params.k_CTLA4_ADCC.Units = '1/day';
params.k_CTLA4_ADCC.Notes = 'Anti-CTLA4 ADCC (antibody-dependent cellular cytotoxicity) rate of Treg (Richards 2008, PMID: 18723496)';
% IgG1-FcRIIIa V158 kon
% params.kon_IgG1_FcR.Value =  8.2e3;
% params.kon_IgG1_FcR.Units = '1/(molarity*second)';
% params.kon_IgG1_FcR.Notes = '(Li 2007, PMID: 17202140)';
% IgG1-FcRIIIa V158 Kd
% params.kd_IgG1_FcR.Value = 0.28; % 0.28 - 0.5 uM
% params.kd_IgG1_FcR.Units = 'micromolarity';
% params.kd_IgG1_FcR.Notes = '(Bruhns 2009, PMID: 19018092; Richards 2008, PMID: 18723496)';
% RcRIIIa Expression on Tregs
% params.FcRIIIa_tot.Value = 2.4e5; % 2.51e6 - 3.55e6 - 3.98e6
% params.FcRIIIa_tot.Units = 'molecule';
% params.FcRIIIa_tot.Notes = '(Richards 2008, PMID: 18723496; Robinett 2018, PMID: 29960887)';


%% Macrophage Module
% Macrophage recruitment rate
params.k_Mac_mig.Value = 1.7e5;
params.k_Mac_mig.Units = 'cell/(milliliter*day)';
params.k_Mac_mig.Notes = '(Yang 2018, PMID: 30026826)';
% Death rate
params.k_Mac_death.Value = 0.02;
params.k_Mac_death.Units = '1/day';
params.k_Mac_death.Notes = '(Ginhoux 2016, PMID: 26982352)';
% Secretion rate of TGFb by M2 macrophage
params.k_TGFb_Msec.Value = 2e-11;
params.k_TGFb_Msec.Units = 'nanomole/cell/day';
params.k_TGFb_Msec.Notes = '(Matsuo 2000, PMID: 10652049)';
% Secretion rate of angiogenic factor by M2 macrophage
params.k_vas_Msec.Value = 1.7e-5;
params.k_vas_Msec.Units = 'picogram/cell/day';
params.k_vas_Msec.Notes = '(Wu 2010, PMID: 20692534)';
% Secretion rate of IL12 by mAPC
params.k_IL12_sec.Value = 8.5e-12;
params.k_IL12_sec.Units = 'nanomole/cell/day';
params.k_IL12_sec.Notes = '(Heufler 1996, PMID: 8605935)';
% Secretion rate of IL12 by M2 macrophage
params.k_IL12_Msec.Value = 5e-13;
params.k_IL12_Msec.Units = 'nanomole/cell/day';
params.k_IL12_Msec.Notes = '(Clough 2007, PMID: 17597818; Sica 2000, PMID: 10623821)';
% Degradation rate of IL12
params.k_IL12_deg.Value = 0.0231;
params.k_IL12_deg.Units = '1/hour';
params.k_IL12_deg.Notes = '(Carreño 2000, PMID: 10707873)';
% Secretion rate of IL10
params.k_IL10_sec.Value = 3e-12;
params.k_IL10_sec.Units = 'nanomole/cell/day';
params.k_IL10_sec.Notes = '(Autenshlyus 2021, PMID: 34059727)';
% Degradation rate of IL10
params.k_IL10_deg.Value = 4;
params.k_IL10_deg.Units = '1/day';
params.k_IL10_deg.Notes = '(Saxena 2015, PMID: 25481648)';
% Secretion rate of IL13
params.k_IL13_sec.Value = 5.35e-11;
params.k_IL13_sec.Units = 'nanomole/cell/day';
params.k_IL13_sec.Notes = '(Autenshlyus 2021, PMID: 34059727)';
% Degradation rate of IL13
params.k_IL13_deg.Value = 0.05; % 20mins; 72 day^-1
params.k_IL13_deg.Units = '1/minute';
params.k_IL13_deg.Notes = '(Khodoun 2007, PMID: 17982031)';
% M1 to M2 polarization rate
params.k_M2_pol.Value = 0.25;
params.k_M2_pol.Units = '1/day';
params.k_M2_pol.Notes = '(estimated assuming that polarization takes roughly 2-3 days)';
% M2 to M1 polarization rate
params.k_M1_pol.Value = 0.045;
params.k_M1_pol.Units = '1/day';
params.k_M1_pol.Notes = '(estimated based on the M1/M2 ratio)';
% Half-maximal IL10 level for M1 to M2 polarization / maintaining Treg function / mAPC inhibition (*STAT6 related)
params.IL10_50.Value = 8;
params.IL10_50.Units = 'picomolarity';
params.IL10_50.Notes = '(Altan-Bonnet 2019, PMID: 30770905)';
% Half-maximal IL13 level for M1 to M2 polarization
params.IL13_50.Value = 120;
params.IL13_50.Units = 'picomolarity';
params.IL13_50.Notes = '(Altan-Bonnet 2019, PMID: 30770905)';
% Half-maximal IL-12 level for M2 to M1 polarization
params.IL12_50.Value = 0.14; % << 2g/L = 2/7e4 pg/mL (2.86e-5 M)
params.IL12_50.Units = 'picomolarity'; % pM
params.IL12_50.Notes = '(Altan-Bonnet 2019, PMID: 30770905)';
% Half-maximal IFNg level for M2 to M1 polarization
params.IFNg_50.Value = 2.9;
params.IFNg_50.Units = 'picomolarity'; % pM
params.IFNg_50.Notes = '(Altan-Bonnet 2019, PMID: 30770905)';
% M1-mediated cancer killing rate
params.k_M1_phago.Value = 0.33;
params.k_M1_phago.Units = '1/day';
params.k_M1_phago.Notes = '(Pan 2020, PMID: 33365025)';
% Half-maximal IL-10 level for phagocytosis
params.IL10_50_phago.Value = 270;
params.IL10_50_phago.Units = 'picomolarity';
params.IL10_50_phago.Notes = '(Bian 2016, PMID: 27578867)';
% Dependence of phagocytosis rate on M1/C ratio
params.K_Mac_C.Value = 2;
params.K_Mac_C.Units = 'dimensionless';
params.K_Mac_C.Notes = '(estimated)'; % Gardner 2020, PMID: 33209523
% Dependence of Teff killing rate on Teff/C ratio
params.K_T_C.Value = 1.2;
params.K_T_C.Units = 'dimensionless';
params.K_T_C.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Dependence of Teff killing rate on Teff/Treg ratio
params.K_T_Treg.Value = 11;
params.K_T_Treg.Units = 'dimensionless';
params.K_T_Treg.Notes = '(Robertson-Tessi 2012, PMID: 22051568)';
% Macrophage Diameter
params.D_Mcell.Value = 21;
params.D_Mcell.Units = 'micrometer';
params.D_Mcell.Notes = '(Krombach 1997, PMID: 9400735)';
% Surface Area of a Macrophages
params.A_Mcell.Value = [];
params.A_Mcell.Units = 'micrometer^2';
params.A_Mcell.Notes = ['calculated based on the average macrophage diameter ' params.D_Mcell.Notes];
params.A_Mcell.Factors = ["D_Mcell"];
params.A_Mcell.Equation = '4*pi*(p(1)/2)^2';
% Volume of a Macrophage
params.vol_Mcell.Value = [];
params.vol_Mcell.Units = 'micrometer^3/cell';
params.vol_Mcell.Notes = ['Volume of a macrophage calculated based on the average macrophage diameter ' params.D_Mcell.Notes];
params.vol_Mcell.Factors = ["D_Mcell", "cell"];
params.vol_Mcell.Equation = '4/3*pi*(p(1)/2)^3/p(2)';


%% Phagocytosis Module
% Binding affinity of CD47-SIRPa
params.kd_CD47_SIRPa.Value = 1.6;
params.kd_CD47_SIRPa.Units = 'micromolarity';
params.kd_CD47_SIRPa.Notes = '(Hayes 2020, PMID: 31964705)';
% Dissociation rate of CD47-SIRPa
params.koff_CD47_SIRPa.Value = 0.3;
params.koff_CD47_SIRPa.Units = '1/minute';
params.koff_CD47_SIRPa.Notes = '(Hayes 2020, PMID: 31964705)';
% Association rate of CD47-SIRPa
params.kon_CD47_SIRPa.Value = [];
params.kon_CD47_SIRPa.Units = '';
params.kon_CD47_SIRPa.Notes = '(Hayes 2020, PMID: 31964705)';
params.kon_CD47_SIRPa.Factors = ["koff_CD47_SIRPa","kd_CD47_SIRPa","d_syn"];
params.kon_CD47_SIRPa.Equation = 'p(1)/p(2)/p(3)';

% Binding affinity of CD47-aCD47
params.kd_CD47_aCD47.Value = 1; % 10-1000 nM PMID: 30133535
params.kd_CD47_aCD47.Units = 'nanomolarity';
params.kd_CD47_aCD47.Notes = '(Jalil 2020, PMID: 32421049)';
% Dissociation rate of CD47-aCD47
params.koff_CD47_aCD47.Value = 1/20; % 1e-4 - 1e-3 1/s
params.koff_CD47_aCD47.Units = '1/minute';
params.koff_CD47_aCD47.Notes = '(Kauder 2018, PMID: 30133535)';
% Association rate of CD47-aCD47
params.kon_CD47_aCD47.Value = [];
params.kon_CD47_aCD47.Units = '';
params.kon_CD47_aCD47.Notes = '(Kauder 2018, PMID: 30133535)';
params.kon_CD47_aCD47.Factors = ["koff_CD47_aCD47","kd_CD47_aCD47"];
params.kon_CD47_aCD47.Equation = 'p(1)/p(2)';

% CD47-aCD47 Chi (antibody cross-arm binding strength)
params.Chi_CD47_aCD47_3D.Value = 100;
params.Chi_CD47_aCD47_3D.Units = 'dimensionless';
params.Chi_CD47_aCD47_3D.Notes = '(estimated)';
% CD47-aCD47 Chi for 2D
params.Chi_CD47_aCD47.Value = [];
params.Chi_CD47_aCD47.Units = '';
params.Chi_CD47_aCD47.Notes = '(estimated)';
params.Chi_CD47_aCD47.Factors = ["Chi_CD47_aCD47_3D","d_syn","N_avg"];
params.Chi_CD47_aCD47.Equation = 'p(1)/(p(2)*p(3))';
% Half-maximal SIRPa binding for phagocytosis inhibition
params.SIRPa_50.Value = 37;
params.SIRPa_50.Units = 'molecule/micrometer^2';
params.SIRPa_50.Notes = '(Willingham 2012, PMID: 22451913)';
% Hill coefficient for H_SIRPa
params.n_SIRPa.Value = 2;
params.n_SIRPa.Units = 'dimensionless';
params.n_SIRPa.Notes = '(estimated)';
% CD47 expression on cancer cell
params.C_CD47.Value = 400;
params.C_CD47.Units = 'molecule/micrometer^2';
params.C_CD47.Notes = '(Morrissey 2020, PMID: 32768386)';
% SIRPa expression on macrophage
params.M_SIRPa.Value = 100;
params.M_SIRPa.Units = 'molecule/micrometer^2';
params.M_SIRPa.Notes = '(Subramanian 2006, PMID: 16291597)';
% PD-1 expression on macrophage
params.M_PD1.Value = 1.25e4;
params.M_PD1.Units = 'molecule';
params.M_PD1.Notes = 'estimated based on PD-1 expression on T cell (Cheng 2014, PMID: 23417675; Mkrtichyan 2012, PMID: 22837483) and adjusted to fit its inhibitory effect on phagocytosis (Gordon 2017, PMID: 28514441)';
  
