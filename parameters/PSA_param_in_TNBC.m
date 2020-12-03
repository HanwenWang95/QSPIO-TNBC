% Function to generate object with model parameters to include in parameter
% sensitivity analysis
%
% Output: params -- object containing parameters
%                   -> for each parameter:
%                       - Adds the name of the parameter to the list
%                       - defines upper and lower bounds for uniform and
%                       loguniform
%                       - defines Median and Sigma for normal and lognormal
%                       - specifies the Sampling technique choose from:
%                           - uniform
%                           - loguniform
%                           - normal
%                           - lognormal
%         Examples:
%             % k1
%             params.names = [params.names; 'k1'];
%             params.k1.UpperBound = 1;
%             params.k1.LowerBound = 0;
%             params.k1.Sampling   = 'uniform';
%             params.k1.ScreenName = 'k1 binding rate';
%
%             % k2
%             params.names = [params.names; 'k2'];
%             params.k2.UpperBound = 1;
%             params.k2.LowerBound = 0;
%             params.k2.Sampling   = 'loguniform';
%             params.k2.ScreenName = 'k2 binding rate';
%
%             % k3
%             params.names = [params.names; 'k3'];
%             params.k3.Median     = 1;
%             params.k3.Sigma      = 1;
%             params.k3.Sampling   = 'normal';
%             params.k3.ScreenName = 'k3 binding rate';
%
%             % k4
%             params.names = [params.names; 'k4'];
%             params.k4.Median     = 1;
%             params.k4.Sigma      = 1;
%             params.k4.Sampling   = 'lognormal';
%             params.k4.ScreenName = 'k4 binding rate';

function params = PSA_param_in_TNBC

params.names = {};

% C1 growth rate 0.0035 for HER2-negative BC: 197.56 +/- 30 days
% TNBC: 103 +/- 43 days (46-205 days; median, 93 days); ER+: 241 +/- 166 days
% params.names = [params.names; 'k_C1_growth'];
% params.k_C1_growth.Median = log(0.00745); % 0.015
% params.k_C1_growth.Sigma = 1; % 0.00336
% params.k_C1_growth.Sampling   = 'lognormal';
% params.k_C1_growth.ScreenName = 'Cancer Growth Rate';

% C1 growth rate
params.names = [params.names; 'k_C1_growth'];
params.k_C1_growth.Median = log(0.0087);
params.k_C1_growth.Sigma = 1;
params.k_C1_growth.Sampling   = 'lognormal';
params.k_C1_growth.ScreenName = 'Rate of tumor growth';

% C1 basal death rate
params.names = [params.names; 'k_C1_death'];
params.k_C1_death.UpperBound = 0.001;
params.k_C1_death.LowerBound = 0.00001;
params.k_C1_death.Sampling   = 'loguniform';
params.k_C1_death.ScreenName = 'Rate of apoptotic cancer cell death';

% T cell exhaustion by cancer cell
params.names = [params.names; 'k_T1'];
params.k_T1.UpperBound = 1;
params.k_T1.LowerBound = 0.01;
params.k_T1.Sampling   = 'loguniform';
params.k_T1.ScreenName = 'Rate of T cell exhaustion by cancer cells';

% T cell killing of cancer cell
params.names = [params.names; 'k_C_T1'];
params.k_C_T1.Median = log(0.9);
params.k_C_T1.Sigma = 1;
params.k_C_T1.Sampling   = 'lognormal';
params.k_C_T1.ScreenName = 'Rate of cancer cell death by Teffs';

% Treg inhibition of Teff
params.names = [params.names; 'k_Treg'];
params.k_Treg.UpperBound = 1;
params.k_Treg.LowerBound = 0.01;
params.k_Treg.Sampling   = 'loguniform';
params.k_Treg.ScreenName = 'Rate of Teff inhibition by Treg';

% Kd of antigen P1
params.names = [params.names; 'k_P1_d1'];
params.k_P1_d1.Median = log(27e-9);
params.k_P1_d1.Sigma = 1;
params.k_P1_d1.Sampling   = 'lognormal';
params.k_P1_d1.ScreenName = 'Kd of tumor neoantigens';

% Kd of antigen P0
% params.names = [params.names; 'k_P0_d1'];
% params.k_P0_d1.Median = log(100e-9);
% params.k_P0_d1.Sigma = 1;
% params.k_P0_d1.Sampling   = 'lognormal';
% params.k_P0_d1.ScreenName = 'Kd of self antigen';

% Tumor mutational Burden
% PMID: 27240091, PMID: 30832597
params.names = [params.names; 'n_T1_clones'];
params.n_T1_clones.Median = log(63);
params.n_T1_clones.Sigma = 0.7;
params.n_T1_clones.Sampling   = 'lognormal';
params.n_T1_clones.ScreenName = 'Tumor-specific T cell clone';

params.names = [params.names; 'n_T0_clones'];
params.n_T0_clones.Median = log(63);
params.n_T0_clones.Sigma = 0.7;
params.n_T0_clones.Sampling   = 'lognormal';
params.n_T0_clones.ScreenName = 'Self-antigen-specific T cell clone';

% Initial tumor diameter
params.names = [params.names; 'initial_tumour_diameter'];
params.initial_tumour_diameter.Median = log(2.5);
params.initial_tumour_diameter.Sigma = 0.3;
params.initial_tumour_diameter.Sampling   = 'lognormal';
params.initial_tumour_diameter.ScreenName = 'Initial tumor diameter';

% MDSC Module
params.names = [params.names; 'MDSC_max'];
params.MDSC_max.Median = log(1.637e5); % 5e4-5e6
params.MDSC_max.Sigma = 1;
params.MDSC_max.Sampling   = 'lognormal';
params.MDSC_max.ScreenName = 'Steady-state MDSC density in tumor';

% params.names = [params.names; 'T_PD1_total'];
% params.T_PD1_total.UpperBound = 3096*20;
% params.T_PD1_total.LowerBound = 3096;
% params.T_PD1_total.Sampling   = 'loguniform';
% params.T_PD1_total.ScreenName = 'PD1 expression on T cells ';

%params.names = [params.names; 'T_PDL1_total'];
%params.T_PDL1_total.UpperBound = 9.3e3*20; 
%params.T_PDL1_total.LowerBound = 9.3e3;
%params.T_PDL1_total.Sampling   = 'loguniform';
%params.T_PDL1_total.ScreenName = 'PD1 expression on T cells ';

params.names = [params.names; 'C1_PDL1_base'];
params.C1_PDL1_base.UpperBound = 5.4e4*20/6;
params.C1_PDL1_base.LowerBound = 9e3;
params.C1_PDL1_base.Sampling   = 'loguniform';
params.C1_PDL1_base.ScreenName = 'Baseline number of PDL1 per cell in the tumor';

% params.names = [params.names; 'r_PDL2C1'];
% params.r_PDL2C1.UpperBound = 0.07;
% params.r_PDL2C1.LowerBound = 1e-4;
% params.r_PDL2C1.Sampling   = 'loguniform';
% params.r_PDL2C1.ScreenName = 'PDL2/PDL1 ratio on cancer cell';

params.names = [params.names; 'APC_PDL1_base'];
params.APC_PDL1_base.UpperBound = 8e4*20/6;
params.APC_PDL1_base.LowerBound = 1.3e4;
params.APC_PDL1_base.Sampling   = 'loguniform';
params.APC_PDL1_base.ScreenName = 'Baseline number of PDL1 per mAPC';

% params.names = [params.names; 'r_PDL2APC'];
% params.r_PDL2APC.UpperBound = 0.07;
% params.r_PDL2APC.LowerBound = 1e-4;
% params.r_PDL2APC.Sampling   = 'loguniform';
% params.r_PDL2APC.ScreenName = 'PDL2/PDL1 ratio on APCs';

% params.names = [params.names; 'PD1_50'];
% params.PD1_50.Median = log(50);
% params.PD1_50.Sigma = 1;
% params.PD1_50.Sampling   = 'lognormal';
% params.PD1_50.ScreenName = 'PD1 50';

params.names = [params.names; 'k_reg'];
params.k_reg.Median = log(0.022);
params.k_reg.Sigma = 1;
params.k_reg.Sampling   = 'lognormal';
params.k_reg.ScreenName = 'Rate of Th to Treg differentiation';

%params.names = [params.names; 'IC50_vas'];
%params.IC50_vas.Median = log(150000);
%params.IC50_vas.Sigma = 1;
%params.IC50_vas.Sampling   = 'lognormal';
%params.IC50_vas.ScreenName = 'Angiogenic factor level for half-maximal tumor vasculature growth ';

params.names = [params.names; 'k_K_g'];
params.k_K_g.UpperBound = 6.9; % 7.5
params.k_K_g.LowerBound = 2.9; % 4
params.k_K_g.Sampling   = 'uniform';
params.k_K_g.ScreenName = 'Rate of tumor vasculature growth';

% nabPaclitaxel PK Parameters
params.names = [params.names; 'Vmcl'];
params.Vmcl.UpperBound = 9836;
params.Vmcl.LowerBound = 6500;
params.Vmcl.Sampling   = 'uniform';
params.Vmcl.ScreenName = 'Max clearance rate from V1 compartment ';

params.names = [params.names; 'Kcl'];
params.Kcl.UpperBound = 58.9;
params.Kcl.LowerBound = 24.9;
params.Kcl.Sampling   = 'uniform';
params.Kcl.ScreenName = 'Half-max conc. of Nab-P for V1 clearance ';

params.names = [params.names; 'Vmt'];
params.Vmt.UpperBound = 540445;
params.Vmt.LowerBound = 190694;
params.Vmt.Sampling   = 'uniform';
params.Vmt.ScreenName = 'Max clearance rate from V1 to V2 distribution ';

params.names = [params.names; 'Kt'];
params.Kt.UpperBound = 7910;
params.Kt.LowerBound = 2210;
params.Kt.Sampling   = 'uniform';
params.Kt.ScreenName = 'Half-max conc. of Nab-P for V1 to V2 distribution ';

params.names = [params.names; 'BSA'];
params.BSA.UpperBound = 2.4;
params.BSA.LowerBound = 1.3;
params.BSA.Sampling   = 'uniform';
params.BSA.ScreenName = 'Body surface area ';

params.names = [params.names; 'V1'];
params.V1.UpperBound = 17.85;
params.V1.LowerBound = 13.71;
params.V1.Sampling   = 'uniform';
params.V1.ScreenName = 'Peripheral compartment (V1) for Nab-P PK';

params.names = [params.names; 'V2'];
params.V2.UpperBound = 1935;
params.V2.LowerBound = 1396;
params.V2.Sampling   = 'uniform';
params.V2.ScreenName = 'Peripheral compartment (V2) for Nab-P PK ';

params.names = [params.names; 'V3'];
params.V3.UpperBound = 99.1;
params.V3.LowerBound = 59.8;
params.V3.Sampling   = 'uniform';
params.V3.ScreenName = 'Peripheral compartment (V3) for Nab-P PK ';

% PD parameters
params.names = [params.names; 'r_nabp'];
params.r_nabp.UpperBound = 2;
params.r_nabp.LowerBound = 1;
params.r_nabp.Sampling   = 'uniform';
params.r_nabp.ScreenName = 'Tumor to plasma conc. ratio of Nab-P ';

params.names = [params.names; 'IC50_nabp'];
params.IC50_nabp.Median = log(4.7e-8); % 4.7e-8
params.IC50_nabp.Sigma = 1.1; % 1.1
params.IC50_nabp.Sampling   = 'lognormal';
params.IC50_nabp.ScreenName = 'Half-max conc. of Nab-P for cancer killing ';

params.names = [params.names; 'k_C_resist'];
params.k_C_resist.Median = log(1e-4);
params.k_C_resist.Sigma = 1;
params.k_C_resist.Sampling   = 'lognormal';
params.k_C_resist.ScreenName = 'Rate of chemo-resistance development ';

params.names = [params.names; 'k_c_nabp'];
params.k_c_nabp.Median = log(0.017);
params.k_c_nabp.Sigma = 1;
params.k_c_nabp.Sampling   = 'lognormal';
params.k_c_nabp.ScreenName = 'Rate of angiogenic factor induction by Nab-P ';

%% Additional Parameters

% CTLA-4
%params.names = [params.names; 'CD28_CD8X_50'];
%params.CD28_CD8X_50.UpperBound = 200*5;
%params.CD28_CD8X_50.LowerBound = 200/2;
%params.CD28_CD8X_50.Sampling   = 'loguniform';
%params.CD28_CD8X_50.ScreenName = 'CD28-CD80/86 conc. for half-maximal T activation';

% MDSC
% params.names = [params.names; 'EC50_ArgI_Treg'];
% params.EC50_ArgI_Treg.UpperBound = 22.1*10;
% params.EC50_ArgI_Treg.LowerBound = 22.1/10;
% params.EC50_ArgI_Treg.Sampling   = 'loguniform';
% params.EC50_ArgI_Treg.ScreenName = 'ArgI half-maximal conc. on Treg expansion';

% params.names = [params.names; 'IC50_ArgI_CTL'];
% params.IC50_ArgI_CTL.UpperBound = 61.7*10;
% params.IC50_ArgI_CTL.LowerBound = 61.7/10;
% params.IC50_ArgI_CTL.Sampling   = 'loguniform';
% params.IC50_ArgI_CTL.ScreenName = 'ArgI half-maximal conc. on Teff inhibition';

% params.names = [params.names; 'IC50_NO_CTL'];
% params.IC50_NO_CTL.UpperBound = .75e-9*10;
% params.IC50_NO_CTL.LowerBound = .75e-9/10;
% params.IC50_NO_CTL.Sampling   = 'loguniform';
% params.IC50_NO_CTL.ScreenName = 'NO half-maximal conc. on Teff inhibition';

%params.names = [params.names; 'k_sec_NO'];
%params.k_sec_NO.UpperBound = 5.65e-7;
%params.k_sec_NO.LowerBound = 3.95e-7;
%params.k_sec_NO.Sampling   = 'uniform';
%params.k_sec_NO.ScreenName = 'Secretion rate of NO by MDSC';

%params.names = [params.names; 'k_sec_ArgI'];
%params.k_sec_ArgI.UpperBound = 1.6e-2;
%params.k_sec_ArgI.LowerBound = 1.2e-2;
%params.k_sec_ArgI.Sampling   = 'uniform';
%params.k_sec_ArgI.ScreenName = 'Secretion rate of Arg I by MDSC';

%params.names = [params.names; 'k_sec_CCL2'];
%params.k_sec_CCL2.UpperBound = 20.2e-11;
%params.k_sec_CCL2.LowerBound = 8.2e-11;
%params.k_sec_CCL2.Sampling   = 'uniform';
%params.k_sec_CCL2.ScreenName = 'Secretion rate of CCL2 by MDSC';

%params.names = [params.names; 'r_resist'];
%params.r_resist.UpperBound = 1.1;
%params.r_resist.LowerBound = .9;
%params.r_resist.Sampling   = 'uniform';
%params.r_resist.ScreenName = 'Ratio of proliferation rate ';

% params.names = [params.names; 'MTT'];
% params.MTT.UpperBound = 125;
% params.MTT.LowerBound = 109;
% params.MTT.Sampling   = 'uniform';
% params.MTT.ScreenName = 'Mean transit time of neutrophils ';
%
% params.names = [params.names; 'slope'];
% params.slope.UpperBound = 0.00290;
% params.slope.LowerBound = 0.00216;
% params.slope.Sampling   = 'uniform';
% params.slope.ScreenName = 'Nab-paclitaxel effect slope ';
%
% params.names = [params.names; 'age'];
% params.age.Prob = 0.255;
% params.age.Scale = 1;
% params.age.Sampling   = 'binary';
% params.age.ScreenName = 'Age > 65 for nab-paclitaxel effect slope ';
%
% params.names = [params.names; 'age_factor'];
% params.age_factor.UpperBound = 0.830;
% params.age_factor.LowerBound = 0.172;
% params.age_factor.Sampling   = 'uniform';
% params.age_factor.ScreenName = 'Age factor for nab-paclitaxel effect slope ';
%
% params.names = [params.names; 'Neu_base'];
% params.Neu_base.UpperBound = 4.62e9;
% params.Neu_base.LowerBound = 3.94e9;
% params.Neu_base.Sampling   = 'uniform';
% params.Neu_base.ScreenName = 'Baseline circulating neutrophil level ';
%
% params.names = [params.names; 'gamma'];
% params.gamma.UpperBound = 0.203;
% params.gamma.LowerBound = 0.171;
% params.gamma.Sampling   = 'uniform';
% params.gamma.ScreenName = 'Feedback parameter for nab-paclitaxel pd ';
%
%params.names = [params.names; 'IC50_nabp_endo'];
%params.IC50_nabp_endo.Median = log(51e-12);
%params.IC50_nabp_endo.Sigma = 1.1;
%params.IC50_nabp_endo.Sampling   = 'lognormal';
%params.IC50_nabp_endo.ScreenName = 'Half-maximal nab-paclitaxel concentration for endothelial cell killing ';
