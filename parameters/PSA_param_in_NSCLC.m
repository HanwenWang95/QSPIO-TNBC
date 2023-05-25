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

function params = PSA_param_in_NSCLC_VPop

params.names = {};

% Cancer cell growth rate
params.names = [params.names; 'k_C1_growth'];
params.k_C1_growth.Median = log(.007);
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

% Treg inhibition of Teff
params.names = [params.names; 'k_Treg'];
params.k_Treg.UpperBound = 1;
params.k_Treg.LowerBound = 0.01;
params.k_Treg.Sampling   = 'loguniform';
params.k_Treg.ScreenName = 'Rate of Teff inhibition by Treg';

% Cellular clearance rate from tumor
params.names = [params.names; 'k_cell_clear'];
params.k_cell_clear.UpperBound = 1;
params.k_cell_clear.LowerBound = 0.01;
params.k_cell_clear.Sampling   = 'loguniform';
params.k_cell_clear.ScreenName = 'Rate of apoptotic cell clearance';

% T cell killing of cancer cell
params.names = [params.names; 'k_C_T1'];
params.k_C_T1.Median = log(0.9);
params.k_C_T1.Sigma = 1;
params.k_C_T1.Sampling   = 'lognormal';
params.k_C_T1.ScreenName = 'Rate of cancer cell death by Teffs';

% Kd of antigen P1
params.names = [params.names; 'k_P1_d1']; % estimated by Jafarnejad et al. based on Forde 2018, PMID: 30157404
params.k_P1_d1.Median = log(40e-9);
params.k_P1_d1.Sigma = 1;
params.k_P1_d1.Sampling   = 'lognormal';
params.k_P1_d1.ScreenName = 'Kd of tumor neoantigens';

% Tumor mutational Burden
params.names = [params.names; 'n_T1_clones']; % estimated by Jafarnejad et al. based on Forde 2018, PMID: 30157404
params.n_T1_clones.Median = log(92);
params.n_T1_clones.Sigma = 1;
params.n_T1_clones.Sampling   = 'lognormal';
params.n_T1_clones.ScreenName = 'Tumor-specific T cell clone';

params.names = [params.names; 'n_T0_clones'];
params.n_T0_clones.Median = log(100);
params.n_T0_clones.Sigma = 1;
params.n_T0_clones.Sampling   = 'lognormal';
params.n_T0_clones.ScreenName = 'Self-antigen-specific T cell clone';

% Number of T Cell Divisions during Activation by IL-2
params.names = [params.names; 'N_IL2_CD8'];
params.N_IL2_CD8.LowerBound = 9;
params.N_IL2_CD8.UpperBound = 13;
params.N_IL2_CD8.Sampling   = 'uniform';
params.N_IL2_CD8.ScreenName = '# of CD8 T cell division by IL-2';

params.names = [params.names; 'N_IL2_CD4'];
params.N_IL2_CD4.LowerBound = 5;
params.N_IL2_CD4.UpperBound = 11;
params.N_IL2_CD4.Sampling   = 'uniform';
params.N_IL2_CD4.ScreenName = '# of CD4 T cell division by IL-2';

% Activated CD8/4 T Cell Transport to Tumor
params.names = [params.names; 'q_T0_T_in'];
params.q_T0_T_in.LowerBound = 0.000145/2;
params.q_T0_T_in.UpperBound = 0.000145*2;
params.q_T0_T_in.Sampling   = 'uniform';
params.q_T0_T_in.ScreenName = 'Rate of CD4 T cell transport to tumor';

params.names = [params.names; 'q_T1_T_in'];
params.q_T1_T_in.LowerBound = 0.000145/2;
params.q_T1_T_in.UpperBound = 0.000145*2;
params.q_T1_T_in.Sampling   = 'uniform';
params.q_T1_T_in.ScreenName = 'Rate of CD8 T cell transport to tumor';

% Initial tumor diameter
params.names = [params.names; 'initial_tumour_diameter'];
params.initial_tumour_diameter.Median = log(3.7);
params.initial_tumour_diameter.Sigma = .4;
params.initial_tumour_diameter.Sampling   = 'lognormal';
params.initial_tumour_diameter.ScreenName = 'Initial tumor diameter';

% Checkpoint Module
% Average PD1 Expression on T cells
params.names = [params.names; 'T_PD1_total'];
params.T_PD1_total.UpperBound = 3.1e3*20;
params.T_PD1_total.LowerBound = 1.2e4;
params.T_PD1_total.Sampling   = 'loguniform';
params.T_PD1_total.ScreenName = 'Average number of PD1 per T cell';
% Average Baseline PDL1 Expression on Tumor/Immune cells in tumor
% params.names = [params.names; 'C1_PDL1_base'];
% params.C1_PDL1_base.UpperBound = 8e4*20/6;
% params.C1_PDL1_base.LowerBound = 5e4;
% params.C1_PDL1_base.Sampling   = 'loguniform';
% params.C1_PDL1_base.ScreenName = 'Baseline number of PDL1 per cell in the tumor';
% Average Baseline PDL1 Expression on mAPCs in TDLN
% params.names = [params.names; 'APC_PDL1_base'];
% params.APC_PDL1_base.UpperBound = 8e4*20/6;
% params.APC_PDL1_base.LowerBound = 5e4;
% params.APC_PDL1_base.Sampling   = 'loguniform';
% params.APC_PDL1_base.ScreenName = 'Baseline number of PDL1 per mAPC';
% Effective Density of Bound PD-1 on Teff and Phagocytosis Inhibition
params.names = [params.names; 'PD1_50'];
params.PD1_50.Median = log(6);
params.PD1_50.Sigma = 1;
params.PD1_50.Sampling   = 'lognormal';
params.PD1_50.ScreenName = 'Half-maximal PD1-PDL1 for Teff inhibition';

% Ratio between PDL1 and PDL2 expression in tumor
params.names = [params.names; 'r_PDL2C1'];
params.r_PDL2C1.UpperBound = 0.1;
params.r_PDL2C1.LowerBound = 0.001;
params.r_PDL2C1.Sampling   = 'loguniform';
params.r_PDL2C1.ScreenName = 'PDL2/PDL1 ratio on cancer cell';
% Ratio between PDL1 and PDL2 expression on APC
params.names = [params.names; 'r_PDL2APC'];
params.r_PDL2APC.UpperBound = 0.1;
params.r_PDL2APC.LowerBound = 0.001;
params.r_PDL2APC.Sampling   = 'loguniform';
params.r_PDL2APC.ScreenName = 'PDL2/PDL1 ratio on APCs';
% Maximal rate of Th-to-Treg transdifferentiation
params.names = [params.names; 'k_Th_Treg'];
params.k_Th_Treg.Median = log(0.022);
params.k_Th_Treg.Sigma = 1;
params.k_Th_Treg.Sampling   = 'lognormal';
params.k_Th_Treg.ScreenName = 'Rate of Th to Treg differentiation';

%% Macrophages/MDSC Module
% MDSC Recruitment
params.names = [params.names; 'k_MDSC_rec'];
params.k_MDSC_rec.Median = log(2e3);
params.k_MDSC_rec.Sigma = 0.6;
params.k_MDSC_rec.Sampling   = 'lognormal';
params.k_MDSC_rec.ScreenName = 'Recruitment rate of MDSC ';
% Macrophage Recruitment
params.names = [params.names; 'k_Mac_rec'];
params.k_Mac_rec.Median = log(2e5);
params.k_Mac_rec.Sigma = 0.6;
params.k_Mac_rec.Sampling   = 'lognormal';
params.k_Mac_rec.ScreenName = 'Recruitment rate of macrophage ';
% Macrophage Polarization
params.names = [params.names; 'k_M1_pol'];
params.k_M1_pol.Median = log(0.02);
params.k_M1_pol.Sigma = 1;
params.k_M1_pol.Sampling   = 'lognormal';
params.k_M1_pol.ScreenName = 'Rate of M2 to M1 polarization ';
% Macrophage Polarization
params.names = [params.names; 'k_M2_pol'];
params.k_M2_pol.Median = log(0.25);
params.k_M2_pol.Sigma = 1;
params.k_M2_pol.Sampling   = 'lognormal';
params.k_M2_pol.ScreenName = 'Rate of M1 to M2 polarization ';
% CCL2 Secretion
params.names = [params.names; 'k_CCL2_sec'];
params.k_CCL2_sec.Median = log(1.7e-12);
params.k_CCL2_sec.Sigma = 0.6;
params.k_CCL2_sec.Sampling   = 'lognormal';
params.k_CCL2_sec.ScreenName = 'Secretion rate of CCL2 ';
% CD47 Expression
params.names = [params.names; 'C_CD47']; % Morrissey 2020, PMID: 32768386; Zhang 2019, PMID: 31829270
params.C_CD47.UpperBound = 150;
params.C_CD47.LowerBound = 50;
params.C_CD47.Sampling = 'uniform';
params.C_CD47.ScreenName = 'Mean CD47 density on cancer cell ';
% SIRPa Expression
params.names = [params.names; 'M_SIRPa'];
params.M_SIRPa.UpperBound = 180;
params.M_SIRPa.LowerBound = 20;
params.M_SIRPa.Sampling   = 'uniform';
params.M_SIRPa.ScreenName = 'Mean SIRPa density on macrophage ';
% TAM PD1 Expression
params.names = [params.names; 'M_PD1_total'];
params.M_PD1_total.UpperBound = 3.1e3*20;
params.M_PD1_total.LowerBound = 1.5e3;
params.M_PD1_total.Sampling   = 'loguniform';
params.M_PD1_total.ScreenName = 'Mean PD1 density on macrophage ';
% Maximal rate macrophage-mediated phagocytosis
params.names = [params.names; 'k_M1_phago'];
params.k_M1_phago.Median = log(.33);
params.k_M1_phago.Sigma = 1;
params.k_M1_phago.Sampling   = 'lognormal';
params.k_M1_phago.ScreenName = 'Maximal rate of macrophage-mediated phagocytosis ';

params.names = [params.names; 'SIRPa_50'];
params.SIRPa_50.Median = log(37);
params.SIRPa_50.Sigma = .5;
params.SIRPa_50.Sampling   = 'lognormal';
params.SIRPa_50.ScreenName = 'Half-maximal SIRPa-CD47 for phagocytosis inhibition';

params.names = [params.names; 'K_Mac_C'];
params.K_Mac_C.Median = log(2);
params.K_Mac_C.Sigma = 1;
params.K_Mac_C.Sampling   = 'lognormal';
params.K_Mac_C.ScreenName = 'Dependence of phagocytosis on M1/C ratio';

end
