% Function to generate OX40 module parameters from default parameters
%
% Inputs: params_in  -- object containing the physical parameters
%                       - kd_OX40l--kd of OX40-OX40l binding
%                       - kd_aOX40--kd of OX40-aOX40 binding
%                       - OX40--number of OX40 molecules on T cells
%                       - OX40l--number of OX40l molecules on cancer cells
%                       (apc?)
%
% Output: params_out -- object containing parameters
%                       - kd_OX40l--kd of OX40-OX40l binding
%                       - kd_aOX40--kd of OX40-aOX40 binding
%                       - OX40--number of OX40 molecules on T cells
%                       - OX40l--number of OX40l molecules on cancer cells
%
% Created: Jun 14, 2019 (Mohammad Jafarnejad)
% Last Modified: Jun 14, 2019 (MJ)

function params_out = OX40_parameters(params_in)

% kon Values
params_out.kon_OX40_aOX40   = params_in.kon_OX40_aOX40 ;
params_out.kon_OX40_aOX40.Notes   = ['kon for OX40-aOX40 ' params_in.kon_OX40_aOX40.Notes];
% koff Values
params_out.koff_OX40_aOX40   = params_in.kon_OX40_aOX40 * params_in.kd_OX40_aOX40;
params_out.koff_OX40_aOX40.Notes   = ['calculated based on the measured kd and kon ' params_in.kd_OX40_aOX40.Notes];
% Bivalent antibody - antibody cross-arm binding efficiency 
params_out.Chi_OX40_aOX40 = params_in.Chi_OX40_aOX40/(params_in.d_syn*params_in.N_avg);
params_out.Chi_OX40_aOX40.Notes = [' that also includes the conversion of kon from 3D to 2D' params_in.Chi_OX40_aOX40.Notes];

% kon Values
params_out.kon_OX40_OX40l   = params_in.kon_OX40_OX40l / params_in.d_syn ;
params_out.kon_OX40_OX40l.Notes   = ['2D kon for OX40-OX40l ' params_in.kon_OX40_OX40l.Notes];
% koff Values
params_out.koff_OX40_OX40l   = params_in.kon_OX40_OX40l * params_in.kd_OX40_OX40l;
params_out.koff_OX40_OX40l.Notes   = ['calculated based on the measured kd and kon ' params_in.kd_OX40_OX40l.Notes];

% Surface Area of naive T cell
params_out.A_Tcell = 4*pi*(params_in.D_Tcell/2)^2;
params_out.A_Tcell.Notes = ['calculated based on the average T cell diameter ' params_in.D_Tcell.Notes];

params_out.A_syn   = params_in.A_syn;

% OX40 and OX40l Expressions
params_out.nT_OX40_tot    = params_in.nT_OX40_tot;
params_out.APC_OX40l_tot  = params_in.APC_OX40l_tot;

% Hill Function Parameters
params_out.nT_OX40_50 = params_in.nT_OX40_50;
params_out.n_nT_OX40  = params_in.n_nT_OX40;
