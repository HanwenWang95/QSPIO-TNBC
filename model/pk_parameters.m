% PK parameter database
%
% Inputs: name -- string containing drug name <optional>:
%                 - nivolumab
%                 - entinostat
%                 - atezolizumab
%                 - ipilimumab
%                 - pembrolizumab
%                 - tremelimumab
%                 - aCD47
%
% Output: params_out -- object containing model parameters
%                       - q_P--rate of diffusive transport C<->P
%                       - q_T--rate of diffusive transport C<->T
%                       - q_LN--rate of diffusive transport C<->LN
%                       - q_LD--rate of convective transport T->LN
%                       - k_cl--clearence rate from central
%                       - gamma_C--volume fraction in C
%                       - gamma_P--volume fraction in P
%                       - gamma_T--volume fraction in T
%                       - gamma_LN--volume fraction in LN


function params_out = pk_parameters(name)

% Optional Inputs
if ~exist('name','var')
    name = '';
end

% Drug Database
% - Transport Rates (q_P,q_T,q_LN,q_LD)
% - Clearence Rate (k_cl)
% - Volume Fraction (gamma_C,gamma_P,gamma_T,gamma_LN)
switch (name)

    case 'pembrolizumab'
        % Central to Peripheral
        q_P.Value = 9.16e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Ahamadi 2016, PMID: 27863186)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.315;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Ahamadi 2016, PMID: 27863186)';
        % Volume Fractions
        gamma_C.Value = 0.698;
        gamma_C.Notes = '(Ahamadi 2016, PMID: 27863186)';
        gamma_P.Value = 0.0789;
        gamma_P.Notes = '(Ahamadi 2016, PMID: 27863186)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'nivolumab'
        % Central to Peripheral
        q_P.Value = 8.8e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Bajaj 2017, PMID: 28019091)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.37;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Bajaj 2017, PMID: 28019091)';
        % Volume Fractions
        gamma_C.Value = 0.61;
        gamma_C.Notes = '(Bajaj 2017, PMID: 28019091)';
        gamma_P.Value = 0.0452;
        gamma_P.Notes = '(Bajaj 2017, PMID: 28019091)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'entinostat'
        MW_ENT = 3.764085e5; % mg/mole
        % Central to Peripheral
        q_P.Value = 320;
        q_P.Units = 'milliliter/second';
        q_P.Notes = '(optimized)';
        % Central to Tumour
        q_T.Value = 3200;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(estimated)';
        % Central to LN
        q_LN.Value = 320;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(optimized)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = .152;
        k_cl.Units = '1/hour';
        k_cl.Notes = '(optimized)';
        % Volume Fractions
        gamma_C.Value = 0.55;
        gamma_C.Notes = '(estimated)';
        gamma_P.Value = 0.062;
        gamma_P.Notes = '(optimized)';
        gamma_T.Value = 0.611;
        gamma_T.Notes = '(Stefanini 2008)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';
        % Buccal absorption rate
        params_out.k_a1.Value = 1.9;
        params_out.k_a1.Units = '1/hour';
        params_out.k_a1.Notes = '(optimized)';
        % GI absorption
        params_out.k_a2.Value = 2.5;
        params_out.k_a2.Units = '1/hour';
        params_out.k_a2.Notes = '(optimized)';
        % Non-linear clearance rate Vmax
        params_out.k_cln.Value = 46/MW_ENT;
        params_out.k_cln.Units = 'mole/(liter*hour)';
        params_out.k_cln.Notes = '(optimized)';
        % Drug conc. at .5*Vmax
        params_out.Kc.Value = 3.53/MW_ENT;
        params_out.Kc.Units = 'mole/liter';
        params_out.Kc.Notes = '(optimized)';

    case 'ipilimumab'
        % Central to Peripheral
        q_P.Value = 1.15e-5;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Feng 2014, PMID: 24433434)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.524;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Feng 2014, PMID: 24433434)';
        % Volume Fractions
        gamma_C.Value = 0.686;
        gamma_C.Notes = '(Feng 2014, PMID: 24433434)';
        gamma_P.Value = 0.0496;
        gamma_P.Notes = '(Feng 2014, PMID: 24433434)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'atezolizumab'
        % Central to Peripheral
        q_P.Value = 8.7e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Stroh 2017, PMID: 27981577)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.324;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Stroh 2017, PMID: 27981577)';
        % Volume Fractions
        gamma_C.Value = 0.61;
        gamma_C.Notes = '(Stroh 2017, PMID: 27981577)';
        gamma_P.Value = 0.068;
        gamma_P.Notes = '(Stroh 2017, PMID: 27981577)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'durvalumab'
        % Central to Peripheral
        q_P.Value = 5.48e-06;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Baverel 2018, PMID: 29243223)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.335;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Baverel 2018, PMID: 29243223)';
        % Volume Fractions
        gamma_C.Value = 0.7;
        gamma_C.Notes = '(Baverel 2018, PMID: 29243223)';
        gamma_P.Value = 0.061;
        gamma_P.Notes = '(Baverel 2018, PMID: 29243223)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'tremelimumab'
        % Central to Peripheral
        q_P.Value = 5.331e-06;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Wang 2014, PMID: 24737343)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.3276;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Wang 2014, PMID: 24737343)';
        % Volume Fractions
        gamma_C.Value = 0.794;
        gamma_C.Notes = '(Wang 2014, PMID: 24737343)';
        gamma_P.Value = 0.062;
        gamma_P.Notes = '(Wang 2014, PMID: 24737343)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';
end

% Diffusive Transport: C<->P
params_out.q_P = q_P;
% Diffusive Transport: C<->T
params_out.q_T = q_T;
% Diffusive Transport: C<->LN
params_out.q_LN = q_LN;
% Convective Transport: T->LN || LN->C
params_out.q_LD = q_LD;
% Clearence from Central
params_out.k_cl = k_cl;
% Volume Fractions

gamma_C.Units = 'dimensionless';
gamma_P.Units = 'dimensionless';
gamma_T.Units = 'dimensionless';
gamma_LN.Units = 'dimensionless';
params_out.gamma_C = gamma_C;
params_out.gamma_P = gamma_P;
params_out.gamma_T = gamma_T;
params_out.gamma_LN = gamma_LN;
