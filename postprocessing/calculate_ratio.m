% This function calculates the ratios of T cell sunsets
%
% Inputs: simDataPSA     -- Object containing the whole simbiology model
%                           outputs for all batch simulations
%         simDataPSApost -- Object containing the postprocessed simbiology
%                           model outputs for all batch simulations
%         params_out     -- object containing model outputs to be organized
%                           for future sensitivity analysis
%
% Outputs: simDataPSAout  -- Updated object containing postprocessed
%                            outputs


function simDataPSAout = calculate_ratio(simDataPSA,simDataPSApost,params_out)

n_PSA = length(params_out.iPatient);
index = params_out.iPatient;
simDataPSAout = simDataPSApost;

% Need to define these based on the literature or individual patient data
div_CD4     = 1.16e6;      % diversity of naive T cells
div_CD8     = 1.11e6;      % diversity of naive T cells

for i = 1:n_PSA
    nT_C_total = 0; % total naive CD8 T cells in the blood
    T_C_total = 0;  % total Teff cells in the blood
    nT_T_total = 0;
    T_T_total = 0;
    T1exh_total = 0;

    nT0_C_total = 0; % total naive CD4 T cells needs to be determined
    T0_C_total = 0;
    nT0_T_total = 0;
    T0_T_total = 0;
    Th_T_total = 0;
    Thexh_total = 0;

    % Calculates the number of clones for that simulation
    numClones = howManyClones(simDataPSA(index(i)).simData);
    % Calculates total CD8 T cells
    for j=1:numClones
        [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_C.nT',num2str(j)]);
        nT_C_total = nT_C_total + Ti*div_CD8;
        [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_C.T',num2str(j)]);
        T_C_total = T_C_total + Ti;

        [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.nT',num2str(j)]);
        nT_T_total = nT_T_total + Ti*div_CD8;
        [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.T',num2str(j)]);
        T_T_total = T_T_total + Ti;
    end
    % Calculates exhausted CD8 T cells
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.T1_exh']);
    T1exh_total = T1exh_total + Ti;

    % Calculate total CD4 T cells
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_C.nT0']);
    nT0_C_total = nT0_C_total + Ti*div_CD4;
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_C.T0']);
    T0_C_total = T0_C_total + Ti;

    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.nT0']);
    nT0_T_total = nT0_T_total + Ti*div_CD4;
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.T0']);
    T0_T_total = T0_T_total + Ti;
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.Th']);
    Th_T_total = Th_T_total + Ti;
    [~,Ti,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T.Th_exh']);
    Thexh_total = Thexh_total + Ti;

    % Calculate total MDSC in Tumor
    [~,MDSC_total,~] = selectbyname(simDataPSA(index(i)).simData, 'V_T.MDSC');

    % Calculate ratios of T cell subsets
    CD8FoxP3ratio_C = (nT_C_total + T_C_total)./T0_C_total;
    CD8FoxP3ratio_T = (nT_T_total + T_T_total + T1exh_total)./T0_T_total;
    CD4FoxP3ratio_C = (nT0_C_total + T0_C_total)./T0_C_total;
    CD4FoxP3ratio_T = (nT0_T_total + T0_T_total + Th_T_total + Thexh_total)./T0_T_total;
    Teff_Treg_C = T_C_total./T0_C_total;
    Teff_Treg_T = T_T_total./T0_T_total;

    % Calculate T cell density
    [~,V_T,~] = selectbyname(simDataPSA(index(i)).simData, 'V_T');
    Teff_density = T_T_total./V_T;
    Treg_density = T0_T_total./V_T;
    CD8_density = (nT_T_total + T_T_total + T1exh_total)./V_T;
    CD4_density = (nT0_T_total + T0_T_total + Th_T_total + Thexh_total)./V_T;
    % Calculate MDSC density
    MDSC_density = MDSC_total./V_T;

    % Add calculated densities and ratios to postprocess structure
    %simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD8FoxP3ratio_C'}];
    %simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD8FoxP3ratio_C];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD8FoxP3ratio_T'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD8FoxP3ratio_T];
    %simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD4FoxP3ratio_C'}];
    %simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD4FoxP3ratio_C];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD4FoxP3ratio_T'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD4FoxP3ratio_T];

    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'Teff_Treg_C'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , Teff_Treg_C];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'Teff_Treg_T'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , Teff_Treg_T];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD8_density'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD8_density];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'CD4_density'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , CD4_density];
    % Add calculated density to postprocess structure
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'Teff_density'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , Teff_density];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'Treg_density'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , Treg_density];
    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'MDSC_density'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , MDSC_density];
end
