%% Function that save necessary data from an in silico clinical trial
function sprint_data(simDataPSA, simDataPSApost, params_in, params_out, name, varargin)

% Optional Inputs
in = inputParser;
% Frequency of tumor imaging during the trial (in weeks)
addParameter(in,'freq',8);
% Simulation time
addParameter(in,'time',[0:400]');
% Parse Inputs
parse(in,varargin{:});
time = in.Results.time;
f_measure_wk = in.Results.freq;

new_file = false;
try
    load('VCT_data.mat')
catch
    disp('No VCT_data.mat was found')
    disp('A new VCT_data.mat will be generated.')
    disp('---------------------------------------------')
    new_file = true;
end

if ~(new_file)
    temp = load('VCT_data.mat');
    DataNames = fieldnames(temp);
    for i = 1:length(DataNames)
        if isfield(temp.(DataNames{i}), name)
            species = eval(DataNames{i});
            species = rmfield(species, name);
            eval([DataNames{i} ' = species;'])
        end
    end
    disp('New trial data will be added to VCT_data.mat.')
    disp('Data from the pre-existing trial with the same name will be replaced.')
    disp('---------------------------------------------')
end

n_PSA.(name) = length(params_out.iPatientPlaus);
index.(name) = params_out.iPatientPlaus;

%% Time-dependent Variables
for i = 1:n_PSA.(name)

    [~,temp1,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T');
    [~,temp2,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_C.nT1');
    [~,temp3,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_C.nT0');
    [~,temp4,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_C.T1');
    [~,temp5,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_C.T0');
    [~,temp6,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_C.Th');
    [~,temp7,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'A_s.M1p1');
    [~,temp8,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'A_s.M1p0');
    [~,temp9,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.APC');
    [~,temp10,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.mAPC');

    [~,temp11,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.T1_exh');
    [~,temp12,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.Th');
    [~,temp13,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.Th_exh');

    [~,temp14,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.nT1');
    [~,temp15,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.nT0');
    [~,temp16,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.aT1');
    [~,temp17,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.aT1');
    [~,temp18,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.C1');
    % [~,temp19,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.C2');

    [~,temp20,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.aTh');
    [~,temp21,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.T1');
    [~,temp22,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.T0');
    [~,temp23,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.Th');
    [~,temp24,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_P.T1');
    [~,temp25,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_P.T0');
    [~,temp26,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_P.Th');

    [~,temp27,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_C1.PDL1_total');
    [~,temp28,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_APC.PDL1_total');
    [~,temp29,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_C1.PDL2_total');
    [~,temp30,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_APC.PDL2_total');

    [~,temp31,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.Mac_M1');
    [~,temp32,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.Mac_M2');
    [~,temp33,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.C_x');
    [~,temp34,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'H_PD1_M');
    [~,temp35,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'H_Mac_C');
    [~,temp36,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_C1.PD1_PDL1');
    [~,temp37,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'syn_T_C1.PD1_PDL2');
    [~,temp38,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'H_PD1_C1');

    [~,temp39,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_LN.IL2');
    [~,temp40,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.IL10');
    [~,temp41,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.IL12');
    [~,temp42,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.IFNg');
    [~,temp43,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.TGFb');
    [~,temp44,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.CCL2');
    % [~,temp45,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.c_vas');

    V_T.(name)(i,:) = temp1;

    nCD8c.(name)(i,:) = temp2.*1.11e6./5000;
    nCD4c.(name)(i,:) = temp3.*1.16e6./5000;
    CD8c.(name)(i,:) = temp4./5000;
    Tregc.(name)(i,:) = temp5./5000;
    Thc.(name)(i,:) = temp6./5000;
    Ag1.(name)(i,:) = temp7;
    Ag0.(name)(i,:) = temp8;
    APC.(name)(i,:) = temp9;
    mAPC.(name)(i,:) = temp10;
    T1exh.(name)(i,:) = temp11./temp1;
    Th.(name)(i,:) = temp12./temp1;
    Thexh.(name)(i,:) = temp13./temp1;
    nT1ln.(name)(i,:) = temp14;
    nT0ln.(name)(i,:) = temp15;
    aT1ln.(name)(i,:) = temp16;
    aT0ln.(name)(i,:) = temp17;
    C1.(name)(i,:) = temp18;
    % C2.(name)(i,:) = temp19;

    aThln.(name)(i,:) = temp20;
    T1ln.(name)(i,:) = temp21;
    T0ln.(name)(i,:) = temp22;
    Thln.(name)(i,:) = temp23;
    T1p.(name)(i,:) = temp24;
    T0p.(name)(i,:) = temp25;
    Thp.(name)(i,:) = temp26;

    A_cell = 907.92; A_APC = 900; % micrometer^2
    PDL1_tum.(name)(i,:) = temp27.*A_cell;
    PDL1_APC.(name)(i,:) = temp28.*A_APC;
    PDL2_tum.(name)(i,:) = temp29.*A_cell;
    PDL2_APC.(name)(i,:) = temp30.*A_APC;

    M1.(name)(i,:) = temp31./temp1;
    M2.(name)(i,:) = temp32./temp1;
    M_ratio.(name)(i,:) = temp31./temp32;
    Cx.(name)(i,:) = temp33;
    H_PD1_M.(name)(i,:) = temp34;
    H_Mac_C.(name)(i,:) = temp35;
    PD1_PDL1.(name)(i,:) = temp36;
    PD1_PDL2.(name)(i,:) = temp37;
    H_PD1_C.(name)(i,:) = temp38;

    IL2.(name)(i,:) = temp39;
    IL10.(name)(i,:) = temp40;
    IL12.(name)(i,:) = temp41;
    IFN.(name)(i,:) = temp42;
    TGFb.(name)(i,:) = temp43;
    CCL2.(name)(i,:) = temp44;
    % VEGF.(name)(i,:) = temp45;

    param_Obs = 'Teff_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    Teff.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'Treg_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    Treg.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'MDSC_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    MDSC.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'M_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    Mac.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'M_ratio';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    M_ratio.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'CD8_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    CD8.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'CD4_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    CD4.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    D_T.(name)(i,:) = 2*(3/(4*pi).*V_T.(name)(i,:)).^(1/3);
    D_T_perc.(name)(i,:) = (D_T.(name)(i,:) - D_T.(name)(i,1))./D_T.(name)(i,1).*100;

end


%% Outcome Predictions
time = time(1:f_measure_wk*7:end);
for i = 1:n_PSA.(name)
    dt = D_T_perc.(name)(i,1:f_measure_wk*7:end);
    di = D_T.(name)(i,1:f_measure_wk*7:end);
    idx_min = find(dt==min(dt), 1);

    if min(di) <= .2 % PMID: 28678153; 16612588
        CR_(i,1) = 1; % CR
        PR_CR_(i,1) = 1; % PR/CR
        ORR_(i,1) = 1; % ORR: 1=PR/CR; .5=SD; 0=PD
        PR_(i,1) = 0; % PR
        SD_(i,1) = 0; % SD
        PD_(i,1) = 0; % PD
        R_status_(i,1) = {'R'};

        idx = find(dt<=-30, 1);
        if (max(di(idx_min:end)) >= 1.2*min(di(idx_min))) && (max(di(idx_min:end)) - min(di) >= 0.5)
            dor_(i) = time(find(di(idx_min:end) >= 1.2*min(di(idx_min)), 1)+idx_min-1)-time(idx);
            ttp_(i) = time(find(di(idx_min:end) >= 1.2*min(di(idx_min)), 1)+idx_min-1);
        else
            dor_(i) = time(end)-time(idx);
            ttp_(i) = time(end) + 1;
        end

    elseif min(dt) <= -30
        CR_(i,1) = 0; % CR
        PR_CR_(i,1) = 1; % PR/CR
        ORR_(i,1) = 1; % ORR: 1=PR/CR; .5=SD; 0=PD
        PR_(i,1) = 1; % PR
        SD_(i,1) = 0; % SD
        PD_(i,1) = 0; % PD
        R_status_(i,1) = {'R'};

        idx = find(dt<=-30, 1);
        if (max(di(idx_min:end)) >= 1.2*min(di(idx_min))) && (max(di(idx_min:end)) - min(di) >= 0.5)
            dor_(i) = time(find(di(idx_min:end) >= 1.2*min(di(idx_min)), 1)+idx_min-1)-time(idx);
            ttp_(i) = time(find(di(idx_min:end) >= 1.2*min(di(idx_min)), 1)+idx_min-1);
        else
            dor_(i) = time(end)-time(idx);
            ttp_(i) = time(end) + 1;
        end

    elseif (max(di(idx_min:end)) >= 1.2*min(di(idx_min))) && (max(di(idx_min:end)) - min(di) >= 0.5) ...
            && (time(find(di(idx_min:end) >= 1.2*min(di(idx_min)), 1)+idx_min-1) - time(1) <= f_measure_wk*7)

        CR_(i,1) = 0; % CR
        PR_CR_(i,1) = 0; % PR/CR
        ORR_(i,1) = 0; % ORR: 1=PR/CR; .5=SD; 0=PD
        PR_(i,1) = 0; % PR
        SD_(i,1) = 0; % SD
        PD_(i,1) = 1; % PD
        R_status_(i,1) = {'NR'};

        dor_(i) = 0;
        ttp_(i) = time(find(dt >= 1.2*(100+dt(1))-100, 1)) - time(1);

    else
        CR_(i,1) = 0; % CR
        PR_CR_(i,1) = 0; % PR/CR
        ORR_(i,1) = .5; % ORR: 1=PR/CR; .5=SD; 0=PD
        PR_(i,1) = 0; % PR
        SD_(i,1) = 1; % SD
        PD_(i,1) = 0; % PD
        R_status_(i,1) = {'NR'};

        dor_(i) = 0;
        ttp_(i) = time(end) + 1;
    end

end

dor_mo = dor_./30;


%% Save Results
CR.(name) = CR_;
PR_CR.(name) = PR_CR_;
PR.(name) = PR_;
SD.(name) = SD_;
PD.(name) = PD_;
ORR.(name) = ORR_;
R_status.(name) = R_status_;
dor.(name) = dor_;

params_in_.(name)  = params_in;
params_out_.(name) = params_out;

filename = 'VCT_data.mat';
save(filename, 'V_T', 'nT1ln', 'nT0ln', 'aT1ln', 'aT0ln', 'aThln',...
    'T1ln', 'T0ln', 'Thln', 'T1p', 'T0p', 'Thp',...
    'Teff', 'Treg', 'MDSC', 'CD8', 'D_T', 'D_T_perc',...
    'nCD8c', 'nCD4c', 'CD8c', 'Tregc', 'Thc', 'Ag0',...
    'Ag1', 'APC', 'mAPC', 'T1exh', 'Thexh', 'Th',...
    'CD4', 'C1', 'n_PSA', 'index', 'CR', 'PR_CR',...
    'PR', 'SD', 'PD', 'R_status', 'ORR', 'dor', 'PDL1_tum', 'PDL1_APC',...
    'PDL2_tum', 'PDL2_APC', 'M1', 'M2', 'Mac', 'M_ratio',...
    'H_PD1_M','H_Mac_C','Cx','IL2','IL10','IL12','IFN','TGFb','CCL2',...
    'PD1_PDL1', 'PD1_PDL2', 'H_PD1_C', 'params_in_', 'params_out_')


%% Print Results
CD8_tot_t = Teff.(name) + T1exh.(name);
CD4_tot_t = Treg.(name) + Th.(name) + Thexh.(name);

disp(['Objective Response Rate: ' num2str(round(sum(PR_CR_)/n_PSA.(name), 3, 'significant'))])
disp(['Complete Response Rate: ' num2str(round(sum(CR_)/n_PSA.(name), 3, 'significant'))])

disp(['mean CD8 in tumor: ', num2str(round(mean(CD8_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.7e7
disp(['mean CD4 in tumor: ', num2str(round(mean(CD4_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.9e7
disp(['mean Tregs in tumor: ' num2str(round(mean(Treg.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 3.5e6
disp(['median CD8/Treg in tumor: ' num2str(round(median(CD8_tot_t(:,1)./(Treg.(name)(:,1))), 3, 'significant'))])
disp(['median CD4/Treg in tumor: ' num2str(round(median(CD4_tot_t(:,1)./(Treg.(name)(:,1))), 3, 'significant'))])

disp(['median CD8 in tumor: ', num2str(round(median(CD8_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 8.4e7
disp(['median CD4 in tumor: ', num2str(round(median(CD4_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.4e8
disp(['median Teffs in tumor: ', num2str(round(median(Teff.(name)(:,1)), 3, 'significant')), ' cell/mL'])
disp(['median Tregs in tumor: ', num2str(round(median(Treg.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 5.5e7
disp(['median MDSC in tumor: ', num2str(round(median(MDSC.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 8.4e7
disp(['median macrophage in tumor: ', num2str(round(median(Mac.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 8.4e7
disp(['median M1/M2 in tumor: ', num2str(round(median(M_ratio.(name)(:,1)), 3, 'significant'))]) % 8.4e7

disp(['median IL2 (15.3 kDa): ', num2str(round(median(IL2.(name)(:,1).*15300), 3, 'significant')), ' pg/mL'])
disp(['IL2 range: ', num2str(round(min(IL2.(name)(:,1).*15300), 3, 'significant')), ' - ', num2str(round(max(IL2.(name)(:,1).*15300), 3, 'significant')), ' pg/mL'])
disp(['median IL10 (18 kDa) in tumor: ', num2str(round(median(IL10.(name)(:,1).*18000), 3, 'significant')), ' pg/mL'])
disp(['IL10 range in tumor: ', num2str(round(min(IL10.(name)(:,1).*18000), 3, 'significant')), ' - ', num2str(round(max(IL10.(name)(:,1).*18000), 3, 'significant')), ' pg/mL'])
disp(['median IL12 (70 kDa): ', num2str(round(median(IL12.(name)(:,1).*70000), 3, 'significant')), ' pg/mL'])
disp(['IL12 range: ', num2str(round(min(IL12.(name)(:,1).*70000), 3, 'significant')), ' - ', num2str(round(max(IL12.(name)(:,1).*70000), 3, 'significant')), ' pg/mL'])
disp(['median IFN-gamma (17 kDa) in tumor: ', num2str(round(median(IFN.(name)(:,1).*17000), 3, 'significant')), ' pg/mL'])
disp(['IFN-gamma range in tumor: ', num2str(round(min(IFN.(name)(:,1).*17000), 3, 'significant')), ' - ', num2str(round(max(IFN.(name)(:,1).*17000), 3, 'significant')), ' pg/mL'])
disp(['median TGF-beta (25 kDa) in tumor: ', num2str(round(median(TGFb.(name)(:,1).*25000), 3, 'significant')), ' pg/mL'])
disp(['TGF-beta range in tumor: ', num2str(round(min(TGFb.(name)(:,1).*25000), 3, 'significant')), ' - ', num2str(round(max(TGFb.(name)(:,1).*25000), 3, 'significant')), ' pg/mL'])
disp(['median CCL2 (14 kDa) in tumor: ', num2str(round(median(CCL2.(name)(:,1).*14000), 3, 'significant')), ' pg/mL'])
disp(['CCL2 range in tumor: ', num2str(round(min(CCL2.(name)(:,1).*14000), 3, 'significant')), ' - ', num2str(round(max(CCL2.(name)(:,1).*14000), 3, 'significant')), ' pg/mL'])

% for i = 1:n_PSA
%     TVDT(i,1) = 100*log(2)/(log(V_T.(name)(i,100)) - log(V_T.(name)(i,1)));
% end
% disp(['median tumor doubling time: ' num2str(median(TVDT))])
disp(['median duration of response: ' num2str(round(median(dor_mo(dor_mo~=0)), 3, 'significant'))])

% max(CD8_tot_t(:,1)./Treg.(name)(:,1))
% min(CD8_tot_t(:,1)./Treg.(name)(:,1))
% max(CD4_tot_t(:,1)./Treg.(name)(:,1))
% min(CD4_tot_t(:,1)./Treg.(name)(:,1))

end
