%% Function that save necessary data from an in silico clinical trial
function sprint_data(simDataPSA, simDataPSApost, params_out, name)
clear n_PSA
try
    load('VCT_data.mat')
catch
    disp('No VCT_data.mat was found or a loading error occurred.')
end

try
    if isfield(n_PSA, name)
        disp('There is a pre-existing trial with the same code name.')
        V_T = rmfield(V_T, name); nT1ln = rmfield(nT1ln, name);
        nT0ln = rmfield(nT0ln, name); aT1ln = rmfield(aT1ln, name);
        aT0ln = rmfield(aT0ln, name); aThln = rmfield(aThln, name);
        T1ln = rmfield(T1ln, name); T0ln = rmfield(T0ln, name);
        Thln = rmfield(Thln, name); T1p = rmfield(T1p, name);
        T0p = rmfield(T0p, name); Thp = rmfield(Thp, name);
        Teff = rmfield(Teff, name); Treg = rmfield(Treg, name);
        MDSC = rmfield(MDSC, name); CD8 = rmfield(CD8, name);
        D_T = rmfield(D_T, name); D_T_perc = rmfield(D_T_perc, name);
        nCD8c = rmfield(nCD8c, name); nCD4c = rmfield(nCD4c, name);

        CD8c = rmfield(CD8c, name); Tregc = rmfield(Tregc, name);
        Thc = rmfield(Thc, name); Ag0 = rmfield(Ag0, name);
        Ag1 = rmfield(Ag1, name); APC = rmfield(APC, name);
        mAPC = rmfield(mAPC, name); T1exh = rmfield(T1exh, name);
        Thexh = rmfield(Thexh, name); Th = rmfield(Th, name);
        CD4 = rmfield(CD4, name); C1 = rmfield(C1, name);
        C2 = rmfield(C2, name); n_PSA = rmfield(n_PSA, name);
        index = rmfield(index, name); ResObs = rmfield(ResObs, name);
        ResObs_R = rmfield(ResObs_R, name); ORR = rmfield(ORR, name);
        dor = rmfield(dor, name);
        disp('Data from the pre-existing trial has been replaced.')
        disp('---------------------------------------------')
    end
catch
    disp('New VCT_data.mat file is generated.')
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
    [~,temp19,~] = selectbyname(simDataPSA(index.(name)(i)).simData, 'V_T.C2');

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
    C1.(name)(i,:) = temp18./temp1;
    C2.(name)(i,:) = temp19./temp1;

    aThln.(name)(i,:) = temp20;
    T1ln.(name)(i,:) = temp21;
    T0ln.(name)(i,:) = temp22;
    Thln.(name)(i,:) = temp23;
    T1p.(name)(i,:) = temp24;
    T0p.(name)(i,:) = temp25;
    Thp.(name)(i,:) = temp26;

    A_syn = 37.8; % micrometer^2
    PDL1_tum.(name)(i,:) = temp27.*A_syn;
    PDL1_APC.(name)(i,:) = temp28.*A_syn;
    PDL2_tum.(name)(i,:) = temp29.*A_syn;
    PDL2_APC.(name)(i,:) = temp30.*A_syn;

    param_Obs = 'Teff_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    Teff.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'Treg_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    Treg.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

    param_Obs = 'MDSC_density';
    k = find(strcmp(simDataPSApost(index.(name)(i)).simData.DataNames, param_Obs));
    MDSC.(name)(i,:) = simDataPSApost(index.(name)(i)).simData.Data(:,k);

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
time = [0:1:400]';

for i = 1:n_PSA.(name)
    dt = D_T_perc.(name)(i,:);
    idx_min = find(dt==min(dt), 1);

    if min(dt) <= -30
        % ResObs is for easy calculation of ORR
        ResObs_(i,1) = 1; % PR/CR
        % ResObs_R is for generating data for RStudio
        ResObs_R_(i,1) = {'R'}; % PR/CR
        % ORR is for easy calculation and bootstrap of PD, PR/CR, and SD
        ORR_(i,1) = 1; % PR/CR
        idx = find(dt<=-30, 1);
        if (max(D_T.(name)(i,idx_min:end)) >= 1.2*min(D_T.(name)(i,idx_min))) && (max(D_T.(name)(i,idx_min:end)) - min(D_T.(name)(i,:)) >= 0.5)
            dor_(i) = time(find(D_T.(name)(i,idx_min:end) >= 1.2*min(D_T.(name)(i,idx_min)), 1)+idx_min-1)-time(idx);
            ttp_(i) = time(find(D_T.(name)(i,idx_min:end) >= 1.2*min(D_T.(name)(i,idx_min)), 1)+idx_min-1);
        else
            dor_(i) = time(end)-time(idx);
            ttp_(i) = time(end) + 1;
        end

    elseif (max(D_T.(name)(i,idx_min:end)) >= 1.2*min(D_T.(name)(i,idx_min))) && (max(D_T.(name)(i,idx_min:end)) - min(D_T.(name)(i,:)) >= 0.5) ...
            && (time(find(D_T.(name)(i,idx_min:end) >= 1.2*min(D_T.(name)(i,idx_min)), 1)+idx_min-1) - time(1) <= 56)
        ResObs_(i,1) = 0; % PD
        ResObs_R_(i,1) = {'NR'}; % PD
        ORR_(i,1) = 0; % PD
        dor_(i) = 0;
        ttp_(i) = time(find(dt >= 1.2*(100+dt(1))-100, 1)) - time(1);

    else
        ResObs_(i,1) = 0; % SD
        ResObs_R_(i,1) = {'NR'}; % SD
        ORR_(i,1) = 0.5; % SD
        dor_(i) = 0;
        ttp_(i) = time(end) + 1;
    end

end

dor_cor = 56.*ceil(dor_./56)./30;
% median(dor_cor(dor_cor~=0))

% plot(time(1:56:401), D_T_perc.(name)(:,1:56:401), '-.^', 'LineWidth', 2)
% ylim([-100 100])

% Complete Response
len = size(D_T.(name),1);
temp = zeros(len,1);
for i = 1:len
    if D_T.(name)(i,end) <= .2 % PMID: 28678153
        temp(i) = 1;
    else
        temp(i) = 0;
    end
end

% Overall Response
len = size(D_T.(name),1);
temp5 = zeros(len,1);

for i = 1:len
    if D_T_perc.(name)(i,end) <= -30 && 1.2*min(D_T.(name)(i,:)) > min(D_T.(name)(i,end))
        temp5(i) = 1;
    else
        temp5(i) = 0;
    end
end


%% Save Results
ResObs.(name) = ResObs_;
ResObs_R.(name) = ResObs_R_;
ORR.(name) = ORR_;
dor.(name) = dor_;

filename = 'VCT_data.mat';
save(filename, 'V_T', 'nT1ln', 'nT0ln', 'aT1ln', 'aT0ln', 'aThln',...
    'T1ln', 'T0ln', 'Thln', 'T1p', 'T0p', 'Thp',...
    'Teff', 'Treg', 'MDSC', 'CD8', 'D_T', 'D_T_perc',...
    'nCD8c', 'nCD4c', 'CD8c', 'Tregc', 'Thc', 'Ag0',...
    'Ag1', 'APC', 'mAPC', 'T1exh', 'Thexh', 'Th',...
    'CD4', 'C1', 'C2', 'n_PSA', 'index', 'ResObs', ...
    'ResObs_R', 'ORR', 'dor', 'PDL1_tum', 'PDL1_APC', 'PDL2_tum', 'PDL2_APC')


%% Print Results
CD8_tot_t = Teff.(name) + T1exh.(name);
CD4_tot_t = Treg.(name) + Th.(name) + Thexh.(name);

disp(['Objective Response Rate: ' num2str(round(sum(ResObs_)/n_PSA.(name), 3, 'significant'))])
disp(['Complete Response Rate: ' num2str(round(sum(temp)/n_PSA.(name), 3, 'significant'))])
disp(['Response Rate at data cutt-off date: ' num2str(round(sum(temp5)/n_PSA.(name), 3, 'significant'))])

disp(['mean CD8 in tumor: ', num2str(round(mean(CD8_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.7e7
disp(['mean CD4 in tumor: ', num2str(round(mean(CD4_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.9e7
disp(['mean Tregs in tumor: ' num2str(round(mean(Treg.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 3.5e6
disp(['median CD8/Treg in tumor: ' num2str(round(median(CD8_tot_t(:,1)./(Treg.(name)(:,1))), 3, 'significant'))])
disp(['median CD4/Treg in tumor: ' num2str(round(median(CD4_tot_t(:,1)./(Treg.(name)(:,1))), 3, 'significant'))])

disp(['median CD8 in tumor: ', num2str(round(median(CD8_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 8.4e7
disp(['median CD4 in tumor: ', num2str(round(median(CD4_tot_t(:,1)), 3, 'significant')), ' cell/mL']) % 1.4e8
disp(['median Teffs in tumor: ', num2str(round(median(Teff.(name)(:,1)), 3, 'significant')), ' cell/mL'])
disp(['median Tregs in tumor: ', num2str(round(median(Treg.(name)(:,1)), 3, 'significant')), ' cell/mL']) % 5.5e7

% for i = 1:n_PSA
%     TVDT(i,1) = 100*log(2)/(log(V_T.(name)(i,100)) - log(V_T.(name)(i,1)));
% end
% disp(['median tumor doubling time: ' num2str(median(TVDT))])
disp(['median duration of response: ' num2str(round(median(dor_cor(dor_cor~=0)), 3, 'significant'))])

% max(CD8_tot_t(:,1)./Treg.(name)(:,1))
% min(CD8_tot_t(:,1)./Treg.(name)(:,1))
% max(CD4_tot_t(:,1)./Treg.(name)(:,1))
% min(CD4_tot_t(:,1)./Treg.(name)(:,1))

end
