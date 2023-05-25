clear
close all
sbioreset

%% Step 1. Generate plausible patients and model prediction
% LHS sampling
immune_oncology_model_NSCLC
n_PSA = 30000;
params_in  = PSA_param_in_NSCLC;
params_in  = PSA_setup(model,params_in,n_PSA);

trial_name = 'durva';

%% Add PK variability

load('PK_pars.mat')
names.pars = {};
names.pars = [names.pars; 'q_P_aPDL1'; 'k_cl_aPDL1'; 'gamma_P_aPDL1'; 'V_C'; 'k_cln_aPDL1'; 'Kc_aPDL1'];
names.labels = {};
names.labels = [names.labels; 'Capillary filtration rate of aPDL1'; 'Linear clearance rate of aPDL1';...
                              'Volume fraction of interstitium in Vp'; 'Total blood volume';...
                              'Non-linear clearance rate of aPDL1'; 'Michaelis–Menten constant for non-linear CL'];
params_in = PSA_setup_PK(params_in,n_PSA,names,psi,V,opt_avg);

%% Simulate for initial condition (pre-treatment)
endtime = 8000+1; output = zeros(n_PSA,3); dt = 56;
ub = [6 2.3 2.0 50]'; lb = [.8 0 0 .06]';
Vt = zeros(n_PSA,endtime); Treg = zeros(n_PSA,endtime);
CD8 = zeros(n_PSA,endtime); CD4 = zeros(n_PSA,endtime);
VDT = zeros(n_PSA,1); success = zeros(n_PSA,1);

for i = 1:n_PSA
    display(['Sample ',num2str(i),'/',num2str(n_PSA)]);
    % Set the new parameters
    if i > 1
        delete(model_PSA)
    end
    model_PSA = copyobj(model);
    variantObj = addvariant(model_PSA, ['v',num2str(i,'%5.5i')]);
    for j = 1:length(params_in.names)
        if ~isempty(sbioselect (model, 'Type', 'parameter', 'Name', params_in.names{j}))
            addcontent(variantObj, {'parameter', params_in.names{j}, 'Value', params_in.(params_in.names{j}).LHS(i)});
        elseif ~isempty(sbioselect (model, 'Type', 'compartment', 'Name', params_in.names{j}))
            addcontent(variantObj, {'compartment', params_in.names{j}, 'Capacity', params_in.(params_in.names{j}).LHS(i)});
        else
            disp(['Unable to identify parameter/compartment named ', params_in.names{j}, ' in the model'])
        end
    end
    % Set Initial Conditions
    for k = 1:length(variantObj.content)
        variantParam = variantObj.content{k};
        if strcmp(variantParam(2),'initial_tumour_diameter')
            % Get Initial Tumour Diameter from 'variant'
            D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
            tumour_diameter.Value = variantParam{4};
            tumour_diameter.Units = D_Ti.ValueUnits;
            tumour_diameter.Notes = D_Ti.Notes;
        end
    end
    % Calculate Target Tumour Volume
    tumour_volume = 4/3*pi*(tumour_diameter/2)^3;
    % Reset Output Times for IC Simulation (assumes time in days)
    config = getconfigset(model);
    set(config,'StopTime',8000);
    set(config.SolverOptions,'OutputTimes',0:8000);
    set(config.SolverOptions,'MaxStep',1);
    
    try
        % IC Simulation
        simData = sbiosimulate(model,variantObj);
        % Get Tumour Volume from Simulation
        [~,V_T,~] = selectbyname(simData, 'V_T'); % milliliter
        Vt(i,:) = V_T';
        target_V_T = tumour_volume.Value; % tumour volume in mL
        % Get Time of Target Tumour Size
        idx = find((target_V_T-V_T)<0,1);
        if isempty(idx) % tumour did not reach target size
            disp(['Initial conditions: ' num2str(target_V_T) 'mL not reached (max: ' num2str(max(V_T)) 'mL)']);
        else
            success(i) = 1;
            [~,M1,~] = selectbyname(simData, 'V_T.Mac_M1');
            [~,M2,~] = selectbyname(simData, 'V_T.Mac_M2');
            [~,Teff,~] = selectbyname(simData, 'V_T.T1');
            [~,T1exh,~] = selectbyname(simData, 'V_T.T1_exh');
            [~,T0,~] = selectbyname(simData, 'V_T.T0');
            [~,Th,~] = selectbyname(simData, 'V_T.Th');
            [~,Thexh,~] = selectbyname(simData, 'V_T.Th_exh');
            [~,C_total,~] = selectbyname(simData, 'V_T.C1');
            CD8(i,:) = Teff + T1exh;
            CD4(i,:) = T0 + Th + Thexh;
            Treg(i,:) = T0;
            Vt(i,:) = V_T';
            C_diag(i) = C_total(idx);
            C_time(i,:) = C_total;
            
            output(i,1) = M1(idx)./M2(idx);
            output(i,2) = Treg(i,idx)./CD8(i,idx);
            output(i,3) = CD8(i,idx)./CD4(i,idx);
            
            t1 = find(V_T>=target_V_T,1); t2 = t1+dt;
            V1 = V_T(t1); V2 = V_T(t2);
            VDT(i) = dt.*log(2)./(log(V2)-log(V1));
            % filter by lower and upper boundaries of selected model variables
            %             if output(i,1)>ub(1) || output(i,1)<lb(1) || output(i,2)>ub(2) || output(i,2)<lb(2) ||...
            %                     output(i,3)>ub(3) || output(i,3)<lb(3)
            %                 success(i) = 0;
            %                 disp('Does not fall within physiologically reasonable ranges');
            %                 disp(output(i,:))
            %             end
        end
    catch
        disp('There was an error while finding the initial conditions');
    end
end
idx_plausible = find(VDT~=0);
ini_diam      = params_in.initial_tumour_diameter.LHS(idx_plausible);
VDT_plausible = VDT(VDT~=0);
median(VDT_plausible)
mean(VDT_plausible)
% geomean(VDT_plausible)
% idx1 = find(ini_diam <= 3); VDT_subgroups(1) = median(VDT_plausible(idx1));
% idx2 = find(ini_diam > 3 & ini_diam <= 4); VDT_subgroups(2) = median(VDT_plausible(idx2));
% idx3 = find(ini_diam > 4 & ini_diam <= 5); VDT_subgroups(3) = median(VDT_plausible(idx3));
% idx4 = find(ini_diam > 5); VDT_subgroups(4) = median(VDT_plausible(idx4));
% bar(VDT_subgroups)

% save plausible patient predictions
N_PP = sum(success);
pred_PP = output(success==1,:);

save(['plausible_patients_' trial_name '.mat'], 'success', 'output', 'params_in', 'Vt', 'C_diag', 'C_time')


%% Step 2. Select virtual patients
% optimize relative prevalence beta
data = csvread('observed_variables.csv',1,0);
observed_variables = data(:,1:3);
pred_PP = output(success==1,1:3);

Mdl = KDTreeSearcher(observed_variables,'BucketSize',1);
NumPoints=10;
[~,D] = knnsearch(Mdl,log(pred_PP),'K',NumPoints);
BallVolume = nsphereVolume(size(pred_PP,2),max(D,[],2));
Density = NumPoints./BallVolume;
Density = Density/numel(observed_variables); % normalize

Npts=5;
VPtree = KDTreeSearcher(log(pred_PP),'BucketSize',1);
[IDX,DVPs] = knnsearch(VPtree,log(pred_PP),'K',Npts);
BallVolume2 = nsphereVolume(size(pred_PP,2),max(DVPs,[],2));
Density2 = Npts./BallVolume2;
ProbabilityofInclusion = Density./Density2;

runs=100;
fopt = @(p) OptimizeVPgeneration(p,ProbabilityofInclusion,log(pred_PP),runs,observed_variables);
%NumOfVPs
maxsf=1/max(ProbabilityofInclusion);

lower2 = 0;
upper2 = log10(maxsf);
options = optimoptions(@simulannealbnd,'TolFun',1e-15,'ReannealInterval',50000,'InitialTemperature',0.5,'MaxIter',1000,'Display','Iter','TemperatureFcn',@temperatureboltz,'AnnealingFcn', @annealingboltz,'AcceptanceFcn',@acceptancesa);
[k,Score_V] = simulannealbnd(fopt,log10(maxsf),lower2,upper2,options);
disp(['optimal beta: ' num2str(10^k) ' (max: ' num2str(maxsf) ')'])

%% generate virtual population with optimal beta
logbeta = k;
[gof,newselection,pvalue]= OptimizeVPgeneration(logbeta,ProbabilityofInclusion,log(pred_PP),1,observed_variables);
disp(['goodness of fit: ' num2str(gof)])
disp(['p value: ' num2str(pvalue')])

%% prepare final virtual patient cohort
n_PSA = sum(newselection);
idx = find(success==1);
idx = idx(newselection==1);

params_in.all = params_in.all(idx,:);
for i = 1:length(params_in.names)
    params_in.(params_in.names{i}).LHS = params_in.all(:,i);
end
params_out.names = {};
params_in.namesObs = {};

save(['virtual_patients_' trial_name '.mat'], 'newselection', 'idx', 'params_in', 'params_out')


%% PD-L1 inhibition Simulations
config = getconfigset(model);
set(config,'StopTime',400);
set(config.SolverOptions,'OutputTimes',0:400);

dose_schedule = schedule_dosing({'durvalumab'});

warning('off','all')
sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

iPatientPlaus = zeros(n_PSA,1);
% Prepare the data for the rest of the analysis
for i = 1:length(simDataPSApost) 
    iPatientPlaus(i) = ~isempty(simDataPSApost(i).simData);
end
temp = 1:n_PSA;
params_out.iPatientPlaus = temp(iPatientPlaus==1);

% Save and print data of interest (by assigning a unique code name for the trial)
sprint_data(simDataPSA, simDataPSApost, params_in, params_out, trial_name)

