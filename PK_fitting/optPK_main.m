%% Clinical data/model
clear
close all
sbioreset

% popPK model
Two_Compartment_PK

step = 0.1;
time = 0:step:365; % time in day
configsetObj = getconfigset(data_model);
set(configsetObj, 'SolverType', 'ode15s')
set(configsetObj.SolverOptions, 'OutputTimes', time)

% dosing schedule
patient_weight = 75; % kg
dose_durv = 10; % mg/kg
dose_schedule = sbiodose('A','Amount',patient_weight*dose_durv,'AmountUnits','milligram','TargetName','V_C.A');
dose_schedule.StartTime = 0;
dose_schedule.Interval = 14;
dose_schedule.TimeUnits = 'day';
dose_schedule.RepeatCount = 13;
dose_schedule.Active = true;

% Simulation
simData = sbiosimulate(data_model,[],[],dose_schedule);

figure
semilogy(time./7,simData.Data(:,1))
xlim([0  52])
ylim([1 1e4])
hold on

tindex = 1:length(time);
expData = simData.Data(tindex,1);

% Generate pseudo-data
WT = lognrnd(log(69.8),.2,[1,400]);
ALB = lognrnd(log(38),.15,[1,400]);
CRCL = lognrnd(log(85.65),.3,[1,400]);
ECOG = binornd(1,.65,[1,400]);
SEX = binornd(1,.433,[1,400]);
tumor = lognrnd(log(74.8),.5,[1,400]);
sPDL1 = lognrnd(log(124.8),.2,[1,400]); % ?

k_cl = .232 * (1-.063.*ECOG) .* (1-.143.*SEX) .* (1-.035.*(ALB-38))...
    .* (1+.0015.*(CRCL-85.65)) .* (1+.00178.*(tumor-74.8)) .* (WT./69.8).^.389;
V_1  = 3.51 * (1-.165.*SEX) .* (WT./69.8).^.406;
V_2  = 3.45 * (1-.205.*SEX);
V_max= .824 .* (1+.00336.*(sPDL1-124.8));

conc = zeros(length(time), 400);
for i = 1:400
    data_model.parameter(3).Value = k_cl(i);
    data_model.compartment(1).capacity = V_1(i);
    data_model.compartment(2).capacity = V_2(i);
    data_model.parameter(4).Value = V_max(i);
    
    simData   = sbiosimulate(data_model,[],[],dose_schedule);
    conc(:,i) = simData.Data(:,1);
end

plot(time./7, prctile(conc, 5, 2))
plot(time./7, prctile(conc,95, 2))


%% QSP PK module

qsp_pk_module

step = 0.1;
time = 0:step:365; % time in day
configsetObj = getconfigset(model);
set(configsetObj, 'SolverType', 'ode15s')
set(configsetObj.SolverOptions, 'OutputTimes', time)

MW = 1.463e8; % mg/mole durvalumab

%% Fitting to average model
dose_schedule_qsp = schedule_dosing_qsp({'durvalumab'});

%             Q_P,   kcl,   fp,   Vc,  kcln,   Kcl
p0    = log([6e-6;   .33;  .06;    6;   5.6;   2.1]);
upper = log([6e-3      1    .2    13    100    100]);
lower = log([6e-9    .01   .01     2     .1     .1]);
fobjective = @(p) costPK(p, model, dose_schedule_qsp, tindex, median(conc,2));
options = optimoptions(@fmincon, 'Display','notify');
% opt = patternsearch(fobjective,p0,[],[],[],[],lower,upper,[],options);
% opt = fminsearch(fobjective,p0);
opt = fmincon(fobjective,p0,[],[],[],[],lower,upper,[],options);
format shortG
disp(exp([opt(1) opt(2) opt(3) opt(4) opt(5) (opt(6))]))

%% Visualization
% opt = log([6.1e-6  .38  .052   5.78  5.88  1.5]');

Qp    = opt(1);
k_cl  = opt(2);
K_P   = opt(3);
Vc    = opt(4);
k_cln = opt(5);
Kcl   = opt(6);

model.parameter(1).Value = exp(Qp);
model.parameter(5).Value = exp(k_cl);
model.parameter(9).Value = exp(K_P);
model.compartment(1).Capacity = exp(Vc);
model.parameter(6).Value = exp(k_cln);
model.parameter(7).Value = exp(Kcl);

dose_schedule_qsp = schedule_dosing_qsp({'durvalumab'});
simData = sbiosimulate(model,[],[],dose_schedule_qsp);

K_B = model.parameter(8).Value;
Average_par = simData.Data(:,1) .* MW ./ K_B;

time = 0:step:365;

figure
hold on
plot(time, Average_par, 'LineWidth',2, 'Color',[209,41,89]./255)
plot(time, median(conc,2), 'k--');

hTitle = sgtitle('');

hXLabel = xlabel('Time (day)');
hYLabel = ylabel('Durvalumab Concentration (\mug/mL)');
hLegend = legend('Model prediction','Median clinical observation');

set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde', 'FontSize', 16)
set(hLegend, 'FontSize', 14)
set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')

% Adjust axes properties
ax = gca;
ax.FontSize = 14;

set(gca, 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 2)
hold off


%% Generate deviations and perform SVD

Theta = zeros(6,500);
for i = 1:6
    Theta(i,:) = normrnd(opt(i),.2,[1,500]);
end

Y_hat = zeros(length(time), 500);
for i = 1:500
    
    par = Theta(:,i);
    
    model.parameter(1).Value = exp(par(1));
    model.parameter(5).Value = exp(par(2));
    model.parameter(9).Value = exp(par(3));
    model.compartment(1).Capacity = exp(par(4));
    model.parameter(6).Value = exp(par(5));
    model.parameter(7).Value = exp(par(6));
    
    simData = sbiosimulate(model,[],[],dose_schedule_qsp);
    Y_hat(:,i) = simData.Data(:,1) .* MW ./ K_B - Average_par;
    
end

Cov =  Y_hat * (Theta-opt)';
[U,S,V] = svd(Cov);


%% Fitting to pseudo-data
opt_avg = opt;
p0 = opt_avg;

for i = 1:400
    
    fobjective = @(p) costPK_CLP(p, model, dose_schedule_qsp, tindex, conc(:,i), V, opt_avg);
    opt_all(:,i) = fmincon(fobjective,p0,[],[],[],[],lower,upper,[],options);
    
    disp([num2str(i) '/' num2str(400)])
    
end

save results.mat opt_all V opt_avg

%% Latent space
psi = opt_all'*V - opt_avg'*V;
psi = rmoutliers(psi);
psi = psi';
dose_schedule_qsp = schedule_dosing_qsp({'durvalumab'});

names_pars = {'q_P_aPDL1'; 'k_cl_aPDL1'; 'gamma_P_aPDL1'; 'V_C'; 'k_cln_aPDL1'; 'Kc_aPDL1'};
save PK_pars.mat names_pars V psi opt_avg

%% VP simulation
dose_schedule_qsp = schedule_dosing_qsp({'durvalumab'});
K_B = model.parameter(8).Value;

for i = 1:6
    p(i,:) = min(psi(i,:)) + (max(psi(i,:))-min(psi(i,:))) * rand(1,400);
end

opt_new_all = V*p + opt_avg;
for i = 1:400
    
    opt_new = opt_new_all(:,i);
    
    Qp    = opt_new(1);
    k_cl  = opt_new(2);
    K_P   = opt_new(3);
    Vc    = opt_new(4);
    k_cln = opt_new(5);
    Kcl   = opt_new(6);
    
    model.parameter(1).Value = exp(Qp);
    model.parameter(5).Value = exp(k_cl);
    model.parameter(9).Value = exp(K_P);
    model.compartment(1).Capacity = exp(Vc);
    model.parameter(6).Value = exp(k_cln);
    model.parameter(7).Value = exp(Kcl);
    
    simData = sbiosimulate(model,[],[],dose_schedule_qsp);
    
    plasma_conc(:,i) = simData.Data(tindex,1) .* MW ./ K_B;
    
end
    
time = 0:step:365;
figure
semilogy(time./7, median(plasma_conc,2))
xlim([0  52])
ylim([1 1e4])
hold on
plot(time./7, prctile(plasma_conc, 5, 2))
plot(time./7, prctile(plasma_conc,95, 2))
    
figure
boxplot(opt_new_all')

%% Plot histograms 
figure
subplot(2,3,1)
histogram(exp(opt_new_all(1,:)),12)
title({'Capillary filtration rate' 'between central and peripheral compartment' 'Qp'})
xlabel('L/s')
set(gca,'FontSize',14)

subplot(2,3,2)
histogram(exp(opt_new_all(2,:)),12)
title({'Linear clearance rate from central compartment' 'kcl'})
xlabel('L/day')
set(gca,'FontSize',14)

subplot(2,3,3)
histogram(exp(opt_new_all(3,:)),12)
title({'Volume fraction of interstial space' 'available to durvalumab' '\gamma_p'})
xlabel('dimensionless')
set(gca,'FontSize',14)

subplot(2,3,4)
histogram(exp(opt_new_all(4,:)),12)
title({'Blood volume' 'Vc'})
xlabel('L')
set(gca,'FontSize',14)

subplot(2,3,5)
histogram(exp(opt_new_all(5,:)),12)
title({'Maximum non-linear clearance rate' 'Vmax'})
xlabel('nmol/day')
set(gca,'FontSize',14)

subplot(2,3,6)
histogram(exp(opt_new_all(6,:)),12)
title({'Half-maximal durvalumab concentration' '(non-linear clearance)' 'Km'})
xlabel('nmol/L')
set(gca,'FontSize',14)

%% Drug concentration 
figure
hold on
patch([time flip(time)]./7, [prctile(plasma_conc, 5, 2); flip(prctile(plasma_conc, 95, 2))],[255, 190, 106]./255,'FaceAlpha',.2);
plot(time./7, median(plasma_conc,2), 'Color', [64, 176, 166]./255, 'LineWidth',2)
plot(time./7, prctile(plasma_conc, 5, 2), 'Color', [255, 190, 106]./255, 'LineWidth',2)
plot(time./7, prctile(plasma_conc,95, 2), 'Color', [255, 190, 106]./255, 'LineWidth',2)
xlim([0  52])
ylim([1 1e3])
set(gca, 'YScale', 'log')
xlabel('Time (weeks)')
ylabel('Durvalumab serum concentration (\mug/mL)')
set(gca,'FontSize',14)
