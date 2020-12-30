% Dose Schedule
%
% Generates Donsing Schedule
%
% Inputs: drugName -- drug name or character array of drug names
%         varargin -- Name-Value Pairs
%                     - drugName_dose (mg/kg)
%                     - drugName_schedule [start,interval,repeat]
%                     - patientWeight (kg)
%
%       vaild drugName values: nivolumab, atezolizumab, ipilimumab
%
% Outputs: dosing -- SimBiology model object with new antigen module


function dose_schedule = schedule_dosing(drugName,varargin)

% Check if drugName is cell array
if (iscell(drugName))
    N = length(drugName);
else
    N = 1;
    drugName = {drugName};
end

% Optional Inputs
in = inputParser;
% Nivolumab
addParameter(in,'nivolumab_dose',3); % 3 mg/kg every two weeks
addParameter(in,'nivolumab_schedule',[0,14,30]);
% Atezolizumab
addParameter(in,'atezolizumab_dose',15); % 15 mg/kg every three weeks
addParameter(in,'atezolizumab_schedule',[0,21,30]);
% Ipilimumab
addParameter(in,'ipilimumab_dose',1); % 1 mg/kg every three weeks
addParameter(in,'ipilimumab_schedule',[0,21,30]);
% Entinostat
addParameter(in,'entinostat_dose',3); % 3 mg every week
addParameter(in,'entinostat_schedule',[0,7,55]);
% NabPaclitaxel
addParameter(in,'nabPaclitaxel_dose',100); % 100 mg/m2 Q3/4W
addParameter(in,'nabPaclitaxel_schedule',[0,28,15]);
% Patient Weight
addParameter(in,'patientWeight',80);
% Patient Body Surface area
addParameter(in,'patientBSA',1.9);
% Parse Inputs
parse(in,varargin{:});
% Nivolumab
dose_nivo = in.Results.nivolumab_dose;
schedule_nivo = in.Results.nivolumab_schedule;
% atezolizumab
dose_atezo = in.Results.atezolizumab_dose;
schedule_atezo = in.Results.atezolizumab_schedule;
% Ipilimumab
dose_ipil = in.Results.ipilimumab_dose;
schedule_ipil = in.Results.ipilimumab_schedule;
% Entinostat
dose_ENT = in.Results.entinostat_dose;
schedule_ENT = in.Results.entinostat_schedule;
% NabPaclitaxel
dose_nabp = in.Results.nabPaclitaxel_dose;
schedule_nabp = in.Results.nabPaclitaxel_schedule;
% Patient Weight
patient_weight = in.Results.patientWeight;
% Patient BSA
patient_BSA = in.Results.patientBSA;

% Nivolumab
MW_nivo = 1.436E8; % milligrams per mole
doseObj_nivo = sbiodose('nivo','Amount',patient_weight*dose_nivo/MW_nivo,'AmountUnits','mole','TargetName','V_C.nivo');
doseObj_nivo.StartTime = schedule_nivo(1);
doseObj_nivo.Interval = schedule_nivo(2);
doseObj_nivo.TimeUnits = 'day';
doseObj_nivo.RepeatCount = schedule_nivo(3);
doseObj_nivo.Active = true;

% Entinostat
MW_ENT = 3.764085e5; % mg/mole
Bio = 0.18; % Fraction of dose through buccal absorption

doseObj_ENT_1 = sbiodose('ENT_1','Amount',Bio*dose_ENT/MW_ENT,'AmountUnits','mole','TargetName','V_C.ENT_Buccal');
doseObj_ENT_1.StartTime = schedule_ENT(1);
doseObj_ENT_1.Interval = schedule_ENT(2);
doseObj_ENT_1.TimeUnits = 'day';
doseObj_ENT_1.RepeatCount = schedule_ENT(3);
doseObj_ENT_1.Active = true;
doseObj_ENT_1.DurationParameterName = 'durP';

doseObj_ENT_2 = sbiodose('ENT_2','Amount',(1-Bio)*dose_ENT/MW_ENT,'AmountUnits','mole','TargetName','V_C.Dose2');
doseObj_ENT_2.StartTime = schedule_ENT(1);
doseObj_ENT_2.Interval = schedule_ENT(2);
doseObj_ENT_2.TimeUnits = 'day';
doseObj_ENT_2.RepeatCount = schedule_ENT(3);
doseObj_ENT_2.Active = true;
doseObj_ENT_2.LagParameterName = 'lagP';

% nabPaclitaxel
% MW = 853.9; % g/mol
doseObj_nabp_1 = sbiodose('nabp_1','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_1.Rate = doseObj_nabp_1.Amount/30;
doseObj_nabp_1.RateUnits = 'milligram/minute';
doseObj_nabp_1.StartTime = schedule_nabp(1);
doseObj_nabp_1.Active = true;
doseObj_nabp_1.Interval = schedule_nabp(2);
doseObj_nabp_1.TimeUnits = 'day';
doseObj_nabp_1.RepeatCount = schedule_nabp(3);

doseObj_nabp_2 = sbiodose('nabp_2','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_2.Rate = doseObj_nabp_2.Amount/30;
doseObj_nabp_2.RateUnits = 'milligram/minute';
doseObj_nabp_2.StartTime = schedule_nabp(1)+7;
doseObj_nabp_2.Active = true;
doseObj_nabp_2.Interval = schedule_nabp(2);
doseObj_nabp_2.TimeUnits = 'day';
doseObj_nabp_2.RepeatCount = schedule_nabp(3);

doseObj_nabp_3 = sbiodose('nabp_3','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_3.Rate = doseObj_nabp_3.Amount/30;
doseObj_nabp_3.RateUnits = 'milligram/minute';
doseObj_nabp_3.StartTime = schedule_nabp(1)+14;
doseObj_nabp_3.Active = true;
doseObj_nabp_3.Interval = schedule_nabp(2);
doseObj_nabp_3.TimeUnits = 'day';
doseObj_nabp_3.RepeatCount = schedule_nabp(3);

% atezolizumab
MW_atezo = 1.436E8; % milligrams per mole
doseObj_atezo = sbiodose('atezo','Amount',patient_weight*dose_atezo/MW_atezo,'AmountUnits','mole','TargetName','V_C.atezo');
doseObj_atezo.StartTime = schedule_atezo(1);
doseObj_atezo.Interval = schedule_atezo(2);
doseObj_atezo.TimeUnits = 'day';
doseObj_atezo.RepeatCount = schedule_atezo(3);
doseObj_atezo.Active = true;

% Ipilimumab
MW_ipil = 1.486349E8; % milligrams per mole
doseObj_ipi = sbiodose('ipi','Amount',patient_weight*dose_ipil/MW_ipil,'AmountUnits','mole','TargetName','V_C.ipi');
doseObj_ipi.StartTime = schedule_ipil(1);
doseObj_ipi.Interval = schedule_ipil(2);
doseObj_ipi.TimeUnits = 'day';
doseObj_ipi.RepeatCount = schedule_ipil(3);
doseObj_ipi.Active = true;

% Dose Schedule Array
dose_schedule(N) = sbiodose('empty'); % preallocate array
for i = 1:N
    switch drugName{i}
        case 'nivolumab'
            dose_schedule(i) = doseObj_nivo;
        case 'atezolizumab'
            dose_schedule(i) = doseObj_atezo;
        case 'ipilimumab'
            dose_schedule(i) = doseObj_ipi;
        case 'nabPaclitaxel'
            dose_schedule(i) = doseObj_nabp_1;
        case 'entinostat'
            dose_schedule(i) = doseObj_ENT_1;
        otherwise
            error('No match for drug name');
    end
end

for i = 1:N
    switch drugName{i}
        case 'nabPaclitaxel'
            dose_schedule(end+1) = doseObj_nabp_2;
            dose_schedule(end+1) = doseObj_nabp_3;
        case 'entinostat'
            dose_schedule(end+1) = doseObj_ENT_2;
    end
end
