function dose_schedule = schedule_dosing_qsp(drugName,varargin)

% Check if drugName is cell array
if (iscell(drugName))
    N = length(drugName);
else
    N = 1;
    drugName = {drugName};
end

% Optional Inputs
in = inputParser;
% Anti-CD47
addParameter(in,'aCD47_dose',1); % 1,3,10 mg/kg Q1W, 30 mg/kg Q2W
addParameter(in,'aCD47_schedule',[0,14,0]);
% Pembrolizumab
addParameter(in,'pembrolizumab_dose',2); % 2 mg/kg every 3 weeks
addParameter(in,'pembrolizumab_schedule',[0,21,0]);
% Nivolumab
addParameter(in,'nivolumab_dose',3); % 3 mg/kg every two weeks
addParameter(in,'nivolumab_schedule',[0,14,0]);
% Atezolizumab
addParameter(in,'atezolizumab_dose',15); % 15 mg/kg every three weeks
addParameter(in,'atezolizumab_schedule',[0,21,0]);
% Ipilimumab
addParameter(in,'ipilimumab_dose',1); % 1 mg/kg every three weeks
addParameter(in,'ipilimumab_schedule',[0,21,0]);
% Durvalumab
addParameter(in,'durvalumab_dose',10); % 20 mg/kg every 4 weeks
addParameter(in,'durvalumab_schedule',[0,14,13]);

% Patient Weight
addParameter(in,'patientWeight',75);

% Parse Inputs
parse(in,varargin{:});

% Anti-CD47
dose_aCD47 = in.Results.aCD47_dose;
schedule_aCD47 = in.Results.aCD47_schedule;
% Pembrolizumab
dose_pembro = in.Results.pembrolizumab_dose;
schedule_pembro = in.Results.pembrolizumab_schedule;
% Nivolumab
dose_nivo = in.Results.nivolumab_dose;
schedule_nivo = in.Results.nivolumab_schedule;
% atezolizumab
dose_atezo = in.Results.atezolizumab_dose;
schedule_atezo = in.Results.atezolizumab_schedule;
% Ipilimumab
dose_ipil = in.Results.ipilimumab_dose;
schedule_ipil = in.Results.ipilimumab_schedule;
% Durvalumab
dose_durv = in.Results.durvalumab_dose;
schedule_durv = in.Results.durvalumab_schedule;

% Patient Weight
patient_weight = in.Results.patientWeight;

% Anti-CD47
MW_aCD47 = 7.6E7; % milligrams per mole
doseObj_aCD47 = sbiodose('aCD47','Amount',patient_weight*dose_aCD47/MW_aCD47,'AmountUnits','mole','TargetName','V_C.aCD47');
doseObj_aCD47.Rate = doseObj_aCD47.Amount;
doseObj_aCD47.RateUnits = 'mole/hour';
doseObj_aCD47.StartTime = schedule_aCD47(1);
doseObj_aCD47.Interval = schedule_aCD47(2);
doseObj_aCD47.TimeUnits = 'day';
doseObj_aCD47.RepeatCount = schedule_aCD47(3);
doseObj_aCD47.Active = true;

% Pembrolizumab
MW_pembro = 1.49E8; % milligrams per mole
doseObj_pembro = sbiodose('A','Amount',76.5*dose_pembro/MW_pembro,'AmountUnits','mole','TargetName','V_C.A');
doseObj_pembro.StartTime = schedule_pembro(1);
doseObj_pembro.Interval = schedule_pembro(2);
doseObj_pembro.TimeUnits = 'day';
doseObj_pembro.RepeatCount = schedule_pembro(3);
doseObj_pembro.Active = true;

% Nivolumab
MW_nivo = 1.436E8; % milligrams per mole
doseObj_nivo = sbiodose('A','Amount',patient_weight*dose_nivo/MW_nivo,'AmountUnits','mole','TargetName','V_C.A');
doseObj_nivo.StartTime = schedule_nivo(1);
doseObj_nivo.Interval = schedule_nivo(2);
doseObj_nivo.TimeUnits = 'day';
doseObj_nivo.RepeatCount = schedule_nivo(3);
doseObj_nivo.Active = true;

% atezolizumab
MW_atezo = 1.436E8; % milligrams per mole
doseObj_atezo = sbiodose('A','Amount',patient_weight*dose_atezo/MW_atezo,'AmountUnits','mole','TargetName','V_C.A');
doseObj_atezo.StartTime = schedule_atezo(1);
doseObj_atezo.Interval = schedule_atezo(2);
doseObj_atezo.TimeUnits = 'day';
doseObj_atezo.RepeatCount = schedule_atezo(3);
doseObj_atezo.Active = true;

% Ipilimumab
MW_ipil = 1.486349E8; % milligrams per mole
doseObj_ipi = sbiodose('A','Amount',patient_weight*dose_ipil/MW_ipil,'AmountUnits','mole','TargetName','V_C.A');
doseObj_ipi.StartTime = schedule_ipil(1);
doseObj_ipi.Interval = schedule_ipil(2);
doseObj_ipi.TimeUnits = 'day';
doseObj_ipi.RepeatCount = schedule_ipil(3);
doseObj_ipi.Active = true;

% Durvalumab
MW_durv = 1.463E8; % milligrams per mole
doseObj_durv = sbiodose('A','Amount',patient_weight*dose_durv/MW_durv,'AmountUnits','mole','TargetName','V_C.A');
doseObj_durv.StartTime = schedule_durv(1);
doseObj_durv.Interval = schedule_durv(2);
doseObj_durv.TimeUnits = 'day';
doseObj_durv.RepeatCount = schedule_durv(3);
doseObj_durv.Active = true;

% Dose Schedule Array
dose_schedule(N) = sbiodose('empty'); % preallocate array
for i = 1:N
    switch drugName{i}
        case 'aCD47'
            dose_schedule(i) = doseObj_aCD47;
        case 'pembrolizumab'
            dose_schedule(i) = doseObj_pembro;
        case 'nivolumab'
            dose_schedule(i) = doseObj_nivo;
        case 'atezolizumab'
            dose_schedule(i) = doseObj_atezo;
        case 'ipilimumab'
            dose_schedule(i) = doseObj_ipi;
        case 'durvalumab'
            dose_schedule(i) = doseObj_durv;
        otherwise
            error('No match for drug name');
    end
end
