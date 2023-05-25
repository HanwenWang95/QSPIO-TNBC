% Cancer Module
%
% Models cancer cell growth under baseline conditions assuming logistic
% growth and first-order death
%
% Inputs: model        -- simbio model object with tumour compartment, T
%         species_name -- name of cancer cells [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - k_C_growth--cancer growth rate
%                         - C_max--cancer cell capacity
%                         - k_C_death--cancer death rate
%
% Outputs: model -- SimBiology model object with new cancer module


function model = cancer_module(model,species_name,params,varargin)

in = inputParser;
addOptional(in,'initial_amount',1); % initial cancer cell number
addOptional(in,'model_type','Logistic',@(s)ischar(s)); % tumor growth model type
parse(in,varargin{:});
init = in.Results.initial_amount;
model_type = in.Results.model_type;

first_call = true;
% Add Species
C = addspecies(model.Compartment(3),'C',init,'InitialAmountUnits','cell');
    set(C,'Notes','Number of cancer cells in tumour');

% Add Parameters
% Growth
k_C_growth = addparameter(model,'k_C_growth',params.k_C_growth.Value,'ValueUnits',params.k_C_growth.Units,'ConstantValue',false);
    set(k_C_growth,'Notes',['Cancer cell growth rate ' params.k_C_growth.Notes]);
try % only add C_max and vaculature parameter/species once
C_max = addparameter(model,'C_max',params.C_max.Value,'ValueUnits',params.C_max.Units,'ConstantValue',false);
    set(C_max,'Notes',['Cancer cell capacity ' params.C_max.Notes]);
catch
first_call = false;
end
% Death
k_C_death = addparameter(model,'k_C_death',params.k_C_death.Value,'ValueUnits',params.k_C_death.Units,'ConstantValue',false);
    set(k_C_death,'Notes',['Cancer cell death rate from innate immune cells ' params.k_C_death.Notes]);
% Therapy
param = addparameter(model,['k_' species_name '_therapy'],0,'ValueUnits','1/day','ConstantValue',false);
    set(param,'Notes',['Rate of ' species_name ' killing by therapy (see Rules)']);
addrule(model,['k_' species_name '_therapy = 0/day'],'repeatedAssignment');

% Initial Tumour Diameter
try % only add once
    p = addparameter(model,'initial_tumour_diameter',params.initial_tumour_diameter.Value,'ValueUnits',params.initial_tumour_diameter.Units);
        set(p,'Notes','Pre-treatment tumor diameter');
catch
end


% Add Reactions
% Growth
if strcmp(model_type,'Gompertzian')
reaction = addreaction(model,'null -> V_T.C');
    set(reaction,'ReactionRate','k_C_growth*V_T.C*log(max(C_max/(C_total+cell), 1))');
    set(reaction,'Notes','Cancer cell growth');
else
reaction = addreaction(model,'null -> V_T.C');
    set(reaction,'ReactionRate','k_C_growth*V_T.C*(1-C_total/C_max)');
    set(reaction,'Notes','Cancer cell growth');
end
% Death
reaction = addreaction(model,'V_T.C -> V_T.C_x');
    set(reaction,'ReactionRate','k_C_death*V_T.C');
    set(reaction,'Notes','Cancer cell death');

if (first_call) & strcmp(model_type,'Gompertzian')

% Vasculature species
s = addspecies(model.Compartment(3),'K',params.K0.Value,'InitialAmountUnits',params.K0.Units); % 2.4e8
    set(s,'Notes',['Maximal tumor capacity ' params.K0.Notes]);
s = addspecies(model.Compartment(3),'c_vas',0,'InitialAmountUnits','picogram/milliliter');
    set(s,'Notes','Angiogenic factors ');
  % Vasculature parameter
p = addparameter(model,'k_K_g',params.k_K_g.Value,'ValueUnits',params.k_K_g.Units,'ConstantValue',false); % 5.33
    set(p,'Notes',['Tumour vasculature growth rate ' params.k_K_g.Notes]);
p = addparameter(model,'k_K_d',params.k_K_d.Value,'ValueUnits',params.k_K_d.Units,'ConstantValue',false); % 7.9e-3
    set(p,'Notes',['Tumour vasculature inhibition rate ' params.k_K_d.Notes]);
p = addparameter(model,'k_vas_Csec',params.k_vas_Csec.Value,'ValueUnits',params.k_vas_Csec.Units,'ConstantValue',false);
    set(p,'Notes',['Secretion rate of angiogenic factors by cancer cells ' params.k_vas_Csec.Notes]);
p = addparameter(model,'k_vas_deg',params.k_vas_deg.Value,'ValueUnits',params.k_vas_deg.Units,'ConstantValue',false);
    set(p,'Notes',['Degradation rate of angiogenic factors ' params.k_vas_deg.Notes]);
p = addparameter(model,'c_vas_50',params.c_vas_50.Value,'ValueUnits',params.c_vas_50.Units,'ConstantValue',false);
    set(p,'Notes',['Half-maximal conc. of angiogenic factor on tumor capacity growth ' params.c_vas_50.Notes]);

addrule(model,'C_max = V_T.K','repeatedAssignment');

reaction = addreaction(model,'null -> V_T.c_vas');
    set(reaction,'ReactionRate','k_vas_Csec*C_total');
    set(reaction,'Notes','Secretion of angiogenic factors by cancer cells');
reaction = addreaction(model,'V_T.c_vas -> null');
    set(reaction,'ReactionRate','k_vas_deg*V_T.c_vas');
    set(reaction,'Notes','Degradation of tumor angiogenic factors');
reaction = addreaction(model,'null -> V_T.K');
    set(reaction,'ReactionRate','k_K_g*C_total*V_T.c_vas/(V_T.c_vas+c_vas_50)');
    set(reaction,'Notes','Growth of tumor carrying capacity');
reaction = addreaction(model,'V_T.K -> null');
    set(reaction,'ReactionRate','k_K_d*V_T.K*(nthroot(C_total/cell*2.57e-6,3))^2'); % 2.57e-6 mm^3/cell is the volume of a cancer cell
    set(reaction,'Notes','Endogenous inhibition of previously generated vasculature');

addevent(model,'C_total < 0.5*cell','V_T.K = 0.01*cell');
end

% Tumour Eradication
addevent(model,'V_T.C < 0.5*cell','V_T.C = 0.01*cell');

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update Total Number of Cancer Cells (Rule 2)
total_cancer_rule = model_rules(2);
rule = get(total_cancer_rule,'Rule');
set(total_cancer_rule,'Rule',[rule '+V_T.' species_name]);

% Rename Objects with 'species_name'
rename(C,species_name);
rename(k_C_growth,['k_' species_name '_growth']);
rename(k_C_death,['k_' species_name '_death']);
