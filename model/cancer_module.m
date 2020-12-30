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
parse(in,varargin{:});
init = in.Results.initial_amount;

% Add Species
C = addspecies(model.Compartment(3),'C',init,'InitialAmountUnits','cell');
    set(C,'Notes','Number of cancer cells in tumour');

% Add Parameters
% Growth
k_C_growth = addparameter(model,'k_C_growth',params.k_C_growth.Value,'ValueUnits',params.k_C_growth.Units,'ConstantValue',false);
    set(k_C_growth,'Notes',['Cancer cell growth rate ' params.k_C_growth.Notes]);
try % only add C_max once
C_max = addparameter(model,'C_max',params.C_max.Value,'ValueUnits',params.C_max.Units,'ConstantValue',false);
    set(C_max,'Notes',['Cancer cell capacity ' params.C_max.Notes]);
catch
end
% Death
k_C_death = addparameter(model,'k_C_death',params.k_C_death.Value,'ValueUnits',params.k_C_death.Units,'ConstantValue',false);
    set(k_C_death,'Notes',['Cancer cell death rate from innate immune cells ' params.k_C_death.Notes]);
% Therapy
param = addparameter(model,['k_' species_name '_therapy'],0,'ValueUnits','1/day','ConstantValue',false);
    set(param,'Notes',['Rate of ' species_name ' killing by therapy (see Rules)']);

% Initial Tumour Diameter
try % only add once
    p = addparameter(model,'initial_tumour_diameter',params.initial_tumour_diameter.Value,'ValueUnits',params.initial_tumour_diameter.Units);
        set(p,'Notes','Pre-treatment tumor diameter');
catch
end


% Add Reactions
% Growth
reaction = addreaction(model,'null -> V_T.C');
    set(reaction,'ReactionRate','k_C_growth*V_T.C*log(max(C_max/(C_total+cell), 1))');
    set(reaction,'Notes','Cancer cell growth');
% Death
reaction = addreaction(model,'V_T.C -> V_T.C_x');
    set(reaction,'ReactionRate','k_C_death*V_T.C');
    set(reaction,'Notes','Cancer cell death');

% Tumour Erradication
addevent(model,'V_T.C < 0.5*cell','V_T.C = 0.01*cell');

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update Tumour Volume (Rule 1)
%volume_rule = model_rules(1);
%rule = get(volume_rule,'Rule');
%set(volume_rule,'Rule',[rule '+1/rho_cell*V_T.' species_name]);

% Update Total Number of Cancer Cells (Rule 2)
total_cancer_rule = model_rules(2);
rule = get(total_cancer_rule,'Rule');
set(total_cancer_rule,'Rule',[rule '+' species_name]);

% Rename Objects with 'species_name'
rename(C,species_name);
rename(k_C_growth,['k_' species_name '_growth']);
rename(k_C_death,['k_' species_name '_death']);
