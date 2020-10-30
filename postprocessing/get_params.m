% Function to get parameter values from the SimBiology Model
%
% Inputs: simData -- SimBiology data object
%         name    -- Parameter name
%
% Outputs: value  -- Parameter value


function value = get_params(model, name)

for k = 1:length(model.parameters)
    if strcmp(model.parameters(k).Name, name)
        value = model.parameters(k).Value;
    end
end
