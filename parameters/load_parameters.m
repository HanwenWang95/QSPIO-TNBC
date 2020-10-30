% Function to load model parameters from file
%
% Inputs: params_in  -- object containing parameters
%
% Output: params_out -- object containing processed parameters
%                   -> each parameter has Value, Units and Notes properties.

function params_out = load_parameters(params_in)
% Initialize params_out object
names = fieldnames(params_in);
for i = 1:numel(fieldnames(params_in))
    params_out.(names{i}).Value = params_in.(names{i}).Value;
    params_out.(names{i}).Units = params_in.(names{i}).Units;
    params_out.(names{i}).Notes = params_in.(names{i}).Notes;
end

% Calculate parameters
for i = 1:numel(fieldnames(params_in))
    if isempty(params_in.(names{i}).Value)
        n = numel(params_in.(names{i}).Factors);
        for j = 1:n
            p(j) = params_out.(params_in.(names{i}).Factors{j});
        end
        params_out.(names{i}) = eval(params_in.(names{i}).Equation);
        params_out.(names{i}).Notes = params_in.(names{i}).Notes;
    end
end
