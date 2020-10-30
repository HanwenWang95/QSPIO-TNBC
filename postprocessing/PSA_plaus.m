% Function to find plausible simulations
%
% Inputs: simDataPSAprep   -- Array containing the simbiology model outputs
%                             for all batch simulations
%         params_out       -- object containing model outputs to be organized
%                             for future sensitivity analysis
%
% Outputs: iPlaus          -- Index of plausible simulations



function iPlaus = PSA_plaus(simDataPSAprep,params_out)


% check the lower and upperbounds of the outputs and replace outofbound
% values with NaN
for i = 1:length(params_out.names)
    index = find( ...
        simDataPSAprep(:,i) > params_out.(params_out.names{i}).UpperBound | ...
        simDataPSAprep(:,i) < params_out.(params_out.names{i}).LowerBound ...
        );
    if ~isempty(index)
        simDataPSAprep(index,i) = NaN;
    end
end

% find the index of valid outputs and return
iPlaus = find(~any(isnan(simDataPSAprep),2));
