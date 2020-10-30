% Function to find the number of MHC considered
%
% Inputs: simDataPSA     -- Object containing the whole simbiology model
%                           outputs for all batch simulations
%
% Outputs: numClones     -- number of MHC in the simulation


function numMHC = howManyMHC(simData)

% Finds the number of clones in the simulation by parsing 'Ti's
for i = 1:length(simData.DataNames)
    temp = sscanf((simData.DataNames{i}),'M%d');
    if ~isempty(temp)
        MHC(i) = temp;
    else
        MHC(i) = NaN;
    end
end
numMHC = nanmax(MHC);
