% Function to set Initial Conditions
%
% Inputs: model_in -- SimBiology model object
%         ICs      -- vector of IC values
%
% Outputs: model_out -- SimBiology model object with ICs


function model_out = set_ICs(model_in,ICs)

model_out = copyobj(model_in);

for i = 1:length(model_in.Species)
  model_out.Species(i).InitialAmount = ICs(i);
end
