% Function to generate object with model parameters to include in parameter
% sensitivity analysis
%
% Input:  model  -- Simbiology model object
%
% Output: params -- object containing parameters
%                   -> for each parameter:
%                       - Adds the name of the parameter to the list
%                       - defines upper and lower bounds for uniform and
%                       loguniform
%                       - defines median and sigma for normal and log
%                       normal
%                       - specifies the Sampling technique choose from:
%                           - uniform
%                           - loguniform
%                           - normal
%                           - lognormal
%         Examples:
%             % k1
%             params.names = [params.names; 'k1'];
%             params.k1.UpperBound = 1;
%             params.k1.LowerBound = 0;
%             params.k1.Sampling   = 'uniform';
%             params.k1.ScreenName = 'k1 binding rate';
%
%             % k2
%             params.names = [params.names; 'k2'];
%             params.k2.UpperBound = 1;
%             params.k2.LowerBound = 0;
%             params.k2.Sampling   = 'loguniform';
%             params.k2.ScreenName = 'k2 binding rate';
%
%             % k3
%             params.names = [params.names; 'k3'];
%             params.k3.Median     = 1;
%             params.k3.Sigma      = 1;
%             params.k3.Sampling   = 'normal';
%             params.k3.ScreenName = 'k3 binding rate';
%
%             % k4
%             params.names = [params.names; 'k4'];
%             params.k4.Median     = 1;
%             params.k4.Sigma      = 1;
%             params.k4.Sampling   = 'lognormal';
%             params.k4.ScreenName = 'k4 binding rate';


function params = PSA_param_in_all(model)

percentage = 5;
params.names = {};

for i = 1:length(model.Parameters)
    if (model.Parameters(i).ConstantValue)
        if ~((model.Parameters(i).Name(1)=='H')||(model.Parameters(i).Name(1)=='a')||strcmp(model.Parameters(i).Name,'m')||...
              strcmp(model.Parameters(i).Name,'k')||contains(model.Parameters(i).Name,'tot')||strcmp(model.Parameters(i).Name,'cell')||...
              strcmp(model.Parameters(i).Name,'day')||strcmp(model.Parameters(i).Name,'Treg_'))
            params.names = [params.names; model.Parameters(i).Name];
            params.(model.Parameters(i).Name).UpperBound = model.Parameters(i).Value *(100+percentage)/100;
            params.(model.Parameters(i).Name).LowerBound = model.Parameters(i).Value *(100-percentage)/100;
            params.(model.Parameters(i).Name).Sampling   = 'uniform';
            params.(model.Parameters(i).Name).ScreenName = strrep(model.Parameters(i).Name,'_','\_');
        end
    end
end

% % Initial tumor diameter
% params.names = [params.names; 'initial_tumour_diameter'];
% params.initial_tumour_diameter.UpperBound = 5;
% params.initial_tumour_diameter.LowerBound = 2;
% params.initial_tumour_diameter.Sampling   = 'uniform';
% params.initial_tumour_diameter.ScreenName = 'Initial Tumor Diameter';
