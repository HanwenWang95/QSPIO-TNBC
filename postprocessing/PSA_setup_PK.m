% Function to generate object with random variation in the parameters
% of interest to be varied in parameter snesitivity analysis
%
% Inputs: model        -- simbio model object
%         params       -- object containing model parameter Values, Median, LowerBound, UpperBound, type of distribution:
%         n_PSA        -- number of samples
%
% Outputs: params_in -- sampled parameter values for params and n_PSA



function params_in = PSA_setup_PK(params,n_PSA,names,psi,V,opt_avg)

params_in = params;

for i = 1:length(names.pars)
    p(i,:) = min(psi(i,:)) + (max(psi(i,:))-min(psi(i,:))) * rand(1,n_PSA);
end
PK_pars = exp(V*p + opt_avg);

for i = 1:length(names.pars)

    params_in.names = [params_in.names; names.pars{i}];
    params_in.(names.pars{i}).ScreenName = names.labels{i};
    params_in.(names.pars{i}).LHS = PK_pars(i,:)';
    params_in.all(:,end+1) = params_in.(names.pars{i}).LHS;

end
