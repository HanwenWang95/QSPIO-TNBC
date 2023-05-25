%% Load model
immune_oncology_model_NSCLC

%% Random sampling
params_in  = PSA_param_in_NSCLC_VPop;
k = length(params_in.names);
M = 1e3; r = 40;
p = 11;
% delta = p/(2*(p-1));
delta = 1/(p-1);

% Generate r trajectories
rng(123)
R = randi([1 p-2],k,M)./(p-1);
B = [zeros(1,k); tril(ones(k))];
z = 1;
B0_M = nan(k+1,k);
for i = 1:M
    x0 = R(:,i);
    D0 = eye(k) .* (binornd(1,0.5,k,1).*2-1);
    P0 = eye(k);
    P0 = P0(randperm(k),:);
    
    B0_temp = ones(k+1,1)*x0' + (delta/2).*(ones(k+1,k) + (2.*B-ones(k+1,k))*D0) * P0;
    if sum(B0_temp<=0 | B0_temp>=1, 'all')==0
        B0_M(:,:,z) = B0_temp;
        z = z + 1;
    end
end

% Calculate spread 
% temp = nan(k+1); d = zeros(size(B0_M,3),1);
% for z = 2:size(B0_M,3)
%     for i = 1:k+1
%         for j = 1:k+1
%             temp(i,j) = sqrt(sum((B0_M(i,:,z)-B0_M(j,:,1)).^2));
%         end
%     end
%     d(z) = sum(temp,'all');
% end
% [~,idx_sorted] = sort(d,'descend');
% idx_sorted = sort(idx_sorted(1:r));
% B0 = B0_M(:,:,idx_sorted);

B0 = B0_M;

%% Generate model output
dose_schedule = schedule_dosing({'durvalumab'});

EE = nan(size(B0,3),k);
EE0 = nan(size(B0,3),k);
D_T = nan(k+1,size(B0,3));
for j = 1:size(B0,3)
    lhs = B0(:,:,j);
    params  = PSA_param_in_NSCLC_VPop;
    params_in = params;
    % Calculate QSP parameter values
    for i = 1:length(params.names)
        if     strcmp(params.(params.names{i}).Sampling , 'uniform')
            LB = params.(params.names{i}).LowerBound;
            UB = params.(params.names{i}).UpperBound;
            params_in.(params_in.names{i}).LHS = LB + (UB-LB) * lhs(:,i);
            params_in.all(:,i) = params_in.(params_in.names{i}).LHS;
            
        elseif strcmp(params.(params.names{i}).Sampling , 'loguniform')
            LB = params.(params.names{i}).LowerBound;
            UB = params.(params.names{i}).UpperBound;
            params_in.(params.names{i}).LHS = 10.^( log10(LB) + (log10(UB)-log10(LB)) * lhs(:,i) );
            params_in.all(:,i) = params_in.(params.names{i}).LHS;
            
        elseif strcmp(params.(params.names{i}).Sampling , 'normal')
            Median = params.(params.names{i}).Median;
            Sigma  = params.(params.names{i}).Sigma;
            params_in.(params.names{i}).LHS = icdf('Normal',lhs(:,i),Median,Sigma);
            params_in.all(:,i) = params_in.(params.names{i}).LHS;
            
        elseif strcmp(params.(params.names{i}).Sampling , 'lognormal')
            Median = params.(params.names{i}).Median;
            Sigma  = params.(params.names{i}).Sigma;
            params_in.(params.names{i}).LHS = icdf('Lognormal',lhs(:,i),Median,Sigma);
            params_in.all(:,i) = params_in.(params.names{i}).LHS;
        end
    end
    
    % Simulation
    n_PSA = length(params_in.(params_in.names{1}).LHS);
    for i = 1:n_PSA
        display(['Sample ',num2str(i+(j-1)*n_PSA),'/',num2str(n_PSA*size(B0,3))]);
        % Set the new parameters
        if i > 1
            delete(model_PSA)
        end
        model_PSA = copyobj(model);
        variantObj = addvariant(model_PSA, ['v',num2str(i+(j-1)*n_PSA,'%5.5i')]);
        for p = 1:length(params_in.names)
            if ~isempty(sbioselect (model, 'Type', 'parameter', 'Name', params_in.names{p}))
                addcontent(variantObj, {'parameter', params_in.names{p}, 'Value', params_in.(params_in.names{p}).LHS(i)});
            elseif ~isempty(sbioselect (model, 'Type', 'compartment', 'Name', params_in.names{p}))
                addcontent(variantObj, {'compartment', params_in.names{p}, 'Capacity', params_in.(params_in.names{p}).LHS(i)});
            else
                disp(['Unable to identify parameter/compartment named ', params_in.names{p}, ' in the model'])
            end
        end
        % Set Initial Conditions
        [model_PSA,success,~] = initial_conditions(model_PSA,'Variant',variantObj);
        params_out.patient(i+(j-1)*n_PSA) = double(success);
        % Run the model with drugs
        if (success)
            try
                simData = sbiosimulate(model_PSA,[],variantObj,dose_schedule);
                % Remove simulations that reached the time limit
                if size(simData.Data,1) < length(time)
                    simData = [];
                    params_out.patient(i+(j-1)*n_PSA) = -1;
                    disp('Simulation takes longer than the preset time limit');
                end
            catch
                disp('Integration tolerance not met');
                simData = [];
                params_out.patient(i+(j-1)*n_PSA) = -1;
            end
        else
            disp('Initial conditions not reached');
            simData = [];
        end
        % Save model output struture within an array of structures
        if isempty(simData)
            D_T(i,j) = nan;
        else
            [~,temp,~] = selectbyname(simData, 'V_T');
            D_T(i,j) = (6*temp(end)/pi)^(1/3);
        end
        % Calculate EE
        if i>1
            step = B0(i,:,j)-B0(i-1,:,j);
            idx_par = find(step~=0);
            EE(j,idx_par) = (D_T(i,j)-D_T(i-1,j))/step(idx_par);
            EE0(j,idx_par) = abs((D_T(i,j)-D_T(i-1,j))/step(idx_par));
            disp(['Trajectory ' num2str(j) ': EE = ' num2str(EE(j,idx_par)) ' for par #' num2str(idx_par)])
        end
    end
end

%% Visualization
mu = nanmean(EE);
mu0= nanmean(EE0);
variance = nanvar(EE);
[mu0_sorted,idx_sorted] = sort(mu0);
labels_sorted = params_in.names(idx_sorted);
variance_sorted = variance(idx_sorted);

figure
scatter(mu0,variance)
hold on
for i=2:2:length(mu)
    text(mu0_sorted(i)*1.05,variance_sorted(i)*1.03,labels_sorted{i}, 'Interpreter', 'none')
end
for i=1:2:length(mu)
    text(mu0_sorted(i)*(1-length(labels_sorted{i})/30), variance_sorted(i)*1.03,labels_sorted{i}, 'Interpreter', 'none')
end
set(gca,'xscale','log')
set(gca,'yscale','log')

xlabel('\mu*')
ylabel('\sigma^2')
title('Morris Sensitivity Analysis')

set(gca,'FontSize',14,'FontName','Arial')



