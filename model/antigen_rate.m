% Function to generate object with default model parameters
%
% Inputs: model   -- SimBiology model object
%         ID      -- T cell-antigen ID number
%         cancers -- cell array of string containing cancer names
%         P_C     -- parameter array (same length as 'cancers') of antigen concentrations per cancer cell
%         nTcells -- number of cytotoxic T cell clones in model
%
% Output: rate    -- string containing rate of antigen release from dying cancer cells


function rate = antigen_rate(model,ID,cancers,P_C,nTcells)

% Parse Inputs
P = ['P' ID];
C = cancers;

% Add Antigen Concentrations to Model
for i = 1:length(P_C)
    param = addparameter(model,[P '_' C{i}],P_C(i).Value,'ValueUnits',P_C(i).Units);
        set(param,'Notes',['Concentration of ' P ' in ' C{i} ...
        ' (estimated by 2.7e6 proteins/um^3 in a cancer cell (Milo 2013, PMID: 24114984) * cancer cell volume / Avogadro constant)']);
end

rate = '';
for i = 1:length(C)
    kT = '';
    for j = 1:nTcells
        if strcmp(kT,'')
            % kT = kT + k_Tj*Tj
            kT = sprintf('k_C_T%d*V_T.T%d',j,j);
        else
            kT = sprintf('%s+k_C_T%d*V_T.T%d',kT,j,j);
        end
    end
    % rate = rate + P_Ci*(k_Ci_death+k_Ci_therapy+kT/(C_total+T_total)*(1-H_PD1))*Ci
    % rate = sprintf('%s+%s_%s*(k_%s_death+k_%s_therapy+(%s)/(C_total+T_total+cell)*(1-H_TGF_CTL)*(1-H_PD1_C1)*(1-H_MDSC_C1))*%s',...
    %       rate,P,C{i},C{i},C{i},kT,C{i});
    for k = 1:length(model.reaction)
        if strcmp(model.reaction(k).reaction, ['V_T.' C{i} ' -> V_T.C_x']) ...
            && ~isempty(strfind(model.reaction(k).ReactionRate, 'k_C_T'))
            R_T = [model.reaction(k).ReactionRate, '*(1-H_MDSC_C1)'];
        end
    end
    if strcmp(rate,'')
        % rate = rate + P_Ci*(k_Ci_death+k_Ci_therapy)*Ci + P_Ci*R_Ti
        rate = sprintf('%s_%s*(k_%s_death+k_%s_therapy)*%s+%s_%s*%s',...
            P,C{i},C{i},C{i},C{i},P,C{i},R_T);
    else
        rate = sprintf('%s+%s_%s*(k_%s_death+k_%s_therapy)*%s+%s_%s*%s',...
            rate,P,C{i},C{i},C{i},C{i},P,C{i},R_T);
    end
end
rate = ['n_T' ID '_clones*(' rate ')'];
