% Function to export simbiology model components into excel files
%
% Inputs: model -- SimBiology model object



function simbio_export(model)
m1 = copyobj(model);

for i = 1:size(m1.Compartments,1)
    Compartments(i,1) = {m1.Compartments(i).Name};
    Compartments(i,2) = {m1.Compartments(i).Capacity};
    Compartments(i,3) = {m1.Compartments(i).CapacityUnits};
    Compartments(i,5) = {m1.Compartments(i).Notes};
end

for i = 1:size(m1.Species,1)
    Species(i,1) = {m1.Species(i).Name};
    Species(i,2) = {m1.Species(i).InitialAmount};
    Species(i,3) = {m1.Species(i).InitialAmountUnits};
    Species(i,4) = {m1.Species(i).Parent.Name};
    Species(i,5) = {m1.Species(i).Notes};
end

for i = 1:size(m1.Parameters,1)
    Parameters(i,1) = {m1.Parameters(i).Name};
    Parameters(i,2) = {m1.Parameters(i).Value};
    Parameters(i,3) = {m1.Parameters(i).ValueUnits};
    Parameters(i,4) = {m1.Parameters(i).Notes};
end

for i = 1:size(m1.Reaction,1)
    Reactions(i,1) = {i};
    Reactions(i,2) = {m1.reaction(i).Reaction};
    Reactions(i,3) = {m1.reaction(i).ReactionRate};
    Reactions(i,4) = {m1.reaction(i).Notes};
end

for i = 1:size(m1.Rules,1)
    Rules(i,1) = {m1.Rules(i).Name};
    Rules(i,2) = {m1.Rules(i).Rule};
    Rules(i,3) = {m1.Rules(i).RuleType};
    Rules(i,4) = {m1.Rules(i).Notes};
end
k = 1;
for i = 1:size(m1.Events,1)
    for j = 1:size(m1.Events(i).EventFcns , 1)
        if j == 1
            Events(k,1) = {i};
            Events(k,2) = {m1.Events(i).Trigger};
            Events(k,3) =  m1.Events(i).EventFcns(j);
            Events(k,4) = {m1.Events(i).Notes};
        else
            Events(k,3) =  m1.Events(i).EventFcns(j);
        end
        k = k+1;
    end
end


T1 = array2table(Compartments, 'VariableNames',{'Name','Capacity','Unit','na','Description'});
writetable(T1, 'compartment.xlsx')
T2 = array2table(Species, 'VariableNames',{'Name','InitialAmount','Unit','Location','Description'});
writetable(T2, 'Species.xlsx')
T3 = array2table(Reactions, 'VariableNames',{'Number','Reaction','Rate','Description'});
writetable(T3, 'Reactions.xlsx')
T4 = array2table(Parameters, 'VariableNames',{'Name','Value','Unit','Description'});
writetable(T4, 'Parameters.xlsx')
T5 = array2table(Rules, 'VariableNames',{'na','Rule','Type','Des'});
writetable(T5, 'Rules.xlsx')
T6 = array2table(Events, 'VariableNames',{'Number','Trigger','Event','Des'});
writetable(T6, 'Events.xlsx')


%% Write to excel files individually and all-in-one
% filenameAll          = [ name,'_all.xlsx'              ];
% 
% writetable(T1, filenameAll, 'Sheet', 'Compartments');
% writetable(T2, filenameAll, 'Sheet', 'Species'     );
% writetable(T3, filenameAll, 'Sheet', 'Parameters'  );
% writetable(T4, filenameAll, 'Sheet', 'Reactions'   );
% writetable(T5, filenameAll, 'Sheet', 'Rules'       );
% writetable(T6, filenameAll, 'Sheet', 'Events'      );

disp('Model contents exported successfully')
end
