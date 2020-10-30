% Antigen Database
%
% Creates antigen data structure from database of antigens
%
% Inputs: cancers  -- cell array of strings of cancer names
%         P_C      -- vector of antigen concentrations per cancer type
%
% Optional Name-Value Pairs
%         antigenID -- antigen ID number
%         concentrationUnits -- concentration units
%         numberOfMHCs -- number of MHC molecules in model
%
% Outputs: antigen -- antigen structure


function antigen = create_antigen(cancers,P_C,varargin)

% Optional Inputs
par = inputParser;
addParameter(par,'antigenID',-1);
addParameter(par,'concentrationUnits','mole');
addParameter(par,'numberOfMHCs',1);
parse(par,varargin{:});
ID = par.Results.antigenID;
units = [par.Results.concentrationUnits '/cell'];
nMHCs = par.Results.numberOfMHCs;

% List of Antigens
switch(ID)
    case 0
        % Self Antigen
        kd(1).Value = 100e-9;
        kd(1).Units = 'molarity';
        kd(1).Notes = '(Stone 2015, PMID: 25618219)';
    otherwise
        % Default Antigen
        kd.Value = 10e-9;
        kd.Units = 'molarity';
        kd.Notes = '(Stone 2015, PMID: 25618219)';
end

% Antigen kd
if (nMHCs==1) % model has 1 MHC
    % Reduce kds to one value using harmonic mean
    sum.Value = 0; sum.Units = ['1/' kd(1).Units];
    for i = 1:length(kd)
        sum = sum + 1/kd(i);
    end
    antigen.kd = length(kd)/sum;
    antigen.kd.Notes = kd(1).Notes; % assumes all kds are from the same source
elseif (nMHCs==length(kd)) % model has same number of MHCs as number of kds
    antigen.kd = kd;
elseif (nMHCs>length(kd)) % model has more MHCs than number of kds
    antigen.kd(1:length(kd)) = kd;
    % Set remaining kds to inf
    for i = length(kd)+1:nMHCs
        antigen.kd(i).Value = inf;
        antigen.kd(i).Units = 'molarity';
        antigen.kd(i).Notes = '';
    end
else % model has less MHCs than number of kds
    % Sort kds in Ascending Order
    values = [kd.Value];
    [~,inds] = sort(values);
    % Ignore largest kds
    for i = 1:nMHCs
        j = inds(i);
        antigen.kd(i).Value = kd(j).Value;
        antigen.kd(i).Units = kd(j).Units;
        antigen.kd(i).Notes = kd(j).Notes;
    end
end

% Cancers Containing Antigen
antigen.cancers = cancers;

% Concentration per Cancer Cell
% estimated based on Milo 2013, PMID: 24114984
for i = 1:length(P_C)
    antigen.concentration(i).Value = P_C(i);
    antigen.concentration(i).Units = units;
    antigen.concentration(i).Notes = '';
end
