% Macrophage Module
%
% Models macrophage recruitment and polarization in tumor compartment
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing model parameter Values, Units, and Notes:
%
% Outputs: model -- SimBiology model object with new macrophage module


function model = macrophage_module(model,params,cancer_types,varargin)

% Optional Inputs
in = inputParser;
addParameter(in,'aCD47',1);
% Parse Inputs
parse(in,varargin{:});
antiCD47 = in.Results.aCD47;

% Add Species
Mac_M1  = addspecies(model.Compartment(3),'Mac_M1',0,'InitialAmountUnits','cell');
    set(Mac_M1,'Notes',['Number of M1 macrophage in the tumor compartment']);
Mac_M2  = addspecies(model.Compartment(3),'Mac_M2',0,'InitialAmountUnits','cell');
    set(Mac_M2,'Notes',['Number of M2 macrophage in the tumor compartment']);
IL12  = addspecies(model.Compartment(3),'IL12',0,'InitialAmountUnits','nanomolarity');
    set(IL12,'Notes',['Concentration of IL-12 in tumor']);
IL10  = addspecies(model.Compartment(3),'IL10',0,'InitialAmountUnits','nanomolarity');
    set(IL10,'Notes',['Concentration of IL-10 in tumor']);

% Add Parameters
k_Mac_rec = addparameter(model,'k_Mac_rec',params.k_Mac_rec.Value,'ValueUnits',params.k_Mac_rec.Units);
    set(k_Mac_rec,'Notes',['Recruitment rate of macrophage into tumor compartment ',params.k_Mac_rec.Notes]);
k_Mac_death = addparameter(model,'k_Mac_death',params.k_Mac_death.Value,'ValueUnits',params.k_Mac_death.Units);
    set(k_Mac_death,'Notes',['Death rate of macrophages into tumor compartment ',params.k_Mac_death.Notes]);
k_TGFb_Msec = addparameter(model,'k_TGFb_Msec',params.k_TGFb_Msec.Value,'ValueUnits',params.k_TGFb_Msec.Units);
    set(k_TGFb_Msec,'Notes',['Secretion rate of TGFb by macrophage in tumor compartment ',params.k_TGFb_Msec.Notes]);
k_vas_Msec = addparameter(model,'k_vas_Msec',params.k_vas_Msec.Value,'ValueUnits',params.k_vas_Msec.Units);
    set(k_vas_Msec,'Notes',['Secretion rate of angiogenic factor by macrophage in tumor compartment ',params.k_vas_Msec.Notes]);
k_IL12_sec = addparameter(model,'k_IL12_sec',params.k_IL12_sec.Value,'ValueUnits',params.k_IL12_sec.Units);
    set(k_IL12_sec,'Notes',['Secretion rate of IL-12 by mAPC in tumor compartment ' params.k_IL12_sec.Notes]);
k_IL12_Msec = addparameter(model,'k_IL12_Msec',params.k_IL12_Msec.Value,'ValueUnits',params.k_IL12_Msec.Units);
    set(k_IL12_Msec,'Notes',['Secretion rate of IL-12 by macrophage in tumor compartment ',params.k_IL12_Msec.Notes]);
k_IL12_deg = addparameter(model,'k_IL12_deg',params.k_IL12_deg.Value,'ValueUnits',params.k_IL12_deg.Units);
    set(k_IL12_deg,'Notes',['Degradation rate of IL-12 in tumor compartment ',params.k_IL12_deg.Notes]);
k_IL10_sec = addparameter(model,'k_IL10_sec',params.k_IL10_sec.Value,'ValueUnits',params.k_IL10_sec.Units);
    set(k_IL10_sec,'Notes',['Secretion rate of IL-10 by macrophage in tumor compartment ',params.k_IL10_sec.Notes]);
k_IL10_deg = addparameter(model,'k_IL10_deg',params.k_IL10_deg.Value,'ValueUnits',params.k_IL10_deg.Units);
    set(k_IL10_deg,'Notes',['Degradation rate of IL-10 in tumor compartment ',params.k_IL10_deg.Notes]);

k_M2_pol = addparameter(model,'k_M2_pol',params.k_M2_pol.Value,'ValueUnits',params.k_M2_pol.Units);
    set(k_M2_pol,'Notes',['Rate of M1 to M2 macrophage polarization ',params.k_M2_pol.Notes]);
k_M1_pol = addparameter(model,'k_M1_pol',params.k_M1_pol.Value,'ValueUnits',params.k_M1_pol.Units);
    set(k_M1_pol,'Notes',['Rate of M2 to M1 macrophage polarization ',params.k_M1_pol.Notes]);
IL10_50 = addparameter(model,'IL10_50',params.IL10_50.Value,'ValueUnits',params.IL10_50.Units);
    set(IL10_50,'Notes',['Half-maximal IL10 level for M1 to M2 polarization / maintaining Treg function / mAPC inhibition (*STAT6 related) ',params.IL10_50.Notes]);
IL12_50 = addparameter(model,'IL12_50',params.IL12_50.Value,'ValueUnits',params.IL12_50.Units);
    set(IL12_50,'Notes',['Half-maximal IL-12 level for M2 to M1 macrophage polarization ',params.IL12_50.Notes]);
IFNg_50 = addparameter(model,'IFNg_50',params.IFNg_50.Value,'ValueUnits',params.IFNg_50.Units);
    set(IFNg_50,'Notes',['Half-maximal IFNg level for M2 to M1 macrophage polarization ',params.IFNg_50.Notes]);

k_M1_phago = addparameter(model,'k_M1_phago',params.k_M1_phago.Value,'ValueUnits',params.k_M1_phago.Units);
    set(k_M1_phago,'Notes',['Rate of M1 macrophage-mediated killing of cancer cell ',params.k_M1_phago.Notes]);

vol_Mcell = addparameter(model,'vol_Mcell',params.vol_Mcell.Value,'ValueUnits',params.vol_Mcell.Units);
    set(vol_Mcell,'Notes',['Volume of a single macrophage cell ',params.vol_Mcell.Notes]);

first_call = true;
try
    % CCL2 in tumor
    CCL2 = addspecies(model.Compartment(3),'CCL2',0,'InitialAmountUnits','nanomolarity');
        set(CCL2,'Notes','Concentration of CCL2 in the tumor compartment');
    % CCL2 Parameters
    p = addparameter(model,'k_CCL2_sec',params.k_CCL2_sec.Value,'ValueUnits',params.k_CCL2_sec.Units);
        set(p,'Notes',['rate of CCL2 secretion ' params.k_CCL2_sec.Notes]);
    p = addparameter(model,'k_CCL2_deg',params.k_CCL2_deg.Value,'ValueUnits',params.k_CCL2_deg.Units);
        set(p,'Notes',['rate of CCL2 degradation ' params.k_CCL2_deg.Notes]);
    p = addparameter(model,'CCL2_50',params.CCL2_50.Value,'ValueUnits',params.CCL2_50.Units);
        set(p,'Notes',['Half-maximal CCL2 level of MDSC recruitment ' params.CCL2_50.Notes]);
    % CCL2 secretion by MDSCs and Cancer cells
    reaction = addreaction(model,'null -> V_T.CCL2');
        set(reaction,'ReactionRate','k_CCL2_sec*C_total');
        set(reaction,'Notes','CCL2 secretion by MDSCs and cancer cells');
    % CCL2 Degradation
    reaction = addreaction(model,'V_T.CCL2 -> null');
        set(reaction,'ReactionRate','k_CCL2_deg*V_T.CCL2');
        set(reaction,'Notes','CCL2 degradation');
catch
    first_call = false;
end

if ~isempty(sbioselect(model, 'Name', 'c_vas'))
    % Secretion of angiogenic factor
    reaction = addreaction(model,'null -> V_T.c_vas');
        set(reaction,'ReactionRate','k_vas_Msec*V_T.Mac_M2');
        set(reaction,'Notes','Secretion of angiogenic factor by macrophage in tumor');
end

% Add Reactions
% Recruitment of M1 macrophage
reaction = addreaction(model,'null -> V_T.Mac_M1');
    set(reaction,'ReactionRate','k_Mac_rec*V_T*(V_T.CCL2/(V_T.CCL2 + CCL2_50))');
    set(reaction,'Notes','Recruitment of pro-inflammatory M1 macrophage to tumor');
% Macrophage death upon tumor eradication
reaction = addreaction(model,'V_T.Mac_M1 -> null');
    set(reaction,'ReactionRate','k_cell_clear*V_T.Mac_M1*(Kc_rec/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','Macrophage death upon tumor eradication');
reaction = addreaction(model,'V_T.Mac_M2 -> null');
    set(reaction,'ReactionRate','k_cell_clear*V_T.Mac_M2*(Kc_rec/(C_total^2 + Kc_rec))');
    set(reaction,'Notes','Macrophage death upon tumor eradication');

% Death of macrophage
reaction = addreaction(model,'V_T.Mac_M1 -> null');
    set(reaction,'ReactionRate','k_Mac_death*V_T.Mac_M1');
    set(reaction,'Notes','Death of macrophage in tumor');
reaction = addreaction(model,'V_T.Mac_M2 -> null');
    set(reaction,'ReactionRate','k_Mac_death*V_T.Mac_M2');
    set(reaction,'Notes','Death of macrophage in tumor');

% Secretion of IL-12 by mAPC in the tumor
reaction = addreaction(model,'null -> V_T.IL12');
    set(reaction,'ReactionRate','k_IL12_sec*V_T.mAPC');
    set(reaction,'Notes','Secretion of IL-12 by APC in the tumor');
% Secretion of IL-12 by M1 macrophage
reaction = addreaction(model,'null -> V_T.IL12');
    set(reaction,'ReactionRate','k_IL12_Msec*V_T.Mac_M1');
    set(reaction,'Notes','Secretion of IL-12 by macrophage in tumor');
% Degradation of IL-12
reaction = addreaction(model,'V_T.IL12 -> null');
    set(reaction,'ReactionRate','k_IL12_deg*V_T.IL12');
    set(reaction,'Notes','Degradation of IL-12 in tumor');
% Secretion of TGFb
reaction = addreaction(model,'null -> V_T.TGFb');
    set(reaction,'ReactionRate','k_TGFb_Msec*V_T.Mac_M2');
    set(reaction,'Notes','Secretion of TGFb by macrophage in tumor');

% Secretion of IL-10
reaction = addreaction(model,'null -> V_T.IL10');
    set(reaction,'ReactionRate','k_IL10_sec*V_T.Mac_M2');
    set(reaction,'Notes','Secretion of IL-10 by macrophage in tumor');
% Degradation of IL-10
reaction = addreaction(model,'V_T.IL10 -> null');
    set(reaction,'ReactionRate','k_IL10_deg*V_T.IL10');
    set(reaction,'Notes','Degradation of IL-10 in tumor');

% M1 to M2 polarization
reaction = addreaction(model,'V_T.Mac_M1 -> V_T.Mac_M2');
    set(reaction,'ReactionRate','k_M2_pol*V_T.Mac_M1*(V_T.TGFb/(V_T.TGFb+TGFb_50) + V_T.IL10/(V_T.IL10+IL10_50))');
    set(reaction,'Notes','M1 to M2 macrophage polarization in tumor');
% M2 to M1 polarization
reaction = addreaction(model,'V_T.Mac_M2 -> V_T.Mac_M1');
    set(reaction,'ReactionRate','k_M1_pol*V_T.Mac_M2*(V_T.IL12/(V_T.IL12+IL12_50) + V_T.IFNg/(V_T.IFNg+IFNg_50))');
    set(reaction,'Notes','M2 to M1 macrophage polarization in tumor');

model = phagocytosis_module(model,params,'aCD47',antiCD47);

IL10_50_phago = addparameter(model,'IL10_50_phago',params.IL10_50_phago.Value,'ValueUnits',params.IL10_50_phago.Units);
    set(IL10_50_phago,'Notes',['Half-maximal IL-10 level for phagocytosis by macrophage ',params.IL10_50_phago.Notes]);
K_Mac_C = addparameter(model,'K_Mac_C',params.K_Mac_C.Value,'ValueUnits',params.K_Mac_C.Units);
    set(K_Mac_C,'Notes',['Association coefficient for macrophage and cancer cell ',params.K_Mac_C.Notes]);

M_total = addparameter(model,'M_total',0,'ValueUnits','cell','ConstantValue',false);
    set(M_total,'Notes','Total number of macrophages ');
H_IL10 = addparameter(model,'H_IL10',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(H_IL10,'Notes','Hill function for IL10-mediated M2-like macrophage polarization / Treg functionality / inhibition on Ag presentaion ');
H_IL10_phago = addparameter(model,'H_IL10_phago',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(H_IL10_phago,'Notes','Hill function for IL10-mediated inhibition on phagocytosis ');
H_IL12 = addparameter(model,'H_IL12',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(H_IL12,'Notes','Hill function for IL12 ');

addrule(model,'M_total = V_T.Mac_M1+V_T.Mac_M2','repeatedAssignment');
addrule(model,'H_IL10 = V_T.IL10/(IL10_50+V_T.IL10)','repeatedAssignment');
addrule(model,'H_IL10_phago = V_T.IL10/(IL10_50_phago+V_T.IL10)','repeatedAssignment');
addrule(model,'H_IL12 = V_T.IL12/(IL12_50+V_T.IL12)','repeatedAssignment');

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update tumor Volume (Rule 1)
volume_rule = model_rules(1);
rule = get(volume_rule,'Rule');
set(volume_rule,'Rule',[rule '+M_total*vol_Mcell/Ve_T']); % 4990 µm^3/cell

% Update Cancer Killing by Macrophage
if (exist('cancer_types','var'))
    for i = 1:length(cancer_types)
        reaction = addreaction(model, ['V_T.' cancer_types{i} ' -> V_T.C_x']);
            set(reaction,'ReactionRate',['k_M1_phago*V_T.' cancer_types{i} '*V_T.Mac_M1/(V_T.Mac_M1+K_Mac_C*V_T.' cancer_types{i} '+cell)*(1-H_Mac_C)*(1-H_IL10_phago)']);
            set(reaction,'Notes','M1 macrophage-mediated killing of cancer cell');
        % Update antigen release rate (Rule 5)
        rule = get(model_rules(5),'Rule');
            set(model_rules(5),'Rule',[rule,'+k_M1_phago*V_T.' cancer_types{i} '*V_T.Mac_M1/(V_T.Mac_M1+K_Mac_C*V_T.' cancer_types{i} '+cell)*(1-H_Mac_C)*(1-H_IL10_phago)']);
    end
end

% Update Treg-mediated Teff exhaustion and APC maturation
idx = [];
for i = 1:length(model.reaction)
    if strcmp(model.reaction(i).reaction, 'V_T.T1 -> V_T.T1_exh') && ~isempty(strfind(model.reaction(i).ReactionRate, 'k_Treg'))
        model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*H_IL10'];
    elseif strcmp(model.reaction(i).reaction, 'V_T.APC -> V_T.mAPC') && ~isempty(strfind(model.reaction(i).ReactionRate, 'c/(c+c50)'))
        model.reaction(i).ReactionRate = 'k_APC_mat*V_T.APC*H_IL12*(1-H_IL10)';
    elseif strcmp(model.reaction(i).reaction, 'null -> V_T.c')
        idx = [idx i];
    end
end
delete(model.reaction(idx))

idx = [];
for i = 1:length(model.parameters)
    if strcmp(model.parameters(i).Name, 'k_c')
        idx = [idx i];
    elseif strcmp(model.parameters(i).Name, 'DAMPs')
        idx = [idx i];
    elseif strcmp(model.parameters(i).Name, 'c50')
        idx = [idx i];
    elseif strcmp(model.parameters(i).Name, 'c0')
        idx = [idx i];
    end
end
delete(model.parameters(idx))

idx = [];
for i = 1:length(model.species)
    if strcmp(model.species(i).Name, 'c')
        idx = [idx i];
    end
end
delete(model.species(idx))
