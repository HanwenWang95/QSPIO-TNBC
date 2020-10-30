% MDSC Mini Model
%
%         species      - MDSC--Myeloid-derived suppressor cell
%                      - ArgI--Arginase I
%                      - MCP1--Monocyte Chemoattractant Protein-1 (CCL2)
%                      - NO--Nitrite Oxide
%         params       - kg_C--Growth rate of cancer cell
%                      - kd_C--Death rate of cancer cell
%                      - C_max--Maximal cancer cell number
%                      - k_rec_MDSC--MDSC recruitment rate by MCP-1
%                      - k_brec_MDSC--Baseline MDSC recruitment rate
%                      - kd_MDSC--MDSC death rate
%                      - k_deg_MCP1--Degradation rate of MCP-1
%                      - k_deg_NO--Degradation rate of NO
%                      - k_deg_ArgI--Degradation rate of Arg I
%                      - k_sec_MCP1--Secretion rate of MCP-1
%                      - k_sec_NO--Secretion rate of NO
%                      - k_sec_ArgI--Secretion rate of Arg I
%                      - kd_C_ENT--Cancer cell apoptosis induced by entinostat (negligible due to high EC50)
%                      - IC50_ENT_C--Effective concentration of entinostat on inhibition of cancer proliferation
%                      - IC50_ENT_NO--Effective concentration of entinostat on inhibition of NO production
%                      - IC50_ENT_ArgI--Effective concentration of entinostat on inhibition of Arg I production
%                      - IC50_ENT_MCP1--Effective concentration of entinostat on inhibition of MCP-1 production
%                      - k_Treg--Treg expansion rate
%                      - kd_T--Death rate of T cells
%                      - IC50_ArgI_CTL--Effective concentration of Arg I on inhibition of CTL activity
%                      - IC50_NO_CTL--Effective concentration of NO on inhibition of CTL activity
%                      - EC50_MCP1_rec--Effective concentration of MCP-1 on recruitment of MDSC
%                      - EC50_ArgI_Treg--Effective concentration of Arg I on Treg expansion
%                      - MDSC0--Maximal MDSC available to recruit
%                      - Cap--Total cell capacity
%
%
% Created: Jul 8, 2019 (Hanwen Wang)
% Last Modified: Aug 1, 2019 (Hanwen Wang)

if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

if (isempty(sbioshowunits('mU')))
    mU_unit = sbiounit('mU','molarity');
    sbioaddtolibrary(mU_unit);
end
% Symbolic Unit
u = symunit;
try u.mU;
catch
    newUnit('mU',u.molarity);
end

% Growth Rate
% params.kg_C.Value = 1.5; % MCF-7
params.kg_C.Value = 0.465;
params.kg_C.Units = '1/day';
params.kg_C.Notes = '(varied depends on exp setting)';
% Death Rate
params.kd_C.Value = 0.001;
params.kd_C.Units = '1/day';
params.kd_C.Notes = '(DB - Vicini 2013)';
% Maximal cancer cells
% params.C_max.Value = 1.8e6; % MCF-7
params.C_max.Value = 3.66e6;
params.C_max.Units = 'cell';
params.C_max.Notes = '(varied depends on exp setting)';

% MDSC recruitment rate
params.k_rec_MDSC.Value = 1.2;
params.k_rec_MDSC.Units = '1/day';
params.k_rec_MDSC.Notes = '(Lai 2018 PMID: 29735668; Huang 2006, PMID: 17257744)';
% Baseline MDSC migration rate
params.k_brec_MDSC.Value = 0.0021;
params.k_brec_MDSC.Units = '1/day';
params.k_brec_MDSC.Notes = '(Huang 2006, PMID: 17257744)';
% rate of MDSC death
params.kd_MDSC.Value = 0.015;
params.kd_MDSC.Units = '1/day';
params.kd_MDSC.Notes = '(Lai 2018)';
% cancer death rate by entinostat
% params.kd_C_ENT.Value = 2.2;
% params.kd_C_ENT.Units = '1/day';
% params.kd_C_ENT.Notes = '(Lee 2001, PMID: 11221885)';
% half-maximal ENT concentration for anti-proliferative effect on tumor cells
params.IC50_ENT_C.Value = 3.74e-7;
params.IC50_ENT_C.Units = 'molarity';
params.IC50_ENT_C.Notes = '(Bouchain 2002, PMID: 12593661; Lee 2001)';
% half-life of IL-10: 1.1 - 2.6 day Mueller 2009, PMID: 18818669

% rate of MCP-1 degradtion
params.k_deg_MCP1.Value = 0.06;
params.k_deg_MCP1.Units = '1/hour';
params.k_deg_MCP1.Notes = '(Tanimoto 2007, PMID: 18089573)';
% rate of NO degradtion
params.k_deg_NO.Value = 135;
params.k_deg_NO.Units = '1/day';
params.k_deg_NO.Notes = '(Hakim 1996, PMID: 8953625)';
% rate of ArgI degradtion
params.k_deg_ArgI.Value = 0.173;
params.k_deg_ArgI.Units = '1/day';
params.k_deg_ArgI.Notes = '(Schimke 1964, PMID: 14257612)';
% rate of MCP-1 secretion
params.k_sec_MCP1.Value = 14.2e-11; % normal: 1.065e-11; TNBC: (14.2 +/- 6)*10^-11
params.k_sec_MCP1.Units = 'nanomole/cell/day';
params.k_sec_MCP1.Notes = '(Huang 2007, PMID: 17257744; Dutta, PMID: 29594759)';
% rate of NO secretion
params.k_sec_NO.Value = 4.8e-7;
params.k_sec_NO.Units = 'nanomole/cell/day';
params.k_sec_NO.Notes = '(Serafini 2008, PMID: 18593947)';
% rate of ArgI secretion
params.k_sec_ArgI.Value = 1.4e-2;
params.k_sec_ArgI.Units = '(mU*microliter)/cell/day';
params.k_sec_ArgI.Notes = '(Serafini 2008, PMID: 18593947)';
% half-maximal ENT concentration for NO inhibition
params.IC50_ENT_NO.Value = 1e-9; % <1 nM
params.IC50_ENT_NO.Units = 'molarity';
params.IC50_ENT_NO.Notes = '(Choo 2010, doi:10.1093/rheumatology/keq108)';
% half-maximal ENT concentration for MCP-1 inhibition
params.IC50_ENT_MCP1.Value = 2e-9; % <2 nM
params.IC50_ENT_MCP1.Units = 'molarity';
params.IC50_ENT_MCP1.Notes = '(Choo 2013, PMID: 24241152)';
% half-maximal ENT concentration for Arg I inhibition
params.IC50_ENT_ArgI.Value = 5e-7;
params.IC50_ENT_ArgI.Units = 'molarity';
params.IC50_ENT_ArgI.Notes = '(Orillian 2017, PMID: 28698201)';

% rate of ArgI-induced Treg expension
params.k_Treg.Value = 2.7;
params.k_Treg.Units = '1/day';
params.k_Treg.Notes = '(Serafini 2008, PMID: 18593947)';
% rate of ArgI/NO-induced T cell death Park 2018, PMID: 29491381
% rate of T cell death
params.kd_T.Value = 0.01;
params.kd_T.Units = '1/day';
params.kd_T.Notes = '(de Boer 1995)';
% half-maximal ArgI concentration for CTL inhibition
params.IC50_ArgI_CTL.Value = 10.14;
params.IC50_ArgI_CTL.Units = 'mU';
params.IC50_ArgI_CTL.Notes = '(Serafini 2008, PMID: 18593947)';
% half-maximal NO concentration for CTL inhibition
params.IC50_NO_CTL.Value = 1.87e-9;
params.IC50_NO_CTL.Units = 'molarity';
params.IC50_NO_CTL.Notes = '(Serafini 2008, PMID: 18593947)';
% half-maximal MCP1 concentration for MDSC recruitment
params.EC50_MCP1_rec.Value = 5e-10;
params.EC50_MCP1_rec.Units = 'molarity';
params.EC50_MCP1_rec.Notes = '(Ernst 1994, PMID: 8144933)';
% half-maximal ArgI concentration for Treg expension
params.EC50_ArgI_Treg.Value = 22.1;
params.EC50_ArgI_Treg.Units = 'mU';
params.EC50_ArgI_Treg.Notes = '(estimated; Serafini 2008, PMID: 18593947)';
% Maximal MDSC level
params.MDSC0.Value = 4e5;
params.MDSC0.Units = 'cell/milliliter';
params.MDSC0.Notes = '(estimated)';
% Total cell capacity
params.Cap.Value = 1.2e6; % maximal cell number in 6-well cell culture plate
params.Cap.Units = 'cell';
params.Cap.Notes = '(varied depends on exp setting)';

% Define Cell Dimension
params.cell.Value = 1;
params.cell.Units = 'cell';
% Cancer-Free Tumour Volume ***verify assumption***
params.V_T.Value = 210;
params.V_T.Units = 'microliter';
params.V_T.Notes = '(constant to reproduce in vitro setting)';
% Cell Volumes
one_cell.Value = 1;
one_cell.Units = 'cell';

%% Create the SimBiology Model
% Model Settings
model_name = 'MDSC Mini Model';
start_time = 0.0; % [days]
time_step = 1; % [days] 0.01 days ~ 15 mins
end_time = 5; % [days]
tol_abs = 1e-12;
tol_rel = 1e-9;
solver = 'ode15s';
% Model Object
time = start_time:time_step:end_time;

% Model Object
model = sbiomodel(model_name);
config = getconfigset(model);
options = get(config,'CompileOptions');
set(options,'UnitConversion',true);
set(config,'TimeUnits','day');
set(config,'SolverType',solver);
set(config.SolverOptions,'OutputTimes',time);
% set(config.SolverOptions,'MaxStep',0.01);
set(config,'StopTime',time(end));
set(config.SolverOptions,'AbsoluteTolerance',tol_abs);
set(config.SolverOptions,'RelativeTolerance',tol_rel);
comp_T = addcompartment(model,'V_T',params.V_T.Value,'CapacityUnits',params.V_T.Units,'ConstantCapacity',false);
set (comp_T,'Notes',['Tumor compartment (T) ' params.V_T.Notes]);

% Define Cell and Time
addparameter(model,'cell',1.0,'ValueUnits','cell');
addparameter(model,'day',1.0,'ValueUnits','day');

V_T = addparameter(model,'V_T',params.V_T.Value,'ValueUnits',params.V_T.Units);
    set(V_T,'Notes',['Cancer-Free Tumour compartment volume' params.V_T.Notes]);

% Initial Conditions
MCP10 = 0; NO0 = 0; Arg0 = 0;

% Add Species
% MDSC in tumor
MDSC = addspecies(comp_T,'MDSC',0,'InitialAmountUnits','cell');
    set(MDSC,'Notes',['Number of MDSCs in the tumour compartment']);

ENT = addspecies(comp_T,'ENT',1e-9,'InitialAmountUnits','molarity');
    set(ENT,'Notes',['Concentration of entinostat in the tumor compartment']);

C = addspecies(comp_T,'C',5e5,'InitialAmountUnits','cell');
    set(C,'Notes',['Number of cancer cells']);

Teff = addspecies(comp_T,'Teff',0,'InitialAmountUnits','cell');
    set(Teff,'Notes',['Number of effector T cells']);

Treg = addspecies(comp_T,'Treg',0,'InitialAmountUnits','cell');
    set(Treg,'Notes',['Number of Treg cells']);

% MCP1 in tumor
MCP1 = addspecies(comp_T,'MCP1',MCP10,'InitialAmountUnits','molarity');
    set(MCP1,'Notes','Concentration of MCP1 in the tumor compartment');
% add NO
NO = addspecies(comp_T,'NO',NO0,'InitialAmountUnits','molarity');
    set(NO,'Notes','Concentration of NO in the tumor compartment');
% add ArgI
ArgI = addspecies(comp_T,'ArgI',Arg0,'InitialAmountUnits','mU');
    set(ArgI,'Notes','Concentration of Arg I in the tumor compartment');

% Add Parameters
k_rec_MDSC = addparameter(model,'k_rec_MDSC',params.k_rec_MDSC.Value,'ValueUnits',params.k_rec_MDSC.Units);
    set(k_rec_MDSC,'Notes',['Rate of MDSC recruitment into the tumor' params.k_rec_MDSC.Notes]);
kd_MDSC = addparameter(model,'kd_MDSC',params.kd_MDSC.Value,'ValueUnits',params.kd_MDSC.Units);
    set(kd_MDSC,'Notes',['Rate of MDSC death ' params.kd_MDSC.Notes]);
% kd_C_ENT = addparameter(model,'kd_C_ENT',params.kd_C_ENT.Value,'ValueUnits',params.kd_C_ENT.Units);
%     set(kd_C_ENT,'Notes',['Cancer cell death rate ' params.kd_C_ENT.Notes]);

% ENT Hill Function
IC50_ENT_C = addparameter(model,'IC50_ENT_C',params.IC50_ENT_C.Value,'ValueUnits',params.IC50_ENT_C.Units);
    set(IC50_ENT_C,'Notes',['ENT concentration for half-maximal tumor cell death' params.IC50_ENT_C.Notes]);
% MCP1 Parameters
p = addparameter(model,'k_deg_MCP1',params.k_deg_MCP1.Value,'ValueUnits',params.k_deg_MCP1.Units);
    set(p,'Notes',['rate of MCP1 degradation ' params.k_deg_MCP1.Notes]);
p = addparameter(model,'k_deg_NO',params.k_deg_NO.Value,'ValueUnits',params.k_deg_NO.Units);
    set(p,'Notes',['rate of NO degradation ' params.k_deg_NO.Notes]);
p = addparameter(model,'k_deg_ArgI',params.k_deg_ArgI.Value,'ValueUnits',params.k_deg_ArgI.Units);
    set(p,'Notes',['rate of ArgI degradation ' params.k_deg_ArgI.Notes]);

p = addparameter(model,'k_sec_MCP1',params.k_sec_MCP1.Value,'ValueUnits',params.k_sec_MCP1.Units);
    set(p,'Notes',['rate of MCP1 secretion ' params.k_sec_MCP1.Notes]);
p = addparameter(model,'k_sec_NO',params.k_sec_NO.Value,'ValueUnits',params.k_sec_NO.Units);
    set(p,'Notes',['rate of NO secretion from MDSCs ' params.k_sec_NO.Notes]);
p = addparameter(model,'k_sec_ArgI',params.k_sec_ArgI.Value,'ValueUnits',params.k_sec_ArgI.Units);
    set(p,'Notes',['rate of ArgI secretion ' params.k_sec_ArgI.Notes]);

p = addparameter(model,'IC50_ENT_NO',params.IC50_ENT_NO.Value,'ValueUnits',params.IC50_ENT_NO.Units);
    set(p,'Notes',['half-maximal ENT concentration for NO inhibition ' params.IC50_ENT_NO.Notes]);
p = addparameter(model,'k_Treg',params.k_Treg.Value,'ValueUnits',params.k_Treg.Units);
    set(p,'Notes',['rate of ArgI-induced Treg expension ' params.k_Treg.Notes]);
p = addparameter(model,'kd_T',params.kd_T.Value,'ValueUnits',params.kd_T.Units);
    set(p,'Notes',['rate of T cell death ' params.kd_T.Notes]);

p = addparameter(model,'IC50_ArgI_CTL',params.IC50_ArgI_CTL.Value,'ValueUnits',params.IC50_ArgI_CTL.Units);
    set(p,'Notes',['rate of ArgI-induced T cell death ' params.IC50_ArgI_CTL.Notes]);
p = addparameter(model,'IC50_NO_CTL',params.IC50_NO_CTL.Value,'ValueUnits',params.IC50_NO_CTL.Units);
    set(p,'Notes',['rate of NO-induced T cell death ' params.IC50_NO_CTL.Notes]);
p = addparameter(model,'EC50_MCP1_rec',params.EC50_MCP1_rec.Value,'ValueUnits',params.EC50_MCP1_rec.Units);
    set(p,'Notes',['Half-maximal MCP1 level of MDSC recruitment ' params.EC50_MCP1_rec.Notes]);
p = addparameter(model,'EC50_ArgI_Treg',params.EC50_ArgI_Treg.Value,'ValueUnits',params.EC50_ArgI_Treg.Units);
    set(p,'Notes',['Half-maximal ArgI level of Treg expansion ' params.EC50_ArgI_Treg.Notes]);

p = addparameter(model,'MDSC0',params.MDSC0.Value,'ValueUnits',params.MDSC0.Units);
    set(p,'Notes',['MDSC density in the tumour ' params.MDSC0.Notes]);

kg_C = addparameter(model,'kg_C',params.kg_C.Value,'ValueUnits',params.kg_C.Units);
    set(kg_C,'Notes',['Cancer cell growth rate ' params.kg_C.Notes]);
C_max = addparameter(model,'C_max',params.C_max.Value,'ValueUnits',params.C_max.Units);
    set(C_max,'Notes',['Cancer cell capacity ' params.C_max.Notes]);
% Death
kd_C = addparameter(model,'kd_C',params.kd_C.Value,'ValueUnits',params.kd_C.Units);
    set(kd_C,'Notes',['Cancer cell death rate from innate immune cells ' params.kd_C.Notes]);

Cap = addparameter(model,'Cap',params.Cap.Value,'ValueUnits',params.Cap.Units);
    set(p,'Notes',['Maximal FoxP3+ T cells']);

p = addparameter(model,'IC50_ENT_MCP1',params.IC50_ENT_MCP1.Value,'ValueUnits',params.IC50_ENT_MCP1.Units);
    set(p,'Notes',['half-maximal ENT concentration for MCP1 inhibition ' params.IC50_ENT_MCP1.Notes]);
k_brec_MDSC = addparameter(model,'k_brec_MDSC',params.k_brec_MDSC.Value,'ValueUnits',params.k_brec_MDSC.Units);
    set(k_brec_MDSC,'Notes',['Baseline rate of MDSC migration into the tumor' params.k_brec_MDSC.Notes]);
p = addparameter(model,'IC50_ENT_ArgI',params.IC50_ENT_ArgI.Value,'ValueUnits',params.IC50_ENT_ArgI.Units);
    set(p,'Notes',['half-maximal ENT concentration for Arg I inhibition ' params.IC50_ENT_ArgI.Notes]);

% Add Reactions
% Cancer Proliferation
reaction = addreaction(model,'null -> V_T.C');
set(reaction,'ReactionRate','kg_C*V_T.C*(1-V_T.C/C_max)*IC50_ENT_C/(V_T.ENT+IC50_ENT_C)');
set(reaction,'Notes','Cancer cell proliferation in the tumor compartment');
% Cancer Death
reaction = addreaction(model,'V_T.C -> null');
set(reaction,'ReactionRate','kd_C*V_T.C');
set(reaction,'Notes','Cancer cell death in the tumor compartment');
% ENT-induced Cancer Cell Death
% reaction = addreaction(model,['V_T.C -> null']);
% set(reaction,'ReactionRate',['kd_C_ENT*V_T.C*V_T.ENT/(V_T.ENT+IC50_ENT_C)']);
% set(reaction,'Notes','Cancer cell death by entinostat');

% MDSC Recruitment
reaction = addreaction(model,'null -> V_T.MDSC');
set(reaction,'ReactionRate','k_rec_MDSC*(MDSC0*V_T-V_T.MDSC)*(MCP1/(EC50_MCP1_rec+MCP1))');
set(reaction,'Notes','MDSC recruitment into the tumor compartment');
% MDSC Baseline Migration
reaction = addreaction(model,'null -> V_T.MDSC');
set(reaction,'ReactionRate','k_brec_MDSC*(MDSC0*V_T-V_T.MDSC)');
set(reaction,'Notes','Baseline MDSC Migration into the tumor compartment');
% MDSC Death
reaction = addreaction(model,'V_T.MDSC -> null');
set(reaction,'ReactionRate','kd_MDSC*V_T.MDSC');
set(reaction,'Notes','MDSC death in the tumor compartment');

% MCP1 Degradation
reaction = addreaction(model,'V_T.MCP1 -> null');
set(reaction,'ReactionRate','k_deg_MCP1*V_T.MCP1');
set(reaction,'Notes','MCP1 degradation');
% NO Degradation
reaction = addreaction(model,'V_T.NO -> null');
set(reaction,'ReactionRate','k_deg_NO*V_T.NO');
set(reaction,'Notes','NO degradation');
% ArgI Degradation
reaction = addreaction(model,'V_T.ArgI -> null');
set(reaction,'ReactionRate','k_deg_ArgI*V_T.ArgI');
set(reaction,'Notes','ArgI degradation');

% MCP1 secretion by MDSCs and Cancer cells
reaction = addreaction(model,'null -> V_T.MCP1');
set(reaction,'ReactionRate',['k_sec_MCP1*V_T.C*(1-V_T.ENT/(V_T.ENT+IC50_ENT_MCP1))']);
set(reaction,'Notes','MCP1 secretion by MDSCs and cancer cells');
% NO Secretion by MDSC
reaction = addreaction(model,'null -> V_T.NO');
set(reaction,'ReactionRate','k_sec_NO*V_T.MDSC*(1-V_T.ENT/(V_T.ENT+IC50_ENT_NO))');
% set(reaction,'ReactionRate','k_sec_NO*V_T.MDSC*(MDSC_NO_50/(NO+MDSC_NO_50))*(1-V_T.ENT/(V_T.ENT+IC50_ENT_NO))');
set(reaction,'Notes','NO secretion from MDSC');
% ArgI Secretion by MDSC
reaction = addreaction(model,'null -> V_T.ArgI');
set(reaction,'ReactionRate','k_sec_ArgI*V_T.MDSC/V_T*(1-V_T.ENT/(V_T.ENT+IC50_ENT_ArgI))');
% set(reaction,'ReactionRate','k_sec_ArgI*V_T.MDSC*(MDSC_ArgI_50/(ArgI+MDSC_ArgI_50))');
set(reaction,'Notes','ArgI secretion from MDSC');

% ArgI-induced Treg Expension
reaction = addreaction(model,['null -> V_T.Treg']);
% set(reaction,'ReactionRate',['k_Treg*V_T.Treg*(1-V_T.Treg/(Cap-V_T.MDSC-V_T.Treg-V_T.Teff))']);
set(reaction,'ReactionRate',['k_Treg*V_T.Treg*(1-V_T.Treg/(Cap-V_T.MDSC-V_T.Treg-V_T.Teff))*ArgI/(EC50_ArgI_Treg+ArgI)']);
set(reaction,'Notes','ArgI-induced Treg Expension');

% Treg Death
reaction = addreaction(model,['V_T.Treg -> null']);
set(reaction,'ReactionRate',['kd_T*V_T.Treg']);
set(reaction,'Notes','Treg Death');

% T cell exhaustion
reaction = addreaction(model,['V_T.Teff -> null']);
set(reaction,'ReactionRate',['kd_T*V_T.Teff']);
set(reaction,'Notes','T cell Exhaustion');

% Cancer killing by T cell
% reaction = addreaction(model,['V_T.C -> V_T.C_x']);
% set(reaction,'ReactionRate',['k_C_Tcell*V_T.T*V_T.C/(C_total+T_total+cell)*(1-H_PD1_C1)*(1-NO/(IC50_NO_CTL+NO))*(1-ArgI/(IC50_ArgI_CTL+ArgI))']);
% set(reaction,'Notes','Updated cancer killing by effector T cells');
