% Phagocytosis Module
%
% Sub-module used by macrophage module
%
% Inputs: model        -- simbio model object with four compartments
%         params       -- object containing model parameter Values, Units, and Notes
%
% Outputs: model -- simbio model object with new phagocytosis module


function model = phagocytosis_module(model,params,varargin)

% Optional Inputs
in = inputParser;
addParameter(in,'aCD47',1);
% Parse Inputs
parse(in,varargin{:});
antiCD47 = in.Results.aCD47;

compDrug = model.Compartment(3);
gamma = 'gamma_T';

% Add the synapse compartment
comp = addcompartment(model,'syn_M_C',params.A_syn.Value,'CapacityUnits',params.A_syn.Units);
    set(comp,'Notes',['synapse comparment between macrophage and cancer cell ', params.A_syn.Notes]);

kdeg_CD47   = addparameter(model,'kdeg_CD47',params.kdeg_CD47.Value,'ValueUnits',params.kdeg_CD47.Units);
    set(kdeg_CD47,'Notes',['Degradation rate of free CD47 on RBCs ' params.kdeg_CD47.Notes]);

if (antiCD47)
% Add Species
x = addspecies(comp,'CD47_aCD47',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD47-aCD47 complex in synapse');
x = addspecies(comp,'CD47_aCD47_CD47',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD47-aCD47-CD47 complex in synapse');

% Pharmacokinetics of a CD47 Ab
params_aCD47 = pk_parameters('aCD47');
model = pk_module(model,'aCD47',params_aCD47);
end

% Add kon Values
kon = addparameter(model,'kon_CD47_SIRPa',params.kon_CD47_SIRPa.Value,'ValueUnits',params.kon_CD47_SIRPa.Units);
    set(kon,'Notes',['kon of CD47-SIRPa binding ' params.kon_CD47_SIRPa.Notes]);
% Add koff Values
koff = addparameter(model,'koff_CD47_SIRPa',params.koff_CD47_SIRPa.Value,'ValueUnits',params.koff_CD47_SIRPa.Units);
    set(koff,'Notes',['koff of CD47-SIRPa binding ' params.koff_CD47_SIRPa.Notes]);
% SIRPa-related Hill parameters
p = addparameter(model,'SIRPa_50',params.SIRPa_50.Value,'ValueUnits',params.SIRPa_50.Units);
    set(p,'Notes',['SIRPa occupancy for half-maximal CD47/SIRPa inhibition on phagocytosis of cancer cells by macrophages ' params.SIRPa_50.Notes]);
p = addparameter(model,'n_SIRPa',params.n_SIRPa.Value,'ValueUnits',params.n_SIRPa.Units);
    set(p,'Notes',['Hill coefficient for CD47/SIRPa inhibition on phagocytosis of cancer cells by macrophages ' params.n_SIRPa.Notes]);
% Define Hill function for phagocytosis
p = addparameter(model,'H_Mac_C',0,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function for overall checkpoint inhibition on phagocytosis of cancer cell by macrophage ');
p = addparameter(model,'H_PD1_M',0,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function for PD1-mediated inhibition on phagocytosis of cancer cell by macrophage ');
p = addparameter(model,'H_SIRPa',0,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function for SIRPa-mediated inhibition on phagocytosis of cancer cell by macrophage ');

p = addparameter(model,'C_CD47',params.C_CD47.Value,'ValueUnits',params.C_CD47.Units);
    set(p,'Notes',['CD47 expression on cancer cells ' params.C_CD47.Notes]);
p = addparameter(model,'M_PD1_total',params.M_PD1.Value,'ValueUnits',params.M_PD1.Units);
    set(p,'Notes',['PD-1 expression on macrophages ' params.M_PD1.Notes]);
p = addparameter(model,'M_SIRPa',params.M_SIRPa.Value,'ValueUnits',params.M_SIRPa.Units);
    set(p,'Notes',['SIRPa expression on macrophages ' params.M_SIRPa.Notes]);

% Add Species
x = addspecies(comp,'CD47',params.C_CD47.Value,'InitialAmountUnits',params.C_CD47.Units); % https://doi.org/10.1101/752311
    set(x,'Notes','concentration of CD47 on cancer cell in the synapse');
x = addspecies(comp,'SIRPa',params.M_SIRPa.Value,'InitialAmountUnits',params.M_SIRPa.Units); % 16291597
    set(x,'Notes','concentration of SIRPa on macrophage in the synapse');
x = addspecies(comp,'CD47_SIRPa',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD47-SIRPa complex in synapse');

x = addspecies(comp,'PDL1_total',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of bound and unbound PDL1 molecules ');
x = addspecies(comp,'PDL2_total',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of bound and unbound PDL2 molecules ');
x = addspecies(comp,'PD1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-PDL1 complex');
x = addspecies(comp,'PD1_PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-PDL2 complex');
x = addspecies(comp,'PD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1 in synapse');
x = addspecies(comp,'PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1 in synapse');
x = addspecies(comp,'PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL2 in synapse');
x = addspecies(comp,'PD1_aPD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-aPD1 complex');
x = addspecies(comp,'PD1_aPD1_PD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-aPD1-PD1 complex');
x = addspecies(comp,'PDL1_aPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1-aPDL1 complex');
x = addspecies(comp,'PDL1_aPDL1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1-aPDL1-PDL1 complex');
x = addspecies(comp,'PDL1_CD80',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1-CD80 complex');
x = addspecies(comp,'CD80',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD80 in synapse');
x = addspecies(comp,'CD80m',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD80 monomer in synapse');

p = addparameter(model,'A_Mcell' ,params.A_Mcell.Value ,'ValueUnits',params.A_Mcell.Units);
    set(p,'Notes',['Surface area of the macrophage ' params.A_Mcell.Notes]);
Cname = 'C1';
% Initial Conditions
addrule(model,[comp.Name,'.PD1 = M_PD1_total /A_Mcell' ] ,'initialAssignment');
addrule(model,[comp.Name,'.CD80 = ',Cname,'_CD80_total/A_cell'] ,'initialAssignment');
addrule(model,[comp.Name,'.PDL1 = ',Cname,'_PDL1_base /A_cell'] ,'initialAssignment');
addrule(model,[comp.Name,'.PDL2 = ',Cname,'_PDL1_base*r_PDL2',Cname,' /A_cell'] ,'initialAssignment');
addrule(model,[comp.Name,'.PDL1_total = ',comp.Name,'.PDL1+',comp.Name,'.PD1_PDL1+',comp.Name,'.PDL1_aPDL1+2*',...
                                          comp.Name,'.PDL1_aPDL1_PDL1+',comp.Name,'.PDL1_CD80'] ,'repeatedAssignment');
addrule(model,[comp.Name,'.PDL2_total = ',comp.Name,'.PD1_PDL2+',comp.Name,'.PDL2'] ,'repeatedAssignment');

addrule(model,[comp.Name,'.CD47 = C_CD47' ] ,'initialAssignment');
addrule(model,[comp.Name,'.SIRPa = M_SIRPa' ] ,'initialAssignment');

% PDL1 Secretion
R = addreaction(model,['null -> ',comp.Name,'.PDL1']);
    set (R, 'ReactionRate', ['k_out_PDL1*V_T.IFNg/(V_T.IFNg+IFNg_50_ind)*(1-',...
    comp.Name,'.PDL1_total/(',Cname,'_PDL1_base*r_PDL1_IFNg/A_cell))']);
    set (R, 'Notes'       , 'Translocation of PDL1 between cell surface and cytoplasm');
R = addreaction(model,['null -> ',comp.Name,'.PDL2']);
    set (R, 'ReactionRate', ['k_out_PDL1*r_PDL2',Cname,'*V_T.IFNg/(V_T.IFNg+IFNg_50_ind)*(1-',...
    comp.Name,'.PDL2_total/(',Cname,'_PDL1_base*r_PDL1_IFNg/A_cell*r_PDL2',Cname,'))']);
    set (R, 'Notes'       , 'Translocation of PDL2 between cell surface and cytoplasm');
R = addreaction(model,['null -> ' comp.Name,'.PDL1']);
    set (R, 'ReactionRate', ['k_in_PDL1*(',Cname,'_PDL1_base/A_cell-',comp.Name,'.PDL1_total)']);
    set (R, 'Notes'       , 'Translocation of PDL1 between cell surface and cytoplasm');
R = addreaction(model,['null -> ' comp.Name,'.PDL2']);
    set (R, 'ReactionRate', ['k_in_PDL1*(',Cname,'_PDL1_base/A_cell*r_PDL2',Cname,'-',comp.Name,'.PDL2_total)']);
    set (R, 'Notes'       , 'Translocation of PDL2 between cell surface and cytoplasm');

% Binding b/t CD47 and SIRPa in the synapse of cancell cell and macrophage in trans
R = addreaction(model,['null -> ',comp.Name,'.CD47']);
    set (R, 'ReactionRate', 'kdeg_CD47*A_syn*C_CD47');
    set (R, 'Notes'       , 'Translocation of CD47 between cell surface and cytoplasm');
R = addreaction(model,[comp.Name,'.CD47 -> null']);
    set (R, 'ReactionRate', ['kdeg_CD47*',comp.Name,'.CD47']);
    set (R, 'Notes'       , 'Internalization of CD47 on cancer cell');
R = addreaction(model,[comp.Name,'.CD47 + ',comp.Name,'.SIRPa <-> ',comp.Name,'.CD47_SIRPa']);
    set (R, 'ReactionRate', ['kon_CD47_SIRPa*(',comp.Name,'.CD47)*(',comp.Name,'.SIRPa)  -  koff_CD47_SIRPa*',comp.Name,'.CD47_SIRPa']);
    set (R, 'Notes'       , 'binding and unbinding of CD47-SIRPa in synapse');

% Dynamics of PD1/PDL1/PDL2/aPD1/aPDL1
R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PD1_PDL1']);
    set (R, 'ReactionRate', ['kon_PD1_PDL1*(',comp.Name,'.PD1)*(',comp.Name,'.PDL1)  -  koff_PD1_PDL1*',comp.Name,'.PD1_PDL1']);
    set (R, 'Notes'       , 'binding and unbinding of PD1 PDL1 in synapse');
R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL2 <-> ',comp.Name,'.PD1_PDL2']);
    set (R, 'ReactionRate', ['kon_PD1_PDL2*(',comp.Name,'.PD1)*(',comp.Name,'.PDL2)  -  koff_PD1_PDL2*',comp.Name,'.PD1_PDL2']);
    set (R, 'Notes'       , 'binding and unbinding of PD1 PDL2 in synapse');
R = addreaction(model,[comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1']);
    set (R, 'ReactionRate', ['2*kon_PD1_aPD1*(',comp.Name,'.PD1 * ',compDrug.Name,'.aPD1/',gamma,'_aPD1) -  koff_PD1_aPD1*',comp.Name,'.PD1_aPD1']);
    set (R, 'Notes'       , ['binding and unbinding of PD1 to aPD1 on macrophage surface in synapse']);
R = addreaction(model,[comp.Name,'.PD1_aPD1 + ',comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1_PD1']);
    set (R, 'ReactionRate', ['Chi_PD1_aPD1*kon_PD1_aPD1*(',comp.Name,'.PD1 * ',comp.Name,'.PD1_aPD1) -  2*koff_PD1_aPD1*',comp.Name,'.PD1_aPD1_PD1']);
    set (R, 'Notes'       , ['binding and unbinding of PD1 to PD1-aPD1 on macrophage surface in synapse']);
R = addreaction(model,[comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1']);
    set (R, 'ReactionRate', ['2*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 * ',compDrug.Name,'.aPDL1/',gamma,'_aPDL1) -  koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1']);
    set (R, 'Notes'       , ['binding and unbinding of PDL1 to aPDL1 on ',Cname,' surface in synapse']);
R = addreaction(model,[comp.Name,'.PDL1_aPDL1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1_PDL1']);
    set (R, 'ReactionRate', ['Chi_PDL1_aPDL1*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 * ',comp.Name,'.PDL1_aPDL1) -  2*koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1_PDL1']);
    set (R, 'Notes'       , ['binding and unbinding of PDL1 to PDL1-aPDL1 on cancer cell surface in synapse']);

% CD80 dimer dissociation
R = addreaction(model,[comp.Name,'.CD80m + ',comp.Name,'.CD80m <-> ',comp.Name,'.CD80']);
    set (R, 'ReactionRate', ['kon_CD80_CD80*(',comp.Name,'.CD80m)*(',comp.Name,'.CD80m) - koff_CD80_CD80*',comp.Name,'.CD80']);
    set (R, 'Notes'       , 'self-association and dissociation of CD80 monomers in synapse');
% cis PDL1-CD80-CD28
R = addreaction(model,[comp.Name,'.CD80m + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_CD80']);
    set (R, 'ReactionRate', ['kon_CD80_PDL1*(',comp.Name,'.CD80m)*(',comp.Name,'.PDL1)  -  koff_CD80_PDL1*',comp.Name,'.PDL1_CD80']);
    set (R, 'Notes'       , 'binding and unbinding of CD80 and PDL1 in synapse');

if (antiCD47)
% Binding b/t CD47 and aCD47 in the synapse of cancell cell and macrophage in trans
R = addreaction(model,[comp.Name,'.CD47 <-> ',comp.Name,'.CD47_aCD47']);
    set (R, 'ReactionRate', ['2*kon_CD47_aCD47*(',comp.Name,'.CD47 * ',compDrug.Name,'.aCD47/',gamma,'_aCD47) -  koff_CD47_aCD47*',comp.Name,'.CD47_aCD47']);
    set (R, 'Notes'       , 'binding and unbinding of CD47 to aCD47 on cancer cell surface in synapse');
R = addreaction(model,[comp.Name,'.CD47_aCD47 + ',comp.Name,'.CD47 <-> ',comp.Name,'.CD47_aCD47_CD47']);
    set (R, 'ReactionRate', ['Chi_CD47_aCD47_3D*kon_CD47_aCD47/d_syn*(',comp.Name,'.CD47 * ',comp.Name,'.CD47_aCD47) -  2*koff_CD47_aCD47*',comp.Name,'.CD47_aCD47_CD47']);
    set (R, 'Notes'       , 'binding and unbinding of CD47 to CD47-aCD47 on cancer cell surface in synapse');
R = addreaction(model,[comp.Name,'.CD47_aCD47 -> null']);
    set (R,'ReactionRate',['kint_CD47*',comp.Name,'.CD47_aCD47']);
    set (R,'Notes','Internalization of CD47-aCD47 on cancer cell ');
R = addreaction(model,[comp.Name,'.CD47_aCD47_CD47 -> null']);
    set (R,'ReactionRate',['kint_CD47*',comp.Name,'.CD47_aCD47_CD47']);
    set (R,'Notes','Internalization of CD47-aCD47-CD47 on cancer cell ');
end

% Add the Hill function based on SIRPa occupancy
addrule(model,'H_SIRPa = (syn_M_C.CD47_SIRPa/SIRPa_50)^n_SIRPa/((syn_M_C.CD47_SIRPa/SIRPa_50)^n_SIRPa + 1)','repeatedAssignment');
addrule(model,'H_PD1_M = ((syn_M_C.PD1_PDL1+syn_M_C.PD1_PDL2)/PD1_50)^n_PD1/(((syn_M_C.PD1_PDL1+syn_M_C.PD1_PDL2)/PD1_50)^n_PD1 + 1)','repeatedAssignment');
addrule(model,'H_Mac_C = 1 - (1-H_SIRPa)*(1-H_PD1_M)','repeatedAssignment');
