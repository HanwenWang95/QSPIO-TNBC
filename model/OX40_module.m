% OX40 Module
%
% Models OX40 Interactions
%
% Inputs: model        -- simbio model object with four compartments 
%         ICs          -- array of ICs
%                         - concentration of OX40-aOX40 complex
%                         - concentration of OX40-OX40l complex
%                         - concentration of Teff-OX40 complex
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - kd_aOX40--kd of PD1-nivolumab binding
%                         - OX40--concentration of OX40 expression on T cells
%                         - aOX40--concentration of aOX40 expression
%        
% Outputs: model -- simbio model object with new PD1 module
%
% Created: Feb 10, 2019 (Mohammad Jafarnejad & Pamela Chansky)
% Last Modified: Jun 14, 2019 (MJ)

function model = OX40_module(model,params,aOX40_params)

% LN Compartment
comp_LN = model.Compartment(4);

% Add kon Values
kon = addparameter(model,'kon_OX40_aOX40',params.kon_OX40_aOX40.Value,'ValueUnits',params.kon_OX40_aOX40.Units);
    set(kon,'Notes',['kon of OX40-aOX40 binding ' params.kon_OX40_aOX40.Notes]);
kon = addparameter(model,'kon_OX40_OX40l',params.kon_OX40_OX40l.Value,'ValueUnits',params.kon_OX40_OX40l.Units);
    set(kon,'Notes',['kon of OX40-OX40l binding ' params.kon_OX40_OX40l.Notes]);
% Add koff Values
koff = addparameter(model,'koff_OX40_aOX40',params.koff_OX40_aOX40.Value,'ValueUnits',params.koff_OX40_aOX40.Units);
    set(koff,'Notes',['koff of OX40-aOX40 binding ' params.koff_OX40_aOX40.Notes]);
koff = addparameter(model,'koff_OX40_OX40l',params.koff_OX40_OX40l.Value,'ValueUnits',params.koff_OX40_OX40l.Units);
    set(koff,'Notes',['koff of OX40-OX40l binding ' params.koff_OX40_OX40l.Notes]);
% Bivalent anibody parameters  
chi = addparameter(model,'Chi_OX40_aOX40' ,params.Chi_OX40_aOX40.Value ,'ValueUnits',params.Chi_OX40_aOX40.Units);
    set(chi,'Notes',['Antibody cross-arm binding efficiency ' params.Chi_OX40_aOX40.Notes]);

% Pharmacokinetics of a OX40 Ab
model = pk_module(model,'aOX40',aOX40_params);

try 
    p = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
    set(p,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
catch 
end
try 
    p = addparameter(model,'A_syn' ,params.A_syn.Value ,'ValueUnits',params.A_syn.Units);
    set(p,'Notes',['Surface area of the synapse ' params.A_syn.Notes]);
catch 
end

% Total amount of OX40 on naive T cells
p = addparameter(model,'nT_OX40_tot',params.nT_OX40_tot.Value,'ValueUnits',params.nT_OX40_tot.Units,'ConstantValue',false);
    set(p,'Notes',['Total number of OX40 on naive T cells ' params.nT_OX40_tot.Notes]); 
p = addparameter(model,'APC_OX40l_tot',params.APC_OX40l_tot.Value,'ValueUnits',params.APC_OX40l_tot.Units,'ConstantValue',false);
    set(p,'Notes',['Total number of OX40l on APCs ' params.APC_OX40l_tot.Notes]);   
% OX40-related Hill parameters
p = addparameter(model,'nT_OX40_50',params.nT_OX40_50.Value,'ValueUnits',params.nT_OX40_50.Units);
    set(p,'Notes',['OX40 occupancy for half-maximal naive T cell co-estimulation' params.nT_OX40_50.Notes]);
p = addparameter(model,'n_nT_OX40',params.n_nT_OX40.Value,'ValueUnits',params.n_nT_OX40.Units);
    set(p,'Notes',['OX40 occupancy Hill coefficient for naive T cell co-estimulation' params.n_nT_OX40.Notes]);   
% Define Hill function for OX40
addparameter(model,'nT_OX40_occupancy',0.0,'ValueUnits','molecule','ConstantValue',false);  
addparameter(model,'H_nT_OX40',0.0,'ValueUnits','dimensionless','ConstantValue',false);  

% Species for states of OX40 on naive T cells
x = addspecies(comp_LN,'nT_OX40_out',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of free OX40 molecules on naive T cells outside of the synapse');
x = addspecies(comp_LN,'nT_OX40_aOX40_out',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of OX40-aOX40 complex on naive T cells outside of the synapse');
x = addspecies(comp_LN,'nT_OX40_aOX40_OX40_out',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of OX40-aOX40-OX40 complex on naive T cells outside of the synapse');
x = addspecies(comp_LN,'nT_OX40_syn',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of free OX40 molecules on naive T cells inside the synapse');
x = addspecies(comp_LN,'APC_OX40l_syn',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of free OX40l molecules on APC inside the synapse');
x = addspecies(comp_LN,'nT_OX40_OX40l_syn',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of OX40-OX40l complex inside the synapse');
x = addspecies(comp_LN,'nT_OX40_aOX40_syn',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of OX40-aOX40 complex on naive T cells inside the synapse');
x = addspecies(comp_LN,'nT_OX40_aOX40_OX40_syn',0,'InitialAmountUnits','molecule');
    set(x,'Notes','Number of OX40-aOX40-OX40 complex on naive T cells inside the synapse');
  
% Initialize the total CTLA4 on the cells 
addrule(model,'V_LN.nT_OX40_out = nT_OX40_tot *(1-A_syn/A_Tcell)'   ,'initialAssignment');
addrule(model,'V_LN.nT_OX40_syn = nT_OX40_tot *(A_syn/A_Tcell)'     ,'initialAssignment');
addrule(model,'V_LN.APC_OX40l_syn = APC_OX40l_tot *(A_syn/A_Tcell)'   ,'initialAssignment');

% Binding and unbinding of aOX40 to OX40 outside of the synapse of naive T cell and APC in the LN
R = addreaction(model, 'V_LN.nT_OX40_out <-> V_LN.nT_OX40_aOX40_out');
    set (R, 'ReactionRate', 'kon_OX40_aOX40*(V_LN.nT_OX40_out * V_LN.aOX40/gamma_LN_aOX40) -  koff_OX40_aOX40*V_LN.nT_OX40_aOX40_out');
    set (R, 'Notes'       , 'binding and unbinding of OX40 to aOX40 on naive T cells outside of the synapse');
R = addreaction(model, 'V_LN.nT_OX40_aOX40_out + V_LN.nT_OX40_out <-> V_LN.nT_OX40_aOX40_OX40_out');
    set (R, 'ReactionRate', 'Chi_OX40_aOX40*kon_OX40_aOX40*(V_LN.nT_OX40_aOX40_out * V_LN.nT_OX40_out)/(A_Tcell-A_syn) -  koff_OX40_aOX40*V_LN.nT_OX40_aOX40_OX40_out');
    set (R, 'Notes'       , 'binding and unbinding of OX40 to OX40-aOX40 on naive T cells outside of the synapse');    
% Binding and unbinding of aOX40 to OX40 inside the synapse of naive T cell and APC in the LN
R = addreaction(model, 'V_LN.nT_OX40_syn <-> V_LN.nT_OX40_aOX40_syn');
    set (R, 'ReactionRate', 'kon_OX40_aOX40*(V_LN.nT_OX40_syn * V_LN.aOX40/gamma_LN_aOX40) -  koff_OX40_aOX40*V_LN.nT_OX40_aOX40_syn');
    set (R, 'Notes'       , 'binding and unbinding of OX40 to aOX40 on naive T cells outside of the synapse');
R = addreaction(model, 'V_LN.nT_OX40_aOX40_syn + V_LN.nT_OX40_syn <-> V_LN.nT_OX40_aOX40_OX40_syn');
    set (R, 'ReactionRate', 'Chi_OX40_aOX40*kon_OX40_aOX40*(V_LN.nT_OX40_aOX40_syn * V_LN.nT_OX40_syn)/(A_syn) -  koff_OX40_aOX40*V_LN.nT_OX40_aOX40_OX40_syn');
    set (R, 'Notes'       , 'binding and unbinding of OX40 to OX40-aOX40 on naive T cells outside of the synapse');   
% Binding and unbinding of OX40l to OX40 inisde the synapse of naive T cell and APC in the LN
R = addreaction(model, 'V_LN.nT_OX40_syn + V_LN.APC_OX40l_syn <-> V_LN.nT_OX40_OX40l_syn');
    set (R, 'ReactionRate', 'kon_OX40_OX40l*(V_LN.nT_OX40_syn * V_LN.APC_OX40l_syn)/(A_syn) -  koff_OX40_OX40l*V_LN.nT_OX40_OX40l_syn');
    set (R, 'Notes'       , 'binding and unbinding of OX40 to OX40l inside the synapse'); 
    
% Add the Hill function based on OX40 occupancy 
addrule(model,'nT_OX40_occupancy = V_LN.nT_OX40_aOX40_out + 2*V_LN.nT_OX40_aOX40_OX40_out + V_LN.nT_OX40_aOX40_syn + 2*V_LN.nT_OX40_aOX40_OX40_syn + V_LN.nT_OX40_OX40l_syn','repeatedAssignment');
addrule(model,'H_nT_OX40 = (nT_OX40_occupancy/nT_OX40_50)^n_nT_OX40/((nT_OX40_occupancy/nT_OX40_50)^n_nT_OX40 + 1)','repeatedAssignment');
  
