% Nab-Paclitaxel Module
%
% Models nab-paclitaxel PK-PD
%
% Inputs: model        -- SimBiology model object
%         params       -- object containing the default parameters
%
% Outputs: model -- SimBiology model object with new nab-paclitaxel module



function model = nabpaclitaxel_module(model, params)

% Setup Compartments
comp_V1 = addcompartment(model,'V_1',15.8,'CapacityUnits','liter','ConstantCapacity',false);
    set(comp_V1,'Notes',['Central compartment (V1) ']);
comp_V2 = addcompartment(model,'V_2',1650,'CapacityUnits','liter','ConstantCapacity',false);
    set (comp_V2,'Notes',['Peripheral compartment (V2) ']);
comp_V3 = addcompartment(model,'V_3',75.4,'CapacityUnits','liter','ConstantCapacity',false);
    set (comp_V3,'Notes',['Peripheral compartment (V3) ']);

% comp_M0 = addcompartment(model,'M_0',1,'CapacityUnits','liter','ConstantCapacity',false);
%     set(comp_M0,'Notes',['Maturation compartment (M0) ']);
% comp_M1 = addcompartment(model,'M_1',1,'CapacityUnits','liter','ConstantCapacity',false);
%     set(comp_M1,'Notes',['Maturation compartment (M1) ']);
% comp_M2 = addcompartment(model,'M_2',1,'CapacityUnits','liter','ConstantCapacity',false);
%     set(comp_M2,'Notes',['Maturation compartment (M2) ']);
% comp_M3 = addcompartment(model,'M_3',1,'CapacityUnits','liter','ConstantCapacity',false);
%     set(comp_M3,'Notes',['Maturation compartment (M3) ']);
% comp_M4 = addcompartment(model,'M_4',1,'CapacityUnits','liter','ConstantCapacity',false);
%     set(comp_M4,'Notes',['Maturation compartment (M4) ']);

% Model Parameters
Vmcl = addparameter(model,'Vmcl',params.Vmcl.Value,'ValueUnits',params.Vmcl.Units);
    set(Vmcl,'Notes','Maximum elimination rate of nab-paclitaxel from the central compartment V1');
Kcl = addparameter(model,'Kcl',params.Kcl.Value,'ValueUnits',params.Kcl.Units);
    set(Kcl,'Notes','Nab-paclitaxel concentration in the central compartment V1 at 50% Vmcl');
Q2 = addparameter(model,'Q2',params.Q2.Value,'ValueUnits',params.Q2.Units);
    set(Q2,'Notes','Intercompartmental clearance between the central compartment V1 and the second peripheral compartment V3');
Vmt = addparameter(model,'Vmt',params.Vmt.Value,'ValueUnits',params.Vmt.Units);
    set(Vmt,'Notes','Maximum intercompartmental distribution rate between the central compartment V1 and the first peripheral compartment V2');
Kt = addparameter(model,'Kt',params.Kt.Value,'ValueUnits',params.Kt.Units);
    set(Kt,'Notes','Nab-paclitaxel concentration in the central compartment V1 at 50% Vmt');
BSA = addparameter(model,'BSA',params.BSA.Value,'ValueUnits',params.BSA.Units);
    set(BSA,'Notes','Body surface area');

% PD Parameters
% p = addparameter(model,'E',0,'ValueUnits','dimensionless','ConstantValue',false);
%     set(p,'Notes',['Nab-paclitaxel effect on neutrophils ']);
% p = addparameter(model,'gamma',0.187,'ValueUnits','dimensionless');
%     set(p,'Notes',['Feedback parameter for nab-paclitaxel pd ']);
% p = addparameter(model,'MTT',117,'ValueUnits','hour');
%     set(p,'Notes',['Mean transit time of neutrophils ']);
% p = addparameter(model,'k',0.0342,'ValueUnits','1/hour','ConstantValue',false);
%     set(p,'Notes',['Transit rate of neutrophils ']);
% p = addparameter(model,'slope',0.00253,'ValueUnits','1/(nanogram/milliliter)');
%     set(p,'Notes',['Nab-paclitaxel effect slope ']);
% p = addparameter(model,'Neu_base',4.28e9,'ValueUnits','cell/liter');
%     set(p,'Notes',['Baseline circulating neutrophil level ']);

% Model Species
s = addspecies(comp_V1,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the central compartment V1');
s = addspecies(comp_V2,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the peripheral compartment V2');
s = addspecies(comp_V3,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the peripheral compartment V3');

% s = addspecies(comp_M0,'Neu',4.28e9,'InitialAmountUnits','cell/liter');
%     set(s,'Notes','Neutrophil in bone marrow');
% s = addspecies(comp_M1,'Neu',4.28e9,'InitialAmountUnits','cell/liter');
%     set(s,'Notes','Neutrophil in the maturation compartment M1');
% s = addspecies(comp_M2,'Neu',4.28e9,'InitialAmountUnits','cell/liter');
%     set(s,'Notes','Neutrophil in the maturation compartment M2');
% s = addspecies(comp_M3,'Neu',4.28e9,'InitialAmountUnits','cell/liter');
%     set(s,'Notes','Neutrophil in the maturation compartment M3');
% s = addspecies(comp_M4,'Neu',4.28e9,'InitialAmountUnits','cell/liter');
%     set(s,'Notes','Circulating neutrophil ');

% Model Reaction
r = addreaction(model,'V_1.NabP <-> V_2.NabP');
    set(r,'ReactionRate','Vmt/(V_1.NabP+Kt)*V_1.NabP - Vmt/(V_2.NabP+Kt)*V_2.NabP');
    set(r,'Notes','Intercompartmental distribution of nab-paclitaxel between V1 and V2 compartment');
r = addreaction(model,'V_1.NabP <-> V_3.NabP');
    set(r,'ReactionRate','Q2*V_1.NabP - Q2*V_3.NabP');
    set(r,'Notes','Intercompartmental clearance of nab-paclitaxel between V1 and V3 compartment');
r = addreaction(model,'V_1.NabP -> null');
    set(r,'ReactionRate','Vmcl/(V_1.NabP+Kcl)*V_1.NabP');
    set(r,'Notes','Clearance of nab-paclitaxel from V1 compartment');

% r = addreaction(model,'null -> M_0.Neu');
%     set(r,'ReactionRate','k*M_0.Neu*(1-E)*(Neu_base/M_4.Neu)^gamma');
% r = addreaction(model,'M_4.Neu -> null');
%     set(r,'ReactionRate','k*M_4.Neu');
% r = addreaction(model,'M_0.Neu -> M_1.Neu');
%     set(r,'ReactionRate','k*M_0.Neu');
% r = addreaction(model,'M_1.Neu -> M_2.Neu');
%     set(r,'ReactionRate','k*M_1.Neu');
% r = addreaction(model,'M_2.Neu -> M_3.Neu');
%     set(r,'ReactionRate','k*M_2.Neu');
% r = addreaction(model,'M_3.Neu -> M_4.Neu');
%     set(r,'ReactionRate','k*M_3.Neu');


% Model Rule
% p = addparameter(model,'age',0,'ValueUnits','dimensionless','ConstantValue',false);
%     set(p,'Notes',['Age factor to estimate paclitaxel effect slope ']);

% addrule(model,'k = 4/MTT','repeatedAssignment');
% addrule(model,'E = V_1.NabP*slope*(1+age)','repeatedAssignment');


% PD Cytotoxicity
r_nabp = addparameter(model,'r_nabp',params.r_nabp.Value,'ValueUnits',params.r_nabp.Units,'ConstantValue',false);
    set(r_nabp,'Notes',['Tumour to plasma concentration ratio of nab-paclitaxel ' params.r_nabp.Notes]);
k_C_nabp = addparameter(model,'k_C_nabp',params.k_C_nabp.Value,'ValueUnits',params.k_C_nabp.Units,'ConstantValue',false);
    set(k_C_nabp,'Notes',['Cancer cell killing rate by nab-paclitaxel ' params.k_C_nabp.Notes]);
IC50_nabp = addparameter(model,'IC50_nabp',params.IC50_nabp.Value,'ValueUnits',params.IC50_nabp.Units,'ConstantValue',false);
    set(IC50_nabp,'Notes',['Half-maximal nab-paclitaxel concentration for cancer cell killing ' params.IC50_nabp.Notes]);
MW_nabp = addparameter(model,'MW_nabp',params.MW_nabp.Value,'ValueUnits',params.MW_nabp.Units);
    set(MW_nabp,'Notes',['Molecular weight of nab-paclitaxel ' params.MW_nabp.Notes]);
Kc_cyt = addparameter(model,'Kc_cyt',params.Kc_cyt.Value,'ValueUnits',params.Kc_cyt.Units);
    set(Kc_cyt,'Notes',['Half-Maximal cancer cell number for cytotoxic drug diffusion ' params.Kc_cyt.Notes]);

NabP = addspecies(model.Compartment(3),'NabP',0,'InitialAmountUnits','molarity');
    set(NabP,'Notes','Nab-paclitaxel concentration in tumour ');

addrule(model,'V_T.NabP = r_nabp*V_1.NabP/MW_nabp','repeatedAssignment');

reaction = addreaction(model,'V_T.C1 -> V_T.C_x');
    set(reaction,'ReactionRate','k_C_nabp*V_T.C1*(V_T.NabP/(V_T.NabP+IC50_nabp))*min(C_total,Kc_cyt)/C_total'); % *(1-C_total/(C_total+Kc_cyt))
    set(reaction,'Notes','Cancer cell death by nab-paclitaxel ');

addrule(model,'k_C1_therapy = k_C_nabp*(V_T.NabP/(V_T.NabP+IC50_nabp))*min(C_total,Kc_cyt)/C_total','repeatedAssignment');

%% Resistance
k_C_resist = addparameter(model,'k_C_resist',params.k_C_resist.Value,'ValueUnits',params.k_C_resist.Units,'ConstantValue',false);
    set(k_C_resist,'Notes',['Cancer resistance to nab-paclitaxel ' params.k_C_resist.Notes]);
r_resist = addparameter(model,'r_resist',params.r_resist.Value,'ValueUnits',params.r_resist.Units);
    set(r_resist,'Notes',['Number of folds increase of nab-paclitaxel EC50 in resistant cancer clones ' params.r_resist.Notes]);

reaction = addreaction(model,'V_T.C1 -> V_T.C2');
    set(reaction,'ReactionRate','k_C_resist*V_T.C1*H_TGF_CTL'); % *H_TGF_CTL
    set(reaction,'Notes','Cancer cell resistance to nab-paclitaxel ');
reaction = addreaction(model,'V_T.C2 -> V_T.C_x');
    set(reaction,'ReactionRate','k_C_nabp*V_T.C2*(V_T.NabP/(V_T.NabP+IC50_nabp*r_resist))*min(C_total,Kc_cyt)/C_total');
    set(reaction,'Notes','Resistant cancer cell death by nab-paclitaxel ');

addrule(model,'k_C2_therapy = k_C_nabp*(V_T.NabP/(V_T.NabP+IC50_nabp*r_resist))*min(C_total,Kc_cyt)/C_total','repeatedAssignment');

% set tumour growth rate of the resistant the same as the sensitive clone
addrule(model,'k_C2_growth = k_C1_growth','repeatedAssignment');
addrule(model,'k_C2_death = k_C1_death','repeatedAssignment');

%% Vascularization
s = addspecies(model.Compartment(3),'K',params.K0.Value,'InitialAmountUnits',params.K0.Units); % 2.4e8
    set(s,'Notes',['Maximal tumor capacity ' params.K0.Notes]);
s = addspecies(model.Compartment(3),'c_vas',0,'InitialAmountUnits','picogram/milliliter');
    set(s,'Notes','Angiogenic factors ');

p = addparameter(model,'k_K_g',params.k_K_g.Value,'ValueUnits',params.k_K_g.Units,'ConstantValue',false); % 5.33
    set(p,'Notes',['Tumour vasculature growth rate ' params.k_K_g.Notes]);
p = addparameter(model,'k_K_d',params.k_K_d.Value,'ValueUnits',params.k_K_d.Units,'ConstantValue',false); % 7.9e-3
    set(p,'Notes',['Tumour vasculature inhibition rate ' params.k_K_d.Notes]);

p = addparameter(model,'k_c_nabp',params.k_c_nabp.Value,'ValueUnits',params.k_c_nabp.Units,'ConstantValue',false);
    set(p,'Notes',['Secretion rate of angiogenic factors induced by nab-paclitaxel ' params.k_c_nabp.Notes]);
p = addparameter(model,'IC50_nabp_vas',params.IC50_nabp_vas.Value,'ValueUnits',params.IC50_nabp_vas.Units,'ConstantValue',false);
    set(p,'Notes',['Half-maximal conc. of nab-paclitaxel on angiogenic factor induction ' params.IC50_nabp_vas.Notes]);
p = addparameter(model,'k_c_vas',params.k_c_vas.Value,'ValueUnits',params.k_c_vas.Units,'ConstantValue',false);
    set(p,'Notes',['Secretion rate of angiogenic factors by cancer cells ' params.k_c_vas.Notes]);
p = addparameter(model,'k_c_deg',params.k_c_deg.Value,'ValueUnits',params.k_c_deg.Units,'ConstantValue',false);
    set(p,'Notes',['Degradation rate of angiogenic factors ' params.k_c_deg.Notes]);
p = addparameter(model,'IC50_vas',params.IC50_vas.Value,'ValueUnits',params.IC50_vas.Units,'ConstantValue',false);
    set(p,'Notes',['Half-maximal conc. of angiogenic factor on tumor capacity growth ' params.IC50_vas.Notes]);
p = addparameter(model,'k_endo_nabp',params.k_endo_nabp.Value,'ValueUnits',params.k_endo_nabp.Units,'ConstantValue',false);
    set(p,'Notes',['Inhibition rate of maximal tumor capacity by nab-paclitaxel ' params.k_endo_nabp.Notes]);

reaction = addreaction(model,'null -> V_T.c_vas');
    set(reaction,'ReactionRate','k_c_vas*C_total');
    set(reaction,'Notes','Baseline secretion/degradation of angiogenic factors ');
reaction = addreaction(model,'null -> V_T.c_vas');
    set(reaction,'ReactionRate','k_c_nabp*C_total*V_T.NabP/(V_T.NabP+IC50_nabp_vas)');
    set(reaction,'Notes','Angiogenic factor release in response to nab-paclitaxel ');
reaction = addreaction(model,'V_T.c_vas -> null');
    set(reaction,'ReactionRate','k_c_deg*V_T.c_vas');
    set(reaction,'Notes','Baseline secretion/degradation of angiogenic factors ');

reaction = addreaction(model,'null -> V_T.K');
    set(reaction,'ReactionRate','k_K_g*C_total*V_T.c_vas/(V_T.c_vas+IC50_vas)');
    set(reaction,'Notes','Resistant cancer cell death by nab-paclitaxel ');
reaction = addreaction(model,'V_T.K -> null');
    set(reaction,'ReactionRate','k_K_d*V_T.K*(nthroot(C_total/cell*2.57e-6,3))^2');
    set(reaction,'Notes','Resistant cancer cell death by nab-paclitaxel ');
reaction = addreaction(model,'V_T.K -> null');
    set(reaction,'ReactionRate','k_endo_nabp*V_T.K*V_T.NabP');
    set(reaction,'Notes','Endothelial cell death by nab-paclitaxel ');

addrule(model,'C_max = V_T.K','repeatedAssignment');
addevent(model,'V_T.K < 0.5*cell','V_T.K = 0.01*cell');
