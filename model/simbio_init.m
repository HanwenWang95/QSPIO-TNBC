% Function to initialize simbiology object
%
% Inputs: model_name -- string containing model name
%         time       -- array of time points (in days)
%         solver     -- string containing solver (e.g. ode15s)
%         tol_abs    -- absolute tolerance of solver
%         tol_rel    -- relative tolerance of solver
%         params     -- object containing compartment parameters
%                       - V_C--central compartment volume
%                       - V_P--peripheral compartment volume
%                       - V_Tmin--cancer-free tumour compartment volume
%                       - V_LN--lymph node compartment volume
%                       - vol_cell--cancer cell volume
%                       - vol_Tcell--T cell volume
%                       - k_cell_clear--dead cell clearance rate
%
% Outputs: model -- simbiology model object containing user-specified
%                   configuration and four compartments: C,P,T,LN


function model = simbio_init(model_name,time,solver,tol_abs,tol_rel,params)

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

% Setup Compartments
comp_C = addcompartment(model,'V_C',params.V_C.Value,'CapacityUnits',params.V_C.Units);
    set(comp_C,'Notes',['Central compartment (C) ' params.V_C.Notes]);
comp_P = addcompartment(model,'V_P',params.V_P.Value,'CapacityUnits',params.V_P.Units);
    set (comp_P,'Notes',['Peripheral compartment (P) ' params.V_P.Notes]);
comp_T = addcompartment(model,'V_T',params.V_Tmin.Value,'CapacityUnits',params.V_Tmin.Units,'ConstantCapacity',false);
    set (comp_T,'Notes',['Tumor compartment (T) ' params.V_Tmin.Notes]);
comp_LN = addcompartment(model,'V_LN',params.V_LN.Value,'CapacityUnits',params.V_LN.Units);
    set(comp_LN,'Notes',['Lymph node (LN) compartment volume ' params.V_LN.Notes]);

% Dead Cells
P = addparameter(model,'k_cell_clear',params.k_cell_clear.Value,'ValueUnits',params.k_cell_clear.Units);
    set(P,'Notes',['Rate of dead cell clearance from tumour compartment ' params.k_cell_clear.Notes]);
S = addspecies(comp_T,'C_x',0,'InitialAmountUnits','cell');
    set(S,'Notes','Dead cancer cells in the tumour compartment');
S = addspecies(comp_T,'T1_exh',0,'InitialAmountUnits','cell');
    set(S,'Notes','Exhausted CD8 effector T cells');
S = addspecies(comp_T,'Th_exh',0,'InitialAmountUnits','cell');
    set(S,'Notes','Exhausted CD4 effector T cells');
R = addreaction(model,'V_T.C_x -> null');
    set(R,'ReactionRate','k_cell_clear*V_T.C_x');
R = addreaction(model,'V_T.T1_exh -> null');
    set(R,'ReactionRate','k_cell_clear*V_T.T1_exh');
R = addreaction(model,'V_T.Th_exh -> null');
    set(R,'ReactionRate','k_cell_clear*V_T.Th_exh');

% Define Cell and Time
p = addparameter(model,'cell',1.0,'ValueUnits','cell');
    set(p,'ReactionRate','unit parameter for calculation');
p = addparameter(model,'day',1.0,'ValueUnits','day');
    set(p,'ReactionRate','unit parameter for calculation');

% Set Tumour Volume (Rule 1)
vol_cell = addparameter(model,'vol_cell',params.vol_cell.Value,'ValueUnits',params.vol_cell.Units);
    set(vol_cell,'Notes',['Average volume of cancer cell ' params.vol_cell.Notes]);
vol_Tcell = addparameter(model,'vol_Tcell',params.vol_Tcell.Value,'ValueUnits',params.vol_Tcell.Units);
    set(vol_Tcell,'Notes',['Average volume of T cells ' params.vol_Tcell.Notes]);
% V_Tmin = addparameter(model,'V_Tmin',params.V_Tmin.Value,'ValueUnits',params.V_Tmin.Units);
%     set(V_Tmin,'Notes',['Cancer-Free Tumour compartment volume' params.V_Tmin.Notes]);

% rho_cell = addparameter(model,'rho_cell',params.rho_cell.Value,'ValueUnits',params.rho_cell.Units,'ConstantValue',false);
%     set(rho_cell,'Notes',['Tumor density']);
Ve_T = addparameter(model,'Ve_T',params.Ve_T.Value,'ValueUnits',params.Ve_T.Units,'ConstantValue',false);
    set(Ve_T,'Notes',['Void fraction of the tumor']);

% addrule(model,'V_T = V_Tmin+vol_cell*C_x+vol_Tcell*T_exh','repeatedAssignment');
addrule(model,'V_T = ((C_x+C_total)*vol_cell+(T1_exh+Th_exh+T_total+V_T.Th)*vol_Tcell)/Ve_T','repeatedAssignment');
%addrule(model,'V_T = (C_x+C_total+T1_exh+Th_exh+T_total+V_T.Th)/rho_cell','repeatedAssignment');
%addrule(model,'V_T = (C_x+C_total)*vol_cell/Ve_T','repeatedAssignment');

% Set Total Number of Cancer Cells (Rule 2)
p = addparameter(model,'C_total',0,'ValueUnits','cell','ConstantValue',false);
    set(p,'Notes','Total number of cancer cells');
addrule(model,'C_total = 0*cell','repeatedAssignment');

% Set Total Number of T Cells in Tumour (Rule 3)
p = addparameter(model,'T_total',0,'ValueUnits','cell','ConstantValue',false);
    set(p,'Notes','Total number of activated T cells in tumor');
addrule(model,'T_total = 0*cell','repeatedAssignment');

% Set Total Number of T Cells in LN (Rule 4)
p = addparameter(model,'T_total_LN',0,'ValueUnits','cell','ConstantValue',false);
    set(p,'Notes','Total number of activated T cells in TDLNs');
addrule(model,'T_total_LN = 0*cell','repeatedAssignment');

% Set Total Rate of Cancer Death by T Cells (Rule 5)
p = addparameter(model,'R_Tcell','ValueUnits','cell/day','ConstantValue',false);
    set(p,'Notes','Rate of cancer cell death');
addrule(model,'R_Tcell = 0*cell/day','repeatedAssignment');

% Set Default Number of Tregs
p = addparameter(model,'Tregs_',0,'ValueUnits','cell','ConstantValue',false);
    set(p,'Notes','Total number of activated T cells in TDLNs');

% Set Default Hill Function for APCs
p = addparameter(model,'H_APC',0.5,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of mAPC for self-antigen-specific Treg activation');

% Set Default Hill Function for mAPCs
p = addparameter(model,'H_mAPC',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of mAPC for CD8 T cell activation');

% Set Default Hill Function for mAPCs
p = addparameter(model,'H_APCh',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of mAPC for neoantigen-specific helper T cell activation');

% Set Default Hill Function for PD1 Checkpoint
p = addparameter(model,'H_PD1_C1',0.90,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of PDL1 on tumor cells for T cell exhaustion');

%addparameter(model,'H_PD1_C2',0.90,'ValueUnits','dimensionless','ConstantValue',false);
p = addparameter(model,'H_PD1_APC',0.90,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of PDL1 on APC for T cell exhaustion');

% Set Default Hill Function for CTLA4 Checkpoint
p = addparameter(model,'H_CD28_C1',0.1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of CD28 on tumor cells');

%addparameter(model,'H_CD28_C2',0.1,'ValueUnits','dimensionless','ConstantValue',false);
p = addparameter(model,'H_CD28_APC',0.1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of CD28 on APCs for T cell activation');

% Set Default Hill Function for MDSC and ENT
p = addparameter(model,'H_MDSC_C1',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of MDSC on T cell inhibition');
p = addparameter(model,'H_ENT_C1',1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of entinostat on MDSC inhibition');
