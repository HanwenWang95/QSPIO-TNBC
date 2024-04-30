## Previous modifications were made based on QSP model by Jafarnejad et al. (AAPS 2019):

- Added MDSC parameter file and MDSC module.
- Added default parameter file with baseline parameter for TNBC
- Added hill functions for MDSC and ENT in Cancer killing (T cell module) and Cancer cell growth (cancer module); and added hill function parameters in simbio_init.m.
- Updated pk parameter file and pk module to add entinostat.
- Updated schedule_dosing.m to add entinostat dose regimen.

## Modifications made based on the previous version published on Front. Bioeng. Biotechnol. 2020 (doi.org/10.3389/fbioe.2020.00141):

- Added CD4 T helper cell module, and modified nomenclature of Teff and Treg parameters
- Modified T cell module, including naive T cell dynamics, cytokine secretions, cancer cell killing, and T cell trafficking parameters
- Added nab-paclitaxel PK/PD module
- Modified H_APC calculation to avoid overestimation of Treg density at the beginning of the simulation
- Changed tumor growth dynamic to Gompertzian model for dynamic maximal tumor capacity
- Added dynamic PD-L1 expression that is upregulated by IFN-gamma secreted by T helper cells
- Modified Treg expansion mechanism that is induced by arginase I and TGF-beta

## Modifications made since the published version in J Immunother Cancer. 2021 (DOI: 10.1136/jitc-2020-002100):

- Added a new macrophage module and a submodule for checkpoint interactions between cancer cell and macrophage that inhibit phagocytosis (macrophage_module.m, phagocytosis_module.m)
- Modified MDSC module to recalibrated cellular recruitment and cytokine secretion rates to match clinical measurements from TNBC tumors (MDSC_module.m, parameters_TNBC.m)
- Modified parameter nomenclature to be more consistent among the modules (parameters_TNBC.m, PSA_param_in_TNBC.m)
- Added estimated distribution of cytokine secretion rates and macrophage/MDSC recruitment rates for virtual patient generation (PSA_param_in_TNBC.m)
- Added postprocessing steps to calculate macrophage, MDSC, M1/M2 macrophage ratio, etc. (calculate_ratio.m, sprint_data.m)
- Added T cell death upon antigen clearance (Teff_module.m, Treg_module.m, Th_module.m)

## Modifications made since the published version in iScience 2022 (DOI: 10.1016/j.isci.2022.104702):

- Added parameter files for NSCLC and scripts for virtual patient generation (new methodology) & clinical trial simulation of durvalumab in NSCLC.
- Added a script for PK parameter fitting using compressed latent parameterization & updated durvalumab PK parameters & dosing schedule.
- Updated a model assumption so that IFN-gamma is now secreted by both CD8 and CD4 T cells when activated. 
- Modified nomenclature for some parameters & reorganized structure for some modules (to maintain relative independence & customizability of the modules, delete redundant model elements, and allow selection of Gompertzian/logistic growth dynamic in cancer_module.m).
- Added a script for sensitivity analysis using Morris screening method to account for non-monotonicity of the model.

## Modifications made since the published version in NPJ Precision Oncology (DOI: 10.1038/s41698-023-00405-9):

- Added PK parameter for anti-CD47 antibodies.
- Incoporated dynamics of anti-CD47 antibodies to PK module and phagocytosis module.
- Modified simbio_PSA.m so that the parameter values can be changed at the beginning of in silico clinical to simulate treatment effect.