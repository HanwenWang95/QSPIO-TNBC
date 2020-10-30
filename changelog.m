Previous modifications were made based on QSP model by Jafarnejad et al. (AAPS 2019):

- Added MDSC parameter file and MDSC module.
- Added default parameter file with baseline parameter for TNBC
- Added hill functions for MDSC and ENT in Cancer killing (T cell module) and Cancer cell growth (cancer module); and added hill function parameters in simbio_init.m.
- Updated pk parameter file and pk module to add entinostat.
- Updated schedule_dosing.m to add entinostat dose regimen.

Modifications made based on the previous version published on Front. Bioeng. Biotechnol. 2020 (doi.org/10.3389/fbioe.2020.00141):

- Added CD4 T helper cell module, and modified nomenclature of Teff and Treg parameters
- Modified T cell module, including naive T cell dynamics, cytokine secretions, and T cell trafficking parameters
- Added nab-paclitaxel PK/PD module
- Modified H_APC calculation to avoid overestimation of Treg density at the beginning of the simulation
- Changed tumor growth dynamic to Gompertzian model for dynamic maximal tumor capacity
- Added dynamic PD-L1 expression that is upregulated by IFN-gamma secreted by T helper cells
- Modified Treg expansion mechanism that is induced by arginase I and TGF-beta
