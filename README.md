# SCC
Get concentration ranges for metabolites with structurally constrained concentrations (SCC).

folders:
* F2C2 
  - F2C2 toolbox from Larhlimi et al. 2002 used to calculate flux coupling matrix
* kinetic_model_predictions
  - code and data used for analysis based on a large-scale kinetic model of E. coli
  - includes the kinetic model used from Khodayari et al. 2014
  - compare predicted SCC concentration and concentration simulated from kinetic model (Correlation_prediction_simulation.m)
  - distribution of Euclidean distances between predicted and simulated SCC concentration over all samples (Distribution_euclidean_prediction_simulation.m)
  - get FBA based flux distribution for the kinetic model, which can be then used to predict SCC concentrations (FBA_flux_range.m)
  - SCC concentration of knock-out mutants predicted using MOMA (Moma_ko_analysis.m)
  - get samples of simulated steady-state data by simulations from different initial conditions (perturbationAnalysis.m)
  - effect of unknown kcat values on SCC concentration prediction (removal_kcat_knowledge.m)
  - compare prediction of concentration ranges obtained from SCC prediction to range obtained from shadow price (shadow_price_comparison.m)
  - compare simulated and predicted fold changes in SCC concentration upon reaction knock-ou (sim_pred_FoldChanges.m)
  - simulate knock-out mutants (Simulation_ko.m)
* Models
  - genome-scale metabolic models used in this study (format: .xml or .mat)
* Results
  - SCC metabolites for genome-scale models in Models/
  - SCC metabolites for genome-scale models in Models/ upon split into elementary reactions -Michaelis-Menten like format- (Results_MM_split/)

files:
* SCC concentration prediction for Gerosa data set (Gerosa_data_E_coli_metabolites.csv and Gerosa_data_E_coli_specific_rates.csv)
  and E. coli model iJO1366 (concentration_prediction_iJO1366_Gerosa.m)
* SCC concentration prediction for  Ishii data set (Ishii_Quantitative_data_E_coli_metabolites.csv and Ishii_Quantitative_data_E_coli_specific_rates.csv)
  and E. coli model iJO1366 (concentration_prediction_iJO1366_Ishii.m)
* get cytosolic SCC metabolites for E coli model iJO1366 (E_coli_iJO1366_get_cytosolic_SCC.m)
* fractional LP used in SCC constant prediction (fractional_LP.m)
* find SCC metabolites (get_SCC.m)
* find SCC metabolites and calculate constants based on provided kcat values (get_SCC_kcat.m)
* model should be without blocked reactions and reversible reactions should be split into irreversible once
  model preprocessing (model_preprocessing.m) and SCC metabolite calculation (main.m) 
* model preprocessing, MM-like split of reactions (model_preprocessing_MM.m) and SCC metabolite calculation (main_MM.m) 
* MM-like split of reactions into their elementary reactions (MM_split_rxns.m)
* split reversible reactions into two irreversible once (spli_rxns.m)

