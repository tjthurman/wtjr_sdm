This README explains the files in the src folder.

* `compile_slim_*.R` - These scripts compile the results from the slim simulations. They compile the results early in the generation cycle, late in the generation cycle, and the random seeds used to initialize each simulation. These scripts are run as part of the `slim_simulations.smk` snakemake workflow. 

* `conservation_by_pheno.R` - Performs statistical analysis of how predicted color pattern and conservation status are co-distributed. Run as part of the main `Snakefile` workflow. 

* `data_curation.Rmd` - An R markdown document that describes the curation of the white-tailed jackrabbit occurrence records that went into the SDM modelling. Run as part of the main `Snakefile` workflow. 

* `figure_1_maps.R` - A script to create the map elements that are part of figure 1 (elements were further combined in Adobe Illustrator for the figure presented in the publication). Run as part of the main `Snakefile` workflow. 

* `figure_4_elements.R` - A script to create the elements that are part of figure 4 (elements were further combined in Adobe Illustrator for the figure presented in the publication). Also performs some calculations of the proportion of the range that is brown. Run as part of the main `Snakefile` workflow. 

* `other_RCP_scenario_data.R` - An ancillary script, not part of the main workflow, giving information on the provenance and data processing of the projected Bioclim variables under different RCP scenarios (not the RCP8.5 scenario used for most of the analyses). 

* `predict_current_pheno.R` - A script to predict current color phenotype from current environmental data. Run as part of the main `Snakefile` workflow. 

* `predict_future_pheno.R` - A script to predict/project future color phenotype from projected future environmental data. Run as part of the main `Snakefile` workflow. 

* `prep_bg_points.R` - A script to generate a set of background points used the the SDM. Run as part of the main `Snakefile` workflow.

* `prep_bioclim_data.R` - Downloads and processes current bioclim data for use in the SDM. Run as part of the main `Snakefile` workflow.

* `prep_conservation_rasters.R` - Prepares rasters of the extents of the different qualitative conservation categories. Run as part of the main `Snakefile` workflow.

* `prep_current_pheno_predictors.R` - Prepares the data used for predicting current color phenotypes. Run as part of the main `Snakefile` workflow. 

* `prep_future_bioclim_data.R` - Prepares the bioclim data used for predicting future color phenotypes. Run as part of the main `Snakefile` workflow. 

* `prep_unthinned_occur.R` - Prepares the unthinned occurrence data for use in SDM models. Run as part of the main `Snakefile` workflow. 

* `process_ENMeval_results.R` - Processes the results of the different SDM models and outputs some performance metrics. Run as part of the main `Snakefile` workflow. 

* `run_ENMeval.R` - A script to run a single ENM model. Run as part of the main `Snakefile` workflow. 

* `sdm_from_best_mod.R` - Generate range rasters from the best-fit SDM model. Run as part of the main `Snakefile` workflow. 

* `slim_calculate_fitness_function_sd.R` - An ancillary script, not part of the main workflow, used to figure out parameters in the Gaussian fitness function used in the SLiM simulations. 

* `slim_simulations/` - A folder containing SLiM code to perform each of the various SLiM simulations. These scripts are run as part of the `slim_simulations.smk` snakemake workflow. 

* `supp_figs_and_tables.R` - Generates supplementary figures, tables, and some statistical results. Run as part of the main `Snakefile` workflow. 

* `thin_records.R` - Thins occurrence records for input into the SDM. Run as part of the main `Snakefile` workflow. 