## Snakemake workflow for WTJR SDM and phenotypic prediction analysis
## Pipeline order goes from top to bottom.

# Set up some wildcards for various aspects of the analysis
# Spatial thinning
thin_dists=[1,5,10,50]
dataset_dists=[0,1,5,10,50]
thin_dists_test=[50]
thin_res_pattern= "processed_data/thin/{thin_dist}km/"
current_bioclim="processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
# ENMeval
feature_class=["L", "LQ", "H", "LQH", "LQHP", "LQHPT"]
# Conservation
CON_STATUSES=["extirpated","broad_extirp","local_extirp","poss_decline","pres_stable"]


localrules: all, rename_ENMeval_res, dag, filegraph, help


rule all:
    input:
        "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd",
        expand("results/sdm/sdm_rangemap_best_{thresh}.{ext}", 
               thresh = ["specsens", "sens95", "sens99"],
               ext=["grd", "gri"]),
        "results/pheno/current_pheno_glm.RData",
        "results/pheno/current_predicted_probWhite.tif",
        "results/pheno/current_predicted_probWhite_SDMrange.tif",
        expand("processed_data/conservation_rasters/{cons_status}.tif", cons_status=CON_STATUSES)
        
## curate_occur_data   : process WTJR occurrence data
rule curate_occur_data:
    input:
        arctos="raw_data/ArctosData_43FA2173A0.csv",
        gbif= "raw_data/GBIF/verbatim.txt",
        vertnet= "raw_data/vertnet_leptownsendii_allrecords_apr2_2020.txt"
    output: 
        "processed_data/leptown_db_occurrences.csv",
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    resources:
        cpus=1
    shell:
        "Rscript -e \"rmarkdown::render('src/data_curation.Rmd', knit_root_dir = '../', output_dir = 'docs/')\""
        
## thin_occur_data  : spatial thinning of WTJR occurrence data
# 5 arguments
# input file, thinning distance, number of reps, output directory, and random seed
# for number of reps, did 3*number of unique gps records
rule thin_occur_data:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    output:
        directory(thin_res_pattern)
    resources:
        cpus=1
    shell:
      """
      Rscript src/thin_records.R {input} {wildcards.thin_dist} 3000 {output} 462
      """

## check_thin_data  : Double-check that thinning worked properly
rule check_thin_data:
    input:
        "processed_data/thin/"
    resources:
        cpus=1
    shell:
        """
        Rscript src/check_thinned_records.R {input}
        """

## prep_unthinned_data:
rule prep_unthinned_data:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    output:
        "processed_data/thin/0km/wtjr_occ_0km_thin1.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_unthinned_occur.R {input} {output}
        """
# prep_conservation_rasters
rule prep_conservation_rasters:
    input:
        "processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
    output:
        expand("processed_data/conservation_rasters/{cons_status}.tif", cons_status=CON_STATUSES)
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_conservation_rasters.R  {input}
        """
    
## process_bioclim_data: process current bioclim data
rule process_bioclim_data:
    output:
        current_bioclim
    resources:
        cpus=1
    shell:
        "Rscript src/prep_bioclim_data.R"

## process_current_pheno_predictors
rule process_current_pheno_predictors:
    input:
        snow="raw_data/rasterstack4cov.tif",
        bc=current_bioclim
    output:
        "processed_data/pheno_predictors_millsetal2018.tif"
    resources:
        cpus=1
    shell:
        """
        Rscript src/process_current_pheno_predictors.R {input.snow} {input.bc} {output}
        """

## process_future_bioclim: process future bioclim data
# inputs: dirctory containing the raw data
# current bioclim data to use as a template
# output file
rule process_future_bioclim:
    input:
        cmip5_dir="raw_data/worldclim_projections",
        current=current_bioclim
    output:
        "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_future_bioclim_data.R {input.cmip5_dir} {input.current} {output}
        """

## make_bg_points
rule make_bg_points:
    input:
        envir=current_bioclim,
    output:
        "processed_data/bg_points_for_sdm.RData"
    resources:
        cpus=1
    shell:
        """
        Rscript src/bg_points.R {input.envir} 782
        """

## run_ENMeval: create SDM with ENMeval
rule run_ENMeval:
    input:
       occur="processed_data/thin/{dataset}/wtjr_occ_{dataset}_thin1.csv",
       envir=current_bioclim,
       bg= "processed_data/bg_points_for_sdm.RData"
    output:
        "results/enmeval/enmeval_res_{dataset}_thin1_{feature_class}.RData"
    resources:
        cpus=1
    shell:
        """
        Rscript src/run_ENMeval.R {input.occur} {input.envir} {input.bg} {wildcards.feature_class} 782 {resources.cpus}

        """
## rename_ENMeval_res: rename files for 1m and 5km results
## Didn't use thin1 for the, used thin3 and thin2
rule rename_ENMeval_res:
    input: 
        expand("results/enmeval/enmeval_res_{km}km_thin1_{fc}.RData", km=["1", "5"], fc=feature_class)
    output:
        expand("results/enmeval/enmeval_res_1km_thin3_{fc}.Rdata", fc=feature_class),
        expand("results/enmeval/enmeval_res_5km_thin2_{fc}.Rdata", fc=feature_class)
    shell: 
        """
        rename thin1 thin3 results/enmeval/enmeval_res_1km_*
        rename thin1 thin2 results/enmeval/enmeval_res_5km_*
        """

## process_ENMeval_res: process results from ENMeval, use to decide best model
rule process_ENMeval_res:
    input:
        expand("results/enmeval/enmeval_res_{thin_dist}km_thin1_{fc}.RData", fc=feature_class, thin_dist=["0", "10", "50"]),
        expand("results/enmeval/enmeval_res_1km_thin3_{fc}.RData", fc=feature_class),
        expand("results/enmeval/enmeval_res_5km_thin2_{fc}.RData", fc=feature_class),
        res_dir="results/enmeval/"
    output:
        "results/enmeval/enmeval_metrics.csv",
        "results/enmeval/enmeval_best_model_per_thin_AIC.csv",
        expand("results/enmeval/performance_plot_{dist}km.pdf", dist = dataset_dists)
    resources:
        cpus=1
    shell:
        """
        Rscript src/process_ENMeval_results.R {input.res_dir}
        """

## sdm_range_raster: make rasters of the SDM rangemaps from the best model
rule sdm_range_raster: 
    input:
        "results/enmeval/enmeval_best_model_per_thin_AIC.csv",
        enm_res="results/enmeval/enmeval_res_0km_thin1_LQHPT.RData"
    output:
        expand("results/sdm/sdm_rangemap_best_{thresh}.{ext}", 
               thresh = ["specsens", "sens95", "sens99"],
               ext=["grd", "gri"])
    resources:
        cpus=1
    params:
        best_params="LQHPT_1"
    shell:
        """
        Rscript src/sdm_from_best_mod.R {input.enm_res} {params.best_params}
        """

## predict_current_phenotypes
rule predict_current_phenotypes:
    input:
        environment="processed_data/pheno_predictors_millsetal2018.tif",
        pheno_data="raw_data/Ltowsendii_database_FINAL.xlsx",
        gbif="raw_data/GBIF/verbatim.txt",
        range="results/sdm/sdm_rangemap_best_sens95.gri"
    output:
        "results/pheno/current_pheno_glm.RData",
        "results/pheno/current_predicted_probWhite.tif",
        "results/pheno/current_predicted_probWhite_SDMrange.tif"
    resources:
        cpus=1
    shell:
        """
        Rscript src/predict_current_pheno.R {input.pheno_data} {input.environment} {input.gbif} {input.range}
        """

## figure_1_maps
# make map of US and colorado for Figure 1
rule figure_1_maps:
    input:
        "results/pheno/current_predicted_probWhite_SDMrange.tif",
        "raw_data/DMNS_spectrometry_PCs.txt"






# --- Misc --- #
## dag               : makes a DAG graph for this pipeline
rule dag:
    input: "Snakefile"
    resources:
        cpus=1
    shell:
        "snakemake --dag | dot -Tpng > DAG.png"
        
## filegraph          : makes a file graph for this pipeline
rule filegraph:
    input: "Snakefile"
    shell:
        "snakemake --filegraph | dot -Tpng > filegraph.png"

## help               : prints help comments for this snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"
