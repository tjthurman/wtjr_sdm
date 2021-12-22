## Snakemake workflow for WTJR SDM and phenotypic prediction analysis
## Pipeline order goes from top to bottom, after rule all.

# Set up some wildcards for various aspects of the analysis
# Spatial thinning
thin_dists=[1,5,10,50]
dataset_dists=[0,1,5,10,50]
thin_res_pattern= "processed_data/thin/{thin_dist}km/"
current_bioclim="processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
# ENMeval
feature_class=["L", "LQ", "H", "LQH", "LQHP", "LQHPT"]
# Conservation
CON_STATUSES=["extirpated","broad_extirp","local_extirp","poss_decline","pres_stable"]


localrules: all, dag, filegraph, help


subworkflow slim_simulations:
    workdir:
        "/home/tt164677e/tim_beegfs/wtjr_sdm"
    snakefile:
        "/home/tt164677e/tim_beegfs/wtjr_sdm/slim_simulations_final.smk"

# Final desired outputs:
# 1) Performance plots from ENMeval (ensures that all ENMeval stuff runs)
# 2) the statistics reported in the main text
# 3) The elements of figure 1
# 4) The elements of figure 4
# 5) All the supplemental figures and tables. 
rule all:
    input:
        # Performance plots to ensure all thinning and ENMeval runs happen
        # and not just the final dataset
        expand("results/enmeval/performance_plot_{dist}km.pdf", dist = dataset_dists),
        # Slim simulations
        slim_simulations("results/slim_summaries/additive_constantK.csv"),
        slim_simulations("results/slim_summaries/additive_varyK.csv"),
        slim_simulations("results/slim_summaries/recessive_constantK.csv"),
        slim_simulations("results/slim_summaries/SSH_constantK.csv"),
        # Stats in main text
        "results/conservation/broad_chisq_res.csv", # Chisq results as table
        "results/conservation/cons_by_current_color.RData", # Chisq results as R object
        "results/pheno/percent_brown_by_time.csv", # Percent brown at different times
        # Figure 1 elements
        "results/figures/current_pheno_map.pdf", # Figure 1A
        "results/figures/colorado.pdf", # Figure 1C 
        # Figure 4 elements
        "results/figures/pheno_change_map.pdf", # Figure 4A
        "results/figures/horizontal_consv.pdf", # Figure 4B
        "results/figures/density_probBrown_insert.pdf", # Figure 4A insert
        # Supplementary figures, tables, and analysis
        "results/enmeval/performance_plot_best_models.pdf", # Figure S1
        "results/figures/supplemental/pheno_compare_maps.pdf", # Figure S24
        "results/figures/supplemental/model_difference_map.pdf", # Figure S25
        "results/pheno/model_difference_metrics.csv", # stats for Figure S25 caption
        "results/figures/supplemental/percent_brown_change.pdf", # Figure S26
        "results/figures/supplemental/discrete_current_pheno_map.pdf", # Figure S27
        "results/pheno/glm_table_current_snow_cover.csv", # Table S1A
        "results/pheno/glm_table_current_srt.csv", # Table S2B
        "results/pheno/glm_metrics_current_snow_cover.csv", # Table S2A
        "results/pheno/glm_metrics_current_srt.csv" # Table S2B
        
                
## curate_occur_data   : process WTJR occurrence data
# Renders an Rmarkdown report on the data curation process
# inputs and outputs are hard-coded into the data_curation.Rmd script
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
# for number of reps, did roughly 3*number of unique gps records
rule thin_occur_data:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    output:
        expand("processed_data/thin/wtjr_occ_{{thin_dist}}km_thin{index}.csv", index = ["1", "2", "3", "4", "5"])
    params:
        outdir="processed_data/thin/"
    resources:
        cpus=1
    shell:
      """
      Rscript src/thin_records.R {input} {wildcards.thin_dist} 3000 {params.outdir} 462
      """

## prep_unthinned_data: prepare unthinned data
# Takes unthinned data from data_curation
# and reformats and renames it to match the thinned data
rule prep_unthinned_data:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    output:
        "processed_data/thin/wtjr_occ_0km_thin1.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_unthinned_occur.R {input} {output}
        """

## prep_bioclim_data: process current bioclim data
# Download current bioclim data, merge and crop tiles
# into our study area. 
rule prep_bioclim_data:
    output:
        current_bioclim
    resources:
        cpus=1
    shell:
        "Rscript src/prep_bioclim_data.R {output}"

# prep_conservation_rasters: make rasters of conservation statuses
# Output filenames are hard-coded into the R script. 
rule prep_conservation_rasters:
    input:
        current_bioclim
    output:
        expand("processed_data/conservation_rasters/{cons_status}.tif", cons_status=CON_STATUSES)
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_conservation_rasters.R  {input}
        """
    
## prep_current_pheno_predictors: process Mills et al. 2018 predictors
# Crop the Mills et al. 2018 predictors down to area for this study
rule prep_current_pheno_predictors:
    input:
        snow="raw_data/rasterstack4cov.tif",
        bc=current_bioclim
    output:
        "processed_data/pheno_predictors_millsetal2018.tif"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_current_pheno_predictors.R {input.snow} {input.bc} {output}
        """

## prep_future_bioclim: process future bioclim data
# inputs: dirctory containing the raw data
# current bioclim data to use as a template
# output: model-averaged BC2 and BC3 for the 2080s
rule prep_future_bioclim:
    input:
        cmip5_dir="raw_data/worldclim_projections_2080s",
        current=current_bioclim
    output:
        "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_future_bioclim_data.R {input.cmip5_dir} {input.current} {output}
        """

## prep_bg_points: generate bg points for sdms
# Uses dismo R package to generate 10k points, with seed
rule prep_bg_points:
    input:
        envir=current_bioclim,
    output:
        "processed_data/bg_points_for_sdm.RData"
    resources:
        cpus=1
    shell:
        """
        Rscript src/prep_bg_points.R {input.envir} 782 {output}
        """

## run_ENMeval: run SDMs with ENMeval
# Inputs are occurrence datasets, environmental data, and bg points
# output is an .RData object for each dataset/featureclass combination
# same seed for each model
# Outfile template hard-coded into Rscript
rule run_ENMeval:
    input:
       occur="processed_data/thin/wtjr_occ_{dataset}_thin1.csv",
       envir=current_bioclim,
       bg= "processed_data/bg_points_for_sdm.RData"
    output:
        "results/enmeval/enmeval_res_{dataset}_{feature_class}.RData"
    resources:
        cpus=1
    shell:
        """
        Rscript src/run_ENMeval.R {input.occur} {input.envir} {input.bg} {wildcards.feature_class} 782 {resources.cpus}
        """

## process_ENMeval_res: process results from ENMeval, use to decide best model
# Collates per-dataset metrics into a single .csv (output)
# Gets performance metrics of best model for each thinning distance
# Makes performance plots of all models for each thinning distance

rule process_ENMeval_res:
    input:
        expand("results/enmeval/enmeval_res_{thin_dist}km_{fc}.RData", fc=feature_class, thin_dist=dataset_dists)
    output:
        expand("results/enmeval/performance_plot_{dist}km.pdf", dist = dataset_dists),
        all_models_metrics = "results/enmeval/enmeval_metrics.csv",
        best_models_metrics = "results/enmeval/enmeval_best_model_per_thin_AIC.csv",
        best_models_plot = "results/enmeval/performance_plot_best_models.pdf"
    params:
        indir="results/enmeval/"
    resources:
        cpus=1
    shell:
        """
        Rscript src/process_ENMeval_results.R {params.indir} {output.all_models_metrics} {output.best_models_metrics} {output.best_models_plot}
        """

## sdm_range_raster: make rasters of the SDM rangemaps from the best model
# Output name template hard coded into R script
# Best model parameter needed to be added here in params. 
rule sdm_range_raster: 
    input:
        "results/enmeval/enmeval_best_model_per_thin_AIC.csv",
        enm_res="results/enmeval/enmeval_res_0km_LQHPT.RData"        
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

## predict_current_phenotypes: predict winter white/brown for current times
# outputs both full predictions, and restricted to sdm range
rule predict_current_phenotypes:
    input:
        environment="processed_data/pheno_predictors_millsetal2018.tif",
        pheno_data="raw_data/Ltowsendii_database_FINAL.xlsx",
        gbif="raw_data/GBIF/verbatim.txt",
        range="results/sdm/sdm_rangemap_best_sens95.grd"
    output:
        glm = "results/pheno/current_pheno_glm.RData",
        full ="results/pheno/current_predicted_probWhite.tif",
        SDMrange = "results/pheno/current_predicted_probWhite_SDMrange.tif"
    resources:
        cpus=1
    shell:
        """
        Rscript src/predict_current_pheno.R {input.pheno_data} {input.environment} {input.gbif} {input.range} {output.glm} {output.full} {output.SDMrange}
        """

## figure_1_maps: make maps for panels of Fig. 1
# make base maps of US and Colorado for Figure 1
# outputs both pdfs and pngs of the maps for each panel of figure 1. 
rule figure_1_maps:
    input:
        pheno_range = "results/pheno/current_predicted_probWhite_SDMrange.tif",
        sample_coords = "raw_data/sample_coordinates_74individuals.txt"
    output:
        colo_pdf= "results/figures/colorado.pdf",
        colo_png= "results/figures/colorado.png",
        us_pdf = "results/figures/current_pheno_map.pdf",
        us_png = "results/figures/current_pheno_map.png"
    resources:
        cpus=1
    shell:
        """
        Rscript src/figure_1_maps.R {input.pheno_range} {input.sample_coords} {output.us_pdf} {output.us_png} {output.colo_pdf} {output.colo_png}
        """

## conserve_by_pheno: analysis of conservation status by phenotype
# calculate overlap between phenotypic categories and conservation categories
# conservation status rasters are hard-coded into script
rule conserve_by_pheno:
    input:
        expand("processed_data/conservation_rasters/{cons_status}.tif", cons_status=CON_STATUSES), 
        pheno_range = "results/pheno/current_predicted_probWhite_SDMrange.tif",
    output:
        "results/conservation/cons_by_current_color.RData"
    resources:
        cpus=1
    shell:
        """
        Rscript src/conservation_by_pheno.R {input.pheno_range} {output}
        """

## predict_future_pheno: predict winter white/brown for future
# also outputs current pheno estimated using SRT, with associated GLM
rule predict_future_phenotypes:
    input:
        pheno="raw_data/Ltowsendii_database_FINAL.xlsx",
        curr_bc="processed_data/bioclim_30arcsec_for_WTJR_SDM.tif",
        curr_srt="raw_data/SRT/SRT_historical/SRT_historical.tif",
        fut_bc="processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd",
        fut_srt="raw_data/SRT/SRT_RCP85_2080s/SRT_RCP85_2080s.tif",
        range="results/sdm/sdm_rangemap_best_sens95.grd",
        gbif="raw_data/GBIF/verbatim.txt"
    output:
        glm = "results/pheno/current_pheno_glm_SRT.RData",
        current_SRT = "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif",
        future_SRT = "results/pheno/future_predicted_probWhite_SDMrange.tif"
    resources:
        cpus=1
    shell:
        """
        Rscript src/predict_future_pheno.R {input.pheno} {input.curr_bc} {input.curr_srt} {input.fut_bc} {input.fut_srt} {input.range} {input.gbif} {output.glm} {output.current_SRT} {output.future_SRT}
        """
     
## figure_4_elements: make panels for Fig. 4
# Create panels for figure 4: map of phenotypic change as pdf and png
# conservation status by phenotype
# and density insert
# Plus a couple stats for the main text
rule figure_4_elements:
    input:
        current_srt_pheno = "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif",
        future_srt_pheno = "results/pheno/future_predicted_probWhite_SDMrange.tif",
        consv_pheno_overlap = "results/conservation/cons_by_current_color.RData"
    output:
        change_map_pdf= "results/figures/pheno_change_map.pdf",
        change_map_png= "results/figures/pheno_change_map.png",
        cons_plot = "results/figures/horizontal_consv.pdf",
        density_insert = "results/figures/density_probBrown_insert.pdf",
        prec_brown_by_time="results/pheno/percent_brown_by_time.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/figure_4_elements.R {input.current_srt_pheno} {input.future_srt_pheno} {input.consv_pheno_overlap} {output.change_map_pdf} {output.change_map_png} {output.cons_plot} {output.density_insert} {output.prec_brown_by_time}
        """

## Make supplementary figure and table rule
# Creates all supplemental tables and figures
# And outputs some bits of analysis for the main text
rule supplemental_and_analysis:
    input: 
        snowcover_glm_rdata = "results/pheno/current_pheno_glm.RData",
        srt_glm_rdata = "results/pheno/current_pheno_glm_SRT.RData",
        consv_pheno_overlap = "results/conservation/cons_by_current_color.RData",
        current_srt_pheno = "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif",
        future_srt_pheno = "results/pheno/future_predicted_probWhite_SDMrange.tif",
        current_cover_pheno = "results/pheno/current_predicted_probWhite_SDMrange.tif"
    output:
        snow_cover_table = "results/pheno/glm_table_current_snow_cover.csv",
        snow_cover_metrics = "results/pheno/glm_metrics_current_snow_cover.csv",
        srt_table = "results/pheno/glm_table_current_srt.csv",
        srt_metrics = "results/pheno/glm_metrics_current_srt.csv",
        broad_chisq_res = "results/conservation/broad_chisq_res.csv",
        pheno_compare_map = "results/figures/supplemental/pheno_compare_maps.pdf",
        glm_method_compare_map = "results/figures/supplemental/model_difference_map.pdf",
        percent_brown_change = "results/figures/supplemental/percent_brown_change.pdf",
        discrete_current_pheno_map = "results/figures/supplemental/discrete_current_pheno_map.pdf",
        glm_method_compare_metrics = "results/pheno/model_difference_metrics.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/supp_figs_and_tables.R {input.snowcover_glm_rdata} {input.srt_glm_rdata} {input.consv_pheno_overlap} {input.current_srt_pheno} {input.future_srt_pheno} {input.current_cover_pheno} {output.snow_cover_table} {output.snow_cover_metrics} {output.srt_table} {output.srt_metrics} {output.broad_chisq_res} {output.pheno_compare_map} {output.glm_method_compare_map} {output.percent_brown_change} {output.discrete_current_pheno_map} {output.glm_method_compare_metrics}
        """
 



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
