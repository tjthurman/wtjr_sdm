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

# Sliim sim scenarios
slim_scenarios = ["additive_constantK_2locus", "additive_varyK_2locus", "recessive_constantK_1locus", "recessive_constantK_2locus"]
# Slim sim file ext
slim_exts = ["early.csv", "late.csv", "seed.csv"]


localrules: all


# May need to change the working directory and path to the snakefile to work on your machine
subworkflow slim_simulations:
    workdir:
        "."
    snakefile:
        "./slim_simulations.smk"

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
        # Do all output files, to ensure everything (not just late) files compile
        slim_simulations(expand("results/slim_summaries/{scen}_{ext}", scen = slim_scenarios, ext = slim_exts)),
        # Stats in main text
        "results/conservation/broad_chisq_res.csv", # Chisq results as table
        "results/conservation/cons_by_current_color.RData", # Chisq results as R object
        "results/pheno/percent_brown_by_time.csv", # Percent brown at different times
        # Figure 1 elements
        "results/figures/current_pheno_map_120mm.pdf", # Figure 1A
        "results/figures/colorado_120mm.pdf", # Figure 1C 
        # Figure 4 elements
        "results/figures/pheno_change_map_55mm.pdf", # Figure 4A
        "results/figures/sim_pop_trajectories_55mm.pdf", # Figure 4B
        "results/figures/density_probBrown_insert_55mm.pdf", # Figure 4A insert
        # Supplementary figures, tables, and analysis
        "results/pheno/glm_table_current_snow_cover.csv",
        "results/pheno/glm_metrics_current_snow_cover.csv",
        "results/pheno/glm_table_current_srt.csv",
        "results/pheno/glm_metrics_current_srt.csv",
        "results/conservation/broad_chisq_res.csv",
        "results/figures/supplemental/extended_data_SDMs.pdf",
        "results/figures/supplemental/extended_data_SDMs.jpeg",
        "results/pheno/model_difference_metrics.csv",
        "results/figures/supplemental/extended_data_arch_varyK_sim_res.pdf",
        "results/figures/supplemental/extended_data_arch_varyK_sim_res.jpeg",
        "results/figures/supplemental/extended_data_sim_robust.pdf",
        "results/figures/supplemental/extended_data_sim_robust.jpeg",
        "results/figures/supplemental/sampling_map.pdf",
        "results/figures/supplemental/sampling_map.jpeg"

      
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
        best_models_metrics = "results/enmeval/enmeval_best_model_per_thin_AIC.csv"
    params:
        indir="results/enmeval/"
    resources:
        cpus=1
    shell:
        """
        Rscript src/process_ENMeval_results.R {params.indir} {output.all_models_metrics} {output.best_models_metrics}
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
        colo_pdf= "results/figures/colorado_120mm.pdf",
        colo_png= "results/figures/colorado_120mm.png",
        us_pdf = "results/figures/current_pheno_map_120mm.pdf",
        us_png = "results/figures/current_pheno_map_120mm.png"
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
        slim_res_additive = slim_simulations("results/slim_summaries/additive_constantK_2locus_late.csv"),
        slim_res_recessive = slim_simulations("results/slim_summaries/recessive_constantK_2locus_late.csv")
    output:
        change_map_pdf= "results/figures/pheno_change_map_55mm.pdf",
        change_map_png= "results/figures/pheno_change_map_55mm.png",
        density_insert = "results/figures/density_probBrown_insert_55mm.pdf",
        pop_traj_plot_out = "results/figures/sim_pop_trajectories_55mm.pdf",
        prec_brown_by_time= "results/pheno/percent_brown_by_time.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/figure_4_elements.R {input.current_srt_pheno} {input.future_srt_pheno} {input.slim_res_additive} {input.slim_res_recessive} {output.change_map_pdf} {output.change_map_png} {output.density_insert} {output.pop_traj_plot_out} {output.prec_brown_by_time}
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
        current_cover_pheno = "results/pheno/current_predicted_probWhite_SDMrange.tif",
        best_metrics_file = "results/enmeval/enmeval_best_model_per_thin_AIC.csv",
        additive_consK_2locus_late_file = "results/slim_summaries/additive_constantK_2locus_late.csv",
        recessive_consK_2locus_late_file = "results/slim_summaries/recessive_constantK_2locus_late.csv",
        additive_varyK_2locus_late_file = "results/slim_summaries/additive_varyK_2locus_late.csv",
        recessive_consK_1locus_late_file = "results/slim_summaries/recessive_constantK_1locus_late.csv",
        shapefile = "raw_data/redlist_species_data_79e8f518-a14c-4d0e-9640-5bea84d7c1b8/data_0.shp",
        occurrence_csv = "processed_data/thin/wtjr_occ_0km_thin1.csv",
        sdm_rangemap = "results/sdm/sdm_rangemap_best_sens95.grd",
        gen_sample_info = "raw_data/genetics/samples_info.csv",
        gwas_sample_file = "raw_data/genetics/WTJR_74lowcovsamples_code_disamb.xlsx"
    output:
        snow_cover_table = "results/pheno/glm_table_current_snow_cover.csv",
        snow_cover_metrics = "results/pheno/glm_metrics_current_snow_cover.csv",
        srt_table = "results/pheno/glm_table_current_srt.csv",
        srt_metrics = "results/pheno/glm_metrics_current_srt.csv",
        broad_chisq_res = "results/conservation/broad_chisq_res.csv",
        ext_data_sdm_fig = "results/figures/supplemental/extended_data_SDMs.pdf", 
        ext_data_sdm_fig_jpg = "results/figures/supplemental/extended_data_SDMs.jpeg", 
        glm_method_compare_metrics = "results/pheno/model_difference_metrics.csv",
        arch_varyK_fig_pdf = "results/figures/supplemental/extended_data_arch_varyK_sim_res.pdf",
        arch_varyK_fig_jpeg = "results/figures/supplemental/extended_data_arch_varyK_sim_res.jpeg",
        robust_fig_pdf = "results/figures/supplemental/extended_data_sim_robust.pdf",
        robust_fig_jpeg = "results/figures/supplemental/extended_data_sim_robust.jpeg",
        sampling_map_pdf = "results/figures/supplemental/sampling_map.pdf",
        sampling_map_jpeg = "results/figures/supplemental/sampling_map.jpeg"
    resources:
        cpus=1
    shell:
        """
        Rscript src/supp_figs_and_tables.R {input.snowcover_glm_rdata} {input.srt_glm_rdata} {input.consv_pheno_overlap} {input.current_srt_pheno} {input.future_srt_pheno} {input.current_cover_pheno} {input.best_metrics_file} {input.additive_consK_2locus_late_file} {input.recessive_consK_2locus_late_file} {input.additive_varyK_2locus_late_file} {input.recessive_consK_1locus_late_file} {input.shapefile} {input.occurrence_csv} {input.sdm_rangemap} {input.gen_sample_info} {input.gwas_sample_file} {output.snow_cover_table} {output.snow_cover_metrics} {output.srt_table} {output.srt_metrics} {output.broad_chisq_res} {output.ext_data_sdm_fig} {output.ext_data_sdm_fig_jpg} {output.glm_method_compare_metrics} {output.arch_varyK_fig_pdf} {output.arch_varyK_fig_jpeg} {output.robust_fig_pdf} {output.robust_fig_jpeg} {output.sampling_map_pdf} {output.sampling_map_jpeg}
        """
