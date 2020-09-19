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


localrules: all, rename_ENMeval_res, dag, filegraph, help


rule all:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv",
        current_bioclim,
        # expand(thin_res_pattern, thin_dist=thin_dists),
        "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd",
        expand("results/enmeval_res_{thin_dist}km_thin1_{fc}.RData", fc=feature_class, thin_dist=["0", "10", "50"]),
        expand("results/enmeval_res_1km_thin3_{fc}.RData", fc=feature_class),
        expand("results/enmeval_res_5km_thin2_{fc}.RData", fc=feature_class)
        
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
    shell:
      """
      Rscript src/thin_records.R {input} {wildcards.thin_dist} 3000 {output} 462
      """

## check_thin_data  : Double-check that thinning worked properly
rule check_thin_data:
    input:
        "processed_data/thin/"
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
    shell:
        """
        Rscript src/prep_unthinned_occur.R {input} {output}
        """
    
## process_bioclim_data: process current bioclim data
rule process_bioclim_data:
    output:
        current_bioclim
    shell:
        "Rscript src/prep_bioclim_data.R"

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
        "results/enmeval_res_{dataset}_thin1_{feature_class}.RData"
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
        expand("results/enmeval_res_{km}km_thin1_{fc}.RData", km=["1", "5"], fc=feature_class)
    output:
        expand("results/enmeval_res_1km_thin3_{fc}.Rdata", fc=feature_class),
        expand("results/enmeval_res_5km_thin2_{fc}.Rdata", fc=feature_class)
    shell: 
        """
        rename thin1 thin3 results/enmeval_res_1km_*
        rename thin1 thin2 results/enmeval_res_5km_*
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
