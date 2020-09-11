## Snakemake workflow for WTJR SDM and phenotypic prediction analysis
## Pipeline order goes from top to bottom.

# Set up some wildcards for various aspects of the analysis
# Spatial thinning
thin_dists=[1,5,10,50]
thin_dists_test=[50]
thin_res_pattern= "processed_data/thin/{thin_dist}km/"
current_bioclim="processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
# ENMeval
feature_class=["L", "LQ"]



rule all:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv",
        current_bioclim,
        expand(thin_res_pattern, thin_dist=thin_dists),
        "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd",
        expand("results/enmeval_res_{thin_dist}km_thin1_{fc}.RData", fc=feature_class, thin_dist=thin_dists_test), 
        


## curate_occur_data   : process WTJR occurrence data
rule curate_occur_data:
    input:
        arctos="raw_data/ArctosData_43FA2173A0.csv",
        gbif= "raw_data/GBIF/verbatim.txt",
        vertnet= "raw_data/vertnet_leptownsendii_allrecords_apr2_2020.txt"
    output: 
        "processed_data/leptown_db_occurrences.csv",
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    shell:
        "Rscript -e \"rmarkdown::render('src/data_curation.Rmd', knit_root_dir = '../', output_dir = 'docs/')\""
        
## thin_occur_data  : spatial thinning of WTJR occurrence data
# 5 arguments
# input file, thinning distance, number of reps, output directory, and random seed
# for number of reps, did 3*number of unique gps records
checkpoint thin_occur_data:
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


# --- Misc --- #
## dag               : makes a DAG graph for this pipeline
rule dag:
    input: "Snakefile"
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
