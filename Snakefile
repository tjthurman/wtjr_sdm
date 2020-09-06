## Snakemake workflow for WTJR SDM and phenotypic prediction analysis
## Pipeline order goes from top to bottom.

# Set up some wildcards for various aspects of the analysis
# Spatial thinning
thin_dists=[1,5,10,50]
thin_res_pattern= "processed_data/thin/{thin_dist}km/"


rule all:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv",
        "processed_data/bioclim_30arcsec_for_WTJR_SDM.tif",
         expand(thin_res_pattern, 
                thin_dist=thin_dists)


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
# 4 arguments
# inut file, thinning distance, number of reps, and output directory
# for number of reps, did 3*number of unique gps records
rule thin_occur_data:
    input:
        "processed_data/leptown_db_occurrences_unique_gps.csv"
    output:
        directory(thin_res_pattern)
    shell:
      """
      Rscript src/thin_records.R {input} 
      """
    
    
## process_bioclim_data: process bioclim data
rule process_bioclim_data:
    output:
        "processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
    shell:
        "Rscript src/prep_bioclim_data.R"


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
