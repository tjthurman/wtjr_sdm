This README explains the files in the processed_data folder.

* `bc23_CMIP5_RCP85_2080s_5modavg.gr*` - these files are raster files of BioClim variables 2 and 3, as predicted for the year 2080 under the RCP8.5 scenarion. They are an average of five climate models. These files are generated by `src/prep_future_bioclim_data.R`. They are too large to be tracked by Git. 

* `bg_points_for_sdm.Rdata`- This file is an Rdata object of the 10000 points used as background points in out SDM. It was generated by `src/prep_bg_points.R`. We track a copy on Git, so that the exact points we used can be re-used by others (though we also use a set seed for generating these random points, so they should be reproducible). 

* `bioclim_30arcsec_for_WTJR_SDM.tif` - This is a raster file of BioClim variables for our stydu area for the year 2000. It was generated by `src/prep_bioclim_data.R`. This file is too large to be tracked by Git. 

* `conservation_rasters/` - This folder contains raster files representing the areas of the 5 conservation statuses of WTJR. It was generated by `src/prep_conservation_rasters.R`, and is not tracked on Git. 

* `leptown_db_occurrences.csv` and `leptown_db_occurrences_unique.csv`- This files list the WTJR occurrences used in the SDM. The first lists all locations (including multiple records at the same location), while the second provides only the unique GPS coordinates. These files are generated by `src/data_curation.Rmd`, and are tracked on Git. 

* `pheno_predictors_millsetal2018.tif` - This is a raster file containing the variables (Bioclim 2, Bioclim 3, and snow cover) used to predict current phenotype for WTJR, taken from Mills et al. 2018 and processed by `src/prep_current_pheno_predictors.R`. This file is too large to track on Git. 

* `thin/`- This folder contains the distance-thinned records of WTJR occurrence, as generated by `src/prep_unthinned_occur.R` and `src/thin_records.R`. It is unclear to me if the package we used for thinning (spThin) is deterministic, so these files are tracked on Git so that the exact versions we used can be used by others. 