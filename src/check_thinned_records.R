

# Load packages -----------------------------------------------------------
library(fields)
library(stringr)


# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 directory containing thinning results
# Should work in top-level thin, or in the thin1/thin5/thin10 ones
dir <- args[1]

# Check the thinning outputs ----------------------------------------------


for (file in list.files(dir, recursive = T, full.names = T)) {
  if (str_detect(file, pattern = "csv$")) {
    thinned_data <- read.csv(file)
    thin_dist <- as.numeric(str_remove(str_extract(string = file, pattern = "\\d+km"), pattern = "km"))
    dists <- fields::rdist.earth(x1 = thinned_data[,2:3], miles = F) # GOTTA BE IN KM, NOT MILES
    # Set diags (comparisons to self) to something large, so the don't get picked up
    diag(dists) <- 500
    too_close <- sum(dists < thin_dist)
    if (too_close != 0) {
      stop(paste("Some records too close together in file:\n",
                 file))
    }
  }
}
