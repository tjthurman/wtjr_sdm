
library(tidyverse)
# library(maptools)
# library(rgdal)
# library(maps)
library(raster)
library(dismo)
# library(maxnet)
library(ENMeval)

load("results/enmeval_res_0km_thin1_LQHPT.RData")
 
# Extract predictions from best model
best.preds <- subset(enmeval_res@predictions, "LQHPT_1")

# evaluate the best model
eval.best.mod <- dismo::evaluate(p = raster::extract(x = best.preds, y = enmeval_res@occ.pts),
                                 a = raster::extract(x = best.preds, y = enmeval_res@bg.pts))
# Use evaluation to get a threshold
spec_sens_thresh <- dismo::threshold(eval.best.mod, stat = "spec_sens")
sens_95_thresh <- dismo::threshold(eval.best.mod, stat = "sensitivity", sensitivity = 0.95)
sens_99_thresh <- dismo::threshold(eval.best.mod, stat = "sensitivity", sensitivity = 0.99)

# Get sensitivity at spec_sens thresh
eval.best.mod@TPR[eval.best.mod@t == spec_sens_thresh]
# Get specificity at spec_sens thresh
eval.best.mod@TNR[eval.best.mod@t == spec_sens_thresh]
eval.best.mod@TNR[eval.best.mod@t == sens_95_thresh]
eval.best.mod@TNR[eval.best.mod@t == sens_99_thresh]



# Make the range raster and plot it, saving as both a png and a pdf
range.raster.spec_sens <- best.preds > spec_sens_thresh
range.raster.sens95 <- best.preds > sens_95_thresh
range.raster.sens99 <- best.preds > sens_99_thresh

# Save as a binary pres/absence file --------------------------------------
values(range.raster.spec_sens) <- ifelse(values(range.raster.spec_sens) < 0.5, NA, 1)
values(range.raster.sens95) <- ifelse(values(range.raster.sens95) < 0.5, NA, 1)
values(range.raster.sens99) <- ifelse(values(range.raster.sens99) < 0.5, NA, 1)

plot(range.raster.spec_sens)
plot(range.raster.sens95)
plot(range.raster.sens99)


writeRaster(range.raster.spec_sens, filename = "results/sdm_rangemap_best_specsens", overwrite = T)
writeRaster(range.raster.sens95, filename = "results/sdm_rangemap_best_sens95", overwrite = T)
writeRaster(range.raster.sens99, filename = "results/sdm_rangemap_best_sens99", overwrite = T)
