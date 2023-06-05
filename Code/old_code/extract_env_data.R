###############################################################################
                    # Extract bioclim values #
###############################################################################

library(dplyr) # for easy data manipulation
library(stringr)
#library(stringi)
library(tidyr) # ditto
library(magrittr) # for piping
library(raster) # for working with raster data
#library(geosphere) # calculating geographic distance between sites
library(maptools)
library(rgdal)
library(sp)
library(sf)
library(rgeos)
library(foreach)
library(doParallel)
library(sm)
library(terra)

# Clear environment. 
rm(list=ls())

setwd("Manuscript/")

# Functions.
source("../Code/functions.R")


################################################################################
                              #### Data ####



# Get the environmental world clim rasters.
#envar <- getData("worldclim", var = "bio", res=2.5, download = FALSE, path = "..") #download the env rasters

envar <- rast(lapply(list.files("../Data/Chelsa/", full.names = TRUE), rast))
#change order()
envar <- rast(list(envar[[5]], # Mean Annual Temp
               envar[[6]], # Mean Dirunal Range
               envar[[7]], # Temperature seasonalilty.
               envar[[8]], # Max Temp warmest month
               envar[[9]],  # Min temp coldest month
               envar[[10]],  # Temp annual range
               envar[[1]], # Annual precipitation
               envar[[2]], # Precip wettest month.
               envar[[3]], # Precip driest month.
               envar[[4]], # Precip seasonality
               envar[[11]], # Season length
               envar[[12]], # Number of growning degree days above 10 degrees
               envar[[13]])) # Net primary productivity.

test <- rast(list(envar[[1]], envar[[2]]))
library(tictoc)
tic()
envar_be <- project(envar, behr)
toc()

writeRaster(envar_be, "test.tif")

stack("test.tif")

# envar <- stack(envar[[1]], # Mean Annual Temp
#                envar[[2]], # Mean Dirunal Range
#                envar[[4]], # Temperature seasonalilty.
#                envar[[5]], # Max Temp warmest month
#                envar[[6]],  # Min temp coldest month
#                envar[[7]],  # Temp annual range
#                envar[[12]], # Annual precipitation
#                envar[[13]], # Precip wettest month.
#                envar[[14]], # Precip driest month.
#                envar[[15]]) # Precip seasonality

# Read in growing degree days.

#bio1 <- raster("../Data/CHELSA_bio1_1981-2010_V.2.1.tif")
#plot(bio1)
#plot(envar[[1]])
#gdd <- raster("../Data/CHELSA_gdd5_1981-2010_V.2.1.tif")
#gdd <- resample(gdd, envar)
#envar <- stack(envar, gdd)

# Create the Behrmann's CRS.
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
behr <-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
registerDoParallel(cores = 7)
library(tictoc)
tic()
envar_be <- foreach(x = 1:13, .packages = "terra") %dopar% {
  loop_raster <- rast(list.files("../Data/Chelsa/", full.names = TRUE)[x])
  project(loop_raster, behr)
}
toc()
str(envar)
?pack

?project


rast(envar[[1]])
?pack
# Clean up any previous clusters.
registerDoSEQ()

# Set a cluster for raster parallel functions.
raster::beginCluster(13)

# Project to new coordinates.
envar_be <- projectRaster(envar, crs = behr)

raster::endCluster()
registerDoSEQ()

save(envar_be, file = "../Data/behr_raster_layers.RData") 

test_2 <- readRDS("../Data/behr_raster_layers.rds")
# Extract values from the raster.
var_vals <- raster::getValues(envar_be[[1]])

#create an empty normal raster for rasterizing in the loop
r <- raster(ncol = 42542, nrow = 13825, ymn = -60)
empty_be <- raster(ncol = 42542, nrow = 13825, ext = extent(envar_be), crs = crs(envar_be))

# #create an empty Behrmann's projection raster as a template for the loop
# empty_be <- projectRaster(r, crs = behr)
# 
# ?raster()
# 
# empty_be <- envar_be[[1]]
# values(empty_be)
# Load in the range maps.
# load("../../../Phd_data/Spatial/clean_jetz_ranges.RData")

load("../Data/clean_jetz_ranges.RData")

# Load in the trait data.
phylo_data <- read.csv("../Data/sexual_traits.csv")

# Join up the data.
sex_ranges <- left_join(phylo_data, new_jetz_ranges,  by = c("bird_tree_name" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)

# sex_ranges[111,]
# plot(sex_ranges[sex_ranges$bird_tree_name == "Ptilinopus chalcurus","Shape"], axes = TRUE)
# 
# length(na.omit(sex_ranges$Shape))
# 
# 
# plot(crop(envar[[10]], c(-149,-148,-16,-15)))
# plot(sex_ranges[sex_ranges$bird_tree_name == "Ptilinopus chalcurus","Shape"], add=TRUE)
# plot(crop(envar[[11]], c(-149,-148,-16,-15)))
# 
# 
# plot(crop(envar[[11]]))
# 
# extent(envar)
registerDoParallel(cores = 4)
sex_ranges <- sex_ranges[1:4,]


# #split <- test_ranges[1,]
library(fasterize)
#rasterOptions(memfrac=.3)
rasterOptions(maxmemory=1e+08)
names(envar_be)

# Use the model that extracts values first.
results <- foreach(split = isplit(sex_ranges, sex_ranges$bird_tree_name),
                   .combine = rbind, .errorhandling = "remove",
                   .packages = c("raster", "fasterize", "sf", "lwgeom")) %dopar% {
                     rasterOptions(maxmemory=1e+08)
                     
                  
                     
                     raster_i <- fasterize(st_as_sf(split$value$Shape), r) #fasterize polygon onto raster
                     #values <- getValues(envar_be[[1]])[which(getValues(projectRaster(fasterize(st_as_sf(test), r), empty_be)) > 0), ]  #fasterize polygon onto raster
                     
                     
                     #Reproject the polygon onto the empty_Behrman's raster
                     raster_i_be <- projectRaster(raster_i, empty_be) #reproject raster to equalarea Behrmann's
                     
                     #Select the projected cells of the raster and take the values at these cells by comparison to the var_vals array
                     rastercells <- which(getValues(raster_i_be) > 0) #select the cell
                     values <- getValues(envar_be)[rastercells, ]
                     
                     #some only have one rastercell and therefore only one value - take this value
                     #cannot divide by 10 if value is NA therefore ifelse 
                     if(nrow(values) == 1) {
                       ifelse(is.na(values[1]), 
                              bio1 <- NA, 
                              bio1 <-values[1]/10)
                       ifelse(is.na(values[2]), 
                              bio2 <- NA, 
                              bio2 <-values[2]/10)
                       ifelse(is.na(values[3]), 
                              bio4 <- NA, 
                              bio4 <-values[3])
                       ifelse(is.na(values[4]), 
                              bio5 <- NA, 
                              bio5 <-values[4]/10)
                       ifelse(is.na(values[5]), 
                              bio6 <- NA, 
                              bio6 <-values[5]/10)
                       ifelse(is.na(values[6]), 
                              bio7 <- NA, 
                              bio7 <-values[6]/10)
                       ifelse(is.na(values[7]), 
                              bio12 <- NA, 
                              bio12 <-values[7])
                       ifelse(is.na(values[8]), 
                              bio13 <- NA, 
                              bio13 <-values[8])
                       ifelse(is.na(values[9]), 
                              bio14 <- NA, 
                              bio14 <-values[9])
                       ifelse(is.na(values[10]), 
                              bio15 <- NA, 
                              bio15 <-values[10])
                       ifelse(is.na(values[11]), 
                              gsl <- NA, 
                              gsl <-values[11])
                       ifelse(is.na(values[12]), 
                              ngd <- NA, 
                              ngd <-values[12])
                       ifelse(is.na(values[13]), 
                              npp <- NA, 
                              npp <-values[13])
                       
                       
                     } else { 
                       #most have multiple so take the mean of all the values
                       #However cannot take the mean or divide by 10 if all vals are NA hense ifelse()
                       
                       ifelse(all(is.na(values[,1])), 
                              bio1 <- NA,
                              bio1 <-mean(values[,1], na.rm = TRUE)/10) #write them in
                       ifelse(all(is.na(values[,2])),
                              bio2 <- NA,
                              bio2 <-mean(values[,2], na.rm = TRUE)/10) #write them in
                       ifelse(all(is.na(values[,3])), 
                              bio4 <- NA,
                              bio4 <-mean(values[,3], na.rm = TRUE)) #write them in
                       ifelse(all(is.na(values[,4])),
                              bio5 <- NA,
                              bio5 <-mean(values[,4], na.rm = TRUE)/10) #write them in
                       ifelse(all(is.na(values[,5])), 
                              bio6 <- NA,
                              bio6 <-mean(values[,5], na.rm = TRUE)/10) #write them in
                       ifelse(all(is.na(values[,6])), 
                              bio7 <- NA,
                              bio7 <-mean(values[,6], na.rm = TRUE)/10) #write them in
                       ifelse(all(is.na(values[,7])), 
                              bio12 <- NA,
                              bio12 <-mean(values[,7], na.rm = TRUE)) #write them in
                       ifelse(all(is.na(values[,8])), 
                              bio13 <- NA,
                              bio13 <-mean(values[,8], na.rm = TRUE)) #write them in
                       ifelse(all(is.na(values[,9])), 
                              bio14 <- NA,
                              bio14 <-mean(values[,9], na.rm = TRUE)) #write them in
                       ifelse(all(is.na(values[,10])), 
                              bio15 <- NA,
                              bio15 <-mean(values[,10], na.rm = TRUE)) #write them in
                       ifelse(all(is.na(values[,11])), 
                              gsl <- NA,
                              gsl <-mean(values[,11], na.rm = TRUE))
                       ifelse(all(is.na(values[,12])), 
                              ngd <- NA,
                              ngd <-mean(values[,12], na.rm = TRUE))
                       ifelse(all(is.na(values[,13])), 
                              npp <- NA,
                              npp <-mean(values[,13], na.rm = TRUE))
                     }
                     
                     #get the study_results
                     data.frame(species = split$value$bird_tree_name[1],
                                bio1,
                                bio2,
                                bio4,
                                bio5,
                                bio6,
                                bio7,
                                bio12,
                                bio13,
                                bio14,
                                bio15,
                                gsl,
                                ngd,
                                npp)
                   }



##There will be some that have maps but such small ranges that there is no value
##Create a "No Maps" and a "Small Range" csv
small_range <- unique(sex_ranges$bird_tree_name)[-which(results$species %in%  unique(sex_ranges$bird_tree_name))]

write.csv(results, "../Data/jetz_environmental.csv")
write.csv(small_range, "../Data/resolution_misses_range.csv")





