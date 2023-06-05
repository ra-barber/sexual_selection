###############################################################################
                     # Extract chelsa values #
###############################################################################

# Packages.
library(dplyr)
library(raster)
library(sf)
library(foreach)
library(doParallel)
library(terra)

# Clear environment. 
rm(list=ls())

#setwd("Manuscript/")

# For on the cluster.
# setwd("..")

# Functions.
#source("../Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
#array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number <- 1

# Get raster names.
raster_name <- c("bio1_", "bio2", "bio4","bio5", "bio6", "bio7", "bio12", "bio13", "bio14",
"bio15", "gsl", "ngd", "npp")[array_number]



################################################################################
                        #### Data ####

# Read in the phylo data.
phylo_data <- read.csv("../Data/sexual_traits.csv")

# Load in the range maps.
#load("../../../Phd_data/Spatial/clean_jetz_ranges.RData")
load("../Data/clean_jetz_ranges.RData")

# Join up the data.
sex_ranges <- left_join(phylo_data, new_jetz_ranges,  by = c("bird_tree_name" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)


################################################################################
                    #### Prepare loop ####

library(tictoc)
tic()

# Create the Behrmann's CRS.
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

# Register cores.
registerDoParallel(cores = 100)

# Foreach loop to extract climate data.
results <- foreach(split = isplit(sex_ranges, sex_ranges$bird_tree_name),
                   .packages = c("terra", "sf" ,"raster"),
                   .errorhandling = "remove", .verbose = TRUE,
                   .combine = "rbind") %dopar% {
                     
                     # This file is absolutely huge.
                     envar_be <- stack("../Data/envar_behr.tif")
                     
                     # Get specific layer.
                     envar_be <- rast(envar_be[[grep(raster_name, names(envar_be))]])
                     
                     # Pull out species name.
                     species_name <- split$value$bird_tree_name
                     
                     # Project the polygon.
                     poly_be <- project(vect(split$value$Shape), behr)
                     
                     # Extract the data.
                     mean_value <- extract(envar_be, poly_be, mean, na.rm=TRUE)
                     
                     extracted_data <- data.frame(bird_tree_name = species_name,
                                                  variable = mean_value[[2]])
                     
                     # Change raster name now.
                     if (raster_name == "bio1_"){
                       colnames(extracted_data)[2] <- "bio1"
                     } else{
                       colnames(extracted_data)[2] <- raster_name
                     }
                     extracted_data
                   }

# Change name for pathway.
if (raster_name == "bio1_"){
  raster_name <- "bio1"
}

# Save data.
pathway <- paste0("../Data/Extracts/jetz_", raster_name, ".csv")
write.csv(results, pathway, row.names = FALSE)



# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# sex_ranges <- sex_ranges[1:8,]
# #x <- 1
# # Foreach loop to extract climate data.
results <- foreach(x = 1:nrow(sex_ranges),
                   .packages = c("terra", "sf" ,"raster"),
                   .combine = "rbind") %dopar% {

                     # This file is absolutely huge.
                     envar_be <- stack("../Data/envar_behr.tif")

                     # Get specific layer.
                     envar_be <- rast(envar_be[[grep(raster_name, names(envar_be))]])

                     species_name <- sex_ranges$bird_tree_name[x]

                     # Project the polygon.
                     poly_be <- project(vect(sex_ranges$Shape[x]), behr)

                     # Extract the data.
                     mean_value <- extract(envar_be, poly_be, mean, na.rm=TRUE)

                     extracted_data <- data.frame(bird_tree_name = species_name,
                                                  variable = mean_value[[2]])

                     # Change raster name now.
                     if (raster_name == "bio1_"){
                       colnames(extracted_data)[2] <- "bio1"
                     } else{
                       colnames(extracted_data)[2] <- raster_name
                     }

                     extracted_data

                   }


# if (raster_name == "bio1_"){
#   raster_name <- "bio1"
# }
# 
# pathway <- paste0("../Data/Extracts/jetz_", raster_name, ".csv")
# 
# 
# write.csv(results, pathway, row.names = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Foreach loop to extract climate data.
# results <- foreach(split = isplit(sex_ranges, sex_ranges$bird_tree_name),
#                    .packages = c("terra", "sf" ,"raster"), 
#                    .combine = "rbind") %dopar% {
#   
#   # This file is absolutely huge.
#   envar_be <- stack("../Data/envar_behr.tif")
#   
#   # Get specific layer.
#   envar_be <- rast(envar_be[[grep(raster_name, names(envar_be))]])
#   
#   species_name <- split$value$bird_tree_name
#   
#   # Project the polygon.
#   poly_be <- project(vect(split$value$Shape), behr)
#   
#   # Extract the data.
#   mean_value <- extract(envar_be, poly_be, mean, na.rm=TRUE)
#   
#   extracted_data <- data.frame(bird_tree_name = species_name,
#                                variable = mean_value[[2]])
#   
#   # Change raster name now.
#   if (raster_name == "bio1_"){
#     colnames(extracted_data)[2] <- "bio1"
#   } else{
#     colnames(extracted_data)[2] <- raster_name
#   }
#   extracted_data
# }
# 
# 
# if (raster_name == "bio1_"){
#   raster_name <- "bio1"
# }
# 
# pathway <- paste0("../Data/Extracts/jetz_", raster_name, ".csv")
# 
# 
# write.csv(results, pathway, row.names = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Project the polygon.
# poly_be <- project(vect(test$Shape), behr)
# 
# # Extract the data.
# mean_value <- extract(rast(envar_be), poly_be, mean, na.rm=TRUE)
# 
# extracted_data <- data.frame(bird_tree_name = species_name,
#                              variable = mean_value)
# 
# colnames(extracted_data)[2] <- raster_name
# toc()
# registerDoParallel(cores = 4)
# results <- foreach(x = 1:4, .packages = c("terra", "sf" ,"raster"), .combine = "rbind") %dopar% {
#   
#   
#   # This file is absolutely huge.
#   envar_be <- stack("../Data/envar_behr.tif")
#   
#   # Get specific layer.
#   envar_be <- rast(envar_be[[grep(raster_name, names(envar_be))]])
#   
#   species_name <- test$bird_tree_name[x]
#   
#   # Project the polygon.
#   poly_be <- project(vect(test$Shape[x]), behr)
#   
#   # Extract the data.
#   mean_value <- extract(envar_be, poly_be, mean, na.rm=TRUE)
#   
#   extracted_data <- data.frame(bird_tree_name = species_name,
#                                variable = mean_value[[2]])
#   
#   # Change raster name now.
#   if (raster_name == "bio1_"){
#     colnames(extracted_data)[2] <- "bio1"
#   } else{
#     colnames(extracted_data)[2] <- raster_name
#   }
#   
#   
# }
# 
# if (raster_name == "bio1_"){
#   raster_name <- "bio1"
# }
# 
# pathway <- paste0("../Data/Extracts/jetz_", raster_name, ".csv")
# 
# 
# write.csv(results, pathway, row.names = FALSE)
# 
# 
# ##There will be some that have maps but such small ranges that there is no value
# ##Create a "No Maps" and a "Small Range" csv
# small_range <- unique(sex_ranges$bird_tree_name)[-which(results$species %in%  unique(sex_ranges$bird_tree_name))]
# 
# 
# write.csv(results, "../Data/jetz_environmental.csv")
# write.csv(small_range, "../Data/resolution_misses_range.csv")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # This file is absolutely huge.
# envar_be <- stack("../Data/envar_behr.tif")
# 
# # Get specific layer.
# envar_be <- rast(envar_be[[grep(raster_name, names(envar_be))]])
# 
# species_name <- test$bird_tree_name[x]
# 
# # Project the polygon.
# poly_be <- project(vect(test$Shape[x]), behr)
# 
# # Extract the data.
# mean_value <- extract(envar_be, poly_be, mean, na.rm=TRUE)
# 
# extracted_data <- data.frame(bird_tree_name = species_name,
#                              variable = mean_value)
# 
# 
# 
# #empty_be <- raster(ncol = 42542, nrow = 13825, ext = extent(envar_be), crs = crs(envar_be))
# 
# 
# 
# 
# 
# # Extract values from the raster.
# #var_vals <- raster::getValues(envar_be)
# 
# #create an empty normal raster for rasterizing in the loop
# r <- raster(ncol = 42542, nrow = 13825, ymn = -60)
# empty_be <- raster(ncol = 42542, nrow = 13825, ext = extent(envar_be), crs = crs(envar_be))
# 
# 
# registerDoParallel(cores = 4)
# sex_ranges <- sex_ranges[1:4,]
# 
# test <- sex_ranges[1:4,]
# # #split <- test_ranges[1,]
# library(fasterize)
# #rasterOptions(memfrac=.3)
# rasterOptions(maxmemory=1e+08)
# names(envar_be)
# 
# 
# poly_be <- project(vect(test$Shape), behr)
# 
# extract(rast(envar_be), poly_be, mean, na.rm=TRUE)
# 
# 
# vect(test$Shape)
# 
# species_raster <- fasterize(st_as_sf(test$Shape), r)
# 
# project(rast(species_raster), empty_be)
# 
# ?fasterize
# #Reproject the polygon onto the empty_Behrman's raster
# library(tictoc)
# tic()
# raster_i_be <- projectRaster(species_raster, empty_be) 
# toc()
# 
# rastercells <- which(getValues(raster_i_be) > 0) #select the cell
# values <- getValues(envar_be)[rastercells, ]
# 
# if(nrow(values) == 1) {
#   ifelse(is.na(values), 
#          env_value <- NA, 
#          env_value <- values)
#   
# } else { 
#   ifelse(all(is.na(values[,1])), 
#          env_value <- NA,
#          env_value <-mean(values[,1], na.rm = TRUE)/10) 
#   
# }
# 
# #some only have one rastercell and therefore only one value - take this value
# #cannot divide by 10 if value is NA therefore ifelse 
# if(nrow(values) == 1) {
#   ifelse(is.na(values[1]), 
#          bio1 <- NA, 
#          bio1 <-values[1]/10)
#   ifelse(is.na(values[2]), 
#          bio2 <- NA, 
#          bio2 <-values[2]/10)
#   ifelse(is.na(values[3]), 
#          bio4 <- NA, 
#          bio4 <-values[3])
#   ifelse(is.na(values[4]), 
#          bio5 <- NA, 
#          bio5 <-values[4]/10)
#   ifelse(is.na(values[5]), 
#          bio6 <- NA, 
#          bio6 <-values[5]/10)
#   ifelse(is.na(values[6]), 
#          bio7 <- NA, 
#          bio7 <-values[6]/10)
#   ifelse(is.na(values[7]), 
#          bio12 <- NA, 
#          bio12 <-values[7])
#   ifelse(is.na(values[8]), 
#          bio13 <- NA, 
#          bio13 <-values[8])
#   ifelse(is.na(values[9]), 
#          bio14 <- NA, 
#          bio14 <-values[9])
#   ifelse(is.na(values[10]), 
#          bio15 <- NA, 
#          bio15 <-values[10])
#   ifelse(is.na(values[11]), 
#          gsl <- NA, 
#          gsl <-values[11])
#   ifelse(is.na(values[12]), 
#          ngd <- NA, 
#          ngd <-values[12])
#   ifelse(is.na(values[13]), 
#          npp <- NA, 
#          npp <-values[13])
#   
#   
# 
# 
# 
# 
# 
# 
# 
# 
# # Use the model that extracts values first.
# results <- foreach(split = isplit(sex_ranges, sex_ranges$bird_tree_name),
#                    .combine = rbind, .errorhandling = "remove",
#                    .packages = c("raster", "fasterize", "sf", "lwgeom")) %dopar% {
#                      rasterOptions(maxmemory=1e+08)
#                      
# 
#                      
#                      raster_i <- fasterize(st_as_sf(split$value$Shape), r) #fasterize polygon onto raster
#                      #values <- getValues(envar_be[[1]])[which(getValues(projectRaster(fasterize(st_as_sf(test), r), empty_be)) > 0), ]  #fasterize polygon onto raster
#                      
#                      
#                      #Reproject the polygon onto the empty_Behrman's raster
#                      raster_i_be <- projectRaster(raster_i, empty_be) #reproject raster to equalarea Behrmann's
#                      
#                      #Select the projected cells of the raster and take the values at these cells by comparison to the var_vals array
#                      rastercells <- which(getValues(raster_i_be) > 0) #select the cell
#                      values <- getValues(envar_be)[rastercells, ]
#                      
#                      #some only have one rastercell and therefore only one value - take this value
#                      #cannot divide by 10 if value is NA therefore ifelse 
#                      if(nrow(values) == 1) {
#                        ifelse(is.na(values[1]), 
#                               bio1 <- NA, 
#                               bio1 <-values[1]/10)
#                        ifelse(is.na(values[2]), 
#                               bio2 <- NA, 
#                               bio2 <-values[2]/10)
#                        ifelse(is.na(values[3]), 
#                               bio4 <- NA, 
#                               bio4 <-values[3])
#                        ifelse(is.na(values[4]), 
#                               bio5 <- NA, 
#                               bio5 <-values[4]/10)
#                        ifelse(is.na(values[5]), 
#                               bio6 <- NA, 
#                               bio6 <-values[5]/10)
#                        ifelse(is.na(values[6]), 
#                               bio7 <- NA, 
#                               bio7 <-values[6]/10)
#                        ifelse(is.na(values[7]), 
#                               bio12 <- NA, 
#                               bio12 <-values[7])
#                        ifelse(is.na(values[8]), 
#                               bio13 <- NA, 
#                               bio13 <-values[8])
#                        ifelse(is.na(values[9]), 
#                               bio14 <- NA, 
#                               bio14 <-values[9])
#                        ifelse(is.na(values[10]), 
#                               bio15 <- NA, 
#                               bio15 <-values[10])
#                        ifelse(is.na(values[11]), 
#                               gsl <- NA, 
#                               gsl <-values[11])
#                        ifelse(is.na(values[12]), 
#                               ngd <- NA, 
#                               ngd <-values[12])
#                        ifelse(is.na(values[13]), 
#                               npp <- NA, 
#                               npp <-values[13])
#                        
#                        
#                      } else { 
#                        #most have multiple so take the mean of all the values
#                        #However cannot take the mean or divide by 10 if all vals are NA hense ifelse()
#                        
#                        ifelse(all(is.na(values[,1])), 
#                               bio1 <- NA,
#                               bio1 <-mean(values[,1], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,2])),
#                               bio2 <- NA,
#                               bio2 <-mean(values[,2], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,3])), 
#                               bio4 <- NA,
#                               bio4 <-mean(values[,3], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,4])),
#                               bio5 <- NA,
#                               bio5 <-mean(values[,4], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,5])), 
#                               bio6 <- NA,
#                               bio6 <-mean(values[,5], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,6])), 
#                               bio7 <- NA,
#                               bio7 <-mean(values[,6], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,7])), 
#                               bio12 <- NA,
#                               bio12 <-mean(values[,7], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,8])), 
#                               bio13 <- NA,
#                               bio13 <-mean(values[,8], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,9])), 
#                               bio14 <- NA,
#                               bio14 <-mean(values[,9], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,10])), 
#                               bio15 <- NA,
#                               bio15 <-mean(values[,10], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,11])), 
#                               gsl <- NA,
#                               gsl <-mean(values[,11], na.rm = TRUE))
#                        ifelse(all(is.na(values[,12])), 
#                               ngd <- NA,
#                               ngd <-mean(values[,12], na.rm = TRUE))
#                        ifelse(all(is.na(values[,13])), 
#                               npp <- NA,
#                               npp <-mean(values[,13], na.rm = TRUE))
#                      }
#                      
#                      #get the study_results
#                      data.frame(species = split$value$bird_tree_name[1],
#                                 bio1,
#                                 bio2,
#                                 bio4,
#                                 bio5,
#                                 bio6,
#                                 bio7,
#                                 bio12,
#                                 bio13,
#                                 bio14,
#                                 bio15,
#                                 gsl,
#                                 ngd,
#                                 npp)
#                    }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Use the model that extracts values first.
# results <- foreach(split = isplit(sex_ranges, sex_ranges$bird_tree_name),
#                    .combine = rbind, .errorhandling = "remove",
#                    .packages = c("raster", "fasterize", "sf", "lwgeom")) %dopar% {
#                      rasterOptions(maxmemory=1e+08)
#                      
#                      ##Convert the maps into multipolygon format
#                      # for (j in 1:length(split$value$Shape)) { 
#                      #   split$value$Shape[j] <- lwgeom_make_valid(split$value$Shape[j])
#                      #   split$value$Shape[[j]] <- st_cast(split$value$Shape[[j]], 'MULTIPOLYGON')  #changes the shapes in the feature
#                      # }
#                      # split$value$Shape <- st_cast(split$value$Shape, 'MULTIPOLYGON')
#                      # 
#                      #Rasterize the polygons
#                      #polygon <- st_combine(split$value$Shape) #combine to one shape layer (union not required as fasterize sorts out the overlaps)
#                      
#                      
#                      raster_i <- fasterize(st_as_sf(split$value$Shape), r) #fasterize polygon onto raster
#                      #values <- getValues(envar_be[[1]])[which(getValues(projectRaster(fasterize(st_as_sf(test), r), empty_be)) > 0), ]  #fasterize polygon onto raster
#                      
#                      
#                      #Reproject the polygon onto the empty_Behrman's raster
#                      raster_i_be <- projectRaster(raster_i, empty_be) #reproject raster to equalarea Behrmann's
#                      
#                      #Select the projected cells of the raster and take the values at these cells by comparison to the var_vals array
#                      rastercells <- which(getValues(raster_i_be) > 0) #select the cell
#                      values <- getValues(envar_be)[rastercells, ]
#                      
#                      #some only have one rastercell and therefore only one value - take this value
#                      #cannot divide by 10 if value is NA therefore ifelse 
#                      if(nrow(values) == 1) {
#                        ifelse(is.na(values[1]), 
#                               bio1 <- NA, 
#                               bio1 <-values[1]/10)
#                        ifelse(is.na(values[2]), 
#                               bio2 <- NA, 
#                               bio2 <-values[2]/10)
#                        ifelse(is.na(values[3]), 
#                               bio4 <- NA, 
#                               bio4 <-values[3])
#                        ifelse(is.na(values[4]), 
#                               bio5 <- NA, 
#                               bio5 <-values[4]/10)
#                        ifelse(is.na(values[5]), 
#                               bio6 <- NA, 
#                               bio6 <-values[5]/10)
#                        ifelse(is.na(values[6]), 
#                               bio7 <- NA, 
#                               bio7 <-values[6]/10)
#                        ifelse(is.na(values[7]), 
#                               bio12 <- NA, 
#                               bio12 <-values[7])
#                        ifelse(is.na(values[8]), 
#                               bio13 <- NA, 
#                               bio13 <-values[8])
#                        ifelse(is.na(values[9]), 
#                               bio14 <- NA, 
#                               bio14 <-values[9])
#                        ifelse(is.na(values[10]), 
#                               bio15 <- NA, 
#                               bio15 <-values[10])
#                        ifelse(is.na(values[11]), 
#                               gsl <- NA, 
#                               gsl <-values[11])
#                        ifelse(is.na(values[12]), 
#                               ngd <- NA, 
#                               ngd <-values[12])
#                        ifelse(is.na(values[13]), 
#                               npp <- NA, 
#                               npp <-values[13])
#                        
#                        
#                      } else { 
#                        #most have multiple so take the mean of all the values
#                        #However cannot take the mean or divide by 10 if all vals are NA hense ifelse()
#                        
#                        ifelse(all(is.na(values[,1])), 
#                               bio1 <- NA,
#                               bio1 <-mean(values[,1], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,2])),
#                               bio2 <- NA,
#                               bio2 <-mean(values[,2], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,3])), 
#                               bio4 <- NA,
#                               bio4 <-mean(values[,3], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,4])),
#                               bio5 <- NA,
#                               bio5 <-mean(values[,4], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,5])), 
#                               bio6 <- NA,
#                               bio6 <-mean(values[,5], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,6])), 
#                               bio7 <- NA,
#                               bio7 <-mean(values[,6], na.rm = TRUE)/10) #write them in
#                        ifelse(all(is.na(values[,7])), 
#                               bio12 <- NA,
#                               bio12 <-mean(values[,7], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,8])), 
#                               bio13 <- NA,
#                               bio13 <-mean(values[,8], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,9])), 
#                               bio14 <- NA,
#                               bio14 <-mean(values[,9], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,10])), 
#                               bio15 <- NA,
#                               bio15 <-mean(values[,10], na.rm = TRUE)) #write them in
#                        ifelse(all(is.na(values[,11])), 
#                               gsl <- NA,
#                               gsl <-mean(values[,11], na.rm = TRUE))
#                        ifelse(all(is.na(values[,12])), 
#                               ngd <- NA,
#                               ngd <-mean(values[,12], na.rm = TRUE))
#                        ifelse(all(is.na(values[,13])), 
#                               npp <- NA,
#                               npp <-mean(values[,13], na.rm = TRUE))
#                      }
#                      
#                      #get the study_results
#                      data.frame(species = split$value$bird_tree_name[1],
#                                 bio1,
#                                 bio2,
#                                 bio4,
#                                 bio5,
#                                 bio6,
#                                 bio7,
#                                 bio12,
#                                 bio13,
#                                 bio14,
#                                 bio15,
#                                 gsl,
#                                 ngd,
#                                 npp)
#                    }
# 
# 
# 
# ##There will be some that have maps but such small ranges that there is no value
# ##Create a "No Maps" and a "Small Range" csv
# small_range <- unique(sex_ranges$bird_tree_name)[-which(results$species %in%  unique(sex_ranges$bird_tree_name))]
# 
# write.csv(results, "../Data/jetz_environmental.csv")
# write.csv(small_range, "../Data/resolution_misses_range.csv")
# 
# 
