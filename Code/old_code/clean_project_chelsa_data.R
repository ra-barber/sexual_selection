

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


behr <-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
registerDoParallel(cores = 7)
library(tictoc)
tic()
envar_be <- foreach(x = 1:13, .packages = "terra") %dopar% {
  loop_raster <- rast(list.files("../Data/Chelsa/", full.names = TRUE)[x])
  project(loop_raster, behr)
}
toc()

envar_be <- rast(envar_be)

writeRaster(envar_be, "../Data/Chelsa/bergmans_chelsa.tif")