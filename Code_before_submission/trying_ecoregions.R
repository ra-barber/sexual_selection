###############################################################################
                        ##### Script Title #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(raster)
library(sf)
library(exactextractr)


# Read in the functions. 
source("Code/functions.R")


###############################################################################
                       #### Read in the data ####


# Read in ecoregions.
wwf_ecoregions <- read_sf("Data/ecoregions/official/wwf_terr_ecos.shp")

# Behrmann projection.
behrmannCRS <- CRS('+proj=cea +lat_ts=30')

# Aggregate polygons into a single polygon per ecoregion using st_union.
ecor.poly <- wwf_ecoregions[,c("ECO_ID")]
ecor.poly <- aggregate(ecor.poly, by=list(ecoregion=ecor.poly$ECO_ID), FUN=mean, do_union=T); ecor.poly$ECO_ID <- NULL

# Add centroid latitude of multipolygon. (For some reason didn't work with lapply)
library(PBSmapping)
ecor.poly$centroid_lat <- NA
for (x in 1:827){
  polygon_i <- ecor.poly$geometry[x]
  ecor.poly$centroid_lat[x] <-as.numeric(calcCentroid(raptr::convert2PolySet(sf::as_Spatial(polygon_i)), rollup = 1)[3])
}
ecor.poly$abs_lat <- abs(ecor.poly$centroid_lat)

# Use the more simple centroid that is fine for equal area.
library(sf)
ecor.poly$st_lat <- NA
ecor.poly$st_lat <- st_coordinates(st_centroid(ecor.poly))[,"Y"] 
ecor.poly$abs_lat_2 <- abs(ecor.poly$st_lat)
ecor.poly$abs_lat <- ecor.poly$abs_lat_2

# Project into behr.
ecor.poly_behr <- st_transform(ecor.poly, crs = behr_proj)
ecor.poly <- ecor.poly_behr

###############################################################################
            #### Extract sexual selection from rasters ####


# Figure 2.
ecor.poly$sex_cellmean <- exact_extract(sex_raster, ecor.poly, 'mean')
ecor.poly$cert_cellmean <- exact_extract(cert_raster, ecor.poly, 'mean')

# Figure 3.
ecor.poly$p_sex_cellmean <- exact_extract(p_sex_raster, ecor.poly, 'mean')
ecor.poly$s_sex_cellmean <- exact_extract(s_sex_raster, ecor.poly, 'mean')
ecor.poly$fruit_sex_cellmean <- exact_extract(fruit_sex_raster, ecor.poly, 'mean')
ecor.poly$invert_sex_cellmean <- exact_extract(invert_sex_raster, ecor.poly, 'mean')

# Extended data figure 3.
ecor.poly$h_sex_cellmean <- exact_extract(h_sex_raster, ecor.poly, 'mean')
ecor.poly$m_sex_cellmean <- exact_extract(m_sex_raster, ecor.poly, 'mean')

# Extended data figure 4.
ecor.poly$mig_sex_cellmean <- exact_extract(mig_sex_raster, ecor.poly, 'mean')
ecor.poly$no_mig_sex_cellmean <- exact_extract(no_mig_sex_raster, ecor.poly, 'mean')

# Extended data figure 5.
ecor.poly$terr_sex_cellmean <- exact_extract(terr_sex_raster, ecor.poly, 'mean')
ecor.poly$no_terr_sex_cellmean <- exact_extract(no_terr_sex_raster, ecor.poly, 'mean')

# Extended data figure 8.
ecor.poly$p_terr_cellmean <- exact_extract(p_terr_raster, ecor.poly, 'mean')
ecor.poly$s_terr_cellmean <- exact_extract(s_terr_raster, ecor.poly, 'mean')
ecor.poly$p_year_terr_cellmean <- exact_extract(p_year_terr_raster, ecor.poly, 'mean')
ecor.poly$s_year_terr_cellmean <- exact_extract(s_year_terr_raster, ecor.poly, 'mean')


###############################################################################
               ###### Calculate ecoregion filll ######


# Calculate size of each ecoregion.
ecor.poly$region_size <- exact_extract(land, ecor.poly, 'count')

# Figure 2.
ecor.poly$sex_cellcount <- exact_extract(sex_raster, ecor.poly, 'count')
ecor.poly$cert_cellcount <- exact_extract(cert_raster, ecor.poly, 'count')

# Figure 3.
ecor.poly$p_sex_cellcount <- exact_extract(p_sex_raster, ecor.poly, 'count')
ecor.poly$s_sex_cellcount <- exact_extract(s_sex_raster, ecor.poly, 'count')
ecor.poly$fruit_sex_cellcount <- exact_extract(fruit_sex_raster, ecor.poly, 'count')
ecor.poly$invert_sex_cellcount <- exact_extract(invert_sex_raster, ecor.poly, 'count')

# Extended data figure 3.
ecor.poly$h_sex_cellcount <- exact_extract(h_sex_raster, ecor.poly, 'count')
ecor.poly$m_sex_cellcount <- exact_extract(m_sex_raster, ecor.poly, 'count')

# Extended data figure 4.
ecor.poly$mig_sex_cellcount <- exact_extract(mig_sex_raster, ecor.poly, 'count')
ecor.poly$no_mig_sex_cellcount <- exact_extract(no_mig_sex_raster, ecor.poly, 'count')

# Extended data figure 5.
ecor.poly$terr_sex_cellcount <- exact_extract(terr_sex_raster, ecor.poly, 'count')
ecor.poly$no_terr_sex_cellcount <- exact_extract(no_terr_sex_raster, ecor.poly, 'count')

# Extended data figure 8.
ecor.poly$p_terr_cellcount <- exact_extract(p_terr_raster, ecor.poly, 'count')
ecor.poly$s_terr_cellcount <- exact_extract(s_terr_raster, ecor.poly, 'count')
ecor.poly$p_year_terr_cellcount <- exact_extract(p_year_terr_raster, ecor.poly, 'count')
ecor.poly$s_year_terr_cellcount <- exact_extract(s_year_terr_raster, ecor.poly, 'count')





###############################################################################
                       #### Loot at plots #####


# Figure 2.
plot(sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(cert_cellmean ~ abs(centroid_lat), data = ecor.poly)

# Figure 3.
plot(p_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(s_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(fruit_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(invert_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)


# Extended data figure 3.
plot(h_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(m_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)

# Extended data figure 4.
plot(mig_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(no_mig_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)

# Extended data figure 5.
plot(terr_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(no_terr_sex_cellmean ~ abs(centroid_lat), data = ecor.poly)

# Extended data figure 8.
plot(p_terr_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(s_terr_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(p_year_terr_cellmean ~ abs(centroid_lat), data = ecor.poly)
plot(s_year_terr_cellmean ~ abs(centroid_lat), data = ecor.poly)


# Function to run spatial model.
run_spatial_model <- function(data_sf){
  
  # Calculate neighbourhood weights.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize) # + 1)
  queen <- nb2listw(queen, style='W')
  
  # Run model.
  errorsarlm(scale(layer) ~ scale(abs_lat), 
             data=data_sf, listw=queen)
  errorsarlm(layer ~ abs_lat, 
             data=data_sf, listw=queen)
}



###############################################################################
           ##### Run the ecoregion spatial models ######



# Set up look parameters.
fit.vars <- c("sex_cellmean", "cert_cellmean", "p_sex_cellmean", "s_sex_cellmean", 
              "fruit_sex_cellmean", "invert_sex_cellmean", "m_sex_cellmean", 
              "h_sex_cellmean", "mig_sex_cellmean", "no_mig_sex_cellmean", "terr_sex_cellmean",
              "no_terr_sex_cellmean", "p_terr_cellmean", "s_terr_cellmean", "p_year_terr_cellmean",
              "s_year_terr_cellmean")

region.size <- c("sex_cellcount", "cert_cellcount", "p_sex_cellcount", "s_sex_cellcount", 
              "fruit_sex_cellcount", "invert_sex_cellcount", "m_sex_cellcount", 
              "h_sex_cellcount", "mig_sex_cellcount", "no_mig_sex_cellcount", "terr_sex_cellcount",
              "no_terr_sex_cellcount", "p_terr_cellcount", "s_terr_cellcount", "p_year_terr_cellcount",
              "s_year_terr_cellcount")


out <- as.list(rep(NA, length(fit.vars)))
names(out) <- fit.vars


library(foreach)
library(doParallel)

# Register cores.
registerDoParallel(cores = 8)


# Use the iterator package to split the sf dataframe into chunks by species name.
results <- foreach(i = 1:length(fit.vars),
                     .packages=c("sf", "sp", "spatialreg", "spdep", "stats", "dplyr", "tidyr"),
                     .inorder = TRUE,
                     .errorhandling = "stop") %dopar%{
                       
                       print (paste(i, "of", length(fit.vars)))
                       var <- fit.vars[i]
                       size <- region.size[i]
                       
                       # Drop NA ecoregions from model.
                       ecor_poly_no_na <- ecor.poly %>% drop_na(!!!var)
                       
                       # Remove ecoregions with less than 75% full.
                       ecor_poly_no_na$region_fill <-  as.data.frame(ecor_poly_no_na[,size] / ecor_poly_no_na$region_size)[,1]
                       ecor_poly_no_na <- ecor_poly_no_na %>% filter(region_fill > 0.75)
                       
                       colnames(ecor_poly_no_na)[colnames(ecor_poly_no_na)==var] <- "layer"
                       
                       # Create neighbourhood matrix.
                       nb <- poly2nb(ecor_poly_no_na, row.names = ecor_poly_no_na$ecoregion)
                       
                       # fit ols
                       fit.ols <- lm(scale(layer) ~ scale(abs_lat), data = ecor_poly_no_na)
                       fit.ols_raw <- lm(layer ~ abs_lat, data = ecor_poly_no_na)
                       
                       # build SAR models with different combinations of distance and neighbourhood style
                       neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
                       AICvec <- numeric(length = length(neiStyle))
                       for (k in 1:length(neiStyle)) {
                         nlw <- nb2listw(nb, style=neiStyle[k], zero.policy=TRUE)
                         sar_e <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
                         AICvec[k] <- AIC(sar_e)			
                       }
                       
                       bestStyle <- neiStyle[which.min(AICvec)]
                       nlw <- nb2listw(nb, style=bestStyle, zero.policy=TRUE)
                       fit.sar <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
                       fit.sar_raw <- errorsarlm(fit.ols_raw, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
                       
                       nlw_W <- nb2listw(nb, style="W", zero.policy=TRUE)
                       fit.sar_W <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
                       fit.sar_raw_W <- errorsarlm(fit.ols_raw, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
                       
                       
                       mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
                       mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
                       aic.ols <- AIC(fit.ols)
                       aic.sar <- AIC(fit.sar)
                       res <- list(fit.ols = fit.ols, bestStyle = bestStyle, nb = nb,
                                   fit.sar = fit.sar, fit.sar_raw = fit.sar_raw, 
                                   aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
                       out[[i]] <- res
                       
                       res
                     }


results

names(results) <- fit.vars

saveRDS(results, "Results/ecoregion_spatiallm.rds")


x <- numeric(length(fit.vars))
results_df <- data.frame(metric = fit.vars, AIC.ols=x, AIC.sar=x, dAIC.ols=x, dAIC.sar=x,
                         ols.slope=x, ols.pvalue=x,
                         sar.slope=x, sar.pvalue=x, slope.ratio=x,
                         moran.ols = x, moran.ols.pvalue = x, moran.sar = x, moran.sar.pvalue = x, stringsAsFactors=F)


for (i in 1:length(out)) {
  
  fres <- results[[i]]
  
  results_df$AIC.ols[i] <- fres$aic.ols
  results_df$bestStyle[i] <- fres$bestStyle
  results_df$AIC.sar[i] <- fres$aic.sar
  results_df$dAIC.ols[i] <- fres$aic.ols - (min(fres$aic.ols, fres$aic.sar))
  results_df$dAIC.sar[i] <- fres$aic.sar - (min(fres$aic.ols, fres$aic.sar))
  results_df$ols.slope[i] <- coef(fres$fit.ols)[2]
  results_df$ols.pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
  results_df$sar.slope[i] <- coef(fres$fit.sar)[2]
  results_df$sar.pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
  results_df$slope.ratio[i] <- results_df$sar.slope[i] / results_df$ols.slope[i]
  results_df$moran.ols[i] <- fres$moran.ols$estimate[1]
  results_df$moran.ols.pvalue[i] <- fres$moran.ols$p.value
  results_df$moran.sar[i] <- fres$moran.sar$estimate[1]
  results_df$moran.sar.pvalue[i] <- fres$moran.sar$p.value
  
}


# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
splm_lat_side_plot <- function(variable, ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 70, lab_ypos = 1.2, plot_label = "b", 
                               plot_model = sex_model, #plot_raster = sex_raster,
                               sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Pull out models for predicting.
  raw_model <- results[[variable]]$fit.sar_raw
  scaled_model <- results[[variable]]$fit.sar
  
  raw_model$y
  
  # Create dataset for plotting.
  data_set <- data.frame( abs_lat = raw_model$X[,2], layer = raw_model$y)
  
  # # Extract model predictions with new data.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  new_data$preds <- predict(raw_model, newdata = new_data)[,1]
  
  # Estimate
  estimate <- last(coef(summary(scaled_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(plot_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(scaled_model))[,4])
  
  # Change p value to a string, using standard thresholds. 
  if (p_value < 0.001 ){
    p_value <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_value <- "p < 0.01"
  } else if (p_value < 0.05) {
    p_value <- "p < 0.05"
  } else {
    p_value <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)))
  }
  
  # Create a label 
  stats_label <- paste0(estimate, "\n", p_value)
  
  if (p_year_terr){
    # Special case for when primary territory raster is skewed by too many zeroes.
    breaks <- quantile(data_set$layer[data_set$layer>0], seq(0,1,length.out=nbins), na.rm=TRUE)
    breaks <- c(0,breaks) 
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  } else {
    # Create binned colour scale.
    breaks <- quantile(data_set$layer, seq(0,1,length.out=nbins+1), na.rm=TRUE)
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  }
  data_set$binned_layer <- as.factor(as.numeric(cuts))
  
  # Create a label of the bins for plotting. (Only used in map legends)
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  ggplot(data_set, aes(x = abs_lat, y = layer)) +
    geom_point(aes(colour = binned_layer), size = 3, alpha = 0.9) + 
    scale_x_continuous(breaks = seq(from = 0, to = 80, by = 40)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 80), clip = 'off') +
    scale_size_continuous(range = c(1,6))+
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    
    scale_colour_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) +
    annotate("text", x = 40, y = ylimits[2], label = variable, size = 12, fontface = 2) + 
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 1)
}



# Figure 2.
sex_lat_plot <- splm_lat_side_plot(fit.vars[1], ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                   ybreaks =  c(0,0.5,1.0), lab_ypos = 0.12, plot_label = "b")



cert_lat_plot <- splm_lat_side_plot(fit.vars[2], ylabel = "Data certainty", ylimits = c(1,4), 
                                    ybreaks =  c(1,2,3,4), lab_ypos = 1.3, plot_label = "d")


# Primary and secondary consumers.
pri_lat_plot <- splm_lat_side_plot(fit.vars[3], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                   ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                                   plot_label = "b") 

sec_lat_plot <- splm_lat_side_plot(fit.vars[4], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                   ybreaks =  c(0,1.0,2.0), lab_ypos = 1.5,
                                   plot_label = "f") 

# Trophic niche models.
fruit_lat_plot <- splm_lat_side_plot(fit.vars[5], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                                     plot_label = "d") 


invert_lat_plot <- splm_lat_side_plot(fit.vars[6], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                      ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                                      plot_label = "h") 


# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- splm_lat_side_plot(fit.vars[7], ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                       ybreaks =  c(0,0.5,1.0), lab_ypos = 1.11, plot_label = "b",
                                       sex_score = TRUE)

hi_sex_lat_plot <- splm_lat_side_plot(fit.vars[8], ylabel = "Sexual selection", ylimits = c(0,1.8), 
                                      ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.67, plot_label = "d",
                                      sex_score = TRUE)

# Make the side plots.
mig_lat_plot <- splm_lat_side_plot(fit.vars[9], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                   ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39,  plot_label = "b",
                                   sex_score = TRUE)

no_mig_lat_plot <- splm_lat_side_plot(fit.vars[10], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                      ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                                      sex_score = TRUE)

terr_lat_plot <-  splm_lat_side_plot(fit.vars[11], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "b",
                                     sex_score = TRUE)

no_terr_lat_plot <- splm_lat_side_plot(fit.vars[12], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                       ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                                       sex_score = TRUE)

# Primary terr.
pri_terr_lat_plot <- splm_lat_side_plot(fit.vars[13], ylabel = "Proportion of species", ylimits = c(0,1), 
                                        ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "b",
                                        sex_score = FALSE)
# Primary year terr.
pri_yearterr_lat_plot <-  splm_lat_side_plot(fit.vars[14], ylabel = "Proportion of species", ylimits = c(0,1), 
                                             ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, 
                                             plot_label = "d", 
                                             sex_score = FALSE)

# Secondary terr.
sec_terr_lat_plot <- splm_lat_side_plot(fit.vars[15], ylabel = "Proportion of species", ylimits = c(0,1), 
                                        ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "f", 
                                        p_year_terr = TRUE,
                                        sex_score = FALSE)

# Secondary year terr.
sec_yearterr_lat_plot <-  splm_lat_side_plot(fit.vars[16], ylabel = "Proportion of species", ylimits = c(0,1), 
                                             ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "h",
                                             sex_score = FALSE)




################################################################################
                           ##### End ######
################################################################################













out <- as.list(rep(NA, length(fit.vars)))
names(out) <- fit.vars


for (i in 1:length(fit.vars)) {
  print (paste(i, "of", length(fit.vars)))
  var <- fit.vars[i]
  size <- region.size[i]
  
  # Drop NA ecoregions from model.
  ecor_poly_no_na <- ecor.poly %>% drop_na(!!!var)
  
  # Remove ecoregions with less than 75% full.
  ecor_poly_no_na$region_fill <-  as.data.frame(ecor_poly_no_na[,size] / ecor_poly_no_na$region_size)[,1]
  ecor_poly_no_na <- ecor_poly_no_na %>% filter(region_fill > 0.75)
  
  # Create neighbourhood matrix.
  nb <- poly2nb(ecor_poly_no_na, row.names = ecor_poly_no_na$ecoregion)
  
  # fit ols
  fit.ols <- lm(scale(ecor_poly_no_na[,var,drop=T]) ~ scale(ecor_poly_no_na[,"abs_lat",drop=T]))
  fit.ols_raw <- lm(ecor_poly_no_na[,var,drop=T] ~ ecor_poly_no_na[,"abs_lat",drop=T])
  
  # build SAR models with different combinations of distance and neighbourhood style
  neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
  AICvec <- numeric(length = length(neiStyle))
  for (k in 1:length(neiStyle)) {
    nlw <- nb2listw(nb, style=neiStyle[k], zero.policy=TRUE)
    sar_e <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
    AICvec[k] <- AIC(sar_e)			
  }
  
  # Run model again with best style.
  bestStyle <- neiStyle[which.min(AICvec)]
  nlw <- nb2listw(nb, style=bestStyle, zero.policy=TRUE)
  fit.sar <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  fit.sar_raw <- spautolm(fit.ols_raw, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  
  # Run local 
  mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  aic.ols <- AIC(fit.ols)
  aic.sar <- AIC(fit.sar)
  res <- list(fit.ols = fit.ols, fit.sar = fit.sar, fit.sar_raw = fit.sar_raw, aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
  out[[i]] <- res
}

saveRDS(out, "Results/ecoregion_spatiallm.rds")


# Summarise
x <- numeric(length(fit.vars))
dff <- data.frame(metric = fit.vars, AIC.ols=x, AIC.sar=x, dAIC.ols=x, dAIC.sar=x,
                  ols.slope=x, ols.pvalue=x,
                  sar.slope=x, sar.pvalue=x, slope.ratio=x,
                  moran.ols = x, moran.ols.pvalue = x, moran.sar = x, moran.sar.pvalue = x, stringsAsFactors=F)

for (i in 1:length(out)) {
  
  fres <- out[[i]]
  
  dff$AIC.ols[i] <- fres$aic.ols
  dff$AIC.sar[i] <- fres$aic.sar
  dff$dAIC.ols[i] <- fres$aic.ols - (min(fres$aic.ols, fres$aic.sar))
  dff$dAIC.sar[i] <- fres$aic.sar - (min(fres$aic.ols, fres$aic.sar))
  dff$ols.slope[i] <- coef(fres$fit.ols)[2]
  dff$ols.pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
  dff$sar.slope[i] <- coef(fres$fit.sar)[2]
  dff$sar.pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
  dff$slope.ratio[i] <- dff$sar.slope[i] / dff$ols.slope[i]
  dff$moran.ols[i] <- fres$moran.ols$estimate[1]
  dff$moran.ols.pvalue[i] <- fres$moran.ols$p.value
  dff$moran.sar[i] <- fres$moran.sar$estimate[1]
  dff$moran.sar.pvalue[i] <- fres$moran.sar$p.value
  
}



nb
error_model <- errorsarlm(s_year_terr_cellmean ~ abs_lat, data = ecor_poly_no_na, listw = nlw, zero.policy = TRUE)
predict(error_model)
s_year_terr_cellmean

summary(error_model)

dff
out[["sex_cellmean"]]

raw_model <- out[["s_year_terr_cellmean"]]$fit.sar_raw
summary(raw_model)
summary(error_model)
data_set <- data.frame(abs_lat = raw_model$X, layer = raw_model$Y)

p_year_terr <- TRUE
variable <- fit.vars[16]

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
splm_lat_side_plot <- function(variable, ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1), 
                              lab_x_pos = 70, lab_ypos = 1.2, plot_label = "b", 
                              plot_model = sex_model, #plot_raster = sex_raster,
                              sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Pull out models for predicting.
  raw_model <- out[[variable]]$fit.sar_raw
  scaled_model <- out[[variable]]$fit.sar
  
  # Create dataset for plotting.
  data_set <- data.frame(abs_lat = raw_model$X[,2], layer = raw_model$Y, preds = raw_model$fit$fitted.values)


  
  # # Extract model predictions with new data.
  # new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  # new_data$preds <- predict(raw_model, newdata = new_data)[,1]

  
  # Estimate
  estimate <- last(coef(summary(scaled_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(plot_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(scaled_model))[,4])
  
  # Change p value to a string, using standard thresholds. 
  if (p_value < 0.001 ){
    p_value <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_value <- "p < 0.01"
  } else if (p_value < 0.05) {
    p_value <- "p < 0.05"
  } else {
    p_value <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)))
  }
  
  # Create a label 
  stats_label <- paste0(estimate, "\n", p_value)
  
  if (p_year_terr){
    # Special case for when primary territory raster is skewed by too many zeroes.
    breaks <- quantile(data_set$layer[data_set$layer>0], seq(0,1,length.out=nbins), na.rm=TRUE)
    breaks <- c(0,breaks) 
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  } else {
    # Create binned colour scale.
    breaks <- quantile(data_set$layer, seq(0,1,length.out=nbins+1), na.rm=TRUE)
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  }
  data_set$binned_layer <- as.factor(as.numeric(cuts))
  
  # Create a label of the bins for plotting. (Only used in map legends)
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  ggplot(data_set, aes(x = abs_lat, y = layer)) +
    geom_point(aes(colour = binned_layer), size = 3, alpha = 0.9) + 
    scale_x_continuous(breaks = seq(from = 0, to = 80, by = 40)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 80), clip = 'off') +
    scale_size_continuous(range = c(1,6))+
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    
    scale_colour_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) +
    annotate("text", x = 40, y = ylimits[2], label = variable, size = 12, fontface = 2) #+ 
  
   # geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 1)
}


splm_lat_side_plot(fit.vars[1])
splm_lat_side_plot(fit.vars[2])
splm_lat_side_plot(fit.vars[3])
splm_lat_side_plot(fit.vars[4])
splm_lat_side_plot(fit.vars[5])
splm_lat_side_plot(fit.vars[6])
splm_lat_side_plot(fit.vars[7])
splm_lat_side_plot(fit.vars[8])
splm_lat_side_plot(fit.vars[9])
splm_lat_side_plot(fit.vars[10])
splm_lat_side_plot(fit.vars[11])
splm_lat_side_plot(fit.vars[12])
splm_lat_side_plot(fit.vars[13])


# Figure 2.
sex_lat_plot <- splm_lat_side_plot(fit.vars[1], ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                  ybreaks =  c(0,0.5,1.0), lab_ypos = 0.12, plot_label = "b")



cert_lat_plot <- splm_lat_side_plot(fit.vars[2], ylabel = "Data certainty", ylimits = c(1,4), 
                 ybreaks =  c(1,2,3,4), lab_ypos = 1.3, plot_label = "d")


# Primary and secondary consumers.
pri_lat_plot <- splm_lat_side_plot(fit.vars[3], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                  ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                                  plot_label = "b") 

sec_lat_plot <- splm_lat_side_plot(fit.vars[4], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                  ybreaks =  c(0,1.0,2.0), lab_ypos = 1.5,
                                  plot_label = "f") 

# Trophic niche models.
fruit_lat_plot <- splm_lat_side_plot(fit.vars[5], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                    ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                                    plot_label = "d") 


invert_lat_plot <- splm_lat_side_plot(fit.vars[6], ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                                     plot_label = "h") 


# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- splm_lat_side_plot(fit.vars[7], ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                      ybreaks =  c(0,0.5,1.0), lab_ypos = 1.11, plot_label = "b",
                                      sex_score = TRUE)

hi_sex_lat_plot <- splm_lat_side_plot(fit.vars[8], ylabel = "Sexual selection", ylimits = c(0,1.8), 
                                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.67, plot_label = "d",
                                     sex_score = TRUE)

# Make the side plots.
mig_lat_plot <- splm_lat_side_plot(fit.vars[9], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                  ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39,  plot_label = "b",
                                  sex_score = TRUE)

no_mig_lat_plot <- splm_lat_side_plot(fit.vars[10], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                                     sex_score = TRUE)

terr_lat_plot <-  splm_lat_side_plot(fit.vars[11], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                    ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "b",
                                    sex_score = TRUE)

no_terr_lat_plot <- splm_lat_side_plot(fit.vars[12], ylabel = "Sexual selection", ylimits = c(0,1.5), 
                                      ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                                      sex_score = TRUE)

# Primary terr.
pri_terr_lat_plot <- splm_lat_side_plot(fit.vars[13], ylabel = "Proportion of species", ylimits = c(0,1), 
                                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "b",
                                       sex_score = FALSE)
# Primary year terr.
pri_yearterr_lat_plot <-  splm_lat_side_plot(fit.vars[14], ylabel = "Proportion of species", ylimits = c(0,1), 
                                            ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, 
                                            plot_label = "d", 
                                            sex_score = FALSE)

# Secondary terr.
sec_terr_lat_plot <- splm_lat_side_plot(fit.vars[15], ylabel = "Proportion of species", ylimits = c(0,1), 
                                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "f", 
                                       p_year_terr = TRUE,
                                       sex_score = FALSE)

# Secondary year terr.
sec_yearterr_lat_plot <-  splm_lat_side_plot(fit.vars[16], ylabel = "Proportion of species", ylimits = c(0,1), 
                                            ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "h",
                                            sex_score = FALSE)



################################################################################
                 #######  For each loop to run models with predictions #####


out_2 <- as.list(rep(NA, length(fit.vars)))
names(out_2) <- fit.vars



library(foreach)
library(doParallel)

# Register cores.
registerDoParallel(cores = 8)


# Use the iterator package to split the sf dataframe into chunks by species name.
results_2 <- foreach(i = 1:length(fit.vars),
                   .packages=c("sf", "sp", "spatialreg", "spdep", "stats", "dplyr", "tidyr"),
                   .inorder = TRUE,
                   .errorhandling = "stop") %dopar%{

  print (paste(i, "of", length(fit.vars)))
  var <- fit.vars[i]
  size <- region.size[i]
  
  # Drop NA ecoregions from model.
  ecor_poly_no_na <- ecor.poly %>% drop_na(!!!var)
  
  #ecor_poly_no_na <- ecor_poly_no_na %>% filter(!!!size > 0.75)
  
  # Remove ecoregions with less than 75% full.
  ecor_poly_no_na$region_fill <-  as.data.frame(ecor_poly_no_na[,size] / ecor_poly_no_na$region_size)[,1]
  ecor_poly_no_na <- ecor_poly_no_na %>% filter(region_fill > 0.75)
  
  colnames(ecor_poly_no_na)[colnames(ecor_poly_no_na)==var] <- "layer"
  
  # Create neighbourhood matrix.
  nb <- poly2nb(ecor_poly_no_na, row.names = ecor_poly_no_na$ecoregion)
  
  # fit ols
  fit.ols <- lm(scale(layer) ~ scale(abs_lat), data = ecor_poly_no_na)
  fit.ols_raw <- lm(layer ~ abs_lat, data = ecor_poly_no_na)
  
  # build SAR models with different combinations of distance and neighbourhood style
  neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
  AICvec <- numeric(length = length(neiStyle))
  for (k in 1:length(neiStyle)) {
    nlw <- nb2listw(nb, style=neiStyle[k], zero.policy=TRUE)
    sar_e <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
    AICvec[k] <- AIC(sar_e)			
  }
  
  bestStyle <- neiStyle[which.min(AICvec)]
  nlw <- nb2listw(nb, style=bestStyle, zero.policy=TRUE)
  fit.sar <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  fit.sar_raw <- errorsarlm(fit.ols_raw, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  
  nlw_W <- nb2listw(nb, style="W", zero.policy=TRUE)
  fit.sar_W <- errorsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  fit.sar_raw_W <- errorsarlm(fit.ols_raw, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  
  
  mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  aic.ols <- AIC(fit.ols)
  aic.sar <- AIC(fit.sar)
  res <- list(fit.ols = fit.ols, bestStyle = bestStyle, nb = nb,
              fit.sar = fit.sar, fit.sar_raw = fit.sar_raw, 
              aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
  out_2[[i]] <- res
  
  res
}


results_2

names(results_2) <- fit.vars

x <- numeric(length(fit.vars))
results_df <- data.frame(metric = fit.vars, AIC.ols=x, AIC.sar=x, dAIC.ols=x, dAIC.sar=x,
                  ols.slope=x, ols.pvalue=x,
                  sar.slope=x, sar.pvalue=x, slope.ratio=x,
                  moran.ols = x, moran.ols.pvalue = x, moran.sar = x, moran.sar.pvalue = x, stringsAsFactors=F)

fres$bestStyle

for (i in 1:length(out)) {
  
  fres <- results[[i]]
  
  results_df$AIC.ols[i] <- fres$aic.ols
  results_df$bestStyle[i] <- fres$bestStyle
  results_df$AIC.sar[i] <- fres$aic.sar
  results_df$dAIC.ols[i] <- fres$aic.ols - (min(fres$aic.ols, fres$aic.sar))
  results_df$dAIC.sar[i] <- fres$aic.sar - (min(fres$aic.ols, fres$aic.sar))
  results_df$ols.slope[i] <- coef(fres$fit.ols)[2]
  results_df$ols.pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
  results_df$sar.slope[i] <- coef(fres$fit.sar)[2]
  results_df$sar.pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
  results_df$slope.ratio[i] <- results_df$sar.slope[i] / results_df$ols.slope[i]
  results_df$moran.ols[i] <- fres$moran.ols$estimate[1]
  results_df$moran.ols.pvalue[i] <- fres$moran.ols$p.value
  results_df$moran.sar[i] <- fres$moran.sar$estimate[1]
  results_df$moran.sar.pvalue[i] <- fres$moran.sar$p.value
  
}


results_df



variable <- fit.vars[1]

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
splm_lat_side_plot <- function(variable, ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 70, lab_ypos = 1.2, plot_label = "b", 
                               plot_model = sex_model, #plot_raster = sex_raster,
                               sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Pull out models for predicting.
  raw_model <- results_2[[variable]]$fit.sar_raw
  scaled_model <- results_2[[variable]]$fit.sar
  
  raw_model$y
  
  # Create dataset for plotting.
  data_set <- data.frame( abs_lat = raw_model$X[,2], layer = raw_model$y)
  

  # new_data$preds <- predict(plot_model, newdata = new_data)[,1]
  # 
  # 
  # predict(raw_model, newdata)
  
  # # Extract model predictions with new data.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  new_data$preds <- predict(raw_model, newdata = new_data)[,1]
  
  
  # Estimate
  estimate <- last(coef(summary(scaled_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(plot_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(scaled_model))[,4])
  
  # Change p value to a string, using standard thresholds. 
  if (p_value < 0.001 ){
    p_value <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_value <- "p < 0.01"
  } else if (p_value < 0.05) {
    p_value <- "p < 0.05"
  } else {
    p_value <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)))
  }
  
  # Create a label 
  stats_label <- paste0(estimate, "\n", p_value)
  
  if (p_year_terr){
    # Special case for when primary territory raster is skewed by too many zeroes.
    breaks <- quantile(data_set$layer[data_set$layer>0], seq(0,1,length.out=nbins), na.rm=TRUE)
    breaks <- c(0,breaks) 
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  } else {
    # Create binned colour scale.
    breaks <- quantile(data_set$layer, seq(0,1,length.out=nbins+1), na.rm=TRUE)
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  }
  data_set$binned_layer <- as.factor(as.numeric(cuts))
  
  # Create a label of the bins for plotting. (Only used in map legends)
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  ggplot(data_set, aes(x = abs_lat, y = layer)) +
    geom_point(aes(colour = binned_layer), size = 3, alpha = 0.9) + 
    scale_x_continuous(breaks = seq(from = 0, to = 80, by = 40)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 80), clip = 'off') +
    scale_size_continuous(range = c(1,6))+
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    
    scale_colour_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) +
    annotate("text", x = 40, y = ylimits[2], label = variable, size = 12, fontface = 2) + 
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 1)
}






options(scipen = 9)



out$fruit_sex_cellmean$fit.sar$X
options(scipen = 99)

dff




summary(sar_e)

st_area(sex_ranges$)



# Test to see if this works after aggregating.
test_extract <- exact_extract(sex_raster, ecor.poly, 'mean')

ecor.poly$sexual_selection <- test_extract


# Test exact extract first. (doesn't seem to work)
test_extract <- exact_extract(sex_raster, wwf_ecoregions, 'mean')

# This works when we don't use the behrmanns projection.

# Aggregate ecoregions.
ecor.poly <- wwf_ecoregions[,c("ECO_ID")]
ecor.poly <- aggregate(ecor.poly, by=list(ecoregion=ecor.poly$ECO_ID), FUN=mean, do_union=T); ecor.poly$ECO_ID <- NULL

library(PBSmapping)

# Add centroid latitude.
ecor.poly$centroid_lat <- NA
for (x in 1:827){
  polygon_i <- ecor.poly$geometry[x]
  ecor.poly$centroid_lat[x] <-as.numeric(calcCentroid(raptr::convert2PolySet(sf::as_Spatial(polygon_i)), rollup = 1)[3])
}

plot(ecor.poly$sexual_selection ~ abs(ecor.poly$centroid_lat))











test_lat <- lapply(ecor.poly$geometry, extract_centroid)
polygon_i <- ecor.poly$geometry[3]
centroid_lat_PBS <- as.numeric(calcCentroid(raptr::convert2PolySet(sf::as_Spatial(ecor.poly$geometry[1])), rollup = 1)[3])


# Test to see if this works after aggregating.
test_extract <- exact_extract(sex_raster, ecor.poly, 'mean')

ecor.poly$sexual_selection <- test_extract

# Extract latutide from cells.
extract_lat <- function(x){
  nums <- as.numeric(x)
  abs(nums[2])
}
# Add latitude.
list_test <- lapply(ecor.poly$geometry, extract_lat)
ecor.poly$lat <- as.numeric(list_test)
ecor.poly$abs_lat <- abs(ecor.poly$lat)


# Test s
ecor.poly = st_transform(ecor.poly, crs = behrmannCRS)
test_extract <- exact_extract(sex_raster, ecor.poly, 'mean')

# Try adding the shapefiles on top.
plot(sex_raster)
plot(ecor.poly$geometry, add= TRUE)




wwf_ecoregions$geometry[1:2]
?aggregate.sf
crs(sex_raster)
crs(ecor.poly$geometry)
ecor.poly = st_transform(wwf_ecoregions, crs = behr)

plot(sex_raster)
plot(wwf_ecoregions$geometry[1:10000], add = TRUE)

bbox(sex_raster) 

st_cast(wwf_ecoregions$geometry[[1]], 'MULTIPOLYGON')
bbox( SpatialPolygons(wwf_ecoregions$geometry[[1]]))


# Calculate vector of mean December precipitation amount for each municipality
brazil$mean_dec_prec <- exact_extract(prec[[12]], brazil, 'mean')


?raster::extract
str(wwf_ecoregions)

plot(wwf_ecoregions$geometry[2])


colnames(wwf_ecoregions)

CRS(wwf_ecoregions$geometry[1]) <- behr

n_distinct(wwf_ecoregions$ECO_NAME)
janitor::get_dupes(wwf_ecoregions, ECO_NAME) %>% pull(ECO_NAME)

wwf_ecoregions %>% filter(ECO_NAME == "Atlantic mixed forests") %>% pull(geometry) %>% plot()

ecor.poly <- wwf_ecoregions
newproj <- behr

ecor.poly <- ecor.poly[,c("ECO_ID")]
ecor.poly <- aggregate(ecor.poly, by=list(ecoregion=ecor.poly$ECO_ID), FUN=mean, do_union=T); ecor.poly$ECO_ID <- NULL
ecor.poly = st_transform(ecor.poly, crs = newproj)



pam$ecoregion <- ecor.cells$ecoregion_combined
res.cells$ecoregion <- ecor.cells$ecoregion_combined

crs(superhigh_res_template)
crs(ecor.poly)

raster::extract(sex_raster, polygon)

plot(polygon)
polygon <- st_combine(ecor.poly$geometry[2])
raster_i <- fasterize(st_as_sf(polygon), high_res_template) #fasterize polygon onto raster
#Rasterize the polygons
polygon <- st_combine(split$value$Shape) #combine to one shape layer (union not required as fasterize sorts out the overlaps)
raster_i <- fasterize(st_as_sf(ecor.poly[1,]), superhigh_res_template) #fasterize polygon onto raster

#Reproject the polygon onto the empty_Behrman's raster
raster_i_be <- projectRaster(raster_i, empty_be) #reproject raster to equalarea Behrmann's

#Select the projected cells of the raster and take the values at these cells by comparison to the var_vals array
rastercells <- which(getValues(raster_i_be) > 0) #select the cell
values <- var_vals[rastercells, ]


plot(raster_i)


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


