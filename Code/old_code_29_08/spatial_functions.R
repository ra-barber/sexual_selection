###############################################################################
                       ##### Spatial functions #####
###############################################################################


# Extract latutide from cells.
extract_lat <- function(x){
  nums <- as.numeric(x)
  abs(nums[2])
}

# Function to convert raster in sf_data with neighbourhood weights.
create_spatial_data <- function(spatial_raster){
  
  # Create spatial feature dataframe.
  data_sf <- st_as_sf(as(spatial_raster, 'SpatialPixelsDataFrame'))
  
  # Use queen nearest neighbours.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)
  data_sf$card_queen <- card(queen)
  
  # Remove the cells with zero neighbours
  data_sf <- subset(data_sf, card_queen > 0)
  
  # Recalculate neighbours.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)
  data_sf$card_queen <- card(queen)
  
  # Convert to weights.
  queen <- nb2listw(queen, style='W')
  
  # Add latitude.
  list_test <- lapply(data_sf$geometry, extract_lat)
  data_sf$abs_lat <- as.numeric(list_test)/100000
  return(data_sf)
}

# Function to run spatial model.
run_spatial_model <- function(data_sf){
  
  # Calculate neighbourhood weights.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)
  queen <- nb2listw(queen, style='W')
  
  # Run model.
  errorsarlm(layer ~ abs_lat,
             data=data_sf, listw=queen)
}

# Function to run spatial model.
scaled_run_spatial_model <- function(data_sf){
  
  # Calculate neighbourhood weights.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)
  queen <- nb2listw(queen, style='W')
  
  # Run model.
  errorsarlm(scale(layer) ~ scale(abs_lat),
             data=data_sf, listw=queen)
}

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
sar_lat_side_plot <- function(data_set, ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1), 
                              lab_x_pos = 60, lab_ypos = 1.2, plot_label = "b", 
                              plot_model = sex_model, stats_model = scaled_sex_model,
                              sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Extract model predictions with new data.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  new_data$preds <- predict(plot_model, newdata = new_data)[,1]
  
  # Estimate
  estimate <- last(coef(summary(stats_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(stats_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(stats_model))[,4])
  
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
    # Special case for if raster is skewed by too many zeroes.
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
    geom_point(aes(colour = binned_layer), size = 2, alpha = 0.9) + 
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 70), clip = 'off') +
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
    annotate("text", x = 1, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 2)
}



###############################################################################
                             #### END ####
###############################################################################