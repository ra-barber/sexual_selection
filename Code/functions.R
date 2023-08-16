##############################################################################
                        #### Functions #####
###############################################################################

# This script has functions used throughout the sexual selection project.



###############################################################################
                     ##### Data handling ######

# Filter the dataset for high certainty behavioural traits.
filter_high <- function(dataset){
  dataset %<>% filter(incubation_certainty < 3 &
                        sexual_certainty < 3 &
                        nest_certainty < 3 &
                        lifehistory_uncertainty < 3)
  return(dataset)
}

# Quick function for centering factors.
center_categorical <- function(predictor){
  as.numeric(predictor) %>% scale(scale = FALSE)
}


###############################################################################
                  ##### Plotting functions ######

# Function for looking at sexual selection vs predictor boxplots.
sex_boxplot <- function(predictor){
  ggplot(phylo_data, aes(x = predictor, y = sexual_score)) +
    geom_boxplot() + ylab("Sexual Selection") + 
    stat_summary(fun=mean, geom="point", shape=20, size=8) +
    coord_cartesian(ylim=c(0,4)) + theme_classic()
}

# Function for looking at sexual selection vs predictor means and standard error.
sex_meanplot <- function(predictor){
  grouped_data <- phylo_data %>% 
    group_by(!!! syms(predictor)) %>% 
    summarise(trait = first(!!! syms(predictor)),
              sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  ggplot(grouped_data, aes(x = trait, y = sex_mean)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se*1.96, ymax = sex_mean + sex_se*1.96)) + 
    geom_point() + ylab("Sexual Selection") + xlab(predictor) + theme_classic()
}

# Function for looking at sexual selection vs predictor means and standard error.
trait_meanplot <- function(predictor){
  grouped_data <- phylo_data %>% 
    group_by(sexual_binary) %>% 
    summarise(trait = first(!!! syms(predictor)),
              trait_mean = mean(!!! syms(predictor)),
              trait_sd = sd(!!! syms(predictor)),
              trait_se = sd(!!! syms(predictor))/sqrt(length(!!! syms(predictor))))
  
  ggplot(grouped_data, aes(x = as.factor(sexual_binary), y = trait_mean)) +
    scale_x_discrete(labels = c("Monogamous", "Polygamous")) +
    geom_errorbar(aes(ymin = trait_mean - trait_se*1.96, ymax = trait_mean + trait_se*1.96)) + 
    geom_point() + ylab(predictor) + xlab("Sexual selection") + theme_classic()
}

# Custom function for hex scatter plot in ggpairs.
ghex_func <- function(data, mapping, ...) { 
  
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ggplot(data) +
    geom_hex(mapping, bins = 50) + scale_fill_viridis_c() + theme(panel.grid.major = element_line(size = 1))
  
}


################################################################################
                  ##### Making maps ######

# Function to make average raster from variable name.
average_raster <- function(var_name, ranges = sex_ranges, spec_raster = species_raster){
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  return(var_raster)
}

# Function to make average raster using different species richness maps.
average_raster_2 <- function(ranges = sex_ranges, var_name = "sexual_score"){
  
  # Create the species richness raster.
  spec_raster <- fasterize(ranges, raster_template, fun = "sum")
  spec_raster[which(getValues(spec_raster) < 10)] <- NA
  spec_raster <- mask(spec_raster, land)
  
  # Create the average raster.
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  return(var_raster)
}

# Function to make species richness raster.
spec_raster_func <- function(ranges = sex_ranges){
  
  # Create the species richness raster.
  spec_raster <- fasterize(ranges, raster_template, fun = "sum")
  spec_raster[which(getValues(spec_raster) < 10)] <- NA
  spec_raster <- mask(spec_raster, land)

  return(spec_raster)
}


# Function to make average raster but without filtering out less than 10 species.
average_raster_3 <- function(ranges = sex_ranges, var_name = "sexual_score"){
  
  # Create the species richness raster.
  spec_raster <- fasterize(ranges, raster_template, fun = "sum")
  spec_raster <- mask(spec_raster, land)
  
  # Create the average raster.
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  return(var_raster)
}

# Function to make average raster but  filtering out less than 5 species.
average_raster_4 <- function(ranges = sex_ranges, var_name = "sexual_score"){
  
  # Create the species richness raster.
  spec_raster <- fasterize(ranges, raster_template, fun = "sum")
  spec_raster[which(getValues(spec_raster) < 5)] <- NA
  spec_raster <- mask(spec_raster, land)
  
  # Create the average raster.
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  return(var_raster)
}

# Make a function for binning raster.
bin_raster <- function(raster, nbins = 10){
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(raster), breaks = breaks, include.lowest = T)
  
  # Replace cell values with bins.
  values(raster) <- as.numeric(cuts)
  
  # Return 'binned' raster.
  return(raster)
}


# Make a function with the legend as bins.
ggplot_raster <- function(raster, nbins = 6, variable = ""){
  
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)
  
  # Replace cell values with bins.
  raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  # Plot with ggplot.
  ggplot() +
    
    # Add the raster data.
    geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    scale_y_continuous(limits = c(-57, 85), expand = expansion(), ) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.line = element_blank(),
                                          legend.position = c(0.13, 0.2),
                                          legend.key.height = unit(0.2, 'cm'),
                                          legend.spacing.y = unit(0.3, 'cm'),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = NA, colour = "lightgrey"),
                                          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
                                          plot.margin = margin(t = 0, l = 0, b = 0, r = 0)
                                          #  plot.margin = margin(t =0, l = -10, b = 0, r = 0.5, unit = "cm")#unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") 
}


# Make a function to plot a raster using geom_raster.
ggplot_raster <- function(raster, nbins = 6, variable = ""){
  
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)
  
  # Replace cell values with bins.
  raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  # Plot with ggplot.
  ggplot() +
    
    # Add the raster data.
    geom_raster(aes(x=long, y=lat, fill= values), data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    scale_y_continuous(limits = c(-57, 85), expand = expansion(), ) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.line = element_blank(),
                                          legend.position = c(0.13, 0.2),
                                          legend.key.height = unit(0.2, 'cm'),
                                          legend.spacing.y = unit(0.3, 'cm'),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = NA, colour = "lightgrey"),
                                          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
                                          plot.margin = margin(t = 0, l = 0, b = 0, r = 0)
                                          #  plot.margin = margin(t =0, l = -10, b = 0, r = 0.5, unit = "cm")#unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") 
}


# This function takes the ggplot and adds labels so it can be incorporated as a multi panel figure.
panel_ggplot_raster <- function(raster_map, plot_title, plot_label){
  ggplot_raster(raster_map, 6, "Sexual selection") +
    annotate("text", x = 20, y = -48, label  = plot_title, size = 8, fontface = 2) +
    annotate("text", x = -170, y = 80, label  = plot_label, size = 12, fontface = 2)
}


# function for making the side plots.
lat_side_plot <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1),
                          lab_x_pos = 60, lab_ypos = 1, plot_label = "b"){
  linear_model <- lm(trait_mean ~ binned_lat, data = data_set)
  linear_sum <- summary(linear_model)
  r_squared <-  as.character(format(round(linear_sum$r.squared, 2), nsmall = 2))
  p_value <- linear_sum$coefficients[2,4]
  if(p_value < 0.001){
    cor_label <- paste0("p < 0.001", "\nR\u00b2 = ", r_squared)
  } else if (p_value > 0.005) {
    cor_label <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)), "\nR\u00b2 = ", r_squared)
  } else {
    cor_label <- paste0("p = ", as.character(format(round(p_value, 3), nsmall = 3)), "\nR\u00b2 = ", r_squared)
  }
  
  
  ggplot(data_set, aes(x = binned_lat, y = trait_mean)) +
    geom_errorbar(aes(ymin = trait_min, ymax = trait_max), 
                  position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
    geom_point(position = position_dodge(width = 1), col = "black") + 
    geom_smooth(se = FALSE, method = "lm", col = "black") + 
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, clip = 'off') +
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    theme(legend.position = "none", 
          text = element_text(face = "bold"),
          axis.title.y = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) +
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2)
}


# Make a function with the legend as bins.
ggplot_legend_raster <- function(raster, nbins = 10, variable = ""){

# Get the limits for breaks, using quantiles so each bin has an equal number of cells.
breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)

# Cut the cell values into bins using the breaks.
cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)

# Replace cell values with bins.
raster@data@values <- as.numeric(cuts)

# Create raster data.
raster_data <- as.data.frame(raster, xy=TRUE)
colnames(raster_data) <- c("long", "lat", "values")
raster_data$values %<>% as.factor()

# Get rid of areas without land.
raster_data$land <- land_data$layer
raster_data %<>% drop_na(land)

# Create a label of the bins for plotting. Unicode is for en-dash.
labels <- as.character(format(round(breaks, 2), nsmall = 2))
labels <- paste0(labels[1:10], " \u2013 ", labels[2:11])
#labels <- c("", labels)
# Plot with ggplot.
ggplot() +
  
  # Clip to world extent.
  xlim(-180, 180) +  
  
  # Add the raster data.
  geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
  
  # Specify colours, legend labels and legend title.
  scale_fill_manual(values = colorRampPalette(pal)(nbins), breaks = 1:10, labels = labels, na.value = "lightgrey") +     #na.value = "grey"
  guides(fill = guide_legend(variable, byrow = TRUE)) +
  
  # Make map closer to the edge.
  scale_x_discrete(expand=c(0,1))+
  scale_y_discrete(expand=c(0,1))+
  # Theme stuff.
  theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                        axis.ticks = element_blank(),
                                        axis.line = element_blank(),
                                        legend.position = c(0.075, 0.3),
                                        legend.key.height = unit(0.2, 'cm'),
                                        legend.spacing.y = unit(0.3, 'cm'),
                                        legend.title = element_text(face = "bold"),
                                        plot.margin = margin(t = -20, l = -20, b = -50, r = 5)     #unit(c(0,0,0,0), "null")
  ) +
  ylab("") + 
  xlab("") + 
  coord_fixed()


}



# Make a function with the legend as bins.
ggplot_terr_raster <- function(raster, nbins = 6, variable = ""){
  
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster)[values(raster)>0], seq(0,1,length.out=nbins), na.rm=TRUE)
  breaks <- c(0,breaks) 
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)
  #cuts <-cut(values(raster), breaks = c(0,0.05,0.1,0.15,0.2,0.25, 0.66666667), include.lowest = TRUE)
  # Replace cell values with bins.
  raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  #labels <- c("", labels)
  # Plot with ggplot.
  ggplot() +
    
    # Clip to world extent.
    xlim(-180, 180) +  
    ylim(-57, 80) +
    
    # Add the raster data.
    geom_raster(aes(x=long, y=lat, fill= values), data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    #scale_x_discrete(expand=c(0,1))+
    #scale_y_discrete(expand=c(0,1))+
    scale_y_continuous(limits = c(-57, 85), expand = expansion()) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.line = element_blank(),
                                          # legend.position = c(0.075, 0.3),
                                          legend.position = c(0.13, 0.2),
                                          legend.key.height = unit(0.2, 'cm'),
                                          legend.spacing.y = unit(0.3, 'cm'),
                                          #legend.title = element_text(face = "bold"),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = NA, colour = "lightgrey"),
                                          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
                                          #plot.margin = margin(t = 15, l = -20000, b = -150, r = 2)
                                          #plot.margin = margin(t =0, l = 0.5, b = 0, r = 0.5, unit = "cm")#unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") #+ 
  #coord_fixed() + 
  #coord_cartesian(ylim=c(-50, 80))
  
  
}



# Make a function with the legend as bins.
ggplot_cont_raster <- function(raster, nbins = 6, variable = ""){
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  #raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  #labels <- c("", labels)
  # Plot with ggplot.
  ggplot() +
    
    # Clip to world extent.
    xlim(-180, 180) +  
    ylim(-57, 80) +
    
    # Add the raster data.
    geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_gradient(low = pal[1], high = pal[5], na.value = "lightgrey") +
    #scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    #scale_x_discrete(expand=c(0,1))+
    #scale_y_discrete(expand=c(0,1))+
    scale_y_continuous(limits = c(-57, 80), expand = expansion()) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.line = element_blank(),
                                          # legend.position = c(0.075, 0.3),
                                          legend.position = c(0.13, 0.2),
                                          legend.key.height = unit(0.2, 'cm'),
                                          legend.spacing.y = unit(0.3, 'cm'),
                                          #legend.title = element_text(face = "bold"),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = NA, colour = "lightgrey"),
                                          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
                                          #plot.margin = margin(t = 15, l = -20000, b = -150, r = 2)
                                          plot.margin = margin(t =0, l = -10, b = 0, r = 0.5, unit = "cm")#unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") #+ 
  #coord_fixed() + 
  #coord_cartesian(ylim=c(-50, 80))
  
  
}

# Function to make the histogram of the colour scale for inserting in plots.
ggplot_colour_hist <- function(raster, legend_title, x_axis_breaks, 
                               x_axis_lim, scaled = TRUE, nbins = 10){
  
  # First convert female elaboration raster into xy data.
  plot_data <- as.data.frame(raster, xy=TRUE) %>% drop_na()
  colnames(plot_data) <- c("long", "lat", "values")
  
  # Scale the values so that zero is near the middle and looks nicer.
  if (scaled){
    plot_data$values <- scale(plot_data$values)
  }
  
  # Create breaks for the new scaled values.
  breaks <- quantile(plot_data$values, seq(0,1,length.out=(nbins+1)), na.rm=TRUE)
  
  
  # Create a plot that could be used as a legend.
  (colour_plot <- ggplot(plot_data,
                         aes(x = values, 
                             # Get fill to be discrete version of x axis, using cuts from quantile breaks.
                             fill = cut(..x.., breaks = breaks, 
                                        include.lowest = TRUE, na.rm = TRUE))) + 
      # Set the number of breaks.
      geom_histogram(bins = 150) +
      
      # Colour scale.
      scale_fill_manual(values = pal, labels = labels) +
      
      # Add breaks for the x-axis, and then cut off the greyyed out bars until we can fix that.
      scale_x_continuous(breaks=x_axis_breaks) + coord_cartesian(xlim = x_axis_lim) +
      
      # X axis label.
      labs(x = legend_title, y = NULL, fill = NULL, colour = NULL) + 
      
      # Standard theme stuff.
      theme_pubr(base_size = 20, legend = "none") +
      theme(axis.text.x=element_text(size=rel(0.7), face = "bold"),
            axis.text.y=element_text(size=rel(0.7), face = "bold"),
            axis.title.x=element_text(size=rel(0.8), face = "bold"),
            axis.line = element_line(size = 1)))
}


################################################################################

# Summarise data for means and CIs for boxplots.
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


################################################################################
                         #####  BRMS functions   ######

# Function to extract and combine brms models.
combine_brms <- function(model_pattern, file_pathways = lat_files){
  pattern_matched <- file_pathways %>% str_subset(pattern = model_pattern)
  if (length(pattern_matched) == 1){
    brms_models <- readRDS(pattern_matched)
  } else {
    
    for (x in 1:length(pattern_matched)){
      if (x ==1){
        brms_model1 <- readRDS(pattern_matched[x])
        print(summary(brms_model1))
      } else if (x == 2){
        brms_model2 <- readRDS(pattern_matched[x])
        print(summary(brms_model2))
        brms_models <- combine_models(brms_model1, brms_model2, check_data = FALSE)
        rm(brms_model1, brms_model2)
      } else {
        brms_model <- readRDS(pattern_matched[x])
        print(summary(brms_model))
        brms_models <- combine_models(brms_models, brms_model, check_data = FALSE)
        rm(brms_model)
      }
      print(x)
    }
  }
  return(brms_models)
}


# Check p-values.
brms_pmap <- function(model){
  p_map(model, regex_pars = "^b_[a-z]")
}

# Quick function to avoid repeating lines.
brms_forest <- function(model){
  mcmc_areas(model, regex_pars = "^b_[a-z]", prob = 0.95, prob_outer = 0.99)
}


# Extract draws function.
extract_draws <- function(model, column_names){
  # Pull out the draws.
  model_draws <- as_draws_df(model, c("^b"), regex = TRUE)[,-c(1:4)] %>% 
    as.data.frame() %>% dplyr::select(-c(.chain, .iteration, .draw))
  # Change column names.
  colnames(model_draws) <- column_names
  # Pivot longer.
  model_draws %>% tidyr::pivot_longer(cols = column_names)
}


# Bayesian McKelvey-Zavoina R2 ------------------------------------------------
# Bayesian version of McKelvey and Zavoina's pseudo-R2 for binary and ordinal
# brms models (McKelvey and Zavoina, 1975). See also Gelman et al. (2018). This
# pseudo-R2 closely approximates the R2 that would have been obtained if a
# linear model had have been run on observations of the continuous latent
# variable underlying the discrete responses (Veall and Zimmermann, 1992; Hagle
# and Mitchell, 1992; Veall and Zimmermann, 1994).
Bayes_R2_MZ <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
  var_fit <- apply(y_pred, 1, var)
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
    var_res <- sig_res^2
  } 
  R2_MZ <- var_fit / (var_fit + var_res)
  print(
    data.frame(
      Estimate = mean(R2_MZ), 
      Est.Error = sd(R2_MZ), 
      "l-95% CI" = quantile(R2_MZ, 0.025),
      "u-95% CI" = quantile(R2_MZ, 0.975),
      row.names = "Bayes_R2_MZ", 
      check.names = FALSE), 
    digits = 2)
}


# Marginal version
marginal_R2_MZ <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, re_formula = NA, ...)
  var_fit <- apply(y_pred, 1, var)
  
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
    var_res <- sig_res^2
  } 
  R2_MZ <- var_fit / (var_fit + var_res)
  print(
    data.frame(
      Estimate = mean(R2_MZ), 
      Est.Error = sd(R2_MZ), 
      "l-95% CI" = quantile(R2_MZ, 0.025),
      "u-95% CI" = quantile(R2_MZ, 0.975),
      row.names = "Bayes_R2_MZ", 
      check.names = FALSE), 
    digits = 2)
}


