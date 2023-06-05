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

# Custom function for hex scatter plot in ggpairs.
ghex_func <- function(data, mapping, ...) { 
  
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ggplot(data) +
    geom_hex(mapping, bins = 50) + scale_fill_viridis_c() + theme(panel.grid.major = element_line(size = 1))
  
}



# Function to plot phyr models.
phyr_plot <- function(phyr_model){
  center <- phyr_model$B[1:length(labels)+1,1]
  upper <- center + (phyr_model$B.se[1:length(labels)+1]*1.96)
  lower <- center - (phyr_model$B.se[1:length(labels)+1]*1.96)
  p_value <- phyr_model$B.pvalue[1:length(labels)+1]
  forest_data <- data.frame(labels, center, lower, upper, p_value)
  
  forest_data$signif <- "Y"
  forest_data[forest_data$p_value > 0.05,"signif"] <- "N"
  forest_data$labels <- factor(forest_data$labels, levels=rev(forest_data$labels))
  forest_data$signif <- factor(forest_data$signif)
  ggplot(data=forest_data, aes(x=labels, y=center, ymin=lower, ymax=upper, colour=signif)) +
    geom_pointrange(size=0.8, shape=18) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    ylab("Estimate") + expand_limits(y=c(-0.05,0.05)) + #scale_color_gradient(low="orange", high="blue") +
    xlab("")  +
    ggtitle(plot_title) +
    theme_pubr(legend = "none", base_size = 12) +
    theme(axis.text.x=element_text(size=rel(0.7), face = "bold"), 
          axis.text.y=element_text(size=rel(0.7), face = "bold"),
          axis.title.x=element_text(size=rel(0.7), face = "bold"),
          plot.title = element_text(hjust = 1))
}








################################################################################
                  ##### Making maps ######

# Function to make average raster from variable name.
average_raster <- function(var_name, ranges = sex_ranges, spec_raster = species_raster){
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, species_raster)
  var_raster <- var_raster / species_raster
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


# Make a basic function for making maps with ggplot.
ggplot_raster <- function(raster, nbins = 10, variable = ""){
  
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
  labels <- rev(as.character(round(breaks, 2)))
  labels <- paste0(labels[2:11], " \u2013 ", labels[1:10])
  
  # Plot with ggplot.
  ggplot() +
    
    # Clip to world extent.
    xlim(-180, 180) +  
    
    # Add the raster data.
    geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = pal, labels = labels, na.value = "grey") +
    guides(fill = guide_legend(variable)) +
    
    # Make map closer to the edge.
    scale_x_discrete(expand=c(0,1))+
    scale_y_discrete(expand=c(0,1))+
    # Theme stuff.
    theme_classic() + theme(axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"
                            #plot.margin = unit(c(0,0,0,0), "null"),
                            #panel.margin = unit(c(0, 0, 0, 0), "null")
    ) +
    ylab("") + 
    xlab("") + 
    coord_fixed()
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
      theme_pubr(base_size = 12, legend = "none") +
      theme(axis.text.x=element_text(size=rel(0.7), face = "bold"),
            axis.text.y=element_text(size=rel(0.7), face = "bold"),
            axis.title.x=element_text(size=rel(0.8), face = "bold"),
            axis.line = element_line(size = 1)))
}

