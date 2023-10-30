###############################################################################
                        ##### Extra functions #####
###############################################################################

# Functions removed from main functions script which might need to be added in again.






###############################################################################
                 ##### Plotting functions ######


# Function for looking at sexual selection vs predictor boxplots.
sex_boxplot <- function(predictor){
  ggplot(phylo_data, aes(x = predictor, y = sexual_score)) +
    geom_boxplot() + ylab("Sexual Selection") + 
    stat_summary(fun=mean, geom="point", shape=20, size=8) +
    coord_cartesian(ylim=c(0,4)) + theme_classic()
}

# Custom function for hex scatter plot in ggpairs.
ghex_func <- function(data, mapping, ...) { 
  
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ggplot(data) +
    geom_hex(mapping, bins = 50) + scale_fill_viridis_c() + theme(panel.grid.major = element_line(size = 1))
  
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



################################################################################
                        ##### Making maps ######



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


# Function to make the histogram of the colour scale for inserting in plots.
ggplot_colour_hist <- function(raster, legend_title, x_axis_breaks, hist_pal = pal,
                               hist_bins = 150,
                               x_axis_lim, scaled = TRUE, nbins = 10){
  
  # Make a pal that matches number of bins.
  hist_pal <- colorRampPalette(hist_pal)(nbins)
  
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
                                        include.lowest = TRUE, na.rm = TRUE),
                             colour = cut(..x.., breaks = breaks, 
                                          include.lowest = TRUE, na.rm = TRUE)
                         )) + 
      # Set the number of breaks.
      geom_histogram(bins = hist_bins) +
      
      # Colour scale.
      scale_fill_manual(values = hist_pal, labels = labels) +
      scale_colour_manual(values = hist_pal, labels = labels) +
      # Add breaks for the x-axis, and then cut off the greyyed out bars until we can fix that.
      scale_x_continuous(breaks=x_axis_breaks) + coord_cartesian(xlim = x_axis_lim) +
      
      # X axis label.
      labs(x = legend_title, y = NULL, fill = NULL, colour = NULL) + 
      
      # Standard theme stuff.
      theme_pubr(base_size = 30, legend = "none") +
      theme(axis.text.x=element_text(size=rel(0.7)),
            #axis.text.y=element_text(size=rel(0.7), face = "bold"),
            axis.title.x=element_text(size=rel(0.8), face = "bold"),
            axis.line = element_line(size = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA)))
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

################################################################################
                   ##### Brms functions #####


# Check p-values.
brms_pmap <- function(model){
  p_map(model, regex_pars = "^b_[a-z]")
}

# Quick function to avoid repeating lines.
brms_forest <- function(model){
  mcmc_areas(model, regex_pars = "^b_[a-z]", prob = 0.95, prob_outer = 0.99)
}


