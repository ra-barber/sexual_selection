##############################################################################
                        #### Functions #####
###############################################################################

# This script has functions used throughout the sexual selection project.



###############################################################################
                     ##### Data handling ######


# Read in the data without having to load extra packages.
read_ss_data <- function(pathway = "Data/supplementary_dataset_1.xlsx"){
  # Read in the data.
  model_data <- readxl::read_excel(pathway, sheet = 2, na = "NA")
  # Remove upper-case column names.
  model_data <- janitor::clean_names(model_data)
  model_data %<>% dplyr::select(-source)
  # Add in tree tips used in phylo models.
  model_data$tree_tip <- gsub(" ", "_", model_data$scientific_name_bird_tree)
  # Add in dummy territoriality variables used in Figure S9.
  model_data %<>% mutate(
    terr_dummy = 0,
    year_terr_dummy = 0,
    terr_dummy = replace(terr_dummy, territoriality_binary == "Territorial", 1),
    year_terr_dummy = replace(year_terr_dummy, territoriality == "Strong", 1),
    # Absolute latitude.
    abs_lat = abs(latitude))
  # Export model data.
  return(model_data)
}

# Quick function for centering factors.
center_categorical <- function(predictor){
  as.numeric(predictor) %>% scale(scale = FALSE)
}

# Calculate Standard error.
standard_error <- function(vector){
  sd(vector)/sqrt(length(vector))
}


################################################################################
                  ##### Making maps ######


# Function to make species richness raster.
spec_raster_func <- function(ranges = sex_ranges){
  spec_raster <- fasterize(ranges, raster_template, fun = "sum") 
  spec_raster[which(getValues(spec_raster) < 10)] <- NA 
  spec_raster <- mask(spec_raster, land)
  return(spec_raster)
}

# Function to make average raster from variable name.
average_raster <- function(ranges = sex_ranges, var_name = "sexual_selection"){
  
  # Create the species richness raster.
  spec_raster <- spec_raster_func(ranges)
  
  # Create an average raster for the variable and data in question.
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  return(var_raster)
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
    geom_raster(aes(x=long, y=lat, fill= values), data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), 
                      breaks = rev(1:nbins), labels = rev(labels), 
                      na.value = "lightgrey") +
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    scale_y_continuous(limits = c(-57, 85), expand = expansion(), ) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = c(0.13, 0.2),
          legend.key.height = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.3, 'cm'),
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA, colour = "lightgrey"),
          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
          plot.margin = margin(t = 0, l = 0, b = 0, r = 0)) +
    ylab("") + xlab("") 
}


# This function takes the ggplot and adds labels so it can be incorporated as a multi panel figure.
panel_ggplot_raster <- function(raster_map, plot_title, plot_label){
  panel_raster <- ggplot_raster(raster_map, 6, "Sexual selection") +
    annotate("text", x = 20, y = -48, label  = plot_title, size = 8, fontface = 2) +
    annotate("text", x = -170, y = 80, label  = plot_label, size = 12, fontface = 2)
  return(panel_raster)
}

# A custom function used for plotting primary consumer territoriality.
ggplot_terr_raster <- function(raster, nbins = 6, variable = ""){
  breaks <- quantile(values(raster)[values(raster)>0], seq(0,1,length.out=nbins), na.rm=TRUE)
  breaks <- c(0,breaks) 
  cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)
  raster@data@values <- as.numeric(cuts)
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])

  # Plot with ggplot.
  ggplot() +
    xlim(-180, 180) +  
    ylim(-57, 80) +
    geom_raster(aes(x=long, y=lat, fill= values), data=raster_data) +
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    scale_y_continuous(limits = c(-57, 85), expand = expansion()) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = c(0.13, 0.2),
          legend.key.height = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.3, 'cm'),
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA, colour = "lightgrey"),
          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm")) +
    ylab("") + xlab("") 
}


################################################################################
                        ##### Additional plots #######


# Function for making the side plots.
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


################################################################################
                      ##### brms functions ######


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


###############################################################################
                        #### SAR model functions #####


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
    
    scale_colour_manual(values = rev(colorRampPalette(pal)(nbins)), 
                        breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +  
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 1, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 2)
}


################################################################################
                          ##### End ######
################################################################################