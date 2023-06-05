###############################################################################
                  # Phylo models for quick exploration  #
###############################################################################


# Packages to load.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(phyr)

# Clear the workspace.
rm(list=ls())

# Functions.
source("../Code/functions.R")


##############################################################################
#### Functions #####



#phyr_model <- pglmm_model

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
#### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create responses for hpc jobs.
responses <- c("sexual_score")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(responses, data_type)

# Response.
response <- all_combos[array_number, 1] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()




###############################################################################
                  #### Read in the data #####


# Read in the tree.
model_tree <- read.tree("../Data/Trees/model_trees.tre")[[1]]

# Read in the life history traits.
model_data <- read.csv("../Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$bird_tree_name)

# Filter for data type.
if (data_type == "high"){
  model_data %<>% filter_high()
}

# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix.
model_covar <- ape::vcv.phylo(model_tree)

# Make sure the covariance matrix and data are in the same order.
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]

###############################################################################
                ###### Prepare predictor variables ########

# colnames(model_data)

# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  territory_year_round = relevel(as.factor(territory_year_round), ref = "No territory"),
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  
  habitat_binary = relevel(as.factor(habitat_binary), ref = "Open"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_year_c = center_categorical(territory_year_round),
  terr_bi_c = center_categorical(territory_binary),
  habitat_bi_c = center_categorical(habitat_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary)
)


# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  range_z = standardize(range_log, two_sd = TRUE)
)





################################################################################
                      ##### Phyr models ######


# Get the predictors
predictors <- " ~ trophic_level_c + migration_bi_c + habitat_bi_c +
body_mass_log_z + range_z + trophic_level_c*centroid_z"

predictors <- " ~ trophic_binary + migration_binary  + habitat_binary  +
body_mass_log_z + trophic_binary*centroid_z + range_z"

labels <- c("Primary Consumer", "Long Migrants", "Habitat Density", 
             "Body Mass", "Range Size", "Latitude (Secondary Consumers)", 
            "Latitude (Primary Consumers)")

# Make the formula
formula <- paste0(response, predictors)


# # test sample size. 
# x <- 1500
# model_data <- model_data[sample(1:nrow(model_data), x),]
# 
# # Drop tips on the tree.
# model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Run the model.
tic()
pglmm_model <- pglmm_compare(formula, data=model_data, 
                               phy = model_tree, family = "poisson")
toc()

# Save the model.
model_pathway <- paste0("../Results/Models/", data_type, "/", array_number, "_",  response, ".rds")
saveRDS(pglmm_model, model_pathway)

# Get the summary.
summary(pglmm_model)

# Get the adjusted R squared of data vs fitted.
fitted_summary <- lm(pglmm_model$mu ~ pglmm_model$Y) %>% summary()
r_square <- fitted_summary$adj.r.squared %>% round(3)*100

# Pathway to save the model diagnostics.
plot_pathway <- paste0("../Plots/Diagnostics/", data_type, "/", response, ".tiff")

# Make the plot title for the plots.
plot_title <- paste0("Response = ", response, 
                     "      Certainty = ", data_type,
                     "\nModel = ", family,
                     "      Adj-Rsqr = ", r_square, "%")

# Save the base plot.
tiff(plot_pathway, res = 150, width = 1000, height = 2000)
par(mfrow=c(2,1))

# Plot the model fit.
plot(pglmm_model$mu ~ pglmm_model$Y, xlab = "Sexual Score", 
     ylab = "Fitted", main = plot_title)


# Look at residuals vs data.
plot(pglmm_model$H ~ pglmm_model$mu, xlab = "Fitted", ylab = "Residuals")

# Turn device off.
dev.off()

# Save the model results.
plot_pathway <- paste0("../Plots/Results/", data_type, "/", response, ".tiff")

# Plot the model.
phyr_plot(pglmm_model)

# Export the plots.
ggsave(plot_pathway, height = 6, width = 6)

################################################################################
                          ##### End #######
