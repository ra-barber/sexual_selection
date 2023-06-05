###############################################################################
                 # Model sexual selection on the HPC #
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

# For testing on a local machine. 
#setwd("Manuscripts/")



##############################################################################
#### Functions #####


# Filter the dataset for medium certainty.
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
responses <- c("dichro_total", "dichro_vis", "dichro_hid", "female_orn",
               "head_dichro", "upper_dichro", "under_dichro",
               "wing_dichro", "tail_dichro",
               
               # Binary response variables.
               "dichro_bin", "dichro_vis_bin", "dichro_hid_bin", "female_orn_bin",
               "head_bin", "upper_bin", "under_bin", "wing_bin", "tail_bin",
               "bright_mono")

# Set the data types.
data_type <- c("all", "high", "nonbrood", "nointer", "female_elab", "female_noelab")

# Expand the grid.
all_combos <- expand.grid(responses, data_type)

# Response.
response <- all_combos[array_number, 1] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Set the family for the models.
if (response %in% c("dichro_total", "dichro_vis", "dichro_hid", "female_orn",
                    "head_dichro", "upper_dichro", "under_dichro",
                    "wing_dichro", "tail_dichro")){
  family <- "poisson"
  
} else {
  family <- "binomial"
}




###############################################################################
#### Read in the data #####


# Read in the tree.
model_tree <- read.tree("../Data/Trees/dichromatism_trees.tre")[[1]]

# Read in the life history traits.
model_data <- read.csv("../Data/global_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$bird_tree_name)

# Filter for data type.
if (data_type == "high"){
  model_data %<>% filter_high()
}
if (data_type == "nonbrood"){
  model_data %<>% filter(nest_placement != "Brood Parasite")
}
if (data_type == "female_elab"){
  model_data %<>% filter(female_orn_bin == 1)
}
if (data_type == "female_noelab"){
  model_data %<>% filter(female_orn_bin == 0)
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


# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  incubation_binary = relevel(as.factor(incubation_binary), ref = "Both Sexes"),
  placement_binary = relevel(as.factor(placement_binary), ref = "Concealed"),
  predation_disparity = relevel(as.factor(predation_disparity), ref = "Weak"),
  sexual_binary = relevel(as.factor(sexual_binary), ref = "Weak"),
  social_selection = relevel(as.factor(social_selection), ref = "Weak"),
  habitat_binary = relevel(as.factor(habitat_binary), ref = "Open"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  cooperative = relevel(as.factor(cooperative), ref = "Non-Cooperative"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary")
)

# Center categorical predictors.
model_data %<>% mutate(
  incubation_roles_c = center_categorical(incubation_binary),
  nest_placement_c = center_categorical(placement_binary),
  nest_predation_c = center_categorical(predation_disparity),
  sex_bin_c = center_categorical(sexual_binary),
  habitat_bi_c = center_categorical(habitat_binary),
  migration_bi_c = center_categorical(migration_binary),
  cooperative_c = center_categorical(cooperative),
  social_selection_c = center_categorical(social_selection),
  trophic_level_c = center_categorical(trophic_binary)
)

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE)
)





################################################################################
##### Phyr models ######


# Get the predictors
predictors <- " ~ sex_bin_c + social_selection_c + nest_predation_c + 
              trophic_level_c + migration_bi_c + habitat_bi_c +
              cooperative_c + body_mass_log_z + centroid_z"

labels <- c("Sexual Selection", "Social Selection", "Nest Predation", 
            "Primary Consumer", "Long Migrants", "Habitat Density", 
            "Cooperative Breeder", "Body Mass", "Latitude")

if (data_type == "nointer"){
  predictors <- " ~ sex_bin_c + social_selection_c + 
                incubation_roles_c + nest_placement_c + 
                trophic_level_c + migration_bi_c + habitat_bi_c +
                cooperative_c + body_mass_log_z + centroid_z"
  
  labels <- c("Sexual Selection", "Social Selection", "Single Incubators",
              "Exposed Nesters",
              "Primary Consumer", "Long Migrants", "Habitat Density", 
              "Cooperative Breeder", "Body Mass", "Latitude")
}

# Make the formula
formula <- paste0(response, predictors)


# # # test sample size. female orn took 990 seconds to run 3000 species.16 mins.
# x <- 1000
# model_data <- model_data[sample(1:nrow(model_data), x),]
# 
# # Drop tips on the tree.
# model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Run the model.
tic()
if (family == "binomial"){
  pglmm_model <- pglmm_compare(formula, data=model_data, 
                               phy = model_tree, family = family, s2.init = 0.05, #verbose = TRUE, maxit = 10000
  )
} else {
  pglmm_model <- pglmm_compare(formula, data=model_data, 
                               phy = model_tree, family = family #, s2.init = 0.1, #verbose = TRUE, maxit = 10000
  )
}

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
plot(pglmm_model$mu ~ pglmm_model$Y, xlab = "Dichromatism", 
     ylab = "Fitted", main = plot_title)


# # How to calculate a logisitic pseudo r-squared.
# log_test_data <- cbind(names(pglmm_model$Y), pglmm_model$mu, pglmm_model$Y) %>% as.data.frame()
# colnames(log_test_data) <- c("bird_tree_name", "fitted", "data")
# log_test_data$fitted %<>% as.numeric()
# log_test_data$data %<>% as.numeric()
# 
# mean_1 <- log_test_data %>% filter(data == 1) %>% pull(fitted) %>% mean()
# mean_2 <- log_test_data %>% filter(data == 0) %>% pull(fitted) %>% mean()
# 
# mean_1 - mean_2

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
