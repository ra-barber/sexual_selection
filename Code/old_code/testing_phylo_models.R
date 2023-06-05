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
source("Code/functions.R")



###############################################################################
                     #### Read in the data #####


# Read in the tree.
model_tree <- read.tree("Data/Trees/model_trees.tre")[[1]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$bird_tree_name)

# Filter for data type.
# if (data_type == "high"){
#   model_data %<>% filter_high()
# }

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
  range_z = standardize(range_log, two_sd = TRUE),
  temp_stab_z = standardize(log(bio4), two_sd = TRUE),
  rain_stab_z = standardize(bio15, two_sd = TRUE),
  npp_z = standardize(npp, two_sd = TRUE)
)



################################################################################
                 ##### VIF Values ######

vif_lm <- lm(sexual_score ~ trophic_level_c + terr_bi_c + migration_bi_c + habitat_bi_c +
                body_mass_log_z + range_z + trophic_level_c*centroid_z, data = model_data)

vif(vif_lm, type = "predictor")



################################################################################
                     ##### Phyr models ######

colnames(model_data)



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

library(phylolm)

phylolm_model <- phylolm(sexual_score ~  trophic_level_c + migration_bi_c + habitat_bi_c +
          body_mass_log_z + range_z + trophic_level_c*centroid_z + temp_stab_z + rain_stab_z + npp_z,
        data = model_data, phy = model_tree)

summary(phylolm_model)

library(car)

qqPlot(phylolm_model$residuals)
hist(phylolm_model$res, breaks = 200)
model_data$sexual_binary %<>% as.factor() %>% relevel(ref = "Weak")
model_data$sexual_binary %<>% as.numeric()
model_data$sexual_binary <- model_data$sexual_binary -1

ig10_model <- phyloglm(sexual_binary ~  trophic_level_c + migration_bi_c + habitat_bi_c +
           body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + rain_stab_z + npp_z, 
         data = model_data, phy = model_tree, method = "logistic_IG10")

summary(ig10_model)

mple_model <- phyloglm(sexual_binary ~  trophic_level_c + migration_bi_c + habitat_bi_c +
                         body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + rain_stab_z + npp_z, 
                       data = model_data, phy = model_tree, method = "logistic_MPLE")

summary(mple_model)

library(phyr)



#test sample size.
x <- 1000
test_data <- model_data[sample(1:nrow(model_data), x),]

# Drop tips on the tree.
test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Run the model.
tic()
pglmm_model <- pglmm_compare(sexual_binary ~  trophic_level_c + migration_bi_c + habitat_bi_c +
                               body_mass_log_z + trophic_level_c*centroid_z + 
                               temp_stab_z + rain_stab_z + npp_z, data=test_data, s2.init = 0.05, 
                             phy = test_tree, family = "binomial")
toc()

summary(pglmm_model)


tic()
pglmm_model <- pglmm_compare(sexual_binary ~  trophic_level_c + migration_bi_c + habitat_bi_c +
                               body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + 
                               rain_stab_z + npp_z, data=test_data,# s2.init = 0.05, 
                             phy = test_tree, family = "zeroinflated.binomial", bayes = TRUE)
toc()
summary(pglmm_model)

hist(pglmm_model$H, breaks = 100)
hist(test_data$sexual_binary)

test_glm <- glm(sexual_score ~  trophic_level_c + migration_bi_c + habitat_bi_c +
      body_mass_log_z + trophic_level_c*centroid_z, data=test_data, family = poisson)

library(performance)
check_zeroinflation(test_glm)




tic()
pglmm_model <- pglmm_compare(sexual_score ~  trophic_level_c + migration_bi_c + habitat_bi_c +
                               body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + rain_stab_z + npp_z, data=test_data,# s2.init = 0.05, 
                             phy = test_tree, family = "zeroinflated.poisson", bayes = TRUE)
toc()


summary(pglmm_model)


###############################################################################
#### Run brms models ####


library(future)
plan(multisession)



#test sample size.
x <- 500
test_data <- model_data[sample(1:nrow(model_data), x),]

# Drop tips on the tree.
test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix.
test_covar <- ape::vcv.phylo(test_tree)

# Make sure the covariance matrix and data are in the same order.
row.names(test_data) <- test_data$tree_tip
test_data <- test_data[row.names(test_covar),]



#med_data$binary_dichromatism <- ifelse(med_data$dichro_total > 0, 1, 0)

library(brms)

tic()
brms_model <- brm(
  sexual_binary ~  trophic_level_c + migration_bi_c + habitat_bi_c +
    body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + 
    rain_stab_z + npp_z + (1|gr(tree_tip, cov=A)),
  data = test_data,
  data2 = list(A=test_covar),
  family = bernoulli(),
  #prior = normal_priors,
  iter = 300,
  warmup = 100,
  chains = 1,
  cores = 1,
  #control = list(adapt_delta=0.99, max_treedepth = 12), # no control for now to see how it goes.
  #file = "Models/bern_hi_total",
  #backend = "cmdstanr",
  #threads = threading(8),
)
toc()

###############################################################################
#### Run brms models ####


library(future)
plan(multisession)

# Formula
formula <- wing_abs_s ~ sex_bin + territory_bi + hand_wing_z + body_mass_log_z# + (1|gr(jetz_name, cov = A))

# Add on a super super small number so we can fit a gamma model.
hi_data$wing_abs <- abs(hi_data$wing_ratio) + 1e-25
hi_data$wing_abs_s <- scale(hi_data$wing_abs, center = FALSE)

# gamma_priors <- c(prior(normal(0,2),class="Intercept"),
#                   prior(normal(0,2),class="b"),
#                   prior(gamma(1,0.5),class="shape"))




hi_wing_model <- brm(formula,
                     data = hi_data,
                     #data2 = list(A = hi_covar),
                     family = exponential(link="log"),
                     #prior = gamma_priors,
                     iter = 2000,
                     warmup = 500,
                     chains = 4,
                     cores = 4,
                     # no control for now to see how it goes.
                     #control = list(adapt_delta=0.90, max_treedepth = 12), 
                     #file = "wing_hi_model",
                     #backend = "cmdstanr",
                     #threads = threading(4)
)












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
