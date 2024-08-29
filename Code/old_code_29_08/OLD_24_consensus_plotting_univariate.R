###############################################################################
# Assess the fit/plot univariate models  #
###############################################################################

# This code won't run properly yet because i change the object name from "deets" to "data"



# Packages to load.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(bayesplot)
library(ggdist)
library(ggbeeswarm)
library(stringr)
library(bayestestR)

# Clear the workspace.
rm(list=ls())

# Functions.
source("Code/functions.R")

###############################################################################
#### Read in the models #####


# Read in the actual models. (This takes up loads of RAM)
all_temp_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/temp_centered_all_models.rds")
all_mig_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/mig_uncentered_all_models.rds")
all_tro_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/tro_uncentered_all_models.rds")
all_terr_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/terr_centered_all_models.rds")

tic()
temp_p_value <- pd_to_p(last(p_direction(all_temp_models)[,2]))
mig_p_value <- pd_to_p(last(p_direction(all_temp_models)[,2]))
tro_p_value <- pd_to_p(last(p_direction(all_tro_models)[,2]))
terr_p_value <- pd_to_p(last(p_direction(all_terr_models)[,2]))
toc()

### Output ###
# > temp_p_value
# [1] 0
# > mig_p_value
# [1] 0
# > tro_p_value
# [1] 0.00725
# > terr_p_value
# [1] 0

rm(all_temp_models, all_mig_models, all_tro_models, all_terr_models)


###############################################################################
#### Read in the model data for plotting #####


# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/"
first_half <- "Z:/home/sexual_selection/Results/Models/Old_models/Combined_models/Univariate/"

list.files(first_half, pattern = "data")

# Read in the models using all data certainties.
#temp_cen_all_data <- readRDS(paste0(first_half, "temp_centered__data.rds"))
temp_uncen_all_data <- readRDS(paste0(first_half, "temp_uncentered__data.rds"))
#mig_cen_all_data <- readRDS(paste0(first_half, "mig_centered__data.rds"))
mig_uncen_all_data <- readRDS(paste0(first_half, "mig_uncentered__data.rds"))
#tro_cen_all_data <- readRDS(paste0(first_half, "tro_centered__data.rds"))
tro_uncen_all_data <- readRDS(paste0(first_half, "tro_uncentered__data.rds"))
#terr_cen_all_data <- readRDS(paste0(first_half, "terr_centered__data.rds"))
terr_uncen_all_data <- readRDS(paste0(first_half, "terr_uncentered__data.rds"))

temp_uncen_all_data[[1]]
temp_uncen_all_data[[3]]
mig_uncen_all_data[[3]]
tro_uncen_all_data[[3]]
terr_uncen_all_data[[3]]

## Reading in all models takes up too much RAM ## 

# Read in the models using high data certainties.
temp_cen_high_data <- readRDS(paste0(first_half, "temp_centered_high_data.rds"))
# temp_uncen_high_data <- readRDS(paste0(first_half, "temp_uncentered_high_data.rds"))
mig_cen_high_data <- readRDS(paste0(first_half, "mig_centered_high_data.rds"))
# mig_uncen_high_data <- readRDS(paste0(first_half, "mig_uncentered_high_data.rds"))
tro_cen_high_data <- readRDS(paste0(first_half, "tro_centered_high_data.rds"))
# tro_uncen_high_data <- readRDS(paste0(first_half, "tro_uncentered_high_data.rds"))
terr_cen_high_data <- readRDS(paste0(first_half, "terr_centered_high_data.rds"))
# terr_uncen_high_data <- readRDS(paste0(first_half, "terr_uncentered_high_data.rds"))
# 
# 


tro_example_model <- readRDS("Z:/home/sexual_selection/Results/Models/Old_models/Univariate/tro_uncentered_1_1.rds")


fruit_example_model <- readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Tropics/ordinal_fruit_all.rds")

conditional_effects(fruit_example_model, method = "posterior_linpred")
conditional_effects(fruit_example_model, method = "posterior_predict")[[1]]

fruit_brms_predict <- posterior_interval(fruit_example_model)
test_2 <- posterior_predict(fruit_example_model, point_estimate = "mean")

str(test_2)


cat <- conditional_effects(fruit_example_model, categorical = TRUE)
cat[[1]]
cat[[1]] %>% group_by(trop_non_trop) %>% summarise(mean_est = mean(estimate__), 
                                                  est_prod = mean(as.numeric(effect2__)*estimate__),
                                                  se_prod = mean(as.numeric(effect2__)*se__))


cat_2 <- conditional_effects(fruit_example_model)

?conditional_effects
cat_2[[1]]
cat_2$Condition
cond_plot_lin <- conditional_effects(tro_example_model, method = "posterior_linpred")

# Try extracting predictions, to manually create plots.
post_epred <- brms::posterior_epred(tro_example_model)
brms_predict <- posterior_predict(tro_example_model)

# Try with multiple cores to speed it up.
brms_predict <- posterior_predict(tro_example_model, cores = 4)

# Looko at the structure of the predictions.
str(post_epred)
str(brms_predict)

# Pull out the data in the right order.
trophic_model_pred_data <- tro_example_model$data


# Pull out the model predictions.
trophic_model_pred_data$preds <- apply(brms_predict, 2, mean)

# Create a function for standard error.
standard_error <- function(data){
  sd(data)/sqrt(length(data))
}


# Create a function for standard error, converting back to SS score first.
standard_error_2 <- function(data){
  data <- data - 1
  sd(data)/sqrt(length(data))
}


# Average standard error.
trophic_model_pred_data$se <- apply(brms_predict, 2, standard_error)
apply(brms_predict, 2, standard_error_2)

# Create a function for standard error.
lower_ci <- function(data){
  quantile(data, 0.025)
}

# Average low ci (doesn't work)
trophic_model_pred_data$lower_ci <- apply(brms_predict, 2, lower_ci)


# Average the variance.
trophic_model_pred_data$var <- apply(brms_predict, 2, var)


var
# 
# # Try "The product rule" for uncertainties.
# trophic_model_pred_data$se_prod <- trophic_model_pred_data$se/trophic_model_pred_data$


trophic_model_pred_data %>% filter(trophic_level_binary == "Primary") %>% 
  pull(preds) %>% quantile(0.225)



qt(0.999, df = 2777)*(0.474/sqrt(2777))

margin <- qt(0.975,df=n-1)*s/sqrt(n)


# Summmarise.
tro_pred_se <- trophic_model_pred_data %>% group_by(trophic_level_binary) %>% 
  summarise(pred_mean = mean(preds),
            pred_low_ci = quantile(preds, 0.025),
            pred_n = length(preds),
            pred_se = sd(preds)/sqrt(length(preds)),
            se_mean = mean(se),
            pred_var = sum(var),
            pred_sd = sqrt(pred_var),
            se_test = pred_sd/sqrt(pred_n)
           # se_product = sum(se_prod)  # didn't work
            )

# Get the indices for primary and secondary consumer species.
primary_indices <- which(trophic_model_pred_data$trophic_level_binary == "Primary")
secondary_indices <- which(trophic_model_pred_data$trophic_level_binary == "Secondary")


# Pull out primary and secondary means.
mean(brms_predict[,primary_indices])
mean(brms_predict[,secondary_indices])

# test the standard errors.
standard_error(brms_predict[,primary_indices])
standard_error(brms_predict[,secondary_indices])
quantile(brms_predict[,primary_indices], 0.975)
quantile(brms_predict[,primary_indices], 0.025)




rowMeans(brms_predict[1,])

c(brms_predict[,primary_indices])


point_estimate

tro_example_model$data

?posterior_predict

tro_example_model$data
tro_example_model$fit

tro_example_model$fit
extract_draws()

test <- as_draws_df(tro_example_model)
test[,1:10]

# Try with trophic example.
cond_plot <- conditional_effects(tro_example_model)
cond_plot_lin <- conditional_effects(tro_example_model, method = "posterior_linpred")
cond_plot_post <- conditional_effects(tro_example_model, method = "posterior_predict")

cond_plot_null <- conditional_effects(tro_example_model, re_formula = NULL)
cond_plot_lin_null <- conditional_effects(tro_example_model, method = "posterior_linpred",  re_formula = NULL)
cond_plot_post_null <- conditional_effects(tro_example_model, method = "posterior_predict",  re_formula = NULL)

cond_plot_null_robust <- conditional_effects(tro_example_model, re_formula = NULL, robust = FALSE)
cond_plot_lin_null_robust <- conditional_effects(tro_example_model, method = "posterior_linpred",  re_formula = NULL, robust = FALSE)
cond_plot_post_null_robust <- conditional_effects(tro_example_model, method = "posterior_predict",  re_formula = NULL, robust = FALSE)

cond_plot_robust <- conditional_effects(tro_example_model, robust = FALSE)
cond_plot_lin_robust <- conditional_effects(tro_example_model, method = "posterior_linpred", robust = FALSE)
cond_plot_post_robust <- conditional_effects(tro_example_model, method = "posterior_predict", robust = FALSE)


cond_plot_draws <- conditional_effects(tro_example_model, ndraws = 10, cores = 4)

cond_plot_draws_2 <- conditional_effects(tro_example_model, ndraws = 5, cores = 4)


cond_plot_draws_3 <- conditional_effects(tro_example_model, ndraws = 2, cores = 4)

# cond_plot_draws_4 <- conditional_effects(tro_example_model, ndraws = 1, cores = 4) # Got an error.


cond_plot_draws_3[[1]]
cond_plot_draws_2[[1]]

cond_plot[[1]]

cond_plot_draws[[1]]

tro_2_draws <- conditional_effects(tro_uncen_all_data[[1]], ndraws = 2, cores = 2)
terr_2_draws <- conditional_effects(terr_uncen_all_data[[1]], ndraws = 2, cores = 2)
temp_2_draws <- conditional_effects(temp_uncen_all_data[[1]], ndraws = 2, cores = 2)
mig_2_draws <- conditional_effects(mig_uncen_all_data[[1]], ndraws = 2, cores = 2)

terr_uncen_all_data[[1]]
temp_uncen_all_data[[1]]
mig_uncen_all_data[[1]]

tro_2_draws[[1]]


tro_10_draws <- conditional_effects(tro_uncen_all_data[[1]], ndraws = 10, cores = 2)
terr_10_draws <- conditional_effects(terr_uncen_all_data[[1]], ndraws = 10, cores = 2)
temp_10_draws <- conditional_effects(temp_uncen_all_data[[1]], ndraws = 10, cores = 2)
mig_10_draws <- conditional_effects(mig_uncen_all_data[[1]], ndraws = 10, cores = 2)



tro_all_draws <- conditional_effects(tro_uncen_all_data[[1]], cores = 2)
terr_all_draws <- conditional_effects(terr_uncen_all_data[[1]], cores = 2)
temp_all_draws <- conditional_effects(temp_uncen_all_data[[1]], cores = 2)
mig_all_draws <- conditional_effects(mig_uncen_all_data[[1]], cores = 2)


?conditional_effects
# Condtional plot tests:
normal_cond_plot <- conditional_effects(tro_uncen_high_data[[1]])
library(brms)
normal_cond_plot <- conditional_effects(tro_uncen_all_data[[1]])

normal_cond_plot_2 <- conditional_effects(tro_uncen_all_data[[1]], robust = FALSE)

normal_cond_plot_na <- conditional_effects(tro_uncen_all_data[[1]], conditions = NA)
normal_cond_plot_null <- conditional_effects(tro_uncen_all_data[[1]], conditions = NULL)

normal_cond_plot_null_2 <- conditional_effects(tro_uncen_all_data[[1]], re_formula = NULL)

normal_cond_plot_lin <- conditional_effects(tro_uncen_all_data[[1]], method = "posterior_linpred")
normal_cond_plot_lin_2 <- conditional_effects(tro_uncen_all_data[[1]], method = "posterior_linpred", robust = FALSE)


normal_cond_plot_lin_3 <- conditional_effects(tro_uncen_all_data[[1]], method = "posterior_linpred", robust = FALSE, re_formula = NULL)

normal_cond_plot_null_3 <- conditional_effects(tro_uncen_all_data[[1]], robust = FALSE, re_formula = NULL)


method = "posterior_epred" #DEFAULT
method = "posterior_linpred" 


################################################################################
#### Export summary tables ####


# Extract relevant coeffcient information.
tro_all_estimates <- summary(tro_cen_all_data[[1]])$fixed[5,c(1,3,4)]
mig_all_estimates <- summary(mig_cen_all_data[[1]])$fixed[5,c(1,3,4)]
terr_all_estimates <- summary(terr_cen_all_data[[1]])$fixed[5,c(1,3,4)]
temp_all_estimates <- summary(temp_cen_all_data[[1]])$fixed[5,c(1,3,4)]

all_estimates <- rbind(tro_all_estimates, mig_all_estimates, 
                       terr_all_estimates, temp_all_estimates)

# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_univariate_regression.csv", row.names = TRUE)


# High certainty,
tro_high_estimates <- summary(tro_cen_high_data[[1]])$fixed[5,c(1,3,4)]
mig_high_estimates <- summary(mig_cen_high_data[[1]])$fixed[5,c(1,3,4)]
terr_high_estimates <- summary(terr_cen_high_data[[1]])$fixed[5,c(1,3,4)]
temp_high_estimates <- summary(temp_cen_high_data[[1]])$fixed[5,c(1,3,4)]

high_estimates <- rbind(tro_high_estimates, mig_high_estimates, 
                        terr_high_estimates, temp_high_estimates)

# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(high_estimates, "Results/Tables/high_univariate_regression.csv", row.names = TRUE)


###############################################################################
#### Prepare family level data #####

model_data <- read_ss_data()

tropic_se <- model_data %>% group_by(trophic_level_binary) %>% 
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))

terr_se <- model_data %>% group_by(territoriality_binary) %>% 
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))

mig_se <- model_data %>% group_by(migration_binary) %>% 
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))

# Read in model data.
#model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- model_data$scientific_name_bird_tree %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip

# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_prum_clades()

# Average traits for family.
family_average_data_data <- model_data %>% group_by(family_bird_tree) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_selection),
            clade_sum = length(sexual_selection),
            temp_log = mean(log(seasonality)),
            trophic_level_binary = mean(as.numeric(as.factor(trophic_level_binary)))-1,
            migration_binary  = mean(as.numeric(as.factor(migration_binary)))-1,
            territoriality_binary  = mean(as.numeric(as.factor(territoriality_binary)))-1)

# Group by family and trait for binary traits.
average_family <- function(varible){
  var_data <- model_data %>% dplyr::group_by(family_bird_tree, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              sexual_mean = mean(sexual_selection),
              clade_size = length(sexual_selection))
  var_data %>% group_by(family_bird_tree) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}






# Get territory and migration means.
teritory_data <- average_family("territoriality_binary")
migration_data <- average_family("migration_binary")
trophic_data <- average_family("trophic_level_binary")


tropic_se <- trophic_data %>% group_by(trophic_level_binary) %>% 
  summarise(sexual_mean_2 = mean(sexual_mean),
            sexual_se = sd(sexual_mean)/sqrt(length(sexual_mean)))

tropic_se <- trophic_data %>% group_by(trophic_level_binary) %>% 
  summarise(sexual_mean_2 = mean(sexual_mean),
            sexual_se = sd(sexual_mean)/sqrt(length(sexual_mean)))

primary <- trophic_data %>% filter(trophic_level_binary == "Primary")

primary %>% summarise(sexual_mean = mean(sexual_mean), sexual_se = sd(sexual_mean))

trophic_data %>% group_by(trophic_level_binary) %>% 
  summarise(sexual_mean = mean(sexual_mean),
            sexual_se = sd(sexual_mean))

tropic_se <- model_data %>% group_by(trophic_level_binary) %>% 
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))


###############################################################################
                       ##### Extract model estimates #####

terr_estimates <- convert_ordinal_back(terr_uncen_all_data[[1]])
tro_estimates <- convert_ordinal_back(tro_uncen_all_data[[1]])
mig_estimates <- convert_ordinal_back(mig_uncen_all_data[[1]])


###############################################################################
                       #### Prepare plots #####


# Reorder colours for plotting.
prum_clade_colours <- c("#7494EA", # Palaeognathae
                        "#8D0801", # Galloanserae
                        "#77AD78", # Strisores
                        "#5941A9", # Columbaves
                        "#988F2A", # Gruiformes
                        "#136F63", # Aequorlitornithes
                        "#AA4465", # Opisthocomiformes
                        "#C8B8DB", # Accipitriformes
                        "#05299E", # Strigiformes
                        "#EB6424", # Coraciimorphae
                        "#E0CA3C") # Australaves

names(prum_clade_colours) <- c("Palaeognathae", "Galloanserae", "Strisores",
                               "Columbaves", "Gruiformes", "Aequorlitornithes", 
                               "Opisthocomiformes", "Accipitriformes", "Strigiformes",
                               "Coraciimorphae", "Australaves")







# Change predictions back to same scale as original sex scores.
temp_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
mig_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
tro_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
terr_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)

# Estimate
temp_estimate <- as.character(format(round(summary(temp_uncen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
mig_estimate <- as.character(format(round(summary(mig_uncen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
tro_estimate <- as.character(format(round(summary(tro_uncen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
terr_estimate <- as.character(format(round(summary(terr_uncen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))

# Lower CI
temp_lower <- as.character(format(round(summary(temp_uncen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
mig_lower <- as.character(format(round(summary(mig_uncen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
tro_lower <- as.character(format(round(summary(tro_uncen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
terr_lower <- as.character(format(round(summary(terr_uncen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))

# Upper CI
temp_upper <- as.character(format(round(summary(temp_uncen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
mig_upper <- as.character(format(round(summary(mig_uncen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
tro_upper <- as.character(format(round(summary(tro_uncen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
terr_upper <- as.character(format(round(summary(terr_uncen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))

# Plot labels.
temp_cor_label <- paste0("\U03B2 = ", temp_estimate, "\n[", temp_lower, ", ", temp_upper, "]")
mig_cor_label <- paste0("\U03B2 = ", mig_estimate, "\n[", mig_lower, ", ", mig_upper, "]")
tro_cor_label <- paste0("\U03B2 = ", tro_estimate, "\n[", tro_lower, ", ", tro_upper, "]")
terr_cor_label <- paste0("\U03B2 = ", terr_estimate, "\n[", terr_lower, ", ", terr_upper, "]")

# Use calculated p-values instead.
temp_cor_label <- paste0("\U03B2 = ", temp_estimate, "\np < 0.001")
mig_cor_label <- paste0("\U03B2 = ", mig_estimate, "\np < 0.001")
tro_cor_label <- paste0("\U03B2 = ", tro_estimate, "\np < 0.01")
terr_cor_label <- paste0("\U03B2 = ", terr_estimate, "\np < 0.001")


###############################################################################
#### Create continuous side plots ####

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = prum_clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(0,4.5)),
                       theme_classic(base_size = 20),
                       theme(text = element_text(face = "bold"),
                             legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))

# Remove low predictions as family averages don't span as far as single species.
temp_plot_predictions <- temp_uncen_all_data[[2]] %>% filter(temp_log > 3.3)
temp_plot_predictions <- temp_plot_predictions %>% filter(temp_log < 7.2)

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sexual_mean)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1)

max_temp <- temp_plot_predictions$temp_log %>% max()
min_temp <- temp_plot_predictions$temp_log %>% min()
temp_range <- max_temp-min_temp
temp_x_label <- max_temp - (temp_range*0.2)
temp_x_label_2 <- min_temp + (temp_range*0.05)

###############################################################################
#### Create categorical side plots ####

# Plot settings.
categorical_settings <- list(geom_point(aes(group=family_bird_tree, colour = higher_clade, 
                                            size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                        position = position_jitter(0.2)),
                             labs(x = "", y = NULL,  colour = "Clade"),
                             scale_colour_manual(values = prum_clade_colours),
                             scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
                             theme_classic(base_size = 20),
                             theme(text = element_text(face = "bold"),
                                   legend.position = "none",
                                   axis.line = element_line(size = 1),
                                   axis.ticks = element_line(size = 1),
                                   axis.title.x = element_text(size = rel(0.8)),
                                   axis.title.y = element_text(size = rel(0.8))))



# Migration plots.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(sexual_mean))) + 
  categorical_settings + ylim(c(0,2)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + 
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())


migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_pointrange(
    stat = "summary",
    fun.min = standard_error,
    fun.max = standard_error,
    fun = median
  )


mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sexual_mean - sexual_se*1.96, 
                    ymax = sexual_mean + sexual_se*1.96),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data =mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())



# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)


tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = estimate__-(se__*1.96), 
                    ymax = estimate__+(se__*1.96)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(sexual_mean))) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_pred_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(pred_mean - se_mean*1.96), 
                    ymax = sqrt(pred_mean + se_mean*1.96)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_pred_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(pred_mean)), size = 3.6)


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sexual_mean_2 - sexual_se*2, 
                    ymax = sexual_mean_2 + sexual_se*2),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sexual_mean_2), size = 3.6)




trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(sexual_mean))) + 
  categorical_settings + ylim(c(0,2)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(sexual_mean_2 - sexual_se*2), 
                    ymax = sqrt(sexual_mean_2 + sexual_se*2)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(sexual_mean_2)), size = 3.6)



trophic_data %>% ggplot(aes(x=trophic_level_binary, y=log(sexual_mean+1))) + 
  categorical_settings + ylim(c(0,1.5)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = log((sexual_mean_2 - sexual_se)+1), 
                    ymax = log((sexual_mean_2 + sexual_se)+1)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = log(sexual_mean_2+1)), size = 3.6)


# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territoriality_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Add labels.
lab_x_pos <- 2.1
lab_y_pos <- 4.4
temp_plot <- temp_plot + annotate("text", x = 6.5, y = lab_y_pos, label = temp_cor_label, size = 7, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 7, fontface = 2) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 7, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 7, fontface = 2) 



###############################################################################
### Group the plots together ####

side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        tro_plot + ylab("Sexual selection"), 
                        terr_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))
side_plots


###############################################################################
### Non-bold text version ####

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = prum_clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
                       theme_classic(base_size = 24),
                       theme(legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))

# Plot settings.
categorical_settings <- list(geom_point(aes(group=family_bird_tree, colour = higher_clade, 
                                            size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                        position = position_jitter(0.2)),
                             labs(x = "", y = NULL,  colour = "Clade"),
                             scale_colour_manual(values = prum_clade_colours),
                             scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
                             theme_classic(base_size = 24),
                             theme(legend.position = "none",
                                   axis.line = element_line(size = 1),
                                   axis.ticks = element_line(size = 1),
                                   axis.title.x = element_text(size = rel(0.8)),
                                   axis.title.y = element_text(size = rel(0.8))))

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sexual_mean)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


# Migration plots.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), 
                              expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territoriality_binary, y = estimate__), size = 3.6) 


# Non-bold version.
temp_plot <- temp_plot + annotate("text", x = temp_x_label, y = lab_y_pos, label = temp_cor_label, size = 5) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 5) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 5) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 5) 


side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        tro_plot + ylab("Sexual selection"), 
                        terr_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))
side_plots

lab_x_pos_2 <- 0.6
lab_y_pos_2 <- 4.6

# temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
# mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
# tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
# terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 
# 
# side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
#                         mig_plot,
#                         tro_plot + ylab("Sexual selection"), 
#                         terr_plot, 
#                         nrow = 2, ncol = 2, align = "h",
#                         widths = c(1.2,1))


## Rearrange for joes order.

tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 

side_plots <- ggarrange(tro_plot + ylab("Sexual selection"), 
                        mig_plot,
                        terr_plot + ylab("Sexual selection"), 
                        temp_plot,
                        nrow = 2, ncol = 2, align = "h",
                        widths = c(1.2,1))
side_plots



side_plots

######### END #########
###################################################################################






brms_univar_plot <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1), 
                             lab_x_pos = 2.2, lab_ypos = 2.2, plot_label = "b", 
                             plot_model = allbirds_model, sex_score = TRUE,
                             r_include = TRUE){
  
  # Extract predictions from brms model.
  predictions <- conditional_effects(plot_model)[[1]]
  
  # Change sexual score back to original values (have to + 1 for ordinal regression)
  if (sex_score){
    predictions$estimate__ <- predictions$estimate__ - 1
    predictions$lower__ <- predictions$lower__ - 1
    predictions$upper__ <- predictions$upper__ - 1
  }
  
  # Extract p-values.
  p_values <- brms_pmap(plot_model)[2]
  p_value <- last(p_values[,1])
  
  # Extract r-squared.
  r_squared <- Bayes_R2_MZ(plot_model)[[1]]
  
  # Sample size.
  sample_size <- nrow(plot_model$data)
  
  # Estimate
  estimate <- summary(plot_model)$fixed[5,1]
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- summary(plot_model)$fixed[5,1]
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  
  # To include r-squared values in plot.
  if (r_include){
    # With R squared, which is super low without using phylogeny.
    if(p_value < 0.001){
      cor_label <- paste0("p < 0.001", "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    } else if (p_value > 0.005) {
      cor_label <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)), "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    } else {
      cor_label <- paste0("p = ", as.character(format(round(p_value, 3), nsmall = 3)), "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    }
  } else {
    # Without R squared
    if(p_value < 0.001){
      cor_label <- paste0("p < 0.001", "\n\U03B2 = ", estimate)
      cor_label <- paste0("\U03B2 = ", estimate, "\np < 0.001")
    } else if (p_value > 0.005) {
      cor_label <- paste0("\U03B2 = ", estimate, "\np = ", as.character(format(round(p_value, 2), nsmall = 2)))
    } else {
      cor_label <- paste0("\U03B2 = ", estimate, "\np = ", as.character(format(round(p_value, 3), nsmall = 3)))
    }
  }
  # ggplot function for sideplots with annotations.
  ggplot(data_set, aes(x = binned_lat, y = trait_mean)) +
    geom_errorbar(aes(ymin = trait_min, ymax = trait_max), 
                  position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
    geom_point(position = position_dodge(width = 1), col = "black") + 
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    theme(legend.position = "none", 
          text = element_text(face = "bold"),
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) +
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    # geom_ribbon(data = plot_predictions, 
    #             aes(x = complete_latitude,  ymin = (lower__ - 1), ymax = (upper__ - 1)), 
    #             fill = "grey70", colour = NA, alpha = 0.2, inherit.aes = FALSE)  + 
    geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
}


migration_data %>% ggplot(aes(x = migration_binary, y = trait_mean)) +
  geom_errorbar(aes(ymin = trait_min, ymax = trait_max), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
  scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
  coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
  ylab(ylabel) +
  xlab("Latitude") + theme_classic(base_size = 25) + 
  theme(legend.position = "none", 
        text = element_text(face = "bold"),
        axis.title.y = element_text(size = rel(0.85)),
        axis.title.x = element_text(size = rel(0.85)),
        plot.margin = margin(t = 1, l = 0.2, b = 0.2, unit = "cm")) + 
  annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) +
  annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
  # geom_ribbon(data = plot_predictions, 
  #             aes(x = complete_latitude,  ymin = (lower__ - 1), ymax = (upper__ - 1)), 
  #             fill = "grey70", colour = NA, alpha = 0.2, inherit.aes = FALSE)  + 
  geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)

# Settings for binary traits.
plot_settings <- list(geom_point(aes(group=family, colour = higher_clade, 
                                     size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                 position = position_jitter(0.2)),
                      labs(x = "", y = NULL,  colour = "Clade"),
                      scale_colour_manual(values = clade_colours),
                      scale_y_continuous(breaks = c(0, 1, 2, 3)),
                      theme_classic(base_size = 20),
                      theme(text = element_text(face = "bold"),
                            legend.position = "none",
                            axis.line = element_line(size = 1),
                            axis.ticks = element_line(size = 1),
                            axis.title.x = element_text(size = rel(0.8)),
                            axis.title.y = element_text(size = rel(0.8))))
migration_data %>% ggplot(aes(x=migration_binary, y=log(sexual_mean + 1))) + 
  geom_point(aes(group=family, colour = higher_clade, 
                 size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
             position = position_jitter(0.2)) +
  labs(x = "", y = NULL,  colour = "Clade") +
  scale_colour_manual(values = clade_colours) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  theme_classic(base_size = 20) +
  theme(text = element_text(face = "bold") +
          legend.position = "none",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8))) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + theme(axis.text.y = element_blank()) +
  geom_errorbar(data = predictions, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = predictions, inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) +
  annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) 

data_set <- migration_data
