###############################################################################
                   # Assess the fit of brms models  #
###############################################################################


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

# Pull out the file pathways.
model_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

# Read in the new files.
feb_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/Feb/"), full.names = T, include.dirs = FALSE, recursive = FALSE)



# Function to extract and combine brms models.
combine_brms <- function(model_pattern){
  pattern_matched <- feb_files %>% str_subset(pattern = model_pattern)
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

# Check output of pattern first.
model_files %>% str_subset(pattern = "chick_pc")

# Centered files.
feb_files %>% str_subset(pattern = "/npp_centered_all")
feb_files %>% str_subset(pattern = "/npp_centered_high")
feb_files %>% str_subset(pattern = "no_npp_centered_all")
feb_files %>% str_subset(pattern = "no_npp_centered_high")
feb_files %>% str_subset(pattern = "no_inter_centered_all")
feb_files %>% str_subset(pattern = "no_inter_centered_high")

feb_files %>% str_subset(pattern = "/npp_uncentered_all")
feb_files %>% str_subset(pattern = "/npp_uncentered_high")
feb_files %>% str_subset(pattern = "no_npp_uncentered_all")
feb_files %>% str_subset(pattern = "no_npp_uncentered_high")
feb_files %>% str_subset(pattern = "no_inter_uncentered_all")
feb_files %>% str_subset(pattern = "no_inter_uncentered_high")


# Read in centered models first.
npp_cen_all_models <- combine_brms("/npp_centered_all")
nonpp_cen_all_models <- combine_brms("no_npp_centered_all")
no_inter_cen_all_models <- combine_brms("no_inter_centered_all")
no_inter_uncen_all_models <- combine_brms("no_inter_uncentered_all")

tic()
conditional_effects(no_inter_uncen_all_models, "territory_binary", categorical = TRUE, resolution = 100)
toc()
# Main analysis for development.





# Kimball models.
wang_models <- combine_brms("newref_devo")
wang_inter_models <- combine_brms("wang")
wang_center_models <- combine_brms("cen_devo")
wang_seasonal <- combine_brms("edited_season")

# Edited Kimball models. (Changing some species based on cooney / greisser datasets)
edited_models <- combine_brms("/edited_all")
edited_inter_models <- combine_brms("inter_edit")

# Removing generation length and devo from models.
nogen_models <- combine_brms("nogen_devo")
nodevo_models <- combine_brms("nogen_all")        
     
# Kitchen sink.
kitchen_sink_models <- combine_brms("kitchen_sink")        
      
# Small sensitivity analysis for development.

# Basic small models.
chick_models <- combine_brms("/chick")
fed_models <- combine_brms("/fed")
greisser_models <- combine_brms("/greisser")
inter_chick_models <- combine_brms("/inter_chick")
inter_fed_models <- combine_brms("/inter_fed")
inter_greisser_models  <- combine_brms("/inter_greisser")

# All greisser devo models.
max_greisser_models  <- combine_brms("/max")
min_greisser_models  <- combine_brms("/min")
season_greisser_models  <- combine_brms("/season_greisserall")
cen_season_greisser_models  <- combine_brms("/season_greisser_cenall")


# Chick PCA Models.

cp_centered_models <- combine_brms("chick_pc_centeredall")
cp_everything_but_npp_models <- combine_brms("chick_pc_everything_but_nppall")
cp_everything_models <- combine_brms("chick_pc_everythingall")
cp_inter_models <- combine_brms("chick_pc_interall")
cp_season_centered_models <- combine_brms("chick_pc_season_centeredall")
cp_season_inter_centered_models <- combine_brms("chick_pc_season_inter_centeredall")
cpng_centered_models <- combine_brms("cpng_centeredall")


# 
# cp_inter_models looks good (except for including generation length)






#all_devo_models  <- combine_brms("/prum_devo_all")

gc()

all_model_files <- model_files %>% str_subset(pattern = "/all")
hi_model_files <- model_files %>% str_subset(pattern = "/high")
all_model_files <- model_files %>% str_subset(pattern = "densall")
model_files %>% str_subset(pattern = "prum_high")
model_files %>% str_subset(pattern = "prum_devo_high")
model_files %>% str_subset(pattern = "/devo_high")

model_files %>% str_subset(pattern = "prum_devo_all")
model_files %>% str_subset(pattern = "/devo_all")
model_files %>% str_subset(pattern = "/prum_all")


# Read in a single model. (This model uses the chick PC1 scores estimated)
full_chick_model <- readRDS("Z:/home/sexual_selection/Results/Models/testing_models/49_all.rds")



###############################################################################
                    #### Model Diagnostics #####


# Look at the summaries.
summary(nonpp_cen_all_models)


summary(inter_greisser_models)

# # Plot the model.
plot(brms_models)
plot(hi_brms_models)

# Check the posterior distribution.
pp_check(brms_models)
pp_check(hi_brms_models)




###############################################################################
                       #### Forest Plots #####


# Quick function to avoid repeating lines.
brms_forest <- function(model){
  mcmc_areas(model, regex_pars = "^b_[a-z]", prob = 0.95, prob_outer = 0.99)
}

# Feb models.
brms_forest(npp_cen_all_models)
brms_forest(nonpp_cen_all_models)
brms_forest(no_inter_cen_all_models)


plot(conditional_effects(no_inter_cen_all_models), points = TRUE, ca)


# Plot main devo analysis models.
brms_forest(wang_models)
brms_forest(wang_inter_models)
brms_forest(wang_center_models)
brms_forest(wang_seasonal)

brms_forest(edited_models)
brms_forest(edited_inter_models)

brms_forest(nogen_models)
brms_forest(nodevo_models)

brms_forest(kitchen_sink_models)

# Smaller analysis of devo.

brms_forest(chick_models)
brms_forest(fed_models)
brms_forest(greisser_models)

brms_forest(inter_chick_models)
brms_forest(inter_fed_models)
brms_forest(inter_greisser_models)

brms_forest(max_greisser_models)
brms_forest(min_greisser_models)
brms_forest(season_greisser_models)
brms_forest(cen_season_greisser_models)

brms_forest(full_chick_model)

# Chick PCA models.
brms_forest(cp_centered_models)
brms_forest(cp_everything_but_npp_models)
brms_forest(cp_everything_models)
brms_forest(cp_inter_models)
brms_forest(cp_season_centered_models)
brms_forest(cp_season_inter_centered_models)
brms_forest(cpng_centered_models)


###############################################################################
                          #### P mapping #####


# Check p-values.
brms_pmap <- function(model){
  p_map(model, regex_pars = "^b_[a-z]")
}

# Main devo analysis.
brms_pmap(wang_models)
brms_pmap(wang_inter_models)
brms_pmap(wang_center_models)
brms_pmap(wang_seasonal)

brms_pmap(edited_models)
brms_pmap(edited_inter_models)

brms_pmap(nogen_models)
brms_pmap(nodevo_models)

brms_pmap(kitchen_sink_models)

brms_pmap(chick_models)
brms_pmap(inter_chick_models)
brms_pmap(fed_models)
brms_pmap(inter_fed_models)
brms_pmap(greisser_models)
brms_pmap(inter_greisser_models)

brms_pmap(max_greisser_models)
brms_pmap(min_greisser_models)


brms_pmap(cp_centered_models)
brms_pmap(cp_everything_but_npp_models)
brms_pmap(cp_everything_models)
brms_pmap(cp_inter_models)
brms_pmap(cp_season_centered_models)
brms_pmap(cp_season_inter_centered_models)
brms_pmap(cpng_centered_models)


###############################################################################
                      #### Interactions #####


library(interactions)

model_data <- read.csv("Data/sexual_traits.csv")

interact_plot(cp_centered_models, pred = temp_seasonality_z, modx = trophic_level_c, data = cp_centered_models$data)


preds <- predict(cp_centered_models)

model_fits <- fitted(cp_centered_models, scale = "response")


combine_preds <- function(preds){
  colnames(preds) <- paste("intercept", 1:5, sep="_")
  preds %>% as.data.frame() %>%  mutate(
    model_preds = intercept_1 + intercept_2*2 + intercept_3*3 + intercept_4*4 + intercept_5*5,
    test_preds = intercept_1 + intercept_2 + intercept_3 + intercept_4 + intercept_5
  )
}

preditions <- preds %>% combine_preds() %>% dplyr::select(model_preds, test_preds)

test <- cbind(cp_centered_models$data, preditions)


plot(model_preds ~ gen_z, data = test)

test %>% ggplot(aes(x = gen_z, y = test_preds)) + geom_point(alpha = 0.1) + geom_smooth(method = "lm") +
  theme_classic()


hist(model_fits[3,1,])

model_fits[2,,]

class(model_fits)

test_fits <- model_fits[,1,1] +  model_fits[,1,2]*2 +  model_fits[,1,3]*3 + model_fits[,1,4]*4 + model_fits[,1,5]*5

conditional_effects(cp_centered_models, categorical = TRUE)

plot(x = test$gen_z, y = test_fits)



cp_centered_models$data

cp_centered_models$

# Look at lambda values.
hyp <- "sd_tree_tip__Intercept^2 / (sd_tree_tip__Intercept^2 + disc^2) = 0"

(all_lambda <- hypothesis(brms_models, hyp, class = NULL))
(high_lambda <- hypothesis(hi_brms_models, hyp, class = NULL))


###############################################################################
                 #### Preparing plot variables #####


# Do the column names for each model.
summary(cpng_centered_models)
model_colnames <- c("Generation", "Migration" ,"Territory", "Diet", "Temp", 
                    "NPP", "Diet:Temp", "Diet:NPP", "Diet:Terr")

# Colnames as they should be in the plots.
model_colnames <- c("Generation\nLength", "Migration", "Territoriality", 
                    "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nSeasonality", 
               "PC x \nProductivity",
               "PC x \nTerritoriality")

# Including population density.
model_colnames <- c("Generation\nLength", "Migration", "Territoriality",  "Population\nDensity",
                    "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
                    "PC x \nSeasonality", 
                    "PC x \nProductivity",
                    "PC x \nTerritoriality")

# Including developmental mode.
devo_model_colnames <- c("Generation\nLength", "Altricial", "Migration", "Territoriality", 
                    "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
                    "PC x \nSeasonality", 
                    "PC x \nProductivity",
                    "PC x \nTerritoriality")

# Greisser development interaction names.
devo_model_colnames <- c("Generation\nLength", "(Pr) Precocial", "Migration", "Territoriality", 
                         "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
                         "PC x \nSeasonality", 
                         "PC x \nProductivity",
                         "PC x \nTerritoriality",
                         "Pr x \nSeasonality",
                         "Pr x \nProductivity",
                         "Pr x \nTerritoriality")

# Model colnames for chick development without generation length.
cpng_model_colnames <- c("Developmental\nmode", "Migration", "Territoriality", 
                         "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
                         "PC x \nseasonality", 
                         "PC x \nproductivity",
                         "PC x \nterritoriality")

cp_everything__colnames <- c("Developmental\nmode", "Generation\nLength", "Migration", "Territoriality", 
                             "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
                             "DM x \ngeneration length", "DM x Migration", 
                             "PC x \nseasonality", "DM x PC", 
                             "PC x \nproductivity",
                             "PC x \nterritoriality",
                             "DM x \nseasonality",
                             "DM x \nproductivity",
                             "DM x \nterritoriality")

cp_everything_np_colnames <- c("Developmental\nmode", "Generation\nLength", "Migration", "Territoriality", 
                             "(PC) Primary\nConsumer", "Seasonality",
                             "DM x \ngeneration length", "DM x Migration", 
                             "PC x \nseasonality", "DM x PC", 
                             "PC x \nterritoriality",
                             "DM x \nseasonality",
                             "DM x \nterritoriality")


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

# Extract draws.
all_plotdata <- extract_draws(cpng_centered_models, cpng_model_colnames)
everything_data <- extract_draws(cp_everything_models, cp_everything__colnames)
npp_data <- extract_draws(cp_everything_but_npp_models, cp_everything_np_colnames)
# Vector of predictors in reverse order of plotting.
rev_preds <- c("Generation", "Migration", "Diet", 
               "Territory", "Temp", "NPP", 
               "Diet:Terr", "Diet:Temp", "Diet:NPP")
rev_preds <- c("Generation\nLength", "Migration", "Territoriality", 
               "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nTerritoriality",
               "PC x \nSeasonality", 
               "PC x \nProductivity")
# Including population density.
rev_preds <- c("Territoriality", "Generation\nLength", "Migration", "Population\nDensity",
               "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nTerritoriality",
               "PC x \nSeasonality", 
               "PC x \nProductivity"
)
# Including Devlopmental.
rev_preds <- c("Territoriality", "Generation\nLength", "Migration", "Altricial",
               "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nTerritoriality",
               "PC x \nSeasonality", 
               "PC x \nProductivity"
)
# Including Devlopmental interactions.
rev_preds <- c("Territoriality", "Generation\nLength", "Migration", "(Pr) Precocial",
               "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nTerritoriality",
               "PC x \nSeasonality", 
               "PC x \nProductivity",
               "Pr x \nTerritoriality",
               "Pr x \nSeasonality",
               "Pr x \nProductivity")

rev_preds <- c("Territoriality", "Migration", "Developmental\nmode",
               "(PC) Primary\nConsumer", "Seasonality", "Productivity", 
               "PC x \nterritoriality",
               "PC x \nseasonality", 
               "PC x \nproductivity")



cp_everything_preds <- c("Territoriality", "Generation\nLength", "Migration", 
                             "(PC) Primary\nConsumer", "Developmental\nmode", "Seasonality", "Productivity", 
                             "PC x \nseasonality", 
                             "PC x \nproductivity",
                             "PC x \nterritoriality",
                             "DM x PC", 
                             "DM x \ngeneration length", 
                             "DM x Migration", 
                             "DM x \nseasonality",
                             "DM x \nproductivity",
                             "DM x \nterritoriality")


npp_everything_preds <- c("Territoriality", "Generation\nLength", "Migration", 
                         "(PC) Primary\nConsumer", "Developmental\nmode", "Seasonality",
                         "PC x \nseasonality", 
                         "PC x \nterritoriality",
                         "DM x PC", 
                         "DM x \ngeneration length", 
                         "DM x Migration", 
                         "DM x \nseasonality",
                         "DM x \nterritoriality")


# Set data as factors for plotting order.
all_plotdata$name %<>% factor(levels = rev(rev_preds))
everything_data$name %<>% factor(levels = rev(cp_everything_preds))
npp_data$name %<>% factor(levels = rev(npp_everything_preds))

lifehistory <- c("Territoriality", "Generation\nLength", "Migration", 
                 "Population\nDensity", "(PC) Primary\nConsumer", "(Pr) Precocial", "Developmental\nmode")
resource <-   c("Seasonality", "Productivity")
interaction <-  c("PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity")
devo_interaction <- c("Pr x \nTerritoriality", "Pr x \nSeasonality", "Pr x \nProductivity")
devo_interaction <- c("DM x PC", "DM x \ngeneration length", "DM x Migration", "DM x \nseasonality",
                      "DM x \nproductivity", "DM x \nterritoriality")




# Add colours for plotting.
# lifehistory <- c("Generation", "Migration", "Diet")
# resource <- c("Territory", "Temp", "NPP")
# interaction <- c("Diet:Terr", "Diet:Temp", "Diet:NPP")
all_plotdata$trait <- NA
all_plotdata$trait[all_plotdata$name %in% lifehistory] <- "Life history"
all_plotdata$trait[all_plotdata$name %in% resource] <- "Environmental"
all_plotdata$trait[all_plotdata$name %in% interaction] <- "Trophic interaction"
all_plotdata$trait[all_plotdata$name %in% devo_interaction] <- "Devlopment interaction"

everything_data$trait <- NA
everything_data$trait[everything_data$name %in% lifehistory] <- "Life history"
everything_data$trait[everything_data$name %in% resource] <- "Environmental"
everything_data$trait[everything_data$name %in% interaction] <- "Trophic interaction"
everything_data$trait[everything_data$name %in% devo_interaction] <- "Devlopment interaction"

npp_data$trait <- NA
npp_data$trait[npp_data$name %in% lifehistory] <- "Life history"
npp_data$trait[npp_data$name %in% resource] <- "Environmental"
npp_data$trait[npp_data$name %in% interaction] <- "Trophic interaction"
npp_data$trait[npp_data$name %in% devo_interaction] <- "Devlopment interaction"


# Y axis labels for ordering plot.
axes_labs <- c("Territoriality", "Generation\nlength", "Migration", "Population\ndensity", "(PC) Primary\nconsumer", 
                "Seasonality", "Productivity",
               "PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity")
axes_labs <- c("Territoriality", "Generation\nlength", "Migration",  "Altricial", "(PC) Primary\nconsumer", 
               "Seasonality", "Productivity",
               "PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity")
axes_labs <- c("Territoriality", "Generation\nlength", "Migration",  "(Pr) Precocial", "(PC) Primary\nconsumer", 
               "Seasonality", "Productivity",
               "PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity",
               "Pr x \nTerritoriality", "Pr x \nSeasonality", "Pr x \nProductivity")
axes_labs <- c("Territoriality", "Migration",  "Developmental\nmode", "(PC) Primary\nconsumer", 
               "Seasonality", "Productivity",
               "PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity")


axes_labs <- c("Territoriality", "Migration",  "Developmental\nmode", "(PC) Primary\nconsumer", 
               "Seasonality", "Productivity",
               "PC x \nterritoriality", "PC x \nseasonality", "PC x \nproductivity")
axes_labs <- c("Territoriality", "Generation\nLength", "Migration", 
                         "(PC) Primary\nConsumer", "Developmental\nmode", "Seasonality", "Productivity", 
                         "PC x \nseasonality", 
                         "PC x \nproductivity",
                         "PC x \nterritoriality",
                         "DM x PC", 
                         "DM x \ngeneration length", 
                         "DM x Migration", 
                         "DM x \nseasonality",
                         "DM x \nproductivity",
                         "DM x \nterritoriality")

axes_labs <- c("Territoriality", "Generation\nLength", "Migration", 
               "(PC) Primary\nConsumer", "Developmental\nmode", "Seasonality",
               "PC x \nseasonality", 
               "PC x \nterritoriality",
               "DM x PC", 
               "DM x \ngeneration length", 
               "DM x Migration", 
               "DM x \nseasonality",
               "DM x \nterritoriality")



# Group plot settings.
ggplot_settings <- list(
  # geom_point(alpha = .1, stroke = 0.2, colour = "#586994", size = 0.5, 
  #            position = position_jitter(seed = 1993, height = .35)),
  # stat_pointinterval(size = 9, alpha = rep(c(0.95, 0.6), times = 9), 
  #                    colour = "navy"),
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6), # Add lines.
  scale_y_discrete(labels = rev(axes_labs)),
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL),
  theme_classic2(base_size = 12),   # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm')))



###############################################################################
                    #### Basic scatter plots #####

pred_n <- 13

# Standard plot all in navy.
ggplot(npp_data, aes(x = value, y = name)) + 
  geom_point(alpha = .1, stroke = 0.2, colour = "#586994", size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
stat_pointinterval(size = 9, alpha = rep(c(0.95, 0.6), times = pred_n), 
                   colour = "navy") + ggplot_settings
ggsave("Plots/Results/all_model.tiff", dpi = 600, width = 8, height = 4)




library(ggnewscale)

light_colours <- c( "#77AD78","#7494EA", "#AA4465")
light_colours <- c( "#77AD78","#7494EA", "#C98986")



dark_colours <- c("#136F63", "#05299E", "#8D0801")
dark_colours <- c("#214F4B", "#05299E", "#8D0801")
legend_order <- c("Life history", "Environmental", "Trophic interaction")

# cOLOURS WITH PRECOCIAL INTERACTION.
light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")
legend_order <- c("Life history", "Environmental", "Trophic interaction", "Devlopment interaction")


# Plot with  alpha interval bars.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) + 
  geom_point(alpha = .05, stroke = 0.2, size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
  scale_colour_manual(values = light_colours, guide="none", breaks = legend_order) +
  new_scale_colour() +
  stat_pointinterval(aes(colour = trait), size = 9, 
                     alpha = rep(c(0.95, 0.66), times = pred_n),
                     .width = c(0.66, 0.95)) + 
  scale_colour_manual(values = dark_colours, breaks = legend_order) +
  guides(fill = "none") +
  ggplot_settings + theme(legend.position = c(0.175, 0.85),
                          axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
                          axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
                          axis.title.x=element_text(size=rel(0.7), face = "bold"),
                          legend.text = element_text(size=rel(0.8), face = "bold"), 
                          legend.key.height = unit(0.2, 'cm'))

ggsave("Plots/Results/cpng_legend_all_model.tiff", dpi = 600, width = 8, height = 4)


# Plot with coloured width changing interval bars.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) + 
  geom_point(alpha = .05, stroke = 0.2, size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
  scale_colour_manual(values = light_colours, guide="none", breaks = legend_order) +
  new_scale_colour() +
  stat_pointinterval(aes(colour = trait), size = rep(c(9, 3), times = pred_n)) + 
  scale_colour_manual(values = dark_colours, breaks = legend_order) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  ggplot_settings + 
  theme(legend.position = c(0.175, 0.85),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.key.height = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        # legend.margin = margin(b=7, t = 2, l = 2, r = 4),
        # legend.background = element_rect(colour = 'black', 
        #                                  fill = 'white', 
        #                                  linetype='dashed', 
        #                                  linewidth = 0.1)
        )

ggsave("Plots/Results/greisser_legend_all_model.tiff", dpi = 600, width = 8, height = 4)


# With black point intervals that change width.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) + 
  geom_point(alpha = .05, stroke = 0.2, size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 3), times = pred_n)) + 
  scale_colour_manual(values = dark_colours, breaks = legend_order) + xlim(c(-6,6)) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  ggplot_settings + 
  theme(legend.position = c(0.175, 0.85),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.spacing.y = unit(0.2, 'cm'),
        legend.title = element_blank())

#ggsave("Plots/Results/legend_all_model_2.tiff", dpi = 600, width = 8, height = 4)

# Inlcuding the devo variable.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) + 
  geom_point(alpha = .05, stroke = 0.2, size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 3), times = pred_n)) + 
  scale_colour_manual(values = dark_colours, breaks = legend_order) + xlim(c(-6,6)) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  ggplot_settings + 
  theme(legend.position = c(0.85,0.175),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.spacing.y = unit(0.2, 'cm'),
        legend.title = element_blank())


# Inlcuding the devo variable and interactions.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) + 
  geom_point(alpha = .05, stroke = 0.2, size = 0.5, 
             position = position_jitter(seed = 1993, height = .35)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 3), times = pred_n)) + 
  scale_colour_manual(values = dark_colours, breaks = legend_order) + xlim(c(-6,6)) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  ggplot_settings + 
  theme(legend.position = c(0.175,0.65),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.spacing.y = unit(0.2, 'cm'),
        legend.title = element_blank())

ggsave("Plots/Results/greisser_model_2.tiff", dpi = 600, width = 8, height = 4)
###############################################################################
                    #### Quasi random plots #####


# Quasi plots all in navy.
ggplot(everything_data, aes(x = value, y = name)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5, colour = "#586994",
                   stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 0.1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = pred_n),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
      axis.text.x=element_text(size=rel(0.8), face = "bold"), 
      axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
      axis.title.x=element_text(size=rel(0.7), face = "bold"),
      legend.text = element_text(size=rel(0.6), face = "bold"), 
      legend.key.height = unit(0.6, 'cm'))

ggsave("Plots/Results/all_model_2.tiff", dpi = 600, width = 5, height = 5)


# Main publication design, but without facet (facet only works if groups are equal size.)
ggplot(npp_data, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
 # facet_grid(rows = vars(trait), scales = "free_y") +

  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.175, 0.2),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )

ggsave("Plots/Results/all_model_4.tiff", dpi = 600, width = 8, height = 4)



# Main publication design, but without facet and moving legend when altricial is included.
ggplot(everything_data, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) + #xlim(c(-6,6)) +
  # facet_grid(rows = vars(trait), scales = "free_y") +
  
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.15, 0.15),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )

ggsave("Plots/Results/everything_model_4.tiff", dpi = 600, width = 8, height = 7)


# With facet grid, but it makes each facet equal so it doesn't work properly.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
    scale_colour_manual(values = light_colours, breaks = legend_order) +
    guides(colour = 
             guide_legend(byrow = TRUE,
               override.aes = list(alpha = 1, shape = 15, size = 4))) +
    new_scale_colour() +
    stat_pointinterval(size = rep(c(9, 2), times = 10)) + 
    guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  #scale_y_discrete(labels = rev(axes_labs), expand = c(5,5)) +
  facet_grid(rows = vars(trait), scales = "free_y") +
  #scale_y_break(c("Diet","Territory")) +
  #scale_x_break(c(2)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.175, 0.9),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
         legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )
ggsave("Plots/Results/all_model_5.tiff", dpi = 600, width = 8, height = 4)

unique(all_plotdata$name)




library(ggbreak)


ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1) +
  scale_colour_manual(values = light_colours, guide="none", breaks = legend_order) +
  new_scale_colour() +
  stat_pointinterval(colour = "black", size = 5, 
                     alpha = rep(c(0.95, 0.8), times = 10)) +
  scale_colour_manual(values = dark_colours, breaks = legend_order) +
  scale_alpha_continuous(range = c(0.8,1)) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.25, 0.85),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.key.height = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'))



ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5, colour = "#586994",
                   stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 0.1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))


ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_beeswarm(size = 0.5, colour = "#586994", dodge.width = 1, priority = "density",
                stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))


ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5, colour = "#586994",
                   stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))






##############################################################################
                      #### Leftover code #####


ggplot_settings <- list(
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6), # Add lines.
  #coord_cartesian(xlim = c(-0.25, 0.75)),   # Axes.
  scale_y_discrete(labels = rev(axis_labels)),
  # scale_color_manual(values = rev(plot_palette),  # Colours.
  #                    guide = guide_legend(reverse = TRUE, 
  #                                         override.aes = list(linetype = 0))),
  # scale_fill_manual(values = rev(plot_palette), 
  #                    guide = "none"),
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL),
  theme_classic2(base_size = 12),   # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.7), face = "bold"), 
        axis.text.y=element_text(size=rel(0.7), face = "bold"),
        axis.title.x=element_text(size=rel(0.6), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm')))



benchmarking <- read.csv("Z:/home/sexual_selection/benchmarking.csv")
benchmarking %<>% filter(elapsed > 0)
plot(elapsed ~ grain_size, data = benchmarking)
plot(elapsed ~ jitter(thread_number), data = benchmarking)
boxplot(log(elapsed) ~ family, data = benchmarking)
boxplot(log(elapsed) ~ thread_number, data = benchmarking)
boxplot(log(elapsed) ~ grain_size, data = benchmarking)

boxplot(elapsed ~ thread_number, data = benchmarking)

hist(benchmarking$elapsed)

linear_test <- lm(elapsed ~ grain_size*thread_number, data = benchmarking)

summary(linear_test)

threads_16 <- benchmarking %>% filter(thread_number == 16)

plot(elapsed ~ grain_size, data = threads_16)

library(interactions)

interact_plot(linear_test, pred = thread_number, modx = grain_size)


# Most basic models.
model_1_all <- readRDS("Z:/home/sexual_selection/Results/Models/Test/sratio_16_100_1.rds")
cum_model <- readRDS("Z:/home/sexual_selection/Results/Models/Test/cumulative_16_100_1.rds")

plot(model_1_all)

mcmc_areas(cum_model, regex_pars = "^b_[a-z]", prob = 0.95)
