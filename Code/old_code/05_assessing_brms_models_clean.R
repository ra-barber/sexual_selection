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

# Read in the new files.
feb_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/Feb/"), full.names = T, include.dirs = FALSE, recursive = FALSE)
march_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/March/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

simple_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/Simple/"), full.names = T, include.dirs = FALSE, recursive = FALSE)


# Function to extract and combine brms models.
combine_brms <- function(model_pattern){
  pattern_matched <- march_files %>% str_subset(pattern = model_pattern)
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



# Read in centered models first.
npp_cen_all_models <- combine_brms("/npp_centered_all")
nonpp_cen_all_models <- combine_brms("no_npp_centered_all")
no_inter_cen_all_models <- combine_brms("no_inter_centered_all")
npp_cen_high_models <- combine_brms("/npp_centered_high")
nonpp_cen_high_models <- combine_brms("no_npp_centered_high")
no_inter_cen_high_models <- combine_brms("no_inter_centered_high")

# Uncentered models.
npp_uncen_all_models <- combine_brms("/npp_uncentered_all")
nonpp_uncen_all_models <- combine_brms("no_npp_uncentered_all")
no_inter_uncen_all_models <- combine_brms("no_inter_uncentered_all")
npp_uncen_high_models <- combine_brms("/npp_uncentered_high")
nonpp_uncen_high_models <- combine_brms("no_npp_uncentered_high")
no_inter_uncen_high_models <- combine_brms("no_inter_uncentered_high")


# Logistic models.
logistic_all_models <- combine_brms("/logistic_centered_all")

# Logistic chick sqrt.
logistic_all_models <- combine_brms("/logistic_chicksqrt_centered_all")


chicksqrt_all_models <- combine_brms("/chicksqrt_uncentered_all")
equidistant_high_models <- combine_brms("/equidistant_centered_high")
equidistant_all_models <- combine_brms("/equidistant_uncentered_all")

no_terr_all_models <- combine_brms("/no_terr_centered_all")
no_terr_high_models <- combine_brms("/no_terr_centered_high")
equidistant_centered_high

# Read in equidistant all models.
equi_all_models <- combine_brms("/equidistant_centered_all")
norm_all_models <- combine_brms("/normal_centered_all")

terr_equi_all_models <- combine_brms("/terr_equidistant_centered_all")
terr_norm_all_models <- combine_brms("/terr_normal_centered_all")

# High certainty models.
terr_equi_high_models <- combine_brms("/terr_equidistant_centered_high")




tic()
conditional_effects(no_inter_uncen_all_models, "territory_binary", categorical = TRUE, resolution = 100)
toc()
# Main analysis for development.

tic()
npp_plots <- plot(conditional_effects(npp_cen_high_models, categorical = FALSE, conditions = conditions), ask = FALSE)
toc()

tic()
equi_plots <- plot(conditional_effects(equidistant_high_models, categorical = FALSE), ask = FALSE)
toc()

tic()
cat_npp_plots <- plot(conditional_effects(npp_cen_high_models, categorical = FALSE), ask = FALSE)
toc()

tic()
cat_equi_plots <- plot(conditional_effects(equidistant_high_models, categorical = TRUE), ask = FALSE)
toc()
for (x in 1:12){

print(ggarrange(npp_plots[[x]], equi_plots[[x]]))
  
}



length(equi_plots)

?conditional_effects
?as.factor()

make_conditions(trophic)
conditions = data.frame(trophic_level_c = unique(model_data$trophic_level_c))
unique(model_data$trophic_level_c)
###############################################################################
                    #### Model Diagnostics #####


# Look at the summaries.
summary(nonpp_cen_all_models)

summary(npp_uncen_all_models)
summary(nonpp_uncen_all_models)
summary(no_inter_uncen_all_models)
summary(logistic_all_models)


single_model <- readRDS("Z:/home/sexual_selection/Results/Models/Feb/equidistant_centered_high_1.rds")
single_model_2 <- readRDS("Z:/home/sexual_selection/Results/Models/Feb/npp_centered_high_1.rds")
single_model_2 <- readRDS("Z:/home/sexual_selection/Results/Models/Feb/npp_centered_all_10.rds")

summary(single_model)

# # Plot the model.
plot(brms_models)
plot(hi_brms_models)

# Check the posterior distribution.
pp_check(single_model)
pp_check(single_model_2)




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
brms_forest(npp_cen_high_models)
brms_forest(nonpp_cen_high_models)
brms_forest(no_inter_cen_high_models)
brms_forest(npp_uncen_all_models)
brms_forest(nonpp_uncen_all_models)
brms_forest(no_inter_uncen_all_models)
brms_forest(npp_uncen_high_models)
brms_forest(nonpp_uncen_high_models)
brms_forest(no_inter_uncen_high_models)


brms_forest(logistic_all_models)
brms_forest(chicksqrt_all_models)
brms_forest(equidistant_high_models)
brms_forest(no_terr_all_models)
brms_forest(no_terr_high_models)
brms_forest(equidistant_all_models)

brms_forest(equi_all_models)
brms_forest(norm_all_models)
brms_forest(terr_equi_all_models)
brms_forest(terr_norm_all_models)




###############################################################################
                          #### P mapping #####


# Check p-values.
brms_pmap <- function(model){
  p_map(model, regex_pars = "^b_[a-z]")
}

# Feb models.
brms_pmap(npp_cen_all_models)
brms_pmap(nonpp_cen_all_models)
brms_pmap(no_inter_cen_all_models)
brms_pmap(npp_cen_high_models)
brms_pmap(nonpp_cen_high_models)
brms_pmap(no_inter_cen_high_models)
brms_pmap(npp_uncen_all_models)
brms_pmap(nonpp_uncen_all_models)
brms_pmap(no_inter_uncen_all_models)
brms_pmap(npp_uncen_high_models)
brms_pmap(nonpp_uncen_high_models)
brms_pmap(no_inter_uncen_high_models)

brms_pmap(logistic_all_models)
brms_pmap(equidistant_high_models)
brms_pmap(npp_cen_high_models)
brms_pmap()
brms_pmap()


###############################################################################
                      #### Lambda values #####


# Look at lambda values.
hyp <- "sd_tree_tip__Intercept^2 / (sd_tree_tip__Intercept^2 + disc^2) = 0"

(all_lambda <- hypothesis(single_model, hyp, class = NULL))
plot(all_lambda)

###############################################################################
                 #### Making the plot names #####


# No interaction names.
no_inter_names <- c( "Territoriality", "Migration", "(PC) Primary\nconsumer", "Chick maturity",
                     "Seasonality", "Productivity", 
                     "PC x \nSeasonality", "PC x \nTerritoriality", "PC x \nProductivity")

# Npp models.
npp_names <- c("Territoriality", "Migration", "Primary\nconsumer", "Chick maturity",
               "Seasonality", "Productivity", 
               "PC x \nSeasonality", "PC x \nTerritoriality", "PC x \nProductivity",
               "CM x \nSeasonality", "CM x \nTerritoriality", "CM x \nProductivity")

# No Npp models.
nonpp_names <- c("Territoriality", "Migration", "Primary\nconsumer", "Chick maturity",
               "Seasonality", 
               "PC x \nSeasonality", "PC x \nTerritoriality",
               "CM x \nSeasonality", "CM x \nTerritoriality")

# Simple all models.
simple_names <- c("Territoriality", "Migration", "Primary\nconsumer", "Chick maturity",
                 "Seasonality", "PC x \nSeasonality")

# Simple terr models.
terr_simple_names <- c("Territoriality", "Migration", "Primary\nconsumer", "Chick maturity",
                  "Seasonality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nTerritoriality")

###############################################################################
            #### Extract the draws for each model #####

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

# Simple models first.
simple_all_plotdata <- extract_draws(equi_all_models, simple_names)
terr_simple_all_plotdata <- extract_draws(terr_equi_all_models, terr_simple_names)

# Npp and no npp complex models.
npp_all_plotdata <- extract_draws(npp_cen_all_models, npp_names)
nonpp_all_plotdata <- extract_draws(nonpp_cen_all_models, nonpp_names)

# High certainty data.
npp_high_plotdata <- extract_draws(npp_cen_high_models, npp_names)
terr_simple_high_plotdata <- extract_draws(terr_equi_high_models, terr_simple_names)

terr_simple_all_plotdata %<>% filter(name %in% c("Territoriality", "Migration", "Primary\nconsumer",
                                                  "Seasonality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nTerritoriality"))
terr_simple_high_plotdata %<>% filter(name %in% c("Territoriality", "Migration", "Primary\nconsumer",
                       "Seasonality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nTerritoriality"))

# Presentation figure.
all_plotdata <- extract_draws(no_inter_cen_all_models, model_colnames)
all_plotdata %<>% filter(name %in%  c("Territoriality", "Migration", "(PC) Primary\nconsumer", "Chick maturity",
                                      "Seasonality", "PC x \nSeasonality"))



###############################################################################
                #### Add axis label and legend information #####


# Vector of predictors in reverse order of plotting.
simple_preds <- simple_names
simple_terr_preds <- terr_simple_names[-4]

simple_terr_preds <- simple_terr_preds[c(4,2,3,1,5,6)]

npp_preds <- npp_names[c(1:7, 9, 8, 10, 12, 11)]
nonpp_preds <- nonpp_names

# Set data as factors for plotting order.
#all_plotdata$name %<>% factor(levels = rev(rev_preds))

# Simple models first.
simple_all_plotdata$name %<>% factor(levels = rev(simple_preds))
terr_simple_all_plotdata$name %<>% factor(levels = rev(simple_terr_preds))

# Npp and no npp complex models.
npp_all_plotdata$name  %<>% factor(levels = rev(npp_preds))
nonpp_all_plotdata$name  %<>% factor(levels = rev(nonpp_preds))

# High certainty data.
npp_high_plotdata$name  %<>% factor(levels = rev(npp_preds))
terr_simple_high_plotdata$name  %<>% factor(levels = rev(simple_terr_preds))



lifehistory <- c("Territoriality",  "Migration",  "Primary\nconsumer", "(Pr) Precocial", "Chick maturity")
resource <-   c("Seasonality", "PC x \nSeasonality")
resource <-   c("Seasonality", "Productivity")
interaction <-  c("Primary\nconsumer\n x \nTerritoriality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nProductivity")
#devo_interaction <- c("Pr x \nTerritoriality", "Pr x \nSeasonality", "Pr x \nProductivity")
devo_interaction <- c("DM x PC", "DM x \ngeneration length", "DM x Migration", "DM x \nseasonality",
                      "DM x \nproductivity", "DM x \nterritoriality","CM x \nSeasonality",
                      "CM x \nProductivity", "CM x \nTerritoriality")

# Function for adding colour to legends.
add_legend_info <- function(plotdata){
  plotdata$trait <- NA
  plotdata$trait[plotdata$name %in% lifehistory] <- "Life history"
  plotdata$trait[plotdata$name %in% resource] <- "Environmental"
  plotdata$trait[plotdata$name %in% interaction] <- "Diet interaction"
  plotdata$trait[plotdata$name %in% devo_interaction] <- "Devlopment interaction"
  plotdata
}

# Add the legend info.
simple_all_plotdata %<>% add_legend_info()
terr_simple_all_plotdata %<>% add_legend_info()
npp_all_plotdata %<>% add_legend_info()
nonpp_all_plotdata %<>% add_legend_info()

# High certainty.
npp_high_plotdata %<>% add_legend_info()
terr_simple_high_plotdata %<>% add_legend_info()
# Y axis labels for ordering plot.
#axes_labs <- rev_preds


# Group plot settings.
ggplot_settings <- list(
  # geom_point(alpha = .1, stroke = 0.2, colour = "#586994", size = 0.5, 
  #            position = position_jitter(seed = 1993, height = .35)),
  # stat_pointinterval(size = 9, alpha = rep(c(0.95, 0.6), times = 9), 
  #                    colour = "navy"),
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6), # Add lines.
  #scale_y_discrete(labels = rev(axes_labs)),
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL),
  theme_classic2(base_size = 12),   # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm')))



###############################################################################
                         #### Colour pal #####

pred_n <- 6

library(ggnewscale)

light_colours <- c( "#77AD78","#7494EA", "#AA4465")
light_colours <- c( "#77AD78","#7494EA", "#C98986")



dark_colours <- c("#136F63", "#05299E", "#8D0801")
dark_colours <- c("#214F4B", "#05299E", "#8D0801")
legend_order <- c("Life history", "Environmental", "Trophic interaction")

# cOLOURS WITH PRECOCIAL INTERACTION.
light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")
light_colours <- c("#7494EA",  "#77AD78", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
dark_colours <- c( "#05299E","#214F4B", "#8D0801","#A88A05")

legend_order <- c( "Environmental", "Life history","Diet interaction", "Devlopment interaction")

###############################################################################
                    #### Quasi random plots #####

gg_quasi_forst <- function(plotdata, pred_n = 6, legend_pos = c(0.15, 0.15), col_pal = light_colours, alpha_val = 0.1){
  # Main publication design for simple models.
  ggplot(plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
    geom_quasirandom(bandwidth = 0.5, size = 0.5,
                     stroke = 0.1, groupOnX = FALSE, alpha = alpha_val, position = position_dodge(2)) +
    scale_colour_manual(values = col_pal, breaks = legend_order) +
    guides(colour = guide_legend(byrow = TRUE, override.aes = list(alpha = 1, shape = 15, size = 4))) +
    new_scale_colour() +
    stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
    guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
    geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
    #scale_y_discrete(labels = rev(axes_labs)) +
    # facet_grid(rows = vars(trait), scales = "free_y") +
    labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
    theme_classic2(base_size = 12) +  # Theme.
    theme(legend.position = legend_pos,
          axis.text.x=element_text(size=rel(1), face = "bold", colour = "black"), 
          axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
          axis.title.x=element_text(size=rel(0.7), face = "bold"),
          legend.text = element_text(size=rel(0.8), face = "bold"), 
          #legend.key.size = unit(0, 'cm'),
          legend.key.height = unit(0.2, 'cm'),
          # legend.key.width = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title = element_blank(),
          strip.text.y = element_blank())
}

gg_quasi_forst(simple_all_plotdata, legend_pos = c(0.85, 0.9), alpha_val = 0.2)
ggsave("Plots/Results/simple_model.tiff", dpi = 600, width = 8, height = 4)

gg_quasi_forst(terr_simple_all_plotdata, pred_n = 6, legend_pos = c(0.125, 0.125), alpha_val = 0.5)
ggsave("Plots/Results/terr_simple_model.tiff", dpi = 600, width = 8, height = 4)  # Changed for presentation size.
ggsave("Plots/Results/figure_5_model.tiff", dpi = 1000, width = 8, height = 5)  # Changed for presentation size.

#gg_quasi_forst(npp_all_plotdata, pred_n = 12, legend_pos = c(0.15, 0.1))
#ggsave("Plots/Results/npp_model.tiff", dpi = 600, width = 8, height = 6)

#gg_quasi_forst(nonpp_all_plotdata, pred_n = 9, legend_pos = c(0.15, 0.125))
#ggsave("Plots/Results/nonpp_model.tiff", dpi = 600, width = 8, height = 5)


# High certainty data.
#gg_quasi_forst(npp_high_plotdata, pred_n = 12, legend_pos = c(0.125, 0.125))
#ggsave("Plots/Results/npp_high_model.tiff", dpi = 600, width = 8, height = 6)

gg_quasi_forst(terr_simple_high_plotdata, pred_n = 6, legend_pos = c(0.125, 0.125), alpha_val = 0.2)
ggsave("Plots/Results/terr_simple_high_model.tiff", dpi = 1000, width = 8, height = 4)


ggsave("Plots/Results/figure_S2_tall.tiff", dpi = 1000, width = 8, height = 8)


# Main publication design for simple models.
ggplot(simple_all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours[c(2,3, 4)], breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  #scale_y_discrete(labels = rev(axes_labs)) +
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
        strip.text.y = element_blank())

ggsave("Plots/Results/simple_model.tiff", dpi = 600, width = 8, height = 4)

# Main publication design, but without facet (facet only works if groups are equal size.)
ggplot(simple_all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours[c(2,3)], breaks = legend_order) +
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
    stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
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
