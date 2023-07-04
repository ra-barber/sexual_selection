###############################################################################
               ##### Summarise data before modelling #####
###############################################################################

# This script summarises data before subsequent modelling.

# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(effectsize)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

# Read in the life history traits.
full_data <- read.csv("Data/ecotraits_09_05_without_invalids.csv") %>% 
  clean_names()
high_full_data <- full_data %>% filter(sexual_certainty < 3)

# Read in model data.
model_data <- read.csv("Data/sexual_traits.csv")
high_data <- model_data %>% filter(sexual_certainty < 3)

# Read in a tree.
tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Read in clements data.
clements_data <- readxl::read_xlsx("Data/Clements_conversion/clements_ss_scores_12_06.xlsx")
clements_data$sexual_score %<>% as.numeric()
clements_data %<>% tidyr::drop_na(sexual_score)

###############################################################################
                    #### Transform data as usual ####

# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary"),
  devo_binary = relevel(as.factor(devo_mode_wang), ref = "altricial"),
  devo_mode = relevel(as.factor(devo_mode), ref = "altrical"),
  wang_edited = relevel(as.factor(wang_edited), ref = "altricial")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
  devo_bi_c = center_categorical(devo_binary),
  devo_ed_c = center_categorical(wang_edited),
  devo_mode_c = center_categorical(devo_mode)
)

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
  dens_z = standardize(dens_sqrt, two_sd = TRUE),
  fed_z = standardize(fed_sqrt, two_sd = TRUE),
  chick_z = standardize(chick_pc1, two_sd = TRUE),
  chick_sqrt_z = standardize(chick_sqrt, two_sd = TRUE)
)


# Set as factor, then re-level for appropriate reference group.
high_data %<>% mutate(
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary"),
  devo_binary = relevel(as.factor(devo_mode_wang), ref = "altricial"),
  devo_mode = relevel(as.factor(devo_mode), ref = "altrical"),
  wang_edited = relevel(as.factor(wang_edited), ref = "altricial")
)

# Center categorical predictors.
high_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
  devo_bi_c = center_categorical(devo_binary),
  devo_ed_c = center_categorical(wang_edited),
  devo_mode_c = center_categorical(devo_mode)
)

# Scale continuous predictors to two SD.
high_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
  dens_z = standardize(dens_sqrt, two_sd = TRUE),
  fed_z = standardize(fed_sqrt, two_sd = TRUE),
  chick_z = standardize(chick_pc1, two_sd = TRUE),
  chick_sqrt_z = standardize(chick_sqrt, two_sd = TRUE)
)

###############################################################################
                #### Count sexual selection scores ####

# Raw numbers of sexual selection scores.
full_data %>% count(sexual_score) 
high_full_data %>% count(sexual_score)

# Calculate percentages.
counts <- full_data %>% count(sexual_score) 
counts[1,2]/9989     
sum(counts[2:5,2])/9989   

counts[2,2]/9989
counts[3,2]/9989
counts[4,2]/9989
counts[5,2]/9989




################################################################################
    #### Calculate phylogenetic signal sexual selection scores #####

model_data$tree_tip <- model_data$birdtree_name %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip

sex_data <- model_data %>% select(tree_tip, sexual_binary)

# Create a comparative object.
sex_comp <- comparative.data(tree, sex_data, names.col = "tree_tip")

tic()
# Run the phylo D signal analysis for binary traits.
phylo_signal_d <- phylo.d(data = sex_comp, binvar = sexual_binary)
toc()
# All species took 4978 seconds.

# Return the phylo D analysis results.
phylo_signal_d

saveRDS(phylo_signal_d, "Results/Tables/phylo_d.rds")

# Plot the result against expected models.
plot(phylo_signal_d)

###############################################################################
             #### Plot sexual selection scores barchat ####

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

ss_barplot <- function(dataset = full_data){
  dataset %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) +
    geom_bar() +  scale_fill_manual(values = pal) + 
    theme_classic(base_size = 20) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5))
}

fruit_data <- full_data %>% filter(trophic_niche == "Frugivore")
invert_data <- full_data %>% filter(trophic_niche == "Invertivore")

jetz_y <- full_data %>% count(sexual_score) %>% pull(n) %>% max()
clements_y <- clements_data %>% count(sexual_score) %>% pull(n) %>% max()
fruit_y <- fruit_data %>% count(sexual_score) %>% pull(n) %>% max()
invert_y <- invert_data %>% count(sexual_score) %>% pull(n) %>% max()

# Create basic bar plots of varied zoom.
jetz_plot <- full_data %>% ss_barplot() + scale_y_continuous(limits = c(0,jetz_y*1.1)) +
  ylab("Species count")  + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(vjust = 1.5))

clements_plot <- clements_data %>% ss_barplot() + scale_y_continuous(limits = c(0,clements_y*1.1), breaks = c(2500,5000,7500)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank())

fruit_plot <- full_data %>% filter(trophic_niche == "Frugivore") %>%  
  ss_barplot() + xlab("Sexual selection") + ylab("Species count") +
  scale_y_continuous(limits = c(0, fruit_y*1.1), breaks = c(0, 250, 500, 700), labels = c("0", "  250", "  500", "  750")) + 
  theme(axis.title.y = element_text(vjust = 1.5))

invert_plot <- full_data %>% filter(trophic_niche == "Invertivore") %>% 
  ss_barplot() + scale_y_continuous(limits = c(0,invert_y*1.1)) +
  xlab("Sexual selection") + theme(axis.title.y = element_blank())

# Arrange the basic plots together.
combined_barplots <- ggarrange(jetz_plot + annotate("text", x = -0.3, y = jetz_y*1.1, label = "a", size = 10, fontface = 2), 
                               clements_plot  + annotate("text", x = -0.3, y = clements_y*1.1, label = "b", size = 10, fontface = 2),
                               fruit_plot + annotate("text", x = -0.3, y = fruit_y*1.1, label = "c", size = 10, fontface = 2), 
                               invert_plot  +  annotate("text", x = -0.3, y = invert_y*1.1, label = "d", size = 10, fontface = 2), ncol = 2, nrow = 2)

ggsave("Plots/Data/simple_scores.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/simple_scores.tiff", width = 8, height = 8)

# Create zoomed facets for easier displays.
jetz_facet <- full_data %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("") + ylab("") +
  coord_cartesian(ylim = c(0,750)) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  theme_classic(base_size = 16) + 
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), line = element_line(linewidth = 0.5))

clements_factet <- clements_data %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("") + ylab("") +
  coord_cartesian(ylim = c(0,750)) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  theme_classic(base_size = 16) + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), line = element_line(linewidth = 0.5))

fruit_factet <- full_data %>% filter(trophic_niche == "Frugivore") %>%  
  ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("") + ylab("") +
  coord_cartesian(ylim = c(0,100)) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  theme_classic(base_size = 16) + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), line = element_line(linewidth = 0.5))

invert_factet <- full_data %>% filter(trophic_niche == "Invertivore") %>% 
  ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("") + ylab("") +
  coord_cartesian(ylim = c(0,250)) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  theme_classic(base_size = 16) + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), line = element_line(linewidth = 0.5))


jetz_both <- jetz_plot + annotation_custom(ggplotGrob(jetz_facet), xmin = 1, xmax = 4.5, 
                       ymin = jetz_y*0.5, ymax = jetz_y*1.1) +
  annotate("text", x = -0.3, y = jetz_y*1.1, label = "a", size = 10, fontface = 2)

clements_both <- clements_plot + annotation_custom(ggplotGrob(clements_factet), xmin = 1, xmax = 4.5, 
                                                   ymin = clements_y*0.5, ymax = clements_y*1.1) +
  annotate("text", x = -0.3, y = clements_y*1.1, label = "b", size = 10, fontface = 2)

fruit_both <- fruit_plot + annotation_custom(ggplotGrob(fruit_factet), xmin = 1, xmax = 4.5, 
                                             ymin = fruit_y*0.5, ymax = fruit_y*1.1) +
  annotate("text", x = -0.3, y = fruit_y*1.1, label = "c", size = 10, fontface = 2)

invert_both <- invert_plot + annotation_custom(ggplotGrob(invert_factet), xmin = 1, xmax = 4.5, 
                                               ymin = invert_y*0.5, ymax = invert_y*1.1) +
  annotate("text", x = -0.3, y = invert_y*1.1, label = "d", size = 10, fontface = 2)

# Arrange new plots with facets.
ggarrange(jetz_both, clements_both, fruit_both, invert_both, 
          ncol = 2, nrow = 2, widths = c(1.1,1), heights = c(1,1.2))

ggsave("Plots/Data/sexual_selection_scores.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/sexual_selection_scores.tiff", width = 8, height = 8)



###############################################################################
               #### Do it with splits in the y axis ####

library(ggbreak)

jetz_plot <- full_data %>% ss_barplot() + scale_y_continuous(limits = c(0,jetz_y*1.1)) +
  ylab("Species count")  + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(vjust = 1.5)) +
  scale_y_break(c(jetz_y*0.1, jetz_y*0.5), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

clements_plot <- clements_data %>% ss_barplot() + scale_y_continuous(limits = c(0,clements_y*1.1), breaks = c(0, 200, 400, 600, 800)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank()) +
scale_y_break(c(clements_y*0.1, clements_y*0.5), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())


fruit_plot <- full_data %>% filter(trophic_niche == "Frugivore") %>%  
  ss_barplot() + xlab("Sexual selection") + ylab("Species count") +
  scale_y_continuous(limits = c(0, fruit_y*1.1), breaks = c(0, 25, 50, 75), labels = c("0", "  25", "  50", "  75")) + 
  theme(axis.title.y = element_text(vjust = 1.5)) +
scale_y_break(c(fruit_y*0.1, fruit_y*0.5), scales = "free", ticklabels = c(500, 600, 700, 800)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())
fruit_plot

invert_plot <- full_data %>% filter(trophic_niche == "Invertivore") %>% 
  ss_barplot() + scale_y_continuous(limits = c(0,invert_y*1.1)) +
  xlab("Sexual selection") + theme(axis.title.y = element_blank()) +
scale_y_break(c(invert_y*0.1, invert_y*0.49), scales = "free", ticklabels = c(2000, 3000, 4000)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

invert_plot

# Arrange the basic plots together.
combined_barplots <- ggarrange(print(jetz_plot + annotate("text", x = -0.3, y = jetz_y*1.1, label = "a", size = 10, fontface = 2)), 
                                     print(clements_plot  + annotate("text", x = -0.3, y = clements_y*1.1, label = "b", size = 10, fontface = 2)),
                                           print(fruit_plot + annotate("text", x = -0.3, y = fruit_y*1.1, label = "c", size = 10, fontface = 2)), 
                                                 print(invert_plot  +  annotate("text", x = -0.3, y = invert_y*1.1, label = "d", size = 10, fontface = 2)), 
                               ncol = 2, nrow = 2)


ggsave("Plots/Data/example_broken_axis.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/sexual_selection_scores.tiff", width = 8, height = 8)


###############################################################################
                    #### Do it joes way ####


jetz_plot <- full_data %>% ss_barplot() + scale_y_continuous(limits = c(0,jetz_y*1.1)) +
  ylab("Species count")  + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(vjust = 1.5)) +
  scale_y_break(c(1000, 5000), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + ylim(0,10000) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

clements_plot <- clements_data %>% ss_barplot() + scale_y_continuous(limits = c(0,clements_y*1.1), breaks = c(0, 250, 500, 750, 1000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank()) +
  scale_y_break(c(1000, 500), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + ylim(0,10000) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())


fruit_plot <- fruit_data %>%  
  ss_barplot() + xlab("Sexual selection") + ylab("Species count") +
  scale_y_continuous(limits = c(0, 5500)) + #, breaks = c(0, 25, 50, 75), labels = c("0", "  25", "  50", "  75")) + 
  theme(axis.title.y = element_text(vjust = 1.5)) +
  #scale_y_break(c(250, 500), scales = "free", ticklabels = c(1000, 2000, 3000, 4000)) + 
  scale_y_break(c(250, 500), scales = "free", ticklabels = c(500, 1500, 2500, 3500, 4500)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())


invert_plot <- full_data %>% filter(trophic_niche == "Invertivore") %>% 
  ss_barplot() + scale_y_continuous(limits = c(0,5500)) +
  xlab("Sexual selection") + theme(axis.title.y = element_blank()) +
  scale_y_break(c(250, 500), scales = "free", ticklabels = c(500, 1500, 2500, 3500, 4500)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())


# Labels are 95% and 80% low than y lim.
0.95*(5500-500) + 500
0.8*(5500-500) + 500


# Add the annotations. 
jetz_anno <- jetz_plot + annotate("text", x = -0.3, y = 9750, label = "a", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 9000, label = "Taxonomy 1\nn = 9989", size = 7)

clem_anno <- clements_plot  + annotate("text", x = -0.3, y = 9750, label = "b", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 9000, label = "Taxonomy 2\nn = 10671", size = 7)

fruit_anno <- fruit_plot + annotate("text", x = -0.3, y = 5250, label = "c", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Frugivores\nn = 1030", size = 7)

invert_anno <- invert_plot  +  annotate("text", x = -0.3, y = 5250, label = "d", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Invertivores\nn = 4780", size = 7)

# Amaurospiza carrizalensis 
ggarrange(print(jetz_anno), 
          print(clem_anno),
          print(fruit_anno), 
          print(invert_anno), 
          ncol = 2, nrow = 2, widths = c(1.1,1), heights = c(1, 1.1))

# Export.
ggsave("Plots/Data/extended_data_figure_1.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/extended_data_figure_1.tiff", width = 8, height = 8)

?scale_y_break

5000 - 10000
500 - 5000
4750/5000
4000/5000



4500*0.8















jetz_plot <- full_data %>% ss_barplot() + scale_y_log10()
clements_plot <- clements_data %>% ss_barplot() + scale_y_log10()
fruit_plot <- fruit_data %>% ss_barplot() + scale_y_log10()
invert_plot <- invert_data %>% ss_barplot() + scale_y_log10()


ggarrange(print(jetz_plot + annotate("text", x = -0.3, y = jetz_y*1.1, label = "a", size = 10, fontface = 2)), 
          print(clements_plot  + annotate("text", x = -0.3, y = clements_y*1.1, label = "b", size = 10, fontface = 2)),
          print(fruit_plot + annotate("text", x = -0.3, y = fruit_y*1.1, label = "c", size = 10, fontface = 2)), 
          print(invert_plot  +  annotate("text", x = -0.3, y = invert_y*1.1, label = "d", size = 10, fontface = 2)), 
          ncol = 2, nrow = 2)

ggsave("Plots/Data/log_scale.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/log_scale.tiff", width = 8, height = 8)





full_data %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("Sexual selection") + ylab("Species count") +
  ylim(0, 8500) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  scale_y_break(c(1000, 5000), scales = "free", ticklabels = c(5000, 6000, 7000, 8000)) +
  theme_classic() + theme(legend.position = "none",axis.text.y.right = element_blank(),
                                                          axis.ticks.y.right  = element_blank(),
                                                          axis.line.y.right  = element_blank())
                                                          
library(ggforce)


full_data %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) + 
  xlab("Sexual selection") + ylab("Species count") +
  ylim(0, 8500) +
  geom_bar() +  scale_fill_manual(values = pal) + 
  theme_classic() + theme(legend.position = "none") +
  facet_zoom(ylim = c(0, 1000), zoom.size = 1, show.area = FALSE, 
             shrink = FALSE)



################################################################################

model_data %>% group_by(cert_reverse) %>% summarise(mean_ss = mean(sexual_score))

model_data %>% ggplot(aes(x = as.factor(cert_reverse), y = sexual_score)) + geom_boxplot() +
  xlab("Sexual selection") + ylab("Species count") +
  theme_classic() + theme(legend.position = "none")


model_data %>% filter(cert_reverse == 1 & sexual_score > 0)


# Function for looking at sexual selection vs predictor means and standard error.
sex_meanplot <- function(predictor){
  grouped_data <- model_data %>% 
    group_by(!!! syms(predictor)) %>% 
    summarise(trait = first(!!! syms(predictor)),
              sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  ggplot(grouped_data, aes(x = trait, y = sex_mean)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), width =0.5) + 
    geom_point(size =2) + ylab("Sexual selection") + xlab(predictor) + theme_classic()
}


sex_meanplot("cert_reverse") + xlab("Data certainty") +  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))

# Export.
ggsave("Plots/Data/figure_s2.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/figure_s2.png", width = 8, height = 8)


###############################################################################
                    #### Check raw correlations ####

continuous_variables <- model_data %>% dplyr::select(temp_seasonality_z, npp_z, chick_sqrt_z)
library(GGally)

ggcorr(continuous_variables, label = TRUE, method =c("everything","pearson"), label_round = 3)

lat_variables <- model_data %>% dplyr::select(centroid_z, temp_seasonality_z, 
                                              terr_bi_c, migration_bi_c, trophic_level_c)

 
ggcorr(lat_variables, label = TRUE, method =c("everything","spearman"), label_round = 3)

###############################################################################
                       #### Check  VIF values ####

# # Uncentered.
# model_formula <- "sexual_score ~ territory_binary + migration_binary + 
#     trophic_binary + chick_z + trophic_binary*temp_seasonality_z + 
#     trophic_binary*territory_binary + trophic_binary*npp_z + 
#     chick_z*temp_seasonality_z + chick_z*territory_binary + chick_z*npp_z"
# 
# # Centered.
# model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
#     trophic_level_c + chick_sqrt_z + trophic_level_c*temp_seasonality_z + 
#     trophic_level_c*terr_bi_c + trophic_level_c*npp_z + 
#     chick_sqrt_z*temp_seasonality_z + chick_sqrt_z*terr_bi_c + chick_sqrt_z*npp_z"
# 
# # No NPP
# model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
#     trophic_level_c + chick_z + trophic_level_c*temp_seasonality_z + 
#     trophic_level_c*terr_bi_c + 
#     chick_z*temp_seasonality_z + chick_z*terr_bi_c"
# 
# # No NPP centerered.
# model_formula <- "sexual_score ~ territory_binary + migration_binary + 
#     trophic_binary + chick_z + trophic_binary*temp_seasonality_z + 
#     trophic_binary*territory_binary + 
#     chick_z*temp_seasonality_z + chick_z*territory_binary"
# 
# # With chick sqrt.
# model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
#     trophic_level_c + trophic_level_c*temp_seasonality_z  + 
#     chick_sqrt_z*temp_seasonality_z"
# 
# # Simple.
# model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
#     trophic_level_c + chick_z + trophic_level_c*temp_seasonality_z"
# 
# 
# 
# model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
#     trophic_level_c + chick_sqrt_z + trophic_level_c*temp_seasonality_z + 
#     trophic_level_c*terr_bi_c + trophic_level_c*npp_z + 
#     chick_sqrt_z*temp_seasonality_z + chick_sqrt_z*terr_bi_c + chick_sqrt_z*npp_z"


# Model formula for VIF values.
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z + 
    trophic_level_c*terr_bi_c"


library(car)
test_model <- lm(model_formula, data = model_data)

vif(test_model) %>% round(digits = 2)
vif(test_model, type = "predictor")
sqrt(3)


test_model <- lm(model_formula, data = high_data)

vif(test_model) %>% round(digits = 2)

summary(test_model)

model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)
row.names(model_data) <- model_data$tree_tip

phylolm(model_formula, model_data, tree) %>% summary()


###############################################################################
                #### Spatial autocorrelation ####

library(spdep)
library(rgeos)
library(sp)
library(sf)

# Select data.
just_sex_long_lat <- full_data %>% select(sexual_score, centroid_latitude, centroid_longitude)

# Create binary trait.
just_sex_long_lat$sexual_binary <- 0
just_sex_long_lat$sexual_binary[just_sex_long_lat$sexual_score > 0] <- 1
just_sex_long_lat %<>% na.omit()

# Create a SpatialPoints object
long_lat <- data.frame(lat = just_sex_long_lat$centroid_latitude, long = just_sex_long_lat$centroid_longitude)
data_sf <- st_as_sf(long_lat, coords = c("long", "lat"),
                    # Change to your CRS
                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
coords <- st_coordinates(data_sf)
sp_points <- SpatialPoints(coords)

# Define the number of nearest neighbours
k <- 5

# Create a spatial weights object based on k nearest neighbors
nb <- knn2nb(knearneigh(sp_points, k = k, longlat = TRUE))

# Convert the spatial weights object to listw format
w <- nb2listw(nb, style="B")

# Generate random attribute values
x <- just_sex_long_lat$sexual_binary %>% as.factor()
x <- just_sex_long_lat$sexual_score %>% as.factor()
just_sex_long_lat
# Perform Monte Carlo simulation
join_count <- joincount.mc(listw = w, fx = x, nsim = 999, alternative = "two.sided")

# Simulated join count statistics
join_count

joincount.test(x, w)

# Suggests significant autocorrelation.

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####

model_data %>% count(trophic_binary, territory)


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


