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

# Read in the cleaner data for this.

clean_data <- read.csv("Data/sexual_selection_cleaned_01_08.csv") %>% clean_names()  # We have since cleaned this data even more.


clean_data %>% group_by(data_certainty) %>% summarise(mean_ss = mean(sexual_selection))

grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(trait = first(sexual_selection),
            sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

clean_data %>% ggplot(aes(x = as.factor(data_certainty), y = sexual_selection, colour = as.factor(sexual_selection))) +
  xlab("Data certainty") + ylab("Sexual selection") + geom_jitter(alpha = 0.5, position = position_jitternormal(sd_y = 0.05)) +
  theme_classic() + theme(legend.position = "none") +  scale_colour_manual(values = pal) +
  geom_errorbar(data = grouped_data, inherit.aes = FALSE,
                aes(x = data_certainty, ymin = sex_mean - sex_se*1.96, ymax = sex_mean + sex_se*1.96), colour = "black", width = 0.2, linewidth = 1)




library(ggforce)
library(ggdist)
library(ggbeeswarm)


# Function for looking at sexual selection vs predictor means and standard error.
sex_meanplot <- function(predictor){
  grouped_data <- model_data %>% 
    group_by(!!! syms(predictor)) %>% 
    summarise(trait = first(!!! syms(predictor)),
              sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  ggplot(grouped_data, aes(x = trait, y = sex_mean)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se),  width = 0.2, linewidth = 0.5) + 
    geom_point(size =3) + ylab("Sexual selection") + xlab(predictor) + theme_classic()
}


whisker_plot <- sex_meanplot("cert_reverse") + xlab("Data certainty") +  theme_classic(base_size = 30) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) +
  annotate("text", x = 1, y = .15, label = "35", size = 5) +
  annotate("text", x = 2, y = .235, label = "2253", size = 5) +
  annotate("text", x = 3, y = .295, label = "4909", size = 5) +
  annotate("text", x = 4, y = .895, label = "2792", size = 5) 

# Export.
ggsave("Plots/Data/figure_s2.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/figure_s2.png", width = 8, height = 8)


###############################################################################
             ##### Proportion plot #####


grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

merged_data <- clean_data %>% 
  left_join(grouped_data, by = "data_certainty")
merged_data$sex_fact <- factor(merged_data$sexual_selection, levels = 4:0)

prop_plot <- ggplot(merged_data, aes(x = data_certainty, y = sex_mean, col = sex_fact, fill = as.factor(sex_fact))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") + xlab("Data certainty") +
  theme_classic(base_size = 30) + 
  scale_y_continuous(limits = c(0,1),expand = expansion(add = c(0.02, 0.09)))+
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = rel(0.9)))
prop_plot

# Function for looking at sexual selection vs predictor means and standard error.
sex_meanplot <- function(predictor){
  grouped_data <- model_data %>% 
    group_by(!!! syms(predictor)) %>% 
    summarise(trait = first(!!! syms(predictor)),
              sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  ggplot(grouped_data, aes(x = trait, y = sex_mean)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se),  width = 0.2, linewidth = 0.5) + 
    geom_point(size =2) + ylab("Sexual selection") + xlab(predictor) + theme_classic()
}

# Create the labels with italic "n"
one_sample_n <-  expression(~italic(n)~'= 35')
two_sample_n <-  expression(~italic(n)~'= 2253')
three_sample_n <-  expression(~italic(n)~'= 4909')
four_sample_n <-  expression(~italic(n)~'= 2792')

whisker_plot <- sex_meanplot("cert_reverse") + xlab("Data certainty") +  theme_classic(base_size = 30) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = rel(0.9))) +
  scale_y_continuous(limits = c(0,1),expand = expansion(add = c(0.02, 0.09)))+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  coord_cartesian(clip = "off") +
  # annotate("text", x = 1, y = .16, label = "35", size = 5) +
  # annotate("text", x = 2, y = .245, label = "2253", size = 5) +
  # annotate("text", x = 3, y = .305, label = "4909", size = 5) +
  # annotate("text", x = 4, y = .905, label = "2792", size = 5) 
  annotate("text", x = 1, y = 1, label = one_sample_n, size = 6) +
  annotate("text", x = 2, y = 1, label = two_sample_n, size = 6) +
  annotate("text", x = 3, y = 1, label = three_sample_n, size = 6) +
  annotate("text", x = 4, y = 1, label = four_sample_n, size = 6)


ggarrange(whisker_plot, prop_plot, labels = c("a", "b"), align = "v",
          #hjust = c(-8, -7), 
          hjust = c(-1.5, -1.5), 
          font.label = list(size = 28), widths = c(1,1))


ggsave("Plots/Data/figure_s2.png", width = 12, height = 6)

ggarrange(whisker_plot +xlim(c(0.5,4.5)) + rremove("xlab") + rremove("x.text"), prop_plot, ncol = 1, 
          labels = c("a", "b"), heights = c(1,1.25),
          hjust = c(-8, -7), font.label = list(size = 28))


?rremove


################################################################################
              ##### Try bars ######

# Basic bar plot.
grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sexual_selection),
            sex_sd = sd(sexual_selection),
            sex_se = sd(sexual_selection)/sqrt(length(sexual_selection)))
# test_plot <- grouped_data %>% ggplot(aes(x = data_certainty, y = sex_mean, alpha = data_certainty, fill = as.factor(sex_mean))) +
#     geom_bar(stat = "identity") +  scale_fill_manual(values = pal) + 
#     theme_classic(base_size = 20) + 
#     theme(legend.position = "none",
#           line = element_line(linewidth = 0.5))
# 
# after_stat()


# Try merging data.
  merged_data <- clean_data %>% left_join(grouped_data)



merged_data$mean_y <- merged_data$sex_mean

merged_data$mean_y[merged_data$data_certainty == 1] <- mean(merged_data$mean_y[merged_data$data_certainty == 1])/length(merged_data$mean_y[merged_data$data_certainty == 1])
merged_data$mean_y[merged_data$data_certainty == 2] <- mean(merged_data$mean_y[merged_data$data_certainty == 2])/length(merged_data$mean_y[merged_data$data_certainty == 2])
merged_data$mean_y[merged_data$data_certainty == 3] <- mean(merged_data$mean_y[merged_data$data_certainty == 3])/length(merged_data$mean_y[merged_data$data_certainty == 3])
merged_data$mean_y[merged_data$data_certainty == 4] <- mean(merged_data$mean_y[merged_data$data_certainty == 4])/length(merged_data$mean_y[merged_data$data_certainty == 4])

merged_data %>% filter(data_certainty == 4) %>% pull(sex_mean)

merged_data$sex_fact <- factor(merged_data$sexual_selection, levels = 4:0)

merged_data %>% ggplot(aes(x = data_certainty, y = mean_y, 
                                        fill = sex_fact, col = sex_fact)) + 
  geom_col(position = "stack") + 
  scale_fill_manual(values = pal, breaks = 0:4) +
  scale_colour_manual(values = pal, breaks = 0:4) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) + 
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), colour = "black", width = 0.2, linewidth = 0.5) + 
  geom_point(aes(y = sex_mean), size =3, colour = "black", show.legend = FALSE) + ylab("Sexual selection") + xlab("Data certainty")


ggsave("Plots/Data/figure_s2_test.png", width = 8, height = 8)


test_data <- inspect_after_stat(test_plot)

test_plot + aes(y = after_stat(y))

test_data %>% filter( x == 1) %>% pull(y) %>% sum()

test_data %>% ggplot(aes(x = x, y = y/length(y), fill = fill)) + geom_col() 
  

test_data %>% ggplot(aes(x = x, y = y_2, fill = fill, colour = fill)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4) + scale_colour_manual(values = pal, breaks = 0:4) 
ggsave("Plots/Data/test_cert.png")

test_data$y_2 <- test_data$y

test_data$y_2[test_data$x == 1] <- mean(test_data$y_2[test_data$x == 1])/length(test_data$y_2[test_data$x == 1])
test_data$y_2[test_data$x == 2] <- mean(test_data$y_2[test_data$x == 2])/length(test_data$y_2[test_data$x == 2])
test_data$y_2[test_data$x == 3] <- mean(test_data$y_2[test_data$x == 3])/length(test_data$y_2[test_data$x == 3])
test_data$y_2[test_data$x == 4] <- mean(test_data$y_2[test_data$x == 4])/length(test_data$y_2[test_data$x == 4])


grouped_data_both <- clean_data %>% 
  group_by(data_certainty, sexual_selection) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

grouped_data_both %>% filter(data_certainty == 1) %>% ggplot(aes(x = data_certainty, y = sex_mean, fill = as.factor(sexual_selection))) + 
  geom_bar(position = "dodge", stat = "identity")


grouped_data_both %>% filter(data_certainty == 1) %>% ggplot(aes(x = data_certainty, y = mean(sex_mean), fill = as.factor(sexual_selection))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = pal)


p_bar_prop <- ggplot(grouped_data, aes(data_certainty)) +
  geom_bar(
    aes(
      y = sex_mean,
      group = data_certainty == 1 # Assign the 3 bars into 2 groups
    )
  )
p_bar_prop


library(ggplot2)
library(dplyr)

# Assuming you have defined `pal` somewhere for the color palette

grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

merged_data <- clean_data %>% 
  left_join(grouped_data, by = "data_certainty")

ggplot(merged_data, aes(x = data_certainty, y = sex_mean, col = sex_fact, fill = as.factor(sex_fact))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))


library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

# Assuming you have defined `pal` somewhere for the color palette

grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection))) %>%
  mutate(sexual_selection = "All")

ggplot(merged_data, aes(x = data_certainty, y = sex_mean, fill = as.factor(sexual_selection))) +
  geom_bar(position = "fill", stat = "summary", fun = "sum") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))


library(ggplot2)
library(dplyr)

# Assuming you have defined `pal` somewhere for the color palette

grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

# Join the grouped_data back to the clean_data
clean_data_2 <- clean_data %>% left_join(grouped_data, by = "data_certainty")

ggplot(clean_data_2, aes(x = data_certainty, y =sex_mean, fill = as.factor(sexual_selection))) +
  geom_bar(aes(y = sex_mean), position = "stack", stat = "identity") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))

clean_data_2 %>% filter(data_certainty == 4) %>% pull(sex_mean) %>% sum()
clean_data_2 %>% filter(data_certainty == 1) %>% pull(sex_mean) %>% sum()


library(ggplot2)
library(dplyr)

# Assuming you have defined `pal` somewhere for the color palette

grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

# Join the grouped_data back to the clean_data
clean_data_2 <- clean_data %>% left_join(grouped_data, by = "data_certainty")

# Custom position function to stack bars based on sex_mean
position_stackv <- function (width = NULL, preserve = "total", reverse = FALSE) {
  ggproto(NULL, PositionStackv, width = width, preserve = preserve, reverse = reverse)
}

PositionStackv <- ggproto("PositionStackv", PositionStack,
                          setup_params = function(self, data) {
                            if (is.null(data$y)) 
                              stop("This position requires a y aesthetic")
                            y <- data$y
                            y_range <- ggplot2:::position_stack_scales(self, data)$y
                            y_position <- if (self$reverse) {
                              y_range - (1 - y)
                            } else {
                              y
                            }
                            transform(data, y = y_position, ymin = y, ymax = y_position)
                          }
)

ggplot(clean_data_2, aes(x = data_certainty, y = sex_mean, fill = as.factor(sexual_selection))) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))


grouped_data_both <- clean_data %>% 
  group_by(data_certainty, sexual_selection) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

ggplot(grouped_data_both, aes(x = data_certainty, y = sex_mean, fill = as.factor(sexual_selection))) +
  geom_bar(position = "stack", stat = "summary", fun = "mean") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5))





###############################################################################
                 ##### Try average data certainty for each SS #######

grouped_data <- clean_data %>% 
  group_by(sexual_selection) %>% 
  summarise(data_mean = mean(sqrt(data_certainty)),
            data_sd = sd(sqrt(data_certainty)),
            data_se = sd(sqrt(data_certainty))/sqrt(length(data_certainty)))

ggplot(grouped_data, aes(x = sexual_selection, y = data_mean)) +
  geom_errorbar(aes(ymin = data_mean - data_se, ymax = data_mean + data_se),  width = 0.2, linewidth = 0.5) + 
  geom_point(size =3) + ylab("Data certainty") + xlab("Sexual selection") + theme_classic() + 
  geom_jitter(data = clean_data, inherit.aes = FALSE, aes(x = sexual_selection, y = data_certainty),
              position = position_jitternormal(sd_y = 0.05))


ss_barplot <- function(dataset = full_data){
  dataset %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) +
    geom_bar() +  scale_fill_manual(values = pal) + 
    theme_classic(base_size = 20) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5))
}

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


# Model formula for VIF values.
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z + 
    trophic_level_c*terr_bi_c"

model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z + 
    terr_bi_c*temp_seasonality_z"

model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z + 
    terr_bi_c*temp_seasonality_z + migration_bi_c*temp_seasonality_z"


library(car)

# Check full model VIF.
test_model <- lm(model_formula, data = model_data)
vif(test_model) %>% round(digits = 2)

# High certainty VIF.
test_model <- lm(model_formula, data = high_data)
vif(test_model) %>% round(digits = 2)


# See all two-way interactions.
model_formula <- "sexual_score ~ (terr_bi_c + migration_bi_c + 
    trophic_level_c + temp_seasonality_z)^2"
test_model <- lm(model_formula, data = model_data)
vif(test_model) %>% round(digits = 2)




# summary(test_model)
# 
# model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)
# row.names(model_data) <- model_data$tree_tip
# 
# phylolm(model_formula, model_data, tree) %>% summary()


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


