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
library(readxl)


# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

# Read in the life history traits.
data_pathway <- "Data/sexual_selection_dataset_04_09.xlsx"
full_data <- read_excel(data_pathway, sheet = 2, na = "NA") %>% 
  clean_names()
high_data <- full_data %>% filter(data_certainty > 2)

# Read in a tree.
tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Read in Clements data.
clements_data <- read_excel(data_pathway, sheet = 3, na = "NA") %>% 
  clean_names()
clements_data$sexual_selection %<>% as.numeric()
clements_data %<>% tidyr::drop_na(sexual_selection)


###############################################################################
                    #### Transform data as usual ####

# Set as factor, then re-level for appropriate reference group.
full_data %<>% mutate(
  territoriality_binary = relevel(as.factor(territoriality_binary), ref = "Non-territorial"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_level_binary = relevel(as.factor(trophic_level_binary), ref = "Secondary"),

  # Center categorical predictors.
  terr_bi_c = center_categorical(territoriality_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_level_binary),
  
  # Center / scale continuous.
  centroid_z = standardize(sqrt(abs(latitude)), two_sd = TRUE),
  temp_seasonality_z = standardize(log(seasonality), two_sd = TRUE))


# Set as factor, then re-level for appropriate reference group.
high_data %<>% mutate(
  territoriality_binary = relevel(as.factor(territoriality_binary), ref = "Non-territorial"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_level_binary = relevel(as.factor(trophic_level_binary), ref = "Secondary"),
  
  # Center categorical predictors.
  terr_bi_c = center_categorical(territoriality_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_level_binary),
  
  # Center / scale continuous.
  centroid_z = standardize(sqrt(abs(latitude)), two_sd = TRUE),
  temp_seasonality_z = standardize(log(seasonality), two_sd = TRUE))


###############################################################################
                #### Count sexual selection scores ####

# Raw numbers of sexual selection scores.
full_data %>% count(sexual_selection) 
high_data %>% count(sexual_selection)

# Calculate percentages.
counts <- full_data %>% count(sexual_selection) 
high_counts <- full_data %>% filter(data_certainty == 4) %>%  count(sexual_selection) 

round((counts[1,2]/9988)*100, 2)  # % Scored 0
round((counts[2,2]/9988)*100, 2)  # % Scored 1
round((counts[3,2]/9988)*100, 2)  # % Scored 2
round((counts[4,2]/9988)*100, 2)  # % Scored 3
round((counts[5,2]/9988)*100, 2)  # % Scored 4

round((sum(counts[2:5,2])/9988)*100, 2) # % Scored 1 - 4

# High counts.
round((high_counts[1,2]/2793)*100, 2)  # % Scored 0
round((sum(high_counts[2:5,2])/2793)*100, 2) # % Scored 1 - 4


################################################################################
    #### Calculate phylogenetic signal sexual selection scores #####

# This should be a separate script on the HPC.
full_data$tree_tip <- full_data$scientific_name_bird_tree %>% str_replace(" ", "_")
row.names(full_data) <- full_data$tree_tip

full_data %<>% mutate(
  sexual_binary = replace(sexual_selection, sexual_selection > 0, 1)
)

sex_data <- full_data %>% select(tree_tip, sexual_binary)

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
  dataset %>% ggplot(aes(x = sexual_selection, fill = as.factor(sexual_selection))) +
    geom_bar() +  scale_fill_manual(values = pal) + 
    theme_classic(base_size = 20) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5))
}

fruit_data <- full_data %>% filter(trophic_niche == "Frugivore")
invert_data <- full_data %>% filter(trophic_niche == "Invertivore")

jetz_y <- full_data %>% count(sexual_selection) %>% pull(n) %>% max()
clements_y <- clements_data %>% count(sexual_selection) %>% pull(n) %>% max()
fruit_y <- fruit_data %>% count(sexual_selection) %>% pull(n) %>% max()
invert_y <- invert_data %>% count(sexual_selection) %>% pull(n) %>% max()



library(ggbreak)


jetz_plot <- full_data %>% ss_barplot() + scale_y_continuous(limits = c(0,jetz_y*1.1)) +
  ylab("Species count")  + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(vjust = 1.5)) +
  scale_y_break(c(1000, 5000), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + ylim(0,10000) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

clements_plot <- clements_data %>% ss_barplot() + scale_y_continuous(limits = c(0,clements_y*1.1), breaks = c(0, 250, 500, 750, 1000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank()) +
  scale_y_break(c(1000, 5000), scales = "free", ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + ylim(0,10000) +
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
  annotate("text", x = 3, y = 9000, label = "Taxonomy 1\nn = 9988", size = 7)

clem_anno <- clements_plot  + annotate("text", x = -0.3, y = 9750, label = "b", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 9000, label = "Taxonomy 2\nn = 10671", size = 7)

fruit_anno <- fruit_plot + annotate("text", x = -0.3, y = 5250, label = "c", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Frugivores\nn = 1030", size = 7)

invert_anno <- invert_plot  +  annotate("text", x = -0.3, y = 5250, label = "d", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Invertivores\nn = 4779", size = 7)

# Amaurospiza carrizalensis 
ggarrange(print(jetz_anno), 
          print(clem_anno),
          print(fruit_anno), 
          print(invert_anno), 
          ncol = 2, nrow = 2, widths = c(1.1,1), heights = c(1, 1.1))

# Export.
ggsave("Plots/Data/extended_data_figure_1.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/Data/extended_data_figure_1.tiff", width = 8, height = 8)



################################################################################
                   ###### Supplementary figure 1 ######

# Read in the most up to date version of the SS database.
clean_data <- read.csv("Data/sexual_selection_cleaned_01_08.csv") %>% clean_names()  # We have since cleaned this data even more.

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

#pal <- c('#F21A00', '#E1AF00', '#EBCC2A', '#78B7C5', '#3B9AB2')

# Group the data by certainty score, and calculate mean and standard error.
grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

# Merge the averages with full dataset.
merged_data <- clean_data %>% 
  left_join(grouped_data, by = "data_certainty")
merged_data$sex_fact <- factor(merged_data$sexual_selection, levels = 4:0)  # For plotting.

# Create the sample size labels with italic "n"
one_sample_n <-  expression(~italic(n)~'= 35')
two_sample_n <-  expression(~italic(n)~'= 2253')
three_sample_n <-  expression(~italic(n)~'= 4909')
four_sample_n <-  expression(~italic(n)~'= 2792')

# Create the whisker plot for panel a.
whisker_plot <- ggplot(grouped_data, aes(x = data_certainty, y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se),  width = 0.2, linewidth = 0.5) + 
  geom_point(size =2) + ylab("Sexual selection") + xlab("Data certainty") +  theme_classic(base_size = 30) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = rel(0.9))) +
  scale_y_continuous(limits = c(0,1),expand = expansion(add = c(0.02, 0.09)))+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  coord_cartesian(clip = "off") +
  annotate("text", x = 1, y = 1, label = one_sample_n, size = 6) +
  annotate("text", x = 2, y = 1, label = two_sample_n, size = 6) +
  annotate("text", x = 3, y = 1, label = three_sample_n, size = 6) +
  annotate("text", x = 4, y = 1, label = four_sample_n, size = 6)

# Create the plot of proportions for panel b.
prop_plot <- ggplot(merged_data, aes(x = data_certainty, y = sex_mean, col = sex_fact, fill = sex_fact)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4, 
                    labels = c("SS = 0", "SS = 1", "SS = 2", "SS = 3", "SS = 4")) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") + xlab("Data certainty") +
  theme_classic(base_size = 30) + 
  labs(fill = "Sexual\nselection") +
  guides(colour = "none", fill = guide_legend(byrow = TRUE, reverse=TRUE)) +
  scale_x_continuous(breaks = c(1,2,3,4)) +
  scale_y_continuous(limits = c(0,1), expand = expansion(add = c(0.02, 0.09)))+
  theme(legend.position = c(1.2,0.5),
        plot.margin = margin(r = 4, t = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        line = element_line(linewidth = 0.5),
        legend.spacing.y = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        #legend.title = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        axis.title.x = element_text(size = rel(0.9)))
prop_plot

prop_plot <- ggplot(merged_data, aes(x = data_certainty, y = sex_mean, col = sex_fact, fill = sex_fact)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") + xlab("Data certainty") +
  theme_classic(base_size = 30) + 
  labs(fill = "SS") +
  guides(colour = "none", fill = guide_legend(byrow = TRUE, reverse=TRUE)) +
  scale_x_continuous(breaks = c(1,2,3,4)) +
  scale_y_continuous(limits = c(0,1), expand = expansion(add = c(0.02, 0.09)))+
  theme(legend.position = c(1.1,0.5),
        plot.margin = margin(r = 2.5, t = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        line = element_line(linewidth = 0.5),
        legend.spacing.y = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(size = rel(0.8), margin = margin(b = -15)),
        #legend.title = element_blank(),
        axis.title.x = element_text(size = rel(0.9)))

# Put the plots together.
ggarrange(whisker_plot, prop_plot, labels = c("a", "b"),
          hjust = c(-1.5, -1.5), 
          font.label = list(size = 28), widths = c(1,1.2))

# Export the plot.
ggsave("Plots/Data/figure_s1.png", width = 13, height = 6)





###############################################################################
              #### Check raw correlations & VIFs ####

library(GGally)

lat_variables <- sexual_selection %>% dplyr::select(centroid_z, temp_seasonality_z, 
                                              terr_bi_c, migration_bi_c, trophic_level_c)

 
ggcorr(lat_variables, label = TRUE, method =c("everything","spearman"), label_round = 3)

# Model formula for VIF values.
model_formula <- "sexual_selection ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z + 
    terr_bi_c*temp_seasonality_z"


library(car)

# Check full model VIF.
test_model <- lm(model_formula, data = sexual_selection)
vif(test_model) %>% round(digits = 2)

# High certainty VIF.
test_model <- lm(model_formula, data = high_data)
vif(test_model) %>% round(digits = 2)

# See all two-way interactions.
model_formula <- "sexual_selection ~ (terr_bi_c + migration_bi_c + 
    trophic_level_c + temp_seasonality_z)^2"
test_model <- lm(model_formula, data = sexual_selection)
vif(test_model) %>% round(digits = 2)


###############################################################################
                #### Spatial autocorrelation ####

library(spdep)
library(rgeos)
library(sp)
library(sf)

# Select data.
just_sex_long_lat <- full_data %>% select(sexual_selection, centroid_latitude, centroid_longitude)

# Create binary trait.
just_sex_long_lat$sexual_binary <- 0
just_sex_long_lat$sexual_binary[just_sex_long_lat$sexual_selection > 0] <- 1
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
x <- just_sex_long_lat$sexual_selection %>% as.factor()
just_sex_long_lat
# Perform Monte Carlo simulation
join_count <- joincount.mc(listw = w, fx = x, nsim = 999, alternative = "two.sided")

# Simulated join count statistics
join_count

joincount.test(x, w)

# Suggests significant autocorrelation.




###############################################################################
                             #### END ####
###############################################################################



################################################################################
    #### Code to create a bar plot, filled by proportions of category ####


# Basic bar plot.
grouped_data <- clean_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sexual_selection),
            sex_sd = sd(sexual_selection),
            sex_se = sd(sexual_selection)/sqrt(length(sexual_selection)))
# test_plot <- grouped_data %>% ggplot(aes(x = data_certainty, y = sex_mean, alpha = data_certainty, fill = as.factor(sex_mean))) +


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

