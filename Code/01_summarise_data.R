###############################################################################
                 ##### Summarise sexual selection data #####
###############################################################################

# This script summarises the spread of sexual selection across birds. It also 
# generates figures S1 & S2


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(effectsize)
library(readxl)
library(ggbreak)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####


# Read in the life history traits.
data_pathway <- "Data/supplementary_dataset_1.xlsx"
full_data <- read_excel(data_pathway, sheet = 2, na = "NA") %>% 
  clean_names()
high_data <- full_data %>% filter(data_certainty > 2)

# Read in Clements data.
clements_data <- read_excel(data_pathway, sheet = 3, na = "NA") %>% 
  clean_names()
clements_data$sexual_selection %<>% as.numeric()
clements_data %<>% tidyr::drop_na(sexual_selection)


###############################################################################
                #### Count sexual selection scores ####


# Raw numbers of sexual selection scores.
full_data %>% count(sexual_selection) 
high_data %>% count(sexual_selection)

# Calculate percentages.
counts <- full_data %>% count(sexual_selection) 
high_counts <- full_data %>% filter(data_certainty == 4) %>% count(sexual_selection) 

round((counts[1,2]/9988)*100, 2)  # % Scored 0
round((counts[2,2]/9988)*100, 2)  # % Scored 1
round((counts[3,2]/9988)*100, 2)  # % Scored 2
round((counts[4,2]/9988)*100, 2)  # % Scored 3
round((counts[5,2]/9988)*100, 2)  # % Scored 4

round((sum(counts[2:5,2])/9988)*100, 2) # % Scored 1 - 4

# High counts.
round((high_counts[1,2]/2793)*100, 2)  # % Scored 0
round((sum(high_counts[2:5,2])/2793)*100, 2) # % Scored 1 - 4


###############################################################################
                       #### Figure S1 ####


# Colour pal.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


# Function to make barplot.
ss_barplot <- function(dataset = full_data){
  dataset %>% ggplot(aes(x = sexual_selection, 
                         fill = as.factor(sexual_selection))) +
    geom_bar() +  scale_fill_manual(values = pal) + 
    theme_classic(base_size = 20) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5))
}

# Pull out frugivores and invertivores.
fruit_data <- full_data %>% filter(trophic_niche == "Frugivore")
invert_data <- full_data %>% filter(trophic_niche == "Invertivore")

# Get axis limits.
jetz_y <- full_data %>% count(sexual_selection) %>% pull(n) %>% max()
clements_y <- clements_data %>% count(sexual_selection) %>% pull(n) %>% max()
fruit_y <- fruit_data %>% count(sexual_selection) %>% pull(n) %>% max()
invert_y <- invert_data %>% count(sexual_selection) %>% pull(n) %>% max()


# Create plots.
jetz_plot <- full_data %>% ss_barplot() + 
  scale_y_continuous(limits = c(0,jetz_y*1.1)) +
  ylab("Species count")  + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(vjust = 1.5)) +
  scale_y_break(c(1000, 5000), scales = "free", 
                ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + 
  ylim(0,10000) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

clements_plot <- clements_data %>% ss_barplot() + 
  scale_y_continuous(limits = c(0,clements_y*1.1), 
                     breaks = c(0, 250, 500, 750, 1000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_break(c(1000, 5000), scales = "free", 
                ticklabels = c(4000, 5000, 6000, 7000, 8000, 9000)) + 
  ylim(0,10000) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

fruit_plot <- fruit_data %>%  
  ss_barplot() + xlab("Sexual selection") + ylab("Species count") +
  scale_y_continuous(limits = c(0, 5500)) + 
  theme(axis.title.y = element_text(vjust = 1.5)) +
  scale_y_break(c(250, 500), scales = "free", 
                ticklabels = c(500, 1500, 2500, 3500, 4500)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

invert_plot <- full_data %>% filter(trophic_niche == "Invertivore") %>% 
  ss_barplot() + scale_y_continuous(limits = c(0,5500)) +
  xlab("Sexual selection") + theme(axis.title.y = element_blank()) +
  scale_y_break(c(250, 500), scales = "free", 
                ticklabels = c(500, 1500, 2500, 3500, 4500)) + 
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right  = element_blank(),
        axis.line.y.right  = element_blank())

# Add the annotations. 
jetz_anno <- jetz_plot + 
  annotate("text", x = -0.3, y = 9750, label = "a", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 9000, label = "Taxonomy 1\nn = 9988", size = 7)
clem_anno <- clements_plot  + 
  annotate("text", x = -0.3, y = 9750, label = "b", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 9000, label = "Taxonomy 2\nn = 10671", size = 7)
fruit_anno <- fruit_plot + 
  annotate("text", x = -0.3, y = 5250, label = "c", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Frugivores\nn = 1030", size = 7)
invert_anno <- invert_plot  +  
  annotate("text", x = -0.3, y = 5250, label = "d", size = 10, fontface = 2) +
  annotate("text", x = 3, y = 4500, label = "Invertivores\nn = 4779", size = 7)

# Arrange figures.
ggarrange(print(jetz_anno), 
          print(clem_anno),
          print(fruit_anno), 
          print(invert_anno), 
          ncol = 2, nrow = 2, widths = c(1.1, 1), heights = c(1, 1.1))

# Export.
ggsave("Figures/Fig_S1.tiff", width = 8, height = 8, compression = "lzw")


################################################################################
                   ###### Figure S2 ######


# Group the data by certainty score, and calculate mean and standard error.
grouped_data <- full_data %>% 
  group_by(data_certainty) %>% 
  summarise(sex_mean = mean(sqrt(sexual_selection)),
            sex_sd = sd(sqrt(sexual_selection)),
            sex_se = sd(sqrt(sexual_selection))/sqrt(length(sexual_selection)))

# Merge the averages with full dataset.
merged_data <- full_data %>% 
  left_join(grouped_data, by = "data_certainty")
merged_data$sex_fact <- factor(merged_data$sexual_selection, levels = 4:0)

# Create the sample size labels with italic "n"
one_sample_n <-  expression(~italic(n)~'= 35')
two_sample_n <-  expression(~italic(n)~'= 2253')
three_sample_n <-  expression(~italic(n)~'= 4909')
four_sample_n <-  expression(~italic(n)~'= 2792')

# Create the whisker plot for panel a.
whisker_plot <- ggplot(grouped_data, aes(x = data_certainty, y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se),  
                width = 0.2, linewidth = 0.5) + 
  geom_point(size =2) + ylab("Sexual selection") + xlab("Data certainty") +  
  theme_classic(base_size = 30) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = rel(0.9))) +
  scale_y_continuous(limits = c(0,1),expand = expansion(add = c(0.02, 0.09))) +
  scale_x_continuous(expand = expansion(mult = 0.1))+
  coord_cartesian(clip = "off") +
  annotate("text", x = 1, y = 1, label = one_sample_n, size = 6) +
  annotate("text", x = 2, y = 1, label = two_sample_n, size = 6) +
  annotate("text", x = 3, y = 1, label = three_sample_n, size = 6) +
  annotate("text", x = 4, y = 1, label = four_sample_n, size = 6)

# Create the plot of proportions for panel b.
prop_plot <- ggplot(merged_data, aes(x = data_certainty, y = sex_mean, 
                                     col = sex_fact, fill = sex_fact)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4, 
                    labels = c("SS = 0", "SS = 1", "SS = 2", "SS = 3", "SS = 4")) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") + 
  xlab("Data certainty") +
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
        legend.title = element_blank(),
        axis.title.x = element_text(size = rel(0.9)))
prop_plot

prop_plot <- ggplot(merged_data, aes(x = data_certainty, y = sex_mean, 
                                     col = sex_fact, fill = sex_fact)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal, breaks = 0:4) +
  scale_colour_manual(values = pal, breaks = 0:4) + ylab("% of species") +
  xlab("Data certainty") +
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
        axis.title.x = element_text(size = rel(0.9)))

# Put the plots together.
ggarrange(whisker_plot, prop_plot, labels = c("a", "b"),
          hjust = c(-1.5, -1.5), 
          font.label = list(size = 28), widths = c(1,1.2))

# Export the plot.
ggsave("Figures/Fig_S2.tiff", width = 13, height = 6, compression = "lzw")


###############################################################################
                             #### END ####
###############################################################################