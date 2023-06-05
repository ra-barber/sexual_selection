###############################################################################
                        ##### Plotting correlations #####
###############################################################################

# This script does something new. Pretty sick eh.


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

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

# Read in some data.
sexual_traits <- read.csv("Data/sexual_traits.csv")

fast_slow <- read.csv("Data/fast_slot.csv") %>% clean_names()


test$bird_tree_name %<>% gsub(pattern = " ", replacement = "_")
colnames(fast_slow)[1] <- "bird_tree_name"

fast_slow %<>% distinct(bird_tree_name, .keep_all = TRUE)

test <-fast_slow %>% select(bird_tree_name, pc1_f)

test <- left_join(test, sexual_traits, by = c("bird_tree_name"="bird_tree_name"))


ggplot(test, aes(x = as.factor(sexual_score), y = body_mass_log)) + geom_boxplot()

lm(sexual_score ~ pc1_f, data = test) %>% summary()
read.tree()



ggplot(test, aes(x = gen_log, y = pc1_f)) + geom_point() + geom_smooth(method = "lm")
library(GGally)
test %>% select(gen_log, body_mass_log, pc1_f) %>% ggpairs(labels = TRUE)

colnames(test)

library(caper)
# Read in a tree.
tree <- read.tree("Data/Trees/model_trees.tre")[[1]]


sexual_traits$test_score <- sexual_traits$sexual_score
sexual_traits$test_score[sexual_traits$test_score < 3] <- 2


sexual_traits$bird_tree_name %<>% gsub(pattern = " ", replacement = "_")
row.names(sexual_traits) <- sexual_traits$bird_tree_name

library(phylolm)
phylolm(test_score ~ gen_log + body_mass_log, data = sexual_traits, phy = tree) %>% summary()
hist(sexual_traits$test_score)
linear_model <- lm(gen_log ~ body_mass_log, data = sexual_traits)
head(sexual_traits)
sexual_traits$resids <- linear_model$residuals

sexual_traits %>% dplyr::select(gen_log, body_mass_log, resids) %>% ggcorr(label = TRUE)

# Read in a raster.
raster <- raster("")


###############################################################################
                    #### Plot correlations with latitude ####

colnames(sexual_traits)
ggplot(sexual_traits, aes(x = sqrt(abs(centroid_latitude)), y = jitter(sexual_score, amount = 0.5), col = trophic_binary)) + 
  geom_point(alpha = 0.2) + geom_smooth(method = "lm") + theme_classic() + facet_wrap( ~ trophic_binary)

ggplot(sexual_traits, aes(x = npp, y = jitter(sexual_score), col = trophic_binary)) + 
  geom_point(alpha = 0.5) + theme_classic() + facet_wrap( ~ trophic_binary)

ggplot(sexual_traits, aes(x = temp_log, y = jitter(sexual_score, amount = 0.4), col = trophic_binary)) + 
  geom_point(alpha = 0.25) + theme_classic() + facet_wrap( ~ trophic_binary)

ggplot(sexual_traits, aes(x = centroid_latitude, y = npp)) + geom_point()

colnames(sexual_traits)
table(sexual_traits$sexual_score)
###############################################################################
                           #### Section 2 ####


###############################################################################
                           #### Section 3 ####


###############################################################################
                           #### Section 4 ####

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


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


