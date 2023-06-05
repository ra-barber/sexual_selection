###############################################################################
                        ##### Script Title #####
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


parental_care_data <- readLines("Data/parental_care.txt")
parental_care_data <- str_flatten(parental_care_data)

bird_tree_names <- str_split(phylo_data$birdtree_name, pattern = " ", simplify = TRUE) %>% as.data.frame()
colnames(bird_tree_names) <- c("genus", "species")

parental_care_data %>% str_extract_all(bird_tree_names$genus)

length(str_split(parental_care_data, pattern = bird_tree_names$species[10]))

str_split_fixed(parental_care_data, pattern = bird_tree_names$species[10])


str_remove_all(parental_care_data, c("[1-9]", "bbbbbbbba", "--", "Neotropical") )
parental_care_data %<>% str_remove_all("[0-9]")
parental_care_data %<>% str_remove_all(" --")
parental_care_data %<>% str_remove_all("Nearctic")
parental_care_data %<>% str_remove_all("bbdbbbbbbbbbbbbaaba")

parental_care_data %<>% str_remove_all("Appendix A")
parental_care_data %<>% str_remove_all("Neotropical")
parental_care_data %<>% str_remove_all("Australia ")

parental_care_data %<>% str_remove_all("\\[")
parental_care_data %<>% str_remove_all("\\]")
parental_care_data %<>% str_remove_all(",")
parental_care_data %<>% str_remove_all("Africa")
parental_care_data %<>% str_remove_all("Widespread")
parental_care_data %<>% str_remove_all("Page  of  // : AM")

parental_care_data

pc_flat <- parental_care_data %>% str_flatten()

patterns <- c("[A-Z][a-z]+ [a-z]+ Female only Female only", 
              "[A-Z]\\. [a-z]+ Female only Female only",
              "[A-Z]\\. [a-z]+ Unknown Female only")

matches <- pc_flat %>% str_extract_all(patterns)

species <- matches[[2]] %>% str_extract_all("^[a-z]+")

bird_tree_names[bird_tree_names$species %in% species,] 

patterns <- c("[a-z]+ Female only Female only", 
              "[a-z]+ Female only Female only",
              "[a-z]+ Unknown Female only")

pc_flat %>% str_split(": ", simplify = TRUE)

parental_care_data %>% str_extract("[A-Z][a-z]+ [a-z]+ ")



library(stringr)

str_c(parental_care_data, )
parental_care_data %>% str_extract_all("[A-Z][a-z]+ [a-z]+ Female only Female only")
str_flatten(parental_care_data) %>% str_extract_all("[A-Za-z]{1,2}\\. *[A-Za-z]+( +[A-Za-z]+)* Female only Female only")



for (x in 1:10){
  parental_care_data  %>% str_detect(bird_tree_names$genus[1])
}


patterns <- c("[A-Z][a-z]+ [a-z]+", "[A-Za-z]{1,2}\\. *[A-Za-z]+( +[A-Za-z]+)*")
patterns <- c("[A-Z][a-z]+ [a-z]+ Female only Female only", "[A-Za-z]{1,2}\\. *[A-Za-z]+( +[A-Za-z]+)* Female only Female only")
patterns <- c("[A-Z][a-z]+ [a-z]+ Female only Female only", 
              "[A-Z]\\. [a-z]+ Female only Female only",
              "[A-Z]\\. [a-z]+ Unknown Female only")





str_flatten(parental_care_data) %>% str_extract_all(phylo_data$birdtree_name)

parental_care_data[parental_care_data %>% str_detect(pattern = matches[[1]][1])]


?str_extract_all
###############################################################################
                           #### Section 1 ####



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


