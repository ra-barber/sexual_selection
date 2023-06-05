###############################################################################
               ##### Testing the ebird crosswalk #####
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


###############################################################################
                   #### Read in crosswalk ####

# Read in the ebird crosswalk.
ebird_crosswalk <- read.csv("Data/crosswalk_BT_EB.csv") %>% clean_names()

# Read in the ebird tree.
ebird_tree <- read.tree("Data/Trees/ebird_24may.tre")
ebird_tree_species <- ebird_tree$tip.label %>% str_replace("_", " ")


###############################################################################
                  #### Extinct ebird species ####

# Pull out extinct species.
extinct_ebird <- ebird_crosswalk %>% filter(match_notes == "Extinct")

# Check to make sure none are in the tree. (all NA so obvs not in the tree.)
intersect(extinct_ebird$ebird_species, ebird_tree_species)

# Check possib extinct.




ebird_crosswalk_2 <- ebird_crosswalk %>% filter(match_notes != "Extinct")

###############################################################################
                #### Remove 1:1 matches and check ####

ebird_crosswalk %>% count(match_notes)


bt1_eb1 <-  ebird_crosswalk %>% filter(match_notes == "1BT to 1eBird")
nbt1_eb1 <-  ebird_crosswalk %>% filter(match_notes != "1BT to 1eBird")

intersect(bt1_eb1$bird_tree_species, nbt1_eb1$bird_tree_species)

ebird_dupes <- ebird_crosswalk %>% get_dupes(ebird_species)


ebird_dupes %>% count(match_notes)



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


