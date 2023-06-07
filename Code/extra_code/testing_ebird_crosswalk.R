###############################################################################
               ##### Testing the ebird crosswalk #####
###############################################################################


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
ebird_crosswalk <- read.csv("Data/crosswalk_BT_EB_07_06.csv") %>% clean_names()

# Read in the ebird tree.
ebird_tree <- read.tree("Data/Trees/ebird_24may.tre")
ebird_tree_species <- ebird_tree$tip.label %>% str_replace("_", " ")

# Remove redundant rows.
ebird_crosswalk %<>% filter(!is.na(robs_matches))

# Remove species not in the tree using the new rob match column.
ebird_crosswalk %<>% filter(robs_matches != "Not tree species")

# 157 species that should be in the crosswalk but aren't. Only 17 are extant.
setdiff(ebird_tree_species, ebird_crosswalk$ebird_species)

# Print unique list of match types.
unique(ebird_crosswalk$match_notes)
unique(ebird_crosswalk$robs_matches)


################################################################################
                    #### Seperate new tree species ####

# Pull out new species.
new_ebird <- ebird_crosswalk %>% filter(robs_matches == "New tree species")
not_new_ebird <- ebird_crosswalk %>% filter(robs_matches != "New tree species")

# Check to make sure none are in the tree. (all NA so obvs not in the tree.)
intersect(new_ebird$ebird_species, ebird_tree_species)

# Make sure there's no species in both groups.
intersect(new_ebird$ebird_species, not_new_ebird$ebird_species)

# Create a new crosswalk object for next round of filtering.
ebird_crosswalk_2 <- ebird_crosswalk %>% filter(robs_matches != "New tree species")

# Dataset to keep.
new_ebird


################################################################################
                    #### Remove 1 to 1 matches ####

# Print unique list of match types.
unique(ebird_crosswalk_2$match_notes)
unique(ebird_crosswalk_2$robs_matches)

# Pull out one to one ebird species.
one_to_one_ebird <- ebird_crosswalk_2 %>% filter(robs_matches == "1BT to 1eBird")
non_matches <- ebird_crosswalk_2 %>% filter(robs_matches != "1BT to 1eBird")

# Make sure there's no species in both groups.
intersect(one_to_one_ebird$ebird_species, non_matches$ebird_species)

# Create a new crosswalk object for next round of filtering.
ebird_crosswalk_3 <- ebird_crosswalk_2 %>% filter(robs_matches != "1BT to 1eBird")

# Dataset to keep.
one_to_one_ebird


################################################################################
                   #### Remove 1 to many matches ####


# Print unique list of match types.
unique(ebird_crosswalk_3$match_notes)
unique(ebird_crosswalk_3$robs_matches)

# filter for 1bT to many ebird
one_to_many_ebird <-  ebird_crosswalk_3 %>% filter(robs_matches == "1BT to many eBird")
non_matches <-  ebird_crosswalk_3 %>% filter(robs_matches != "1BT to many eBird")

# Make sure there's no species in both groups.
intersect(one_to_many_ebird$ebird_species, non_matches$ebird_species)

# No issues. Filter out these species.
ebird_crosswalk_4 <-  ebird_crosswalk_3 %>% filter(robs_matches != "1BT to many eBird")

# Dataset to keep.
one_to_many_ebird



################################################################################
                #### Remove many to 1 matches ####

# Print unique list of match types.
unique(ebird_crosswalk_4$match_notes)
unique(ebird_crosswalk_4$robs_matches)

# filter for 1bT to many ebird
many_to_one_ebird <-  ebird_crosswalk_4 %>% filter(robs_matches == "Many BT to 1eBird")
non_matches <-  ebird_crosswalk_4 %>% filter(robs_matches != "Many BT to 1eBird")

# Make sure there's no species in both groups.
intersect(many_to_one_ebird$ebird_species, non_matches$ebird_species)

# No issues. Filter out these species.
ebird_crosswalk_5 <-  ebird_crosswalk_4 %>% filter(robs_matches != "Many BT to 1eBird")

# Dataset to keep.
many_to_one_ebird

################################################################################
                   #### Last crosswalk group ####

# Check no other values.
unique(ebird_crosswalk_5$robs_matches)

# Assign name to match others.
many_to_many_ebird <- ebird_crosswalk_5



################################################################################
                    #### Assign easy nominates ####

# No missing nominates from other group.
skim(many_to_one_ebird)

# Check one to many group. (Most are blank)
one_to_many_ebird %>% count(nominate)

# Filter for species where they match exactly. 
matching_binom_data <- one_to_many_ebird %>% filter(bird_tree_species == ebird_species)

# Assign to main one to many dataset the ebird species that are in the matchin binom dataset.
one_to_many_ebird$nominate[one_to_many_ebird$ebird_species %in% matching_binom_data$ebird_species] <- "Yes"

# Now assign all the leftover birdtree species that were part of the same complex.
one_to_many_ebird$nominate[one_to_many_ebird$bird_tree_species %in% matching_binom_data$bird_tree_species & one_to_many_ebird$nominate != "Yes"] <- "No"

# Remove those birdlife species that are matched up.
extra_ebird_data <- one_to_many_ebird %>% filter(!bird_tree_species %in% matching_binom_data$bird_tree_species)

# Create a species bird tree and ebird column.
extra_ebird_data$bt_species <- str_split(extra_ebird_data$bird_tree_species, pattern = " ", simplify = TRUE)[,2]
extra_ebird_data$eb_species <- str_split(extra_ebird_data$ebird_species, pattern = " ", simplify = TRUE)[,2]

# Filter for species where they match exactly. 
matching_species_data <- extra_ebird_data %>% filter(bt_species == eb_species)

# Check for dupes. (Checked each species column indvidually and no dupes)
matching_species_data %>% get_dupes()

# Assign nominates to matching species.
one_to_many_ebird$nominate[one_to_many_ebird$ebird_species %in% matching_species_data$ebird_species]  <- "Yes"

# Assign Nos to birdtree species.
one_to_many_ebird$nominate[one_to_many_ebird$bird_tree_species %in% matching_species_data$bird_tree_species & one_to_many_ebird$nominate != "Yes"] <- "No"

# Only 55 blanks left.
one_to_many_ebird %>% count(nominate)


################################################################################
                #### Check remaining nominates ####

extra_ebird_data_2 <- one_to_many_ebird %>% filter(!bird_tree_species %in% matching_binom_data$bird_tree_species &
                                                   !bird_tree_species %in% matching_species_data$bird_tree_species)

# Easier to do by hand.
clean_ebird_crosswalk <- rbind(new_ebird, one_to_one_ebird, one_to_many_ebird, many_to_one_ebird, many_to_many_ebird)

# Check dupes.
clean_ebird_crosswalk %>% get_dupes()

# Export.  Commented out after making edits.
#write.csv(clean_ebird_crosswalk, "Data/clean_crosswalk_BT_EB_07_06.csv", row.names = FALSE)


################################################################################
            #### add in missing species #####

# 157 species that should be in the crosswalk but aren't. Only 17 are extant.
missing_species <- setdiff(ebird_tree_species, clean_ebird_crosswalk$ebird_species)

write.csv(missing_species, "Data/missing_ebird_species.csv", row.names = FALSE)



setdiff(clean_ebird_crosswalk$ebird_species, ebird_tree_species)


# Check for missing jetz species.

jetz_data <- read.csv("Data/ecotraits_09_05_without_invalids.csv") %>% clean_names()












one_to_many_ebird %>% filter(bird_tree_species %in% matching_binom_data$bird_tree_species)



###############################################################################
                #### Extinct ebird species ####

# Pull out extinct species.
extinct_ebird <- ebird_crosswalk %>% filter(match_notes == "Extinct")

# Check to make sure none are in the tree. (all NA so obvs not in the tree.)
intersect(extinct_ebird$ebird_species, ebird_tree_species)

# Check possib extinct.
poss_extinct_ebird <- ebird_crosswalk %>% filter(match_notes == "Poss. Extinct")
intersect(poss_extinct_ebird$ebird_species, ebird_tree_species)

# Remove species not in the tree using the new rob match column.
ebird_crosswalk %>% filter(robs_matches == "Not tree species")


# Remove extinct species.
ebird_crosswalk_2 <- ebird_crosswalk %>% filter(match_notes != "Extinct")
ebird_crosswalk_2 <- ebird_crosswalk_2 %>% filter(match_notes != "Poss. Extinct")

# There are extra extinct ebird species included in the tree.

intersect(poss_extinct_ebird$ebird_species, ebird_tree_species)


###############################################################################
                #### Look at newly described species ####

# Print unique list of match types.
unique(ebird_crosswalk_2$match_notes)

# Pull out newly described species.
new_ebird <- ebird_crosswalk_2 %>% filter(match_notes == "Newly described species")

# Remove newly described birdlife species.
new_ebird %<>% tidyr::drop_na(ebird_species)

# Remove newly described species.
ebird_crosswalk_3 <- ebird_crosswalk_2 %>% filter(!ebird_species %in% new_ebird$ebird_species)



###############################################################################
                  #### Remove 1:1 matches and check ####

# Print unique list of match types.
unique(ebird_crosswalk_3$match_notes)

# Filter for 1 to 1 matches.
one_bt_to_one_eb <-  ebird_crosswalk_3 %>% filter(match_notes == "1BT to 1eBird")

# Pull out non-matches.
non_one_to_one_matches <-  ebird_crosswalk_3 %>% filter(match_notes != "1BT to 1eBird")

# Check overlap.
intersect(one_bt_to_one_eb$bird_tree_species, non_one_to_one_matches$bird_tree_species)
intersect(one_bt_to_one_eb$ebird_species, non_one_to_one_matches$ebird_species)

# Create new ebird crosswalk list.
ebird_crosswalk_4 <-  ebird_crosswalk_3 %>% filter(match_notes != "1BT to 1eBird")


###############################################################################
              #### Pull out 1 1BT to many eBird splits ####

# Print unique list of match types.
unique(ebird_crosswalk_4$match_notes)

# Filter for 1 to 1 matches.
one_bt_to_many_eb <-  ebird_crosswalk_4 %>% filter(match_notes == "1BT to many eBird")

# Pull out non-matches.
none_matches <-  ebird_crosswalk_4 %>% filter(match_notes != "1BT to many eBird")

# Check overlap.
intersect(one_bt_to_many_eb$bird_tree_species, none_matches$bird_tree_species)
intersect(one_bt_to_many_eb$ebird_species, none_matches$ebird_species)



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


