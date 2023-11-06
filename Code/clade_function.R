# ###############################################################################
#                         ##### Clade function #####
# ###############################################################################
# 
# # A function to add clades for plotting on a tree.
# 
# assign_clades <- function(data){
#  
#   # Orders for PALAEOGNATHAE (in brain paper)
#   Palaeognathae <- c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
#                      "CASUARIIFORMES", "TINAMIFORMES")
#   
#   # Orders for Galloanseres. (in brain paper)
#   Galloanseres <- c("GALLIFORMES", "ANSERIFORMES")
#   
#   # Orders for aequornithes (in brain paper)
#   Aequornithes <- c("GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", "CICONIIFORMES", 
#                     "SULIFORMES", "PELECANIFORMES", "OPISTHOCOMIFORMES")
#   
#   # Telluraves. (Core land birds)
#   Afroaves <- c("ACCIPITRIFORMES", "STRIGIFORMES", "COLIIFORMES", "LEPTOSOMIFORMES",
#                   "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", "PICIFORMES")
#   
#   # Australaves (Austrailian radiation)
#   Australaves <- c("CARIAMIFORMES",  "FALCONIFORMES", "PSITTACIFORMES", "PASSERIFORMES")
#   
#   # Strisores. (in brain paper)
#   Strisores <- c("APODIFORMES", "CAPRIMULGIFORMES")
#   
#   # Columbimorphae (in brain paper)
#   Columbimorphae <- c("COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES")
#   
#   # Mirandornithes
#   Mirandornithes <- c("PHOENICOPTERIFORMES","PODICIPEDIFORMES")
#   
#   # Eurpygimorphae
#   Eurpygimorphae <- c("EURYPYGIFORMES", "PHAETHONTIFORMES")
#   
#   # Otidimorphae (in brain paper)
#   Otidimorphae <- c("MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES")
#   
#   # orphaned groups
#    Opisthocomiformes <- c("OPISTHOCOMIFORMES")
#   # Charadriiformes <- c("CHARADRIIFORMES")
#   # Gruiformes <- c("GRUIFORMES")
#   
#   # Gruimorphae
#  # Gruae <- c("GRUIFORMES", "CHARADRIIFORMES", "OPISTHOCOMIFORMES")
#   Gruimorphae  <- c("GRUIFORMES", "CHARADRIIFORMES")
#   
#   # Add in prum clades.
#   Gruiformes <- c("GRUIFORMES")
#   Charadriifromes <- c("CHARADRIIFORMES")
#   # Create a new empty column.
#   data$higher_clade <- NA
#   
#   # Add the clade info.
#   data$higher_clade[data$order %in% Palaeognathae] <- "Palaeognathae"
#   data$higher_clade[data$order %in% Galloanseres] <- "Galloanseres"
#   data$higher_clade[data$order %in% Aequornithes] <- "Aequornithes"
#   data$higher_clade[data$order %in% Afroaves] <- "Afroaves"
#   data$higher_clade[data$order %in% Australaves] <- "Australaves"
#   data$higher_clade[data$order %in% Strisores] <- "Strisores"
#   data$higher_clade[data$order %in% Columbimorphae] <- "Columbimorphae"
#   data$higher_clade[data$order %in% Mirandornithes] <- "Mirandornithes"
#   data$higher_clade[data$order %in% Eurpygimorphae] <- "Eurpygimorphae"
#   data$higher_clade[data$order %in% Otidimorphae] <- "Otidimorphae"
#   data$higher_clade[data$order %in% Opisthocomiformes] <- "Opisthocomiformes"
#   data$higher_clade[data$order %in% Gruiformes] <- "Gruiformes"
#   data$higher_clade[data$order %in% Charadriifromes] <- "Charadriiformes"
#   
#   # Return the data frame.
#   return(data)
# }
# 



assign_prum_clades <- function(data){
  
  # Orders for PALAEOGNATHAE (in brain paper)
  Palaeognathae <- c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
                     "CASUARIIFORMES", "TINAMIFORMES")
  
  # Orders for Galloanserae (in brain paper)
  Galloanserae <- c("GALLIFORMES", "ANSERIFORMES")
  
  # Strisores. (in brain paper)
  Strisores <- c("APODIFORMES", "CAPRIMULGIFORMES")
  
  
  # Columbaves (in brain paper)
  Columbaves <- c("COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES",
                  "MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES")
  
  # Gruiformes.
  Gruiformes <- c("GRUIFORMES")
  
  # Orders for Aequorlitornithes (in brain paper)
  Aequorlitornithes <- c("GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", 
                         "CICONIIFORMES", "SULIFORMES", "PELECANIFORMES", 
                         "OPISTHOCOMIFORMES", "PHOENICOPTERIFORMES","PODICIPEDIFORMES",
                         "EURYPYGIFORMES", "PHAETHONTIFORMES", "CHARADRIIFORMES")
  
  # orphaned groups
  Opisthocomiformes <- c("OPISTHOCOMIFORMES")
  Accipitriformes <- c("ACCIPITRIFORMES")
  Strigiformes <- c("STRIGIFORMES")
  
  # Coraciimorphae (Core land birds)
  Coraciimorphae <- c("COLIIFORMES", "LEPTOSOMIFORMES",
                "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", "PICIFORMES")
  
  # Australaves (Austrailian radiation)
  Australaves <- c("CARIAMIFORMES",  "FALCONIFORMES", "PSITTACIFORMES", "PASSERIFORMES")

  # Create a new empty column.
  data$higher_clade <- NA
  
  # Add the clade info.
  data$higher_clade[data$order_bird_tree %in% Palaeognathae] <- "Palaeognathae"
  data$higher_clade[data$order_bird_tree %in% Galloanserae] <- "Galloanserae"
  data$higher_clade[data$order_bird_tree %in% Strisores] <- "Strisores"
  data$higher_clade[data$order_bird_tree %in% Columbaves] <- "Columbaves"
  data$higher_clade[data$order_bird_tree %in% Gruiformes] <- "Gruiformes"
  data$higher_clade[data$order_bird_tree %in% Aequorlitornithes] <- "Aequorlitornithes"
  data$higher_clade[data$order_bird_tree %in% Opisthocomiformes] <- "Opisthocomiformes"
  data$higher_clade[data$order_bird_tree %in% Accipitriformes] <- "Accipitriformes"
  data$higher_clade[data$order_bird_tree %in% Strigiformes] <- "Strigiformes"
  data$higher_clade[data$order_bird_tree %in% Coraciimorphae] <- "Coraciimorphae"
  data$higher_clade[data$order_bird_tree %in% Australaves] <- "Australaves"
  
  # Return the data frame.
  return(data)
}
















# 
# 
# 
# 
# 
# 
# 
# 
# grp <- list(# Orders for PALAEOGNATHAE (in brain paper)
#   Palaeognathae = c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
#                      "CASUARIIFORMES", "TINAMIFORMES"),
#   
#   # Orders for Galloanseres. (in brain paper)
#   Galloanseres = c("GALLIFORMES", "ANSERIFORMES"),
#   
#   # Orders for aequornithes (in brain paper)
#   Aequornithes = c("GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", "CICONIIFORMES", 
#                     "SULIFORMES", "PELECANIFORMES"),
#   
#   # Telluraves. (Core land birds)
#   Afroaves = c("ACCIPITRIFORMES", "STRIGIFORMES", "COLIIFORMES", "LEPTOSOMIFORMES",
#                   "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", "PICIFORMES"),
#   
#   # Australaves (Austrailian radiation)
#   Australaves = c("CARIAMIFORMES",  "FALCONIFORMES", "PSITTACIFORMES", "PASSERIFORMES"),
#   
#   # Strisores. (in brain paper)
#   Strisores = c("APODIFORMES", "CAPRIMULGIFORMES"),
#   
#   # Columbimorphae (in brain paper)
#   Columbimorphae = c("COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES"),
#   # Mirandornithes
#   Mirandornithes = c("PHOENICOPTERIFORMES","PODICIPEDIFORMES"),
#   
#   # Eurpygimorphae
#   Eurpygimorphae = c("EURYPYGIFORMES", "PHAETHONTIFORMES"),
#   
#   # Otidimorphae (in brain paper)
#   Otidimorphae = c("MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES"),
#   
#   # Gruimorphae
#   Gruae = c("GRUIFORMES", "CHARADRIIFORMES", "OPISTHOCOMIFORMES"))
# 
# 
# 
# 
# 




###############################################################################
                             #### END ####
###############################################################################


# Remove clades from environment.
 # rm("Palaeognathae", "Galloanseres", "Aequornithes", "Telluraves", "Australaves", 
 #            "Strisores", "Columbimorphae", "Mirandornithes", "Eurpygimorphae", 
 #            "Otidimorphae", "Charadriiformes", "Gruiformes")



# # Leftover neonaves
# c("APODIFORMES", "CAPRIMULGIFORMES", "CHARADRIIFORMES", "COLUMBIFORMES", 
#   "CUCULIFORMES", "EURYPYGIFORMES", "GRUIFORMES", "MESITORNITHIFORMES", 
#   "MUSOPHAGIFORMES", "OTIDIFORMES", "PHAETHONTIFORMES", "PHOENICOPTERIFORMES", 
#   "PODICIPEDIFORMES", "PTEROCLIDIFORMES")

