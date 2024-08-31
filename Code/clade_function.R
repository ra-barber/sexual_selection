###############################################################################
                       ##### Clade function #####
################################################################################

# A function to add clades for plotting on a tree.
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


###############################################################################
                             #### END ####
###############################################################################