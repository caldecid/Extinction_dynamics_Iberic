
# EDGE2 and simulation ----------------------------------------------------

library(tidyverse)
library(phytools)
devtools::install_github("iramosgutierrez/rEDGE")
library(rEDGE)

### Reading mega plant phylogeny
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

######### Peninsula ###########
peninsula_iucn <- read_csv("Data/Processed/peninsula_iucn.csv") 

#unifying species name
peninsula_iucn$species <- str_replace_all(peninsula_iucn$species, " ", "_")

##matching names
common_species_pen <- intersect(plant_phylo$tip.label, peninsula_iucn$species)

##prunning branches
peninsula_phylo <- keep.tip(plant_phylo, tip = common_species_pen)

# saving and reading the phylogeny as ape format
peninsula_phylo <- read.tree(text = write.tree(peninsula_phylo))

#forcing ultrametric
peninsula_phylo <- force.ultrametric(peninsula_phylo)

#fully bifurcating
peninsula_phylo <- multi2di(peninsula_phylo)

#reordering 
peninsula_phylo <- reorder.phylo(peninsula_phylo, order = "cladewise")

##transforming df
peninsula_iucn <- peninsula_iucn %>% select(species, Status) %>% 
         rename(RL.cat = Status) %>% 
          filter(species %in% common_species_pen)

## Regionally extinct as Extinct
peninsula_iucn$RL.cat[peninsula_iucn$RL.cat == "RE"] <- "EX"

peninsula_iucn_clean <- data.frame(
        species = as.character(peninsula_iucn$species),
          RL.cat  = as.character(peninsula_iucn$RL.cat))

#### applying EDGE2
EDGE2_Peninsula <- calculate_EDGE2(tree = peninsula_phylo,
                                   table = peninsula_iucn_clean,
                                   verbose = FALSE)

#saving
write_csv(EDGE2_Peninsula, file = "Data/Processed/EDGE2_Peninsula.csv")


############# Andalusia #################
andalusia_iucn <- read_csv("Data/Processed/andalucia_iucn.csv")

#unifying species name
andalusia_iucn$species <- str_replace_all(andalusia_iucn$species, " ", "_")

##matching names
common_species_and <- intersect(plant_phylo$tip.label, andalusia_iucn$species)

##prunning branches
andalusia_phylo <- keep.tip(plant_phylo, tip = common_species_and)

# saving and reading the phylogeny as ape format
andalusia_phylo <- read.tree(text = write.tree(andalusia_phylo))

#forcing ultrametric
andalusia_phylo <- force.ultrametric(andalusia_phylo)

#fully bifurcating
andalusia_phylo <- multi2di(andalusia_phylo)

#reordering 
andalusia_phylo <- reorder.phylo(andalusia_phylo, order = "cladewise")

##transforming df
andalusia_iucn <- andalusia_iucn %>% select(species, Status) %>% 
  rename(RL.cat = Status) %>% 
  filter(species %in% common_species_and)

#cleaning
andalusia_iucn_clean <- data.frame(
  species = as.character(andalusia_iucn$species),
  RL.cat  = as.character(andalusia_iucn$RL.cat))

#factor
andalusia_iucn_clean$RL.cat <- factor(andalusia_iucn$RL.cat,
                                      levels = unique(andalusia_iucn$RL.cat))

##there are some duplicated species
andalusia_iucn_clean <- andalusia_iucn_clean[
  !duplicated(andalusia_iucn_clean$species), ]

#### applying EDGE2
EDGE2_andalusia <- calculate_EDGE2(tree = andalusia_phylo,
                                   table = andalusia_iucn_clean,
                                   verbose = FALSE)

##saving
write_csv(EDGE2_andalusia, file = "Data/Processed/EDGE2_Andalusia.csv")
