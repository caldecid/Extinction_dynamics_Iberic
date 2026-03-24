# Phylogenetic manipulations Iberian flora --------------------------------

library(phytools)
library(tidyverse)
library(remotes)
library(SpeciesAge) #install_github("thauffe/SpeciesAge")
library(picante)

### Reading mega plant phylogeny#############
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

###collapsing the phylogeny to genera

# obtaining genera name
genera <- sapply(strsplit(plant_phylo$tip.label, "_"), `[`, 1)

# Create a mapping from genus to their corresponding species tips
genus_map <- split(plant_phylo$tip.label, genera)

### Upscaling with the function
source("Scripts/functions.R")
genus_phylo <- collapse_to_genus(plant_phylo, genus_map) # collapse_to_genus inside the function R.file

### Force in ultrametric
genus_phylo <- force.ultrametric(genus_phylo)

##resolving polytomies
genus_phylo <- multi2di(genus_phylo)

### Estimating phylogenetic age without correction
genus_age <- calculateTipAges(genus_phylo)

##saving genus_phylo
save(genus_phylo, file = "Data/Processed/plant_genus_phylo.RData")

# calculating speciation and extinction rates

## Extinction fraction vector
ep.v <- c(0.1, 0.5, 0.9)


## Empty dataframe
div.rates <- data.frame(ep.fr = ep.v,
                        rho.fr = rep(0.85, 3),
                        lambda = c(0,0,0),
                        mu = c(0,0,0))


##for loop estimating diversification rates
n <- nrow(div.rates)
for (i in 1:n) {
  
  div.rates[i, 3:4] <- estimateBD(genus_phylo, epsilon = div.rates[i, 1],
                                  rho = div.rates[i, 2])
  
}

##correcting species ages

#empty list for storing information
list_ages_corrected <- vector(mode = "list", length = n)


##for loop to estimate meanAge for each genus in each list
for (i in 1:n) {
  
  list_ages_corrected[[i]] <- data.frame(genus_age,
                                         mean_age = vector(mode = "numeric",
                                                           length = nrow(genus_age)))  
  for (j in 1:nrow(genus_age)) {
    
    list_ages_corrected[[i]][j,3] <- meanAge(genus_age[j,2],
                                             lambda = div.rates[i, 3],
                                             mu = div.rates[i,4],
                                             rho = div.rates[i,2])
    
  }
  
}

#naming lists according to extinction
names(list_ages_corrected) <- c("low_ex", "int_ex", "high_ex")

##collapsing list into a dataframe
plant_ages_corrected <- bind_rows(list_ages_corrected,
                                  .id = "ext_fraction") %>%
                         rename(Genus = tip)

##saving
write_csv(plant_ages_corrected,
          file = "Data/Processed/plant_ages_corrected.csv")


# Peninsula ---------------------------------------------------------------

### Reading peninsula iucn data
peninsula <- read_csv(file = "Data/Processed/peninsula_iucn.csv")


### Determining the proportion of threatened and non-threatened species
peninsula_genus <- peninsula %>% group_by(Genus) %>%
  summarize(
    richness = n(),
    threatened_species = sum(threat_category == "threatened"),
    proportion_threatened = threatened_species / richness
  )


## filtering in the genus age dataset only the genera present in the Peninsula
peninsula_ages <- plant_ages_corrected %>% 
  filter(Genus %in% peninsula_genus$Genus)

## left joining peninsula ages and proportion data
peninsula_ages_proportion <- left_join(peninsula_ages, peninsula_genus,
                                       by = "Genus")

##saving
write_csv(peninsula_ages_proportion, 
          file = "Data/Processed/peninsula_ages_proportion.csv")


##calculating Bloomberg's K and Pagel's lambda

##filtering data
peninsula_signal <- peninsula_genus[peninsula_genus$Genus %in% genus_phylo$tip.label, ]

# reorder data
peninsula_signal <- peninsula_signal[match(genus_phylo$tip.label,
                                           peninsula_signal$Genus), ]

# extracting the proportion of threatened species
threatened_peninsula <- peninsula_signal$proportion_threatened

names(threatened_peninsula) <- peninsula_signal$Genus

# Calculate Blomberg's K and Pagel's lambda
K_peninsula <- phylosig(genus_phylo, threatened_peninsula,
                        nsim = 100, test = TRUE, method="K")

lambda_peninsula <- phylosig(genus_phylo,
                             threatened_peninsula, nsim = 100,
                             test = TRUE, method="lambda")

# Andalusia ---------------------------------------------------------------

##reading andalucia iucn data
andalucia <- read_csv(file = "Data/Processed/andalucia_iucn.csv")


### Determining the proportion of threatened and non-threatened species
andalucia_genus <- andalucia %>% group_by(Genus) %>%
  summarize(
    richness = n(),
    threatened_species = sum(threat_category == "threatened"),
    proportion_threatened = threatened_species / richness
  )

## filtering in the genus age dataset only the genera present in the Peninsula
andalucia_ages <- plant_ages_corrected %>% 
  filter(Genus %in% andalucia_genus$Genus)

## left joining peninsula ages and proportion data
andalucia_ages_proportion <- left_join(andalucia_ages, andalucia_genus,
                                       by = "Genus")

#saving dataset
write_csv(andalucia_ages_proportion,
          file = "Data/Processed/andalucia_ages_proportion.csv")


##calculating Bloomberg's K and Pagel's lambda

##filtering data
andalucia_signal <- andalucia_genus[andalucia_genus$Genus %in% genus_phylo$tip.label, ]

# reorder data
andalucia_signal <- andalucia_signal[match(genus_phylo$tip.label,
                                           andalucia_signal$Genus), ]

# extracting the proportion of threatened species
threatened_andalucia <- andalucia_signal$proportion_threatened

names(threatened_andalucia) <- andalucia_signal$Genus

# Calculate Blomberg's K and Pagel's lambda
K_andalucia <- phylosig(genus_phylo, threatened_andalucia,
                        nsim = 100, test = TRUE, method="K")

lambda_andalucia <- phylosig(genus_phylo,
                             threatened_andalucia, nsim = 100,
                             test = TRUE, method="lambda")
