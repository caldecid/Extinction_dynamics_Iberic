# Phylogenetic manipulations Iberian Flora 
library(phytools)
library(tidyverse)
library(remotes)
library(SpeciesAge) #install_github("thauffe/SpeciesAge")
library(picante)


## PENINSULA

### Reading peninsula iucn data
peninsula <- read_csv(file = "Data/Processed/peninsula_iucn.csv")

### String from the peninsula species name
peninsula_sp_names <- str_replace(peninsula$species, " ", "_")

### Eliminate subspecies
peninsula_without_sub <- unique(str_replace(peninsula_sp_names, " .*", ""))

### Reading mega plant phylogeny#############
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

### Eliminating species not present in phylogeny
n_pen <- peninsula_without_sub[peninsula_without_sub %in% plant_phylo$tip.label]

### Pruning tree preserving peninsula sp
peninsula_phy <- keep.tip(plant_phylo, n_pen)

### Saving tree
save(peninsula_phy, file = "Data/Processed/peninsula_phy.RData")

### Determining the proportion of threatened and non-threatened species
peninsula_genus <- peninsula %>% group_by(Genus) %>%
  summarize(
    richness = n(),
    threatened_species = sum(threat_category == "threatened"),
    proportion_threatened = threatened_species / richness
  )


## Implementing phylogenetic age analyses 

### Upscaling to the genus level
#### obtaining genera name
genera <- sapply(strsplit(peninsula_phy$tip.label, "_"), `[`, 1)

#### Create a mapping from genus to their corresponding species tips
genus_map <- split(peninsula_phy$tip.label, genera)

### Upscaling with the function
source("Scripts/functions.R")
genus_phylo <- collapse_to_genus(peninsula_phy, genus_map) # collapse_to_genus inside the function R.file

### Estimating phylogenetic age without correction
genus_age <- calculateTipAges(genus_phylo)

### Force in ultrametric
genus_phylo <- force.ultrametric(genus_phylo)

### Calculating evolutionary distinctiviness among genera
ev.distinc.peninsula <- evol.distinct(genus_phylo) %>%
  rename(Genus = Species,
         ev.dist = w)


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

names(list_ages_corrected) <- c("low_ex", "int_ex", "high_ex")

##collapsing list into a dataframe
peninsula_ages_rates <- bind_rows(list_ages_corrected,
                                  .id = "ext_fraction") %>%
                                  rename(Genus = tip)

##merging with peninsula_genus and evo distinctiveness
peninsula_ages_total <- left_join(peninsula_ages_rates, peninsula_genus,
                                  by = "Genus") %>%
  left_join(ev.distinc.peninsula,
            by = "Genus") %>% #calculating speciation rate
  mutate(lambda_a = log(richness)/age,
         lambda_mean = log(richness)/mean_age)


##saving
write_csv(peninsula_ages_total, file = "Data/Processed/peninsula_ages_total.csv")


##calculating Bloomberg's K and Pagel's lambda
##filtering data
data <- peninsula_genus[peninsula_genus$Genus %in% genus_phylo$tip.label, ]

# reorder data
data <- data[match(genus_phylo$tip.label, data$Genus), ]

# extracting the proportion of threatened species
threatened_vector <- data$proportion_threatened

names(threatened_vector) <- data$Genus

# Calculate Blomberg's K or Pagel's lambda
K_result <- phylosig(genus_phylo, threatened_vector,
                     nsim = 100, test = TRUE, method = "K")

lambda_result <- phylosig(genus_phylo,
                          threatened_vector, nsim = 100,
                          test = TRUE, method = "lambda")


#########################ANDALUCIA##############################################


##reading andalucia iucn data
andalucia <- read_csv(file = "Data/Processed/andalucia_iucn.csv")

##string from the andalucia species name
andalucia_sp_names <- str_replace(andalucia$species, " ", "_")

## eliminate subspecies
andalucia_without_sub <- unique(str_replace(andalucia_sp_names, " .*", ""))

##reading mega plant phylogeny#############
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

##eliminating species not present in phylogeny
n_and <- andalucia_without_sub[andalucia_without_sub %in% plant_phylo$tip.label]

##pruning tree preserving andalucia sp
andalucia_phy <- keep.tip(plant_phylo, n_and)

##saving
save(andalucia_phy, file = "Data/Processed/andalucia_phy.RData")


##determining the proportion of threatened and non-threatened species
andalucia_genus <- andalucia %>% group_by(Genus) %>%
  summarize(
    richness = n(),
    threatened_species = sum(threat_category == "threatened"),
    proportion_threatened = threatened_species / richness
  )


# implementing phylogenetic age analyses ----------------------------------

##upscaling to the genus level
#obtaining genera name
genera <- sapply(strsplit(andalucia_phy$tip.label, "_"), `[`, 1)

# Create a mapping from genus to their corresponding species tips
genus_map <- split(andalucia_phy$tip.label, genera)

##upscaling with the function
genus_phylo <- collapse_to_genus(andalucia_phy, genus_map)

##estimating phylogenetic age without correction
genus_age <- calculateTipAges(genus_phylo)

##transforming in ultrametric
genus_phylo <- force.ultrametric(genus_phylo)

##calculating evolutionary distinctiviness among genera
ev.distinc.andalucia <- evol.distinct(genus_phylo) %>%
  rename(Genus = Species,
         ev.dist = w)


########## calculating speciation and extinction rates########

##extinction fraction vector
ep.v <- c(0.1, 0.5, 0.9)


##empty dataframe
div.rates <- data.frame(ep.fr = ep.v,
                        rho.fr = rep(0.85, 3),
                        lambda = c(0,0,0),
                        mu = c(0,0,0))


##for loop estimating diversification rates
for (i in 1:nrow(div.rates)) {
  
  div.rates[i, 3:4] <- estimateBD(genus_phylo, epsilon = div.rates[i, 1],
                                  rho = div.rates[i, 2])
  
}

##correcting species ages

#empty list for storing information
list_ages_corrected <- vector(mode = "list", length = nrow(div.rates))


##for loop to estimate meanAge for each genus in each list
for (i in 1:nrow(div.rates)) {
  
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

names(list_ages_corrected) <- c("low_ex", "int_ex", "high_ex")

##collapsing list into a dataframe
andalucia_ages_rates <- bind_rows(list_ages_corrected,
                                  .id = "ext_fraction") %>%
                                 rename(Genus = tip)

##merging with andalucia_genus and evo distinctiveness
andalucia_ages_total <- left_join(andalucia_ages_rates, andalucia_genus,
                                  by = "Genus") %>%
  left_join(ev.distinc.andalucia,
            by = "Genus") %>% #calculating speciation rate
  mutate(lambda_a = log(richness)/age,
         lambda_mean = log(richness)/mean_age)


##saving
write_csv(andalucia_ages_total, file = "Data/Processed/andalucia_ages_total.csv")


##calculating Bloomberg's K and Pagel's lambda
##filtering data
data <- andalucia_genus[andalucia_genus$Genus %in% genus_phylo$tip.label, ]

# reorder data
data <- data[match(genus_phylo$tip.label, data$Genus), ]

# extracting the proportion of threatened species
threatened_vector <- data$proportion_threatened

names(threatened_vector) <- data$Genus

# Calculate Blomberg's K and Pagel's lambda
K_result_and <- phylosig(genus_phylo, threatened_vector,
                         nsim = 100, test = TRUE, method="K")

lambda_result_and <- phylosig(genus_phylo,
                          threatened_vector, nsim = 100,
                          test = TRUE, method="lambda")
