
# PGLS as sensitivity analysis -------------------------------

library(tidyverse)
library(caper)
library(phytools)
library(writexl)

###############################################################
#########      Prop. of threatened sp  ########################
###############################################################


########################## Peninsula ########################################

#reading dataframe with variables 
peninsula_total <- read_csv("Data/Processed/peninsula_merged_iucn_clads.csv")

# phylogenetic tree
plant_genus_phylo <- read.tree("Data/Processed/plant_genus_phylo.tre")

####### high extinction
peninsula_high <- peninsula_total %>%
               filter(ext_fraction == "high_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
peninsula_high_pgls <- peninsula_high %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_pen_high <- drop.tip(plant_genus_phylo,
                          setdiff(plant_genus_phylo$tip.label,
                                  peninsula_high_pgls$Genus))
                        

# Make sure the order of the data matches the tree
peninsula_high_pgls <- peninsula_high_pgls %>%
                        slice(match(tree_pen_high$tip.label, Genus))

# checking if everything matches
all(peninsula_high_pgls$Genus == tree_pen_high$tip.label)


# checking node labels
head(tree_pen_high$tip.label)

# Check internal node labels
head(tree_pen_high$node.label)

# strings in the nodes
tree_pen_high$node.label <- NULL

# Check
tree_pen_high$node.label

#creating comparative data
comp_pen_high <- comparative.data(phy = tree_pen_high,
                                    data = as.data.frame(peninsula_high_pgls),
                                    names.col = "Genus",
                                    vcv = TRUE,
                                    na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_pen_high_pgls <- pgls(
                        proportion_threatened ~ mean_age + rates,
                        data = comp_pen_high,
                        lambda = "ML")

#calling summary
x_pen_high <- summary(model_pen_high_pgls)

#transforming coefficients in df
coef_pen_high <- as.data.frame(x_pen_high$coefficients)

coef_pen_high$term <- rownames(coef_pen_high)
rownames(coef_pen_high) <- NULL

#table for Peninsula high
coef_pen_high <- coef_pen_high %>%
  dplyr::mutate(
    Region = "Peninsular Spain",
    Extinction_scenario = "High",
    lambda = x_pen_high$param.CI$lambda$opt,
    R_squared = x_pen_high$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared
  )

######## Intermediate

peninsula_int <- peninsula_total %>%
  filter(ext_fraction == "int_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
peninsula_int_pgls <- peninsula_int %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_pen_int <- drop.tip(plant_genus_phylo,
                          setdiff(plant_genus_phylo$tip.label,
                                  peninsula_int_pgls$Genus))


# Make sure the order of the data matches the tree
peninsula_int_pgls <- peninsula_int_pgls %>%
  slice(match(tree_pen_int$tip.label, Genus))

# checking if everything matches
all(peninsula_int_pgls$Genus == tree_pen_int$tip.label)


# checking node labels
head(tree_pen_int$tip.label)

# Check internal node labels
head(tree_pen_int$node.label)

# strings in the nodes
tree_pen_int$node.label <- NULL

# Check
tree_pen_int$node.label

#creating comparative data
comp_pen_int <- comparative.data(phy = tree_pen_int,
                                  data = as.data.frame(peninsula_int_pgls),
                                  names.col = "Genus",
                                  vcv = TRUE,
                                  na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_pen_int_pgls <- pgls(
  proportion_threatened ~ mean_age + rates,
  data = comp_pen_int,
  lambda = "ML")

#calling summary
x_pen_int <- summary(model_pen_int_pgls)

#transforming coefficients in df
coef_pen_int <- as.data.frame(x_pen_int$coefficients)

coef_pen_int$term <- rownames(coef_pen_int)
rownames(coef_pen_int) <- NULL

#table for Peninsula int
coef_pen_int <- coef_pen_int %>%
  dplyr::mutate(
    Region = "Peninsular Spain",
    Extinction_scenario = "Intermediate",
    lambda = x_pen_int$param.CI$lambda$opt,
    R_squared = x_pen_int$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared
  )


######## Low

peninsula_low <- peninsula_total %>%
  filter(ext_fraction == "low_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
peninsula_low_pgls <- peninsula_low %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_pen_low <- drop.tip(plant_genus_phylo,
                         setdiff(plant_genus_phylo$tip.label,
                                 peninsula_low_pgls$Genus))


# Make sure the order of the data matches the tree
peninsula_low_pgls <- peninsula_low_pgls %>%
  slice(match(tree_pen_low$tip.label, Genus))

# checking if everything matches
all(peninsula_low_pgls$Genus == tree_pen_low$tip.label)


# checking node labels
head(tree_pen_low$tip.label)

# Check lowernal node labels
head(tree_pen_low$node.label)

# strings in the nodes
tree_pen_low$node.label <- NULL

# Check
tree_pen_low$node.label

#creating comparative data
comp_pen_low <- comparative.data(phy = tree_pen_low,
                                 data = as.data.frame(peninsula_low_pgls),
                                 names.col = "Genus",
                                 vcv = TRUE,
                                 na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_pen_low_pgls <- pgls(
  proportion_threatened ~ mean_age + rates,
  data = comp_pen_low,
  lambda = "ML")

#calling summary
x_pen_low <- summary(model_pen_low_pgls)

#transforming coefficients in df
coef_pen_low <- as.data.frame(x_pen_low$coefficients)
coef_pen_low$term <- rownames(coef_pen_low)
rownames(coef_pen_low) <- NULL

#table for Peninsula low
coef_pen_low <- coef_pen_low %>%
  dplyr::mutate(
    Region = "Peninsular Spain",
    Extinction_scenario = "Low",
    lambda = x_pen_low$param.CI$lambda$opt,
    R_squared = x_pen_low$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared)
  
### joining all Peninsula coefficients
coef_peninsula_prop_threatened_pgls <- rbind(coef_pen_high, 
                                            coef_pen_int,
                                            coef_pen_low)


########################## Andalusia ########################################

#reading dataframe with variables 
Andalusia_total <- read_csv("Data/Processed/andalucia_merged_iucn_clads.csv")


####### high extinction
Andalusia_high <- Andalusia_total %>%
  filter(ext_fraction == "high_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
Andalusia_high_pgls <- Andalusia_high %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_andalusia_high <- drop.tip(plant_genus_phylo,
                          setdiff(plant_genus_phylo$tip.label,
                                  Andalusia_high_pgls$Genus))


# Make sure the order of the data matches the tree
Andalusia_high_pgls <- Andalusia_high_pgls %>%
  slice(match(tree_andalusia_high$tip.label, Genus))

# checking if everything matches
all(Andalusia_high_pgls$Genus == tree_andalusia_high$tip.label)


# checking node labels
head(tree_andalusia_high$tip.label)

# Check internal node labels
head(tree_andalusia_high$node.label)

# strings in the nodes
tree_andalusia_high$node.label <- NULL

# Check
tree_andalusia_high$node.label

#creating comparative data
comp_andalusia_high <- comparative.data(phy = tree_andalusia_high,
                                  data = as.data.frame(Andalusia_high_pgls),
                                  names.col = "Genus",
                                  vcv = TRUE,
                                  na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_andalusia_high_pgls <- pgls(
  proportion_threatened ~ mean_age + rates,
  data = comp_andalusia_high,
  lambda = "ML")

#calling summary
x_andalusia_high <- summary(model_andalusia_high_pgls)

#transforming coefficients in df
coef_andalusia_high <- as.data.frame(x_andalusia_high$coefficients)

coef_andalusia_high$term <- rownames(coef_andalusia_high)
rownames(coef_andalusia_high) <- NULL

#table for Andalusia high
coef_andalusia_high <- coef_andalusia_high %>%
  dplyr::mutate(
    Region = "Eastern Andalusia",
    Extinction_scenario = "High",
    lambda = x_andalusia_high$param.CI$lambda$opt,
    R_squared = x_andalusia_high$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared
  )

######## Intermediate

Andalusia_int <- Andalusia_total %>%
  filter(ext_fraction == "int_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
Andalusia_int_pgls <- Andalusia_int %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_andalusia_int <- drop.tip(plant_genus_phylo,
                         setdiff(plant_genus_phylo$tip.label,
                                 Andalusia_int_pgls$Genus))


# Make sure the order of the data matches the tree
Andalusia_int_pgls <- Andalusia_int_pgls %>%
  slice(match(tree_andalusia_int$tip.label, Genus))

# checking if everything matches
all(Andalusia_int_pgls$Genus == tree_andalusia_int$tip.label)


# checking node labels
head(tree_andalusia_int$tip.label)

# Check internal node labels
head(tree_andalusia_int$node.label)

# strings in the nodes
tree_andalusia_int$node.label <- NULL

# Check
tree_andalusia_int$node.label

#creating comparative data
comp_andalusia_int <- comparative.data(phy = tree_andalusia_int,
                                 data = as.data.frame(Andalusia_int_pgls),
                                 names.col = "Genus",
                                 vcv = TRUE,
                                 na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_andalusia_int_pgls <- pgls(
  proportion_threatened ~ mean_age + rates,
  data = comp_andalusia_int,
  lambda = "ML")

#calling summary
x_andalusia_int <- summary(model_andalusia_int_pgls)

#transforming coefficients in df
coef_andalusia_int <- as.data.frame(x_andalusia_int$coefficients)

coef_andalusia_int$term <- rownames(coef_andalusia_int)
rownames(coef_andalusia_int) <- NULL

#table for Andalusia int
coef_andalusia_int <- coef_andalusia_int %>%
  dplyr::mutate(
    Region = "Eastern Andalusia",
    Extinction_scenario = "Intermediate",
    lambda = x_andalusia_int$param.CI$lambda$opt,
    R_squared = x_andalusia_int$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared
  )


######## Low

Andalusia_low <- Andalusia_total %>%
  filter(ext_fraction == "low_ex")

# match datasets

# Keep only genera present in both the data and phylogeny
Andalusia_low_pgls <- Andalusia_low %>%
  filter(Genus %in% plant_genus_phylo$tip.label)

# Prune the phylogeny to the genera in the dataset
tree_andalusia_low <- drop.tip(plant_genus_phylo,
                         setdiff(plant_genus_phylo$tip.label,
                                 Andalusia_low_pgls$Genus))


# Make sure the order of the data matches the tree
Andalusia_low_pgls <- Andalusia_low_pgls %>%
  slice(match(tree_andalusia_low$tip.label, Genus))

# checking if everything matches
all(Andalusia_low_pgls$Genus == tree_andalusia_low$tip.label)


# checking node labels
head(tree_andalusia_low$tip.label)

# Check lowernal node labels
head(tree_andalusia_low$node.label)

# strings in the nodes
tree_andalusia_low$node.label <- NULL

# Check
tree_andalusia_low$node.label

#creating comparative data
comp_andalusia_low <- comparative.data(phy = tree_andalusia_low,
                                 data = as.data.frame(Andalusia_low_pgls),
                                 names.col = "Genus",
                                 vcv = TRUE,
                                 na.omit = FALSE)

# model (prop. threated ~ mean age + rates)                                  
model_andalusia_low_pgls <- pgls(
  proportion_threatened ~ mean_age + rates,
  data = comp_andalusia_low,
  lambda = "ML")

#calling summary
x_andalusia_low <- summary(model_andalusia_low_pgls)

#transforming coefficients in df
coef_andalusia_low <- as.data.frame(x_andalusia_low$coefficients)
coef_andalusia_low$term <- rownames(coef_andalusia_low)
rownames(coef_andalusia_low) <- NULL

#table for Andalusia low
coef_andalusia_low <- coef_andalusia_low %>%
  dplyr::mutate(
    Region = "Eastern Andalusia",
    Extinction_scenario = "Low",
    lambda = x_andalusia_low$param.CI$lambda$opt,
    R_squared = x_andalusia_low$r.squared
  ) %>%
  dplyr::select(
    Region,
    Extinction_scenario,
    term,
    Estimate,
    `Std. Error`,
    `t value`,
    `Pr(>|t|)`,
    lambda,
    R_squared)

### joining all Andalusia coefficients
coef_Andalusia_prop_threatened_pgls <- rbind(coef_andalusia_high, 
                                            coef_andalusia_int,
                                            coef_andalusia_low)

#joining pgls results
coef_prop_threatened_pgls <- rbind(coef_peninsula_prop_threatened_pgls,
                                   coef_Andalusia_prop_threatened_pgls)

## saving
write_xlsx(coef_prop_threatened_pgls,
                 "Data/Processed/Sensitivity/PGLS/coef_prop_threat_pgls.xlsx")

###############################################################
#########    Mean extinction probability ######################
###############################################################


########################## Peninsula ###############################

#reading
peninsula_total <- read_csv("Data/Processed/peninsula_merged_iucn_clads.csv")

##calling EDGE 
EDGE2_Peninsula <- read_csv(file = "Data/Processed/EDGE2_Peninsula.csv")

##Adding genus
EDGE2_Peninsula <- EDGE2_Peninsula %>% mutate(genus = str_extract(species, "^[^_]+"))

##grouping genus and obtainig the mean probability of extinction
EDGE2_peninsula_genera <- EDGE2_Peninsula %>% group_by(genus) %>% 
  summarise(mean_pext = mean(pext))

##left join
peninsula_total_pext <- peninsula_total %>% left_join(EDGE2_peninsula_genera,
                                          by = c("Genus" = "genus")) %>% 
                                           drop_na()
