# Phylogenetic diversity loss analyses ------------------------------------

# Load required packages
library(ape)
library(phytools)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyverse)
devtools::install_github("iramosgutierrez/rEDGE")
library(rEDGE)
source("Scripts/functions.R")

### Reading mega plant phylogeny
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

######### Proportion of threatened species #############################

# Peninsula ---------------------------------------------------------------

#iucn data
peninsula_iucn <- read_csv("Data/Processed/peninsula_iucn.csv") 

#proportion data
peninsula_ages_proportion <- read_csv("Data/Processed/peninsula_ages_proportion.csv") %>% 
  filter(ext_fraction == "low_ex") %>% 
  select(Genus, proportion_threatened) 

#joining
peninsula_iucn <- left_join(peninsula_iucn, peninsula_ages_proportion,
                            by = "Genus") %>% 
  drop_na() %>% 
  rename("genus" = "Genus")


#unifying species name
peninsula_iucn$species <- str_replace_all(peninsula_iucn$species, " ", "_")

##matching names
common_species_pen <- intersect(plant_phylo$tip.label, peninsula_iucn$species)

peninsula_iucn <- peninsula_iucn %>% filter(species %in% common_species_pen)

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

##calling the 'calculate_pd_curve_prop' function
pd_curves_peninsula <- replicate(100, calculate_PD_curve_prop(tree = peninsula_phylo,
                                                              df = peninsula_iucn))

##calculating the observed area under the curve (AUC)
mean_pd_curve_peninsula <- rowMeans(pd_curves_peninsula)
observed_auc_peninsula <- mean(colSums(pd_curves_peninsula))

##generating the null test
set.seed(13)

# Generate null PD curves
null_pd_curves_peninsula <- replicate(99, {
  df <- peninsula_iucn
  # Randomize the threatened proportions
  df$proportion_threatened <- sample(df$proportion_threatened)
  # Recalculate PD curve and AUC
  null_pd_curve <- calculate_PD_curve_prop(tree = peninsula_phylo,
                                           df = df)
  
})

#null auc
null_auc_peninsula <- colSums(null_pd_curves_peninsula)

##p value
p_value_peninsula_prob <- mean(null_auc_peninsula <= observed_auc_peninsula)


#mean null auc
mean_null_auc_peninsula <- mean(null_auc_peninsula)

#confidence interval
ci_null_auc_peninsula <- quantile(null_auc_peninsula, probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_peninsula_df <- data.frame(
  step = 1:length(mean_pd_curve_peninsula),
  PD = mean_pd_curve_peninsula)

#saving
write_csv(pd_curves_peninsula_df,
          file = "Data/Processed/pd_curves_df_peninsula_prop.csv")

#read
pd_curves_df_peninsula_prop <- read_csv("Data/Processed/pd_curves_df_peninsula_prop.csv")

# null PD curves for shading 
null_pd_summary_peninsula <- cbind(step = 1:nrow(null_pd_curves_peninsula),
                                   null_pd_curves_peninsula)



colnames(null_pd_summary_peninsula) <- c("step", paste0("Sim", 1:99))  # Name columns

##as dataframe
null_pd_summary_peninsula <- as.data.frame(null_pd_summary_peninsula)

##saving
write_csv(null_pd_summary_peninsula,
          file = "Data/Processed/null_pd_curves_peninsula_prop.csv")

#reading
null_pd_curves_peninsula_prop <- read_csv("Data/Processed/null_pd_curves_peninsula_prop.csv")


# Reshape to long format for ggplot
long_df_peninsula_prop <- null_pd_curves_peninsula_prop %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_Peninsula_PD_prop.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_peninsula_prop_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_peninsula_prop,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_peninsula_prop,
            aes(x = step, y = PD), color = "blue", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = NULL,
    y = "Phylogenetic Diversity (PD)",
    title = "Peninsula Phylogenetic Diversity loss",
    subtitle = "Based on the Proportion of threatened sp."
  ) +
  theme_classic() +
  mynamestheme

pd_peninsula_prop_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_peninsula_df <- data.frame(null = null_auc_peninsula)

# Create the histogram


svg("Figures/Figure_Peninsula_AUC_prop.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_peninsula_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_peninsula), color = "blue",
             linetype = "solid", size = 1.5) +
  labs(x = "AUC", y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


# andalusia ---------------------------------------------------------------

#iucn data
andalusia_iucn <- read_csv("Data/Processed/andalucia_iucn.csv") 

#proportion data
andalusia_ages_proportion <- read_csv("Data/Processed/andalucia_ages_proportion.csv") %>% 
  filter(ext_fraction == "low_ex") %>% 
  select(Genus, proportion_threatened) 

#joining
andalusia_iucn <- left_join(andalusia_iucn, andalusia_ages_proportion,
                            by = "Genus") %>% 
  drop_na() %>% 
  rename("genus" = "Genus")


#unifying species name
andalusia_iucn$species <- str_replace_all(andalusia_iucn$species, " ", "_")

##matching names
common_species_andalusia <- intersect(plant_phylo$tip.label, andalusia_iucn$species)

andalusia_iucn <- andalusia_iucn %>% filter(species %in% common_species_andalusia)

##prunning branches
andalusia_phylo <- keep.tip(plant_phylo, tip = common_species_andalusia)

# saving and reading the phylogeny as ape format
andalusia_phylo <- read.tree(text = write.tree(andalusia_phylo))

#forcing ultrametric
andalusia_phylo <- force.ultrametric(andalusia_phylo)

#fully bifurcating
andalusia_phylo <- multi2di(andalusia_phylo)

#reordering 
andalusia_phylo <- reorder.phylo(andalusia_phylo, order = "cladewise")

##calling the 'calculate_pd_curve_prop' function
pd_curves_andalusia <- replicate(100, calculate_PD_curve_prop(tree = andalusia_phylo,
                                                              df = andalusia_iucn))

##calculating the observed area under the curve (AUC)
mean_pd_curve_andalusia <- rowMeans(pd_curves_andalusia)
observed_auc_andalusia <- mean(colSums(pd_curves_andalusia))

##generating the null test
set.seed(13)

# Generate null PD curves
null_pd_curves_andalusia <- replicate(99, {
  df <- andalusia_iucn
  # Randomize the threatened proportions
  df$proportion_threatened <- sample(df$proportion_threatened)
  # Recalculate PD curve and AUC
  null_pd_curve <- calculate_PD_curve_prop(tree = andalusia_phylo,
                                           df = df)
  
})

#null auc
null_auc_andalusia <- colSums(null_pd_curves_andalusia)

##p value
p_value_andalusia_prob <- mean(null_auc_andalusia <= observed_auc_andalusia)


#mean null auc
mean_null_auc_andalusia <- mean(null_auc_andalusia)

#confidence interval
ci_null_auc_andalusia <- quantile(null_auc_andalusia, probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_andalusia_df <- data.frame(
  step = 1:length(mean_pd_curve_andalusia),
  PD = mean_pd_curve_andalusia)

#saving
write_csv(pd_curves_andalusia_df,
          file = "Data/Processed/pd_curves_df_andalusia_prop.csv")

##reading
pd_curves_df_andalusia_prop <- read_csv("Data/Processed/pd_curves_df_andalusia_prop.csv")

# null PD curves for shading 
null_pd_summary_andalusia <- cbind(step = 1:nrow(null_pd_curves_andalusia),
                                   null_pd_curves_andalusia)



colnames(null_pd_summary_andalusia) <- c("step", paste0("Sim", 1:99))  # Name columns

##as dataframe
null_pd_summary_andalusia <- as.data.frame(null_pd_summary_andalusia)

##saving
write_csv(null_pd_summary_andalusia,
          file = "Data/Processed/null_pd_curves_andalusia_prop.csv")

#reading
null_pd_curves_andalusia_prop <- read_csv("Data/Processed/null_pd_curves_andalusia_prop.csv")


# Reshape to long format for ggplot
long_df_andalusia_prop <- null_pd_curves_andalusia_prop %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_andalusia_PD_prop.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_andalusia_prop_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_andalusia_prop,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_andalusia_prop,
            aes(x = step, y = PD), color = "orange", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = NULL,
    y = "Phylogenetic Diversity (PD)",
    title = "Proportion of threatened sp."
  ) +
  theme_classic() +
  mynamestheme

pd_andalusia_prop_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_andalusia_df <- data.frame(null = null_auc_andalusia)

# Create the histogram

svg("Figures/Figure_andalusia_AUC_prop.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_andalusia_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_andalusia),
             color = "orange",
             linetype = "solid", size = 1.5) +
  labs(x = "AUC", y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


# Using EDGE metrics ------------------------------------------------------

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

#reading
EDGE2_Peninsula <- read_csv(file = "Data/Processed/EDGE2_Peninsula.csv")

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

## reading
EDGE2_andalusia <- read_csv(file = "Data/Processed/EDGE2_Andalusia.csv")


############## Implementing NULL models #######################

# Peninsula ---------------------------------------------------------------

########### probability of extinction ####################

##calling the 'calculate_pd_curve_EDGE2' function, using probability of extinction (pext)
pd_curves_peninsula_pext <- replicate(100, calculate_PD_curve_EDGE2(tree = peninsula_phylo,
                                                                    df = EDGE2_Peninsula))

##calculating the observed area under the curve (AUC)
mean_pd_curve_peninsula_pext <- rowMeans(pd_curves_peninsula_pext)
observed_auc_peninsula_pext <- mean(colSums(pd_curves_peninsula_pext))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_peninsula_pext <- lapply(1:99, function(i) {
  
  df_null <- EDGE2_Peninsula
  
  # Shuffle extinction probabilities
  df_null$pext <- sample(df_null$pext)
  
  # Recalculate PD curve
  pd_curve_null <- calculate_PD_curve_EDGE2(
    tree = peninsula_phylo,
    df = df_null,
    ranking = "sum_pext"
  )
  
  return(pd_curve_null)
})


#null auc
null_auc_peninsula_pext <- sapply(null_pd_curves_peninsula_pext,
                                  function(curve) {
                                    sum(curve, na.rm = TRUE)
                                  })

p_value_pext <- mean(null_auc_peninsula_pext <= observed_auc_peninsula_pext)


#mean null auc
mean_null_auc_peninsula_pext <- mean(null_auc_peninsula_pext)

#confidence interval
ci_null_auc_peninsula_pext <- quantile(null_auc_peninsula_pext,
                                       probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_peninsula_pext <- data.frame(
  step = 1:length(mean_pd_curve_peninsula_pext),
  PD = mean_pd_curve_peninsula_pext)

##saving
write_csv(pd_curves_df_peninsula_pext, file = "Data/Processed/pd_curves_df_peninsula_pext.csv")

#reading
pd_curves_df_peninsula_pext <- read_csv("Data/Processed/pd_curves_df_peninsula_pext.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_peninsula_pext, length))

null_pd_padded_pext <- lapply(null_pd_curves_peninsula_pext, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_peninsula_pext <- as.data.frame(null_pd_padded_pext)

colnames(null_pd_summary_peninsula_pext) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_peninsula_pext$step <- 1:nrow(null_pd_summary_peninsula_pext)

##saving
write_csv(null_pd_summary_peninsula_pext,
          file = "Data/Processed/null_pd_curves_peninsula_pext.csv")

##reading
null_pd_summary_peninsula_pext <- read_csv("Data/Processed/null_pd_curves_peninsula.csv")

# Reshape to long format for ggplot
long_df_peninsula_pext <- null_pd_summary_peninsula_pext %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_Peninsula_PD_pext.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_peninsula_pext_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_peninsula_pext,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_peninsula_pext,
            aes(x = step, y = PD), color = "blue", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = "Removed genera",
    y = "Phylogenetic Diversity (PD)",
    #title = "Peninsula Phylogenetic Diversity loss",
    title = "Based on extinction probability"
  ) +
  theme_classic() +
  mynamestheme

pd_peninsula_pext_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_peninsula_pext_df <- data.frame(null = null_auc_peninsula_pext)

# Create the histogram

svg("Figures/Figure_Peninsula_AUC_pext.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_peninsula_pext_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_peninsula_pext), color = "blue",
             linetype = "solid", size = 1.5) +
  labs(y = "Frequency", x = "AUC") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


################## EDGE ####################

##calling the 'calculate_pd_curve_EDGE2' function, using EDGE
pd_curves_peninsula_edge <- replicate(100, calculate_PD_curve_EDGE2(tree = peninsula_phylo,
                                                                    df = EDGE2_Peninsula,
                                                                    ranking = "sum_EDGE"))

##calculating the observed area under the curve (AUC)
mean_pd_curve_peninsula_edge <- rowMeans(pd_curves_peninsula_edge)
observed_auc_peninsula_edge <- mean(colSums(pd_curves_peninsula_edge))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_peninsula_edge <- lapply(1:99, function(i) {
  
  df_null <- EDGE2_Peninsula
  
  # Shuffle extinction probabilities
  df_null$EDGE <- sample(df_null$EDGE)
  
  # Recalculate PD curve
  pd_curve_null <- calculate_PD_curve_EDGE2(
    tree = peninsula_phylo,
    df = df_null,
    ranking = "sum_EDGE"
  )
  
  return(pd_curve_null)
})


#null auc
null_auc_peninsula_edge <- sapply(null_pd_curves_peninsula_edge, function(curve) {
  sum(curve, na.rm = TRUE)
})

p_value_peninsula_edge <- mean(null_auc_peninsula_edge <= observed_auc_peninsula_edge)



#mean null auc
mean_null_auc_peninsula_edge <- mean(null_auc_peninsula_edge)

#confidence interval
ci_null_auc_peninsula_edge <- quantile(null_auc_peninsula_edge, 
                                       probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_peninsula_edge <- data.frame(
  step = 1:length(mean_pd_curve_peninsula_edge),
  PD = mean_pd_curve_peninsula_edge)

#saving
write_csv(pd_curves_df_peninsula_edge,
          file = "Data/Processed/pd_curves_df_peninsula_edge.csv")

#reading
pd_curves_df_peninsula_edge <- read_csv("Data/Processed/pd_curves_df_peninsula_edge.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_peninsula_edge, length))

null_pd_padded_edge <- lapply(null_pd_curves_peninsula_edge, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_peninsula_edge <- as.data.frame(null_pd_padded_edge)

colnames(null_pd_summary_peninsula_edge) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_peninsula_edge$step <- 1:nrow(null_pd_summary_peninsula_edge)

##saving
write_csv(null_pd_summary_peninsula_edge,
          file = "Data/Processed/null_pd_curves_peninsula_edge.csv")

#reading
null_pd_summary_peninsula_edge <- read_csv("Data/Processed/null_pd_curves_peninsula_edge.csv")

# Reshape to long format for ggplot
long_df_peninsula_edge <- null_pd_summary_peninsula_edge %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_Peninsula_PD_EDGE.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_peninsula_EDGE_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_peninsula_edge,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_peninsula_edge,
            aes(x = step, y = PD), color = "blue", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = NULL,
    y = "Phylogenetic Diversity (PD)",
    #title = "Peninsula Phylogenetic Diversity loss",
    title = "Based on the metric EDGE"
  ) +
  theme_classic() +
  mynamestheme

pd_peninsula_EDGE_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_peninsula_edge_df <- data.frame(null = null_auc_peninsula_edge)

# Create the histogram

svg("Figures/Figure_Peninsula_AUC_EDGE.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_peninsula_edge_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_peninsula_edge), color = "blue",
             linetype = "solid", size = 1.5) +
  labs(y = "Frequency", x = "AUC") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


# andalusia ---------------------------------------------------------------

########### probability of extinction ####################

#creating genus column
EDGE2_andalusia_pext <- EDGE2_andalusia_pext%>%
  mutate(genus = word(species, 1, sep = "_"))

##calling the 'calculate_pd_curve_EDGE2' function, using probability of extinction (pext)
pd_curves_andalusia_pext <- replicate(100, calculate_PD_curve_EDGE2(tree = andalusia_phylo,
                                                                    df = EDGE2_andalusia))

##calculating the observed area under the curve (AUC)
mean_pd_curve_andalusia_pext <- rowMeans(pd_curves_andalusia_pext)
observed_auc_andalusia_pext <- mean(colSums(pd_curves_andalusia_pext))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_andalusia_pext <- lapply(1:99, function(i) {
  
  df_null <- EDGE2_andalusia
  
  # Shuffle extinction probabilities
  df_null$pext <- sample(df_null$pext)
  
  # Recalculate PD curve
  pd_curve_null <- calculate_PD_curve_EDGE2(
    tree = andalusia_phylo,
    df = df_null,
    ranking = "sum_pext"
  )
  
  return(pd_curve_null)
})


#null auc
null_auc_andalusia_pext <- sapply(null_pd_curves_andalusia_pext, function(curve) {
  sum(curve, na.rm = TRUE)
})

p_value_pext <- mean(null_auc_andalusia_pext <= observed_auc_andalusia_pext)

#mean null auc
mean_null_auc_andalusia_pext <- mean(null_auc_andalusia_pext)

#confidence interval
ci_null_auc_andalusia_pext <- quantile(null_auc_andalusia_pext,
                                       probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_andalusia_pext <- data.frame(
  step = 1:length(mean_pd_curve_andalusia_pext),
  PD = mean_pd_curve_andalusia_pext)

#saving
write_csv(pd_curves_df_andalusia_pext, file = "Data/Processed/pd_curves_df_andalusia_pext.csv")

#reading
pd_curves_df_andalusia_pext <- read_csv("Data/Processed/pd_curves_df_andalusia_pext.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_andalusia_pext, length))

null_pd_padded_pext <- lapply(null_pd_curves_andalusia_pext, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_andalusia_pext <- as.data.frame(null_pd_padded_pext)

colnames(null_pd_summary_andalusia_pext) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_andalusia_pext$step <- 1:nrow(null_pd_summary_andalusia_pext)

##saving
write_csv(null_pd_summary_andalusia_pext,
          file = "Data/Processed/null_pd_curves_andalusia_pext.csv")

#reading
null_pd_curves_andalusia_pext <- read_csv("Data/processed/null_pd_curves_andalusia_pext.csv")

# Reshape to long format for ggplot
long_df_andalusia_pext <- null_pd_curves_andalusia_pext %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_andalusia_PD_pext.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_andalusia_pext_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_andalusia_pext,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_andalusia_pext,
            aes(x = step, y = PD), color = "orange", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = "Removed Genera",
    y = NULL,
    #title = "Andalusia Phylogenetic Diversity loss",
    title = "Extinction probability"
  ) +
  theme_classic() +
  mynamestheme

pd_andalusia_pext_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_andalusia_pext_df <- data.frame(null = null_auc_andalusia_pext)

# Create the histogram

svg("Figures/Figure_andalusia_AUC_pext.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_andalusia_pext_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_andalusia_pext),
             color = "orange",
             linetype = "solid", size = 1.5) +
  labs(y = "Frequency", x = "AUC") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


########################### EDGE ####################

##calling the 'calculate_pd_curve_EDGE2' function, using EDGE
pd_curves_andalusia_edge <- replicate(100,
                                      calculate_PD_curve_EDGE2(tree = andalusia_phylo,
                                                               df = EDGE2_andalusia,
                                                               ranking = "sum_EDGE"))

##calculating the observed area under the curve (AUC)
mean_pd_curve_andalusia_edge <- rowMeans(pd_curves_andalusia_edge)
observed_auc_andalusia_edge <- mean(colSums(pd_curves_andalusia_edge))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_andalusia_edge <- lapply(1:99, function(i) {
  
  df_null <- EDGE2_andalusia
  
  # Shuffle extinction probabilities
  df_null$EDGE <- sample(df_null$EDGE)
  
  # Recalculate PD curve
  pd_curve_null <- calculate_PD_curve_EDGE2(
    tree = andalusia_phylo,
    df = df_null,
    ranking = "sum_EDGE"
  )
  
  return(pd_curve_null)
})


#null auc
null_auc_andalusia_edge <- sapply(null_pd_curves_andalusia_edge,
                                  function(curve) {
                                  sum(curve, na.rm = TRUE)
                                                     })

p_value_andalusia_edge <- mean(null_auc_andalusia_edge <= 
                                             observed_auc_andalusia_edge)



#mean null auc
mean_null_auc_andalusia_edge <- mean(null_auc_andalusia_edge)

#confidence interval
ci_null_auc_andalusia_edge <- quantile(null_auc_andalusia_edge, 
                                       probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_andalusia_edge <- data.frame(
  step = 1:length(mean_pd_curve_andalusia_edge),
  PD = mean_pd_curve_andalusia_edge)

#saving
write_csv(pd_curves_df_andalusia_edge,
          file = "Data/Processed/pd_curves_df_andalusia_edge.csv")

#reading
pd_curved_df_andalusia_edge <- read_csv("Data/Processed/pd_curves_df_andalusia_edge.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_andalusia_edge, length))

null_pd_padded_edge <- lapply(null_pd_curves_andalusia_edge, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_andalusia_edge <- as.data.frame(null_pd_padded_edge)

colnames(null_pd_summary_andalusia_edge) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_andalusia_edge$step <- 1:nrow(null_pd_summary_andalusia_edge)

##saving
write_csv(null_pd_summary_andalusia_edge,
          file = "Data/Processed/null_pd_curves_andalusia_edge.csv")

##reading
null_pd_summary_andalusia_edge <- read_csv("Data/Processed/null_pd_curves_andalusia_edge.csv")

# Reshape to long format for ggplot
long_df_andalusia_edge <- null_pd_summary_andalusia_edge %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_andalusia_PD_EDGE.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_andalusia_EDGE_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_andalusia_edge,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_andalusia_edge,
            aes(x = step, y = PD), color = "orange", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = NULL,
    y = NULL,
    title = "EDGE metric"
  ) +
  theme_classic() +
  mynamestheme

pd_andalusia_EDGE_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_andalusia_edge_df <- data.frame(null = null_auc_andalusia_edge)

# Create the histogram

svg("Figures/Figure_andalusia_AUC_EDGE.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_andalusia_edge_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_andalusia_edge), color = "orange",
             linetype = "solid", size = 1.5) +
  labs(y = "Frequency", x = "AUC") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


# Combining plots ---------------------------------------------------------

theme_shared <- theme_classic() +
  mynamestheme +
  theme(
    plot.title = element_text(size = 13, face = "italic"),
    axis.title = element_text(size = 13)
  )

pd_peninsula_prop_plot <- pd_peninsula_prop_plot + theme_shared
pd_peninsula_pext_plot <- pd_peninsula_pext_plot + theme_shared
pd_peninsula_EDGE_plot <- pd_peninsula_EDGE_plot + theme_shared

pd_andalusia_prop_plot <- pd_andalusia_prop_plot + theme_shared
pd_andalusia_pext_plot <- pd_andalusia_pext_plot + theme_shared
pd_andalusia_EDGE_plot <- pd_andalusia_EDGE_plot + theme_shared


# Create row titles
title_peninsula <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Peninsular Spain",
           hjust = 0.5, size = 5, fontface = "bold", family = "serif") +
  theme_void()

title_andalusia <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Eastern Andalusia",
           hjust = 0.5, size = 5, fontface = "bold", family = "serif") +
  theme_void()


combined_plot <-
  (title_peninsula | plot_spacer() | plot_spacer()) /
  (pd_peninsula_prop_plot + pd_peninsula_pext_plot + pd_peninsula_EDGE_plot) /
  (title_andalusia | plot_spacer() | plot_spacer()) /
  (pd_andalusia_prop_plot + pd_andalusia_pext_plot + pd_andalusia_EDGE_plot) +
  plot_layout(heights = c(0.08, 1, 0.08, 1)) 

#saving
svg("Figures/Figure_PD_loss_regions.svg",
    width = 20/2.54,   # convert cm → inches
    height = 16/2.54)

print(combined_plot)

dev.off()

