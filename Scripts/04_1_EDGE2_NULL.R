
# EDGE2 and simulation ----------------------------------------------------

library(tidyverse)
library(phytools)
devtools::install_github("iramosgutierrez/rEDGE")
library(rEDGE)
source("Scripts/functions.R")


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
pd_curves_peninsula <- replicate(100, calculate_PD_curve_EDGE2(tree = peninsula_phylo,
                                                               df = EDGE2_Peninsula))

##calculating the observed area under the curve (AUC)
mean_pd_curve_peninsula <- rowMeans(pd_curves_peninsula)
observed_auc_peninsula <- mean(colSums(pd_curves_peninsula))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_peninsula <- lapply(1:99, function(i) {
  
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
null_auc_peninsula <- sapply(null_pd_curves_peninsula, function(curve) {
  sum(curve, na.rm = TRUE)
})

p_value <- mean(null_auc_peninsula <= observed_auc_peninsula)



#mean null auc
mean_null_auc_peninsula <- mean(null_auc_peninsula)

#confidence interval
ci_null_auc_peninsula <- quantile(null_auc_peninsula, probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_peninsula <- data.frame(
  step = 1:length(mean_pd_curve_peninsula),
  PD = mean_pd_curve_peninsula)

##saving
write_csv(pd_curves_df_peninsula, file = "Data/Processed/pd_curves_df_peninsula.csv")

#reading
pd_curves_df_peninsula <- read_csv("Data/Processed/pd_curves_df_peninsula.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_peninsula, length))

null_pd_padded <- lapply(null_pd_curves_peninsula, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_peninsula <- as.data.frame(null_pd_padded)

colnames(null_pd_summary_peninsula) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_peninsula$step <- 1:nrow(null_pd_summary_peninsula)

##saving
write_csv(null_pd_summary_peninsula,
          file = "Data/Processed/null_pd_curves_peninsula.csv")

##reading
null_pd_summary_peninsula <- read_csv("Data/Processed/null_pd_curves_peninsula.csv")

# Reshape to long format for ggplot
long_df_peninsula <- null_pd_summary_peninsula %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
svg("Figures/Figure_Peninsula_PD_pext.svg",
    width = 14/2.54,
    height = 11/2.54)

pd_peninsula_pext_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_peninsula,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df_peninsula,
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
null_auc_peninsula_df <- data.frame(null = null_auc_peninsula)

# Create the histogram

svg("Figures/Figure_Peninsula_AUC_pext.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_peninsula_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_peninsula), color = "blue",
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
EDGE2_andalusia <- EDGE2_andalusia%>%
                   mutate(genus = word(species, 1, sep = "_"))

##calling the 'calculate_pd_curve_EDGE2' function, using probability of extinction (pext)
pd_curves_andalusia <- replicate(100, calculate_PD_curve_EDGE2(tree = andalusia_phylo,
                                                               df = EDGE2_andalusia))

##calculating the observed area under the curve (AUC)
mean_pd_curve_andalusia <- rowMeans(pd_curves_andalusia)
observed_auc_andalusia <- mean(colSums(pd_curves_andalusia))

##null pd curves

##generating the null test
set.seed(13)

null_pd_curves_andalusia <- lapply(1:99, function(i) {
  
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
null_auc_andalusia <- sapply(null_pd_curves_andalusia, function(curve) {
  sum(curve, na.rm = TRUE)
})

p_value <- mean(null_auc_andalusia <= observed_auc_andalusia)



#mean null auc
mean_null_auc_andalusia <- mean(null_auc_andalusia)

#confidence interval
ci_null_auc_andalusia <- quantile(null_auc_andalusia, probs = c(0.05, 0.975))

# Prepare data for plotting the PD curves

pd_curves_df_andalusia <- data.frame(
  step = 1:length(mean_pd_curve_andalusia),
  PD = mean_pd_curve_andalusia)

#saving
write_csv(pd_curves_df_andalusia, file = "Data/Processed/pd_curves_df_andalusia_pext.csv")

#reading
pd_curves_df_andalusia_pext <- read_csv("Data/Processed/pd_curves_df_andalusia_pext.csv")

# null PD curves for shading 
max_len <- max(sapply(null_pd_curves_andalusia, length))

null_pd_padded <- lapply(null_pd_curves_andalusia, function(curve) {
  length(curve) <- max_len
  return(curve)
})

#as df
null_pd_summary_andalusia <- as.data.frame(null_pd_padded)

colnames(null_pd_summary_andalusia) <- paste0("Sim", 1:99)  # Name columns

null_pd_summary_andalusia$step <- 1:nrow(null_pd_summary_andalusia)

##saving
write_csv(null_pd_summary_andalusia,
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
null_auc_andalusia_df <- data.frame(null = null_auc_andalusia)

# Create the histogram

svg("Figures/Figure_andalusia_AUC_pext.svg",
    width = 12/2.54,
    height = 10/2.54)


ggplot(null_auc_andalusia_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc_andalusia), color = "orange",
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
pd_curves_andalusia_edge <- replicate(100, calculate_PD_curve_EDGE2(tree = andalusia_phylo,
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
null_auc_andalusia_edge <- sapply(null_pd_curves_andalusia_edge, function(curve) {
  sum(curve, na.rm = TRUE)
})

p_value_andalusia_edge <- mean(null_auc_andalusia_edge <= observed_auc_andalusia_edge)



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
