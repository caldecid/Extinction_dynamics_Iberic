calculate_PD_curve_prop <- function(tree, df,
                                    random_ties = TRUE) {
  
  require(dplyr)
  require(ape)
  
  # Copy objects
  tree_sim <- tree
  df_sim <- df
  
  pd_curve <- c()
  
  # Loop until no genera left
  while (nrow(df_sim) > 0) {
    
    # 1. Compute genus-level proportion threatened
    genus_scores <- df_sim %>%
                     group_by(genus) %>%
                    summarise(proportion_threatened = first(proportion_threatened), 
                              .groups = "drop")
    
    # 2. Select genus with highest proportion
    max_val <- max(genus_scores$proportion_threatened)
    
    candidates <- genus_scores$genus[
      genus_scores$proportion_threatened == max_val
    ]
    
    if (length(candidates) > 1 && random_ties) {
      selected_genus <- sample(candidates, 1)
    } else {
      selected_genus <- candidates[1]
    }
    
    # 3. Get species of that genus
    species_to_remove <- df_sim$species[df_sim$genus == selected_genus]
    
    # 4. Prune from tree
    tree_sim <- ape::drop.tip(tree_sim, species_to_remove)
    
    # 5. Calculate PD
    pd_curve <- c(pd_curve, sum(tree_sim$edge.length))
    
    # 6. Update dataframe
    df_sim <- df_sim[df_sim$genus != selected_genus, ]
    
    # Stop if tree empty
    if (length(tree_sim$tip.label) == 0) break
  }
  
  return(pd_curve)
}

# PD loss (proportion of threatened species) ------------------------------

library(tidyverse)
library(phytools)
library(patchwork)
source("Scripts/functions.R")


### Reading mega plant phylogeny
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")


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
null_pd_curves_andalusia_prop <- read_csv( "Data/Processed/null_pd_curves_andalusia_prop.csv")


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
  geom_vline(aes(xintercept = observed_auc_andalusia), color = "orange",
             linetype = "solid", size = 1.5) +
  labs(x = "AUC", y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme


dev.off()


# Summary plot ------------------------------------------------------------

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
