
# phylogenetic diversity loss analyses ------------------------------------

# Load required packages
library(ape)
library(phytools)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyverse)
source("Scripts/functions.R")


# Peninsula ---------------------------------------------------------------


##loading Peninsula tree
load("Data/Processed/peninsula_phy.RData")

##genera manipulation

#obtaining genera name
genera <- sapply(strsplit(peninsula_phy$tip.label, "_"), `[`, 1)

# Create a mapping from genus to their corresponding species tips
genus_map <- split(peninsula_phy$tip.label, genera)

##peninsula genera tree
peninsula_genera_phy <- collapse_to_genus(peninsula_phy, genus_map = genus_map)

##loading dataframe
peninsula_genera <- read_csv("Data/Processed/peninsula_ages_total.csv")

##filtering only one extinction scenario
peninsula_gen_df <- peninsula_genera %>% filter(ext_fraction == "low_ex")

##calling the 'calculate_pd_curve' function
pd_curves <- replicate(100, calculate_PD_curve(tree = peninsula_genera_phy,
                                               df = peninsula_genera))

##calculating the observed area under the curve (AUC)
mean_pd_curve <- rowMeans(pd_curves)
observed_auc <- mean(colSums(pd_curves))

##generating the null test
set.seed(13)

# Generate null PD curves
null_pd_curves <- replicate(99, {
  df <- peninsula_gen_df
  # Randomize the threatened proportions
  df$proportion_threatened <- sample(df$proportion_threatened)
  # Recalculate PD curve and AUC
  null_pd_curve <- calculate_PD_curve(tree = peninsula_genera_phy,
                                      df = df)
  
})

#null auc
null_auc <- colSums(null_pd_curves)

p_value <- sum(abs(null_auc) >= abs(observed_auc)) / length(null_auc)



# Prepare data for plotting the PD curves

pd_curves_df <- data.frame(
  step = 1:length(mean_pd_curve),
  PD = mean_pd_curve)

write_csv(pd_curves_df, file = "Data/Processed/pd_curves_df_peninsula.csv")

# null PD curves for shading 
null_pd_summary_peninsula <- cbind(step = 1:nrow(null_pd_curves),
                                   null_pd_curves)



colnames(null_pd_summary_peninsula) <- c("step", paste0("Sim", 1:99))  # Name columns

##as dataframe
null_pd_summary_peninsula <- as.data.frame(null_pd_summary_peninsula)

##saving
write_csv(null_pd_summary_peninsula,
          file = "Data/Processed/null_pd_curves_peninsula.csv")

# Subset of simulations to plot (25 sim)
subset_columns <- c("step", paste0("Sim", sample(1:99, 99))) 

subset_df_pen <- null_pd_summary_peninsula[, subset_columns]

# Reshape to long format for ggplot
long_df_peninsula <- subset_df_pen %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
png("Figures/Figure_Peninsula_PD.png", 
    width = 17, height = 15,
    units = "cm", pointsize = 8, res = 300)

pd_peninsula_plot <- ggplot() +
  # Shaded area for the range of null PD curves
  geom_line(data = long_df_peninsula,
            aes(x = step, y = PD), color = "gray", 
            size = 0.5, alpha = 0.5) +
  # Observed mean PD curve
  geom_line(data = pd_curves_df,
            aes(x = step, y = PD), color = "blue", size = 1.2) +
  # Overlay some null PD curves for illustration
  labs(
    x = NULL,
    y = "Phylogenetic Diversity (PD)",
    title = "Peninsula Phylogenetic Diversity loss",
  ) +
  theme_classic() +
  mynamestheme
pd_peninsula_plot
dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_df <- data.frame(null = null_auc)

# Create the histogram

png("Figures/Figure_Peninsula_AUC.png", 
    width = 12, height = 10,
    units = "cm", pointsize = 8, res = 300)


ggplot(null_auc_df, aes(x = null)) +
  geom_histogram( fill = "lightgray", color = "gray") +
  geom_vline(aes(xintercept = observed_auc), color = "blue",
             linetype = "solid", size = 1.5) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  mynamestheme
  

dev.off()


# Andalucia ---------------------------------------------------------------

##loading Andalucia tree
load("Data/Processed/andalucia_phy.RData")

#obtaining genera name
genera <- sapply(strsplit(andalucia_phy$tip.label, "_"), `[`, 1)

# Create a mapping from genus to their corresponding species tips
genus_map <- split(andalucia_phy$tip.label, genera)

##peninsula genera tree
andalucia_genera_phy <- collapse_to_genus(andalucia_phy, 
                                          genus_map = genus_map)

##loading dataframe
andalucia_genera <- read_csv("Data/Processed/andalucia_ages_total.csv")

##filtering only one extinction scenario
andalucia_gen_df <- andalucia_genera %>% filter(ext_fraction == "low_ex")

##calling the 'calculate_pd_curve' function
pd_curves_andalucia <- replicate(100,
                                 calculate_PD_curve(tree = andalucia_genera_phy,
                                               df = andalucia_genera))

##calculating the observed area under the curve (AUC)
mean_pd_curve_andalucia <- rowMeans(pd_curves_andalucia)
observed_auc_andalucia <- mean(colSums(pd_curves_andalucia))

##generating the null test
set.seed(14)

# Generate null PD curves
null_pd_curves_andalucia <- replicate(99, {
  df <- andalucia_gen_df
  # Randomize the threatened proportions
  df$proportion_threatened <- sample(df$proportion_threatened)
  # Recalculate PD curve and AUC
  null_pd_curve <- calculate_PD_curve(tree = andalucia_genera_phy,
                                      df = df)
  
})

#null auc
null_auc_andalucia <- colSums(null_pd_curves_andalucia)

p_value_andalucia <- sum(abs(null_auc_andalucia) >= abs(observed_auc_andalucia)) / length(null_auc_andalucia)



# Prepare data for plotting the PD curves

pd_curves_df_andalucia <- data.frame(
  step = 1:length(mean_pd_curve_andalucia),
  PD = mean_pd_curve_andalucia)

##saving
write_csv(pd_curves_df_andalucia, 
          file = "Data/Processed/pd_curves_df_andalucia.csv")

# Calculate the range of null PD curves for shading
null_pd_summary_andalucia <- cbind(step = 1:nrow(null_pd_curves_andalucia),
                                   null_pd_curves_andalucia)

colnames(null_pd_summary_andalucia) <- c("step", paste0("Sim", 1:99))  # Name columns

null_pd_summary_andalucia <- as.data.frame(null_pd_summary_andalucia)

##saving
write_csv(null_pd_summary_andalucia, 
          file = "Data/Processed/null_pd_curves_andalucia.csv")

# Subset of simulations to plot (25 sim)
subset_columns <- c("step", paste0("Sim", sample(1:99, 99))) 

subset_df_and <- null_pd_summary_andalucia[, subset_columns]

# Reshape to long format for ggplot
long_df_and <- as.data.frame(subset_df_and) %>%
  pivot_longer(cols = -step, names_to = "Simulation", values_to = "PD")



# Plot the PD curves
png("Figures/Figure_andalucia_PD.png", 
    width = 17, height = 15,
    units = "cm", pointsize = 8, res = 300)

pd_andalucia_plot <- ggplot() +
  # mean observed pd curve
  geom_line(data = pd_curves_df_andalucia,
            aes(x = step, y = PD), color = "red", alpha = 1, size = 1.2) +
  # null pd curves
  geom_line(data = long_df_and,
            aes(x = step, y = PD), color = "gray", size = 0.5,
            alpha = 0.5) +
  
  labs(
    x = NULL,
    y = "Phylogenetic diversity (PD)",
    title = "Andalusia Phylogenetic Diversity loss",
  ) +
  theme_classic()+
  mynamestheme
pd_andalucia_plot

dev.off()


##plotting AUC

# Convert the vector to a data frame
null_auc_df_andalucia <- data.frame(null = null_auc_andalucia)

# Create the histogram

png("Figures/Figure_andalucia_AUC.png", 
    width = 12, height = 10,
    units = "cm", pointsize = 8, res = 300)


ggplot(null_auc_df_andalucia, aes(x = null)) +
  geom_histogram( fill = "lightgray", 
                  color = "gray") +
  geom_vline(aes(xintercept = observed_auc_andalucia), color = "red",
             linetype = "solid", size = 1.5) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  mynamestheme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

dev.off()

#####uniting PD curves
bottom = text_grob("Genera remove",
                   family = "serif", size = 13, 
                   face = "bold")


png("Figures/Figure_PD_curves.png", 
    width = 30, height = 15,
    units = "cm", pointsize = 8, res = 300)

grid.arrange(pd_peninsula_plot, pd_andalucia_plot,
             ncol = 2)

dev.off()

