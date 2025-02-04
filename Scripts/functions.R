
# Functions ---------------------------------------------------------------


# Phylogenetic manipulation -----------------------------------------------

# Function to keep only one representative per genus
collapse_to_genus <- function(phylo, genus_map) {
  for (genus in names(genus_map)) {
    species_tips <- genus_map[[genus]]
    if (length(species_tips) > 1) {
      # Drop all but one species in this genus
      phylo <- drop.tip(phylo, setdiff(species_tips, species_tips[1]))
    }
    # Rename the remaining species to the genus name
    phylo$tip.label[phylo$tip.label == species_tips[1]] <- genus
  }
  return(phylo)
}


# Phylogenetic diversity loss experiment ----------------------------------


# Function to calculate the PD curve and AUC#########
calculate_PD_curve <- function(tree, df) {
  pd_curve <- numeric(length(tree$tip.label))  # Placeholder for PD values
  for (i in seq_along(tree$tip.label)) {
    # 1. Select genus with the highest proportion of threatened species
    max_genera <- df %>%
      filter(proportion_threatened == max(proportion_threatened))
    
    selected_genus <- sample(max_genera$Genus, 1)
    
    # 2. Prune the selected genus from the tree
    tree <- drop.tip(tree, selected_genus)
    
    # 3. Calculate PD for the pruned tree
    pd_curve[i] <- sum(tree$edge.length)
    
    # Update the dataframe to exclude the pruned genus
    df <- df %>% filter(Genus != selected_genus)
  }
  return(pd_curve)
}


# Define a custom function to extract and format coefficients
tidy_rq_summary <- function(rq_sum) {
  
  require(dplyr)
  require(purrr)
  
  # Extract the coefficients matrix
  coef_mat <- rq_sum$coefficients
  # Convert it to a data frame and add a column for the term names
  coef_df <- as.data.frame(coef_mat)
  coef_df <- tibble::rownames_to_column(coef_df, var = "term")
  # Optionally, rename columns for clarity
  coef_df <- coef_df %>%
    rename(estimate = Value,
           std.error = `Std. Error`,
           t.value = `t value`,
           p.value = `Pr(>|t|)`) %>% 
    filter(term != "(Intercept)")
  return(coef_df)
}
