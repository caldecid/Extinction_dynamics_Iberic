
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

####### adapting the previous function for accepting EDGE2 metrics ########

calculate_PD_curve_EDGE2 <- function(tree, df,
                                     ranking = "sum_pext",
                                     random_ties = TRUE) {
  
  require(dplyr)
  require(ape)
  
  # Copy objects (avoid overwriting originals)
  tree_sim <- tree
  df_sim <- df
  
  pd_curve <- c()
  
  while (nrow(df_sim) > 0) {
    
    # 1. Compute genus-level scores
    genus_scores <- df_sim %>%
      group_by(genus) %>%
      summarise(
        sum_pext = sum(pext, na.rm = TRUE),
        mean_pext = mean(pext, na.rm = TRUE),
        sum_EDGE = sum(EDGE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 2. Select ranking variable
    scores <- genus_scores[[ranking]]
    names(scores) <- genus_scores$genus
    
    max_val <- max(scores)
    candidates <- names(scores)[scores == max_val]
    
    if (length(candidates) > 1 && random_ties) {
      selected_genus <- sample(candidates, 1)
    } else {
      selected_genus <- candidates[1]
    }
    
    # 3. Identify species to remove
    species_to_remove <- df_sim$species[df_sim$genus == selected_genus]
    
    # 4. Prune tree
    tree_sim <- ape::drop.tip(tree_sim, species_to_remove)
    
    # 5. Compute PD
    pd_curve <- c(pd_curve, sum(tree_sim$edge.length))
    
    # 6. Update dataframe
    df_sim <- df_sim[df_sim$genus != selected_genus, ]
    
    # Stop if tree empty
    if (length(tree_sim$tip.label) == 0) break
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


##############Continuous phylogenetic signal
D_scale_cont <- function(ytree, xdata, 
                         max_group = 783,
                         by = 10,
                         min_group = 10) {
  mytree_scale <- ytree
  pos <- match(mytree_scale$tip.label,
               xdata$spp)
  xdata <- xdata[pos,]
  dist_tree <- as.matrix(cophenetic(mytree_scale))
  quantis <- quantile(dist_tree,  seq(0.1, 1, 0.1))
  mytree_scale_clus <-
    as.hclust.phylo(force.ultrametric(mytree_scale))
  grupos <- seq(min_group, max_group, by = by)
  n_groups <- length(grupos)
  groups <- matrix(nrow = Ntip(mytree_scale), ncol = n_groups)
  
  for (i in 1:n_groups) {
    groups[, i] <- paste0("taxa", 
                          cutree(mytree_scale_clus, grupos[i]))
  }
  colnames(groups) <- paste0("group", grupos)
  xdata <- cbind(xdata, groups)
  signals <- matrix(nrow = n_groups, ncol = 3)
  pb = txtProgressBar(min = 0, max = n_groups, initial = 0,style= 3) 
  
  for (i in 1:n_groups) {
    setTxtProgressBar(pb,i)
    
    category <- paste0("group", grupos[i])
    mytree_scale$tip.label <- xdata[, category]
    remove_tip <- which(duplicated(mytree_scale$tip.label))
    higher_tree <- drop.tip(mytree_scale, remove_tip)
    mydata_high <- aggregate_taxa(category = category, x = xdata)
    trait <- setNames(as.numeric(as.character(mydata_high[, 2])),
                      mydata_high[, 1])
    phyloK <- phylosig(higher_tree, trait, "K", nsim = 1, niter=1)
    phylo_lambda <- phylosig(higher_tree, trait, "lambda", nsim = 1, niter=1)
    signals[i,] <- c(grupos[i],
                     phylo_lambda$lambda,
                     phyloK)
  }
  close(pb)
  colnames(signals) <- c("N groups", "lambda", "K")
  return(signals)
}

### Used to aggregate taxa values for the previous function
aggregate_taxa <- function(x, category) {
  taxa <- as.character(unique(x[, category]))
  n <- length(taxa)
  effect_neg <- character(n)
  for (i in 1:n) {
    efeitos <- x[x[, category] == taxa[i], "effect"]
    effect_neg[i] <- mean(efeitos)
  }
  return(as.data.frame(cbind(taxa, effect_neg)))
}

################my theme###############################
## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      plot.title = element_text(family = "serif", size = (15),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (11)),
                      legend.background = element_rect(fill="gray90",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")
