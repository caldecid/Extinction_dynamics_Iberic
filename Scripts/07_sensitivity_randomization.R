# Randomization of DD and NE for Sensitivity analyses --------------------------

## Libraries
library(tidyverse)
library(readxl)
library(writexl)
library(phytools)
library(remotes)
library(SpeciesAge) #install_github("thauffe/SpeciesAge")
library(picante)
library(ggcorrplot)
library(DescTools)
library(gridExtra)
library(purrr)
library(broom)
library(DescTools)
library(glmmTMB) #for using betabinomial model for overdispersion
library(DHARMa)
library(rEDGE)



# Peninsula  ----------------------------

## Iberic Peninsula data
path <- "Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls"
excel_sheets(path)

peninsula <- read_excel(path,
                        sheet = "Lista_Peninsula_Final.txt")

##species age and speciation rates data
pen_ages_rates <- read_csv("Data/Processed/peninsula_merged_iucn_clads.csv")

#extracting variables
pen_ages_rates <- pen_ages_rates %>% select(Genus, ext_fraction, mean_age,
                                            richness, rates)

### Ordering the conservation status for sensitivity analyses
not_threatened_ss <- c("LC", "NT")
threatened_ss <- c("DD", "NE", "VU", "EN", "CR", "EX", "RE", "EW")

#defining status and threat category
peninsula_ss <- peninsula %>%
  mutate(
    Status = factor(
      Status,
      labels = c(not_threatened_ss, threatened_ss),
      ordered = TRUE
    ),
    threat_category = case_when(
      Status %in% not_threatened_ss ~ "not_threatened",
      Status %in% threatened_ss ~ "threatened"
    )
  )

#facts
peninsula_ss$threat_category <- factor(peninsula_ss$threat_category)

peninsula_ss$Status <- factor(peninsula_ss$Status)

### Genus arrangements
peninsula_ss <- peninsula_ss %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")

### Determining the peninsula genus df
peninsula_genus <- peninsula_ss %>% group_by(Genus) %>%
  summarize(
    richness = n(),
    threatened_species = sum(threat_category == "threatened"),
    proportion_threatened = threatened_species / richness
  )

# Filter out NE and DD to calculate proportions
peninsula_assessed <- peninsula_ss %>%
  filter(!Status %in% c("DD", "NE"))

##calculating proportions
prop_threatened <- mean(peninsula_assessed$threat_category == "threatened")
prop_not_threatened <- 1 - prop_threatened


##randomizing the DD and NE species
n_replicates <- 100
sensitivity_list_pen <- vector("list", n_replicates)

set.seed(42)

for (i in 1:n_replicates) {
  sensitivity_list_pen[[i]] <- peninsula_ss %>%
    mutate(threat_category_adjusted = threat_category) %>%
    mutate(threat_category_adjusted = ifelse(
      Status %in% c("DD", "NE"),
      sample(c("threatened", "not_threatened"), size = sum(Status %in% c("DD", "NE")), 
             replace = TRUE, prob = c(prop_threatened, prop_not_threatened)),
      as.character(threat_category_adjusted)
    )) %>% group_by(Genus) %>%
    summarize(
      richness = n(),
      threatened_species = sum(threat_category_adjusted == "threatened"),
      proportion_threatened = threatened_species / richness
    ) %>% 
    left_join(pen_ages_rates, #merging ages and speciation rates
              by = "Genus") %>% 
    drop_na(ext_fraction)
}

#saving
save(sensitivity_list_pen, file = "Data/Processed/Sensitivity/sensitivity_list_peninsula.RData")

#loading
load("Data/Processed/Sensitivity/sensitivity_list_peninsula.RData")

###phylogenetic signal

# Load the genus_phylo obtained from 02_phylo_manipulation 
load("Data/Processed/plant_genus_phylo.RData")

# Create a vector of tip labels once
tip_labels <- genus_phylo$tip.label

# Initialize result lists
K_pen_list <- vector("list", length = 100)
lambda_pen_list <- vector("list", length = 100)

# Start loop
for (i in seq_len(100)) {
  message("Processing replicate ", i)
  
  df <- sensitivity_list_pen[[i]] %>%
    dplyr::filter(ext_fraction == "low_ex") 
  
  ##filtering data
  df_signal <- df[df$Genus %in% genus_phylo$tip.label, ]
  
  # reorder data
  df_signal <- df_signal[match(genus_phylo$tip.label,
                                            df_signal$Genus), ]
  
  # extracting the proportion of threatened species
  threatened_vector <- df_signal$proportion_threatened
  
  names(threatened_vector) <- df_signal$Genus
  
  
  # Try-catch in case phylosig fails
  K_pen_list[[i]] <- tryCatch(
    phylosig(genus_phylo, 
             threatened_vector, nsim = 100, test = TRUE, method = "K"),
    error = function(e) NA
  )
  
  lambda_pen_list[[i]] <- tryCatch(
    phylosig(genus_phylo, 
             threatened_vector, nsim = 100, test = TRUE, method = "lambda"),
    error = function(e) NA
  )
}


##converting to df

###Blomberg
K_pen_df <- purrr::map_dfr(K_pen_list, function(x) {
  if (is.list(x)) {
    tibble(K = x$K, P = x$P)
  } else {
    tibble(K = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(K_pen_df,
          file = "Data/Processed/Sensitivity/Randomizations/K_pen_df.csv")

#reading
K_pen_df <- read_csv("Data/Processed/Sensitivity/Randomizations/K_pen_df.csv")


###lambda
lambda_pen_df <- purrr::map_dfr(lambda_pen_list, function(x) {
  if (is.list(x)) {
    tibble(lambda = x$lambda, logL = x$logL, P = x$P)
  } else {
    tibble(lambda = NA, logL = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(lambda_pen_df,
          file = "Data/Processed/Sensitivity/Randomizations/lambda_pen_df.csv")

#reading
lambda_pen_df <- read_csv("Data/Processed/Sensitivity/Randomizations/lambda_pen_df.csv")


###merging and plotting
phylo_signal_pen <- dplyr::left_join(K_pen_df, lambda_pen_df,
                                     by = "replicate") %>% 
                   rename(p_value_K = P.x,
                          p_value_lambda = P.y)


#saving 
write_csv(phylo_signal_pen, "Results/phylo_signal_peninsula_random.csv")


##summarising
phylo_signal_pen_summary <- phylo_signal_pen %>%
  summarise(
    lambda_mean = mean(lambda),
    lambda_sd = sd(lambda),
    lambda_p_mean = mean(p_value_lambda),
    K_mean = mean(K),
    K_sd = sd(K),
    K_p_mean = mean(p_value_K)
  ) %>%
  tibble::tibble(
    parameter = c("lambda", "K"),
    mean = c(.$lambda_mean, .$K_mean),
    sd = c(.$lambda_sd, .$K_sd),
    p_value_mean = c(.$lambda_p_mean, .$K_p_mean)
  )

phylo_signal_pen_summary <- phylo_signal_pen_summary[,-c(1:6)]

phylo_signal_pen_summary$region <- "Peninsula"

#pivoting for plotting

long_pen_signal <- phylo_signal_pen %>%
  select(K, lambda) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         K = "Blomberg's K",
                         lambda = "Pagel's λ"))

# Define reference values from main analysis (DD and NE as Not Threatened)
ref_main <- data.frame(
  Metric = c("Blomberg's K", "Pagel's λ"),
  ref = c(0.027, 0.195)
)


# Plot
# Create axis limit values
xlims <- tibble(
  Metric = c("Blomberg's K", "Pagel's λ"),
  xmin = c(0, 0),
  xmax = c(0.05, 0.4)
)

# Merge with main data for facet-specific limits
long_pen_signal <- long_pen_signal %>%
  left_join(xlims, by = "Metric")


####plotting
png("Figures/Supplementary/Sensitivity/Phylo_signal_Peninsular.png", 
    width = 20, height = 12,
    units = "cm", pointsize = 8, res = 300)

ggplot(long_pen_signal, aes(x = Value)) +
  geom_histogram(fill = "skyblue", alpha = 0.6, bins = 30) +
  facet_wrap(~Metric, scales = "free") +
  
  # Add invisible points to enforce per-facet x-axis limits
  geom_blank(aes(x = xmin)) +
  geom_blank(aes(x = xmax)) +
  
  # Main analysis vertical line (red, dashed)
  geom_vline(data = ref_main, aes(xintercept = ref), 
             linetype = "dashed", size = 1.2, color = "red") +
  
  
  labs(x = NULL, y = NULL, title = "Phylogenetic Signal in\n Peninsular Spain") +
  theme_minimal(base_size = 14) +
  mynamestheme+
  theme(
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 15)
  ) 

dev.off()


## fitting the glm 

# List of extinction scenarios
ext_scenarios <- c("low_ex", "int_ex", "high_ex")

# Initialize empty list to collect results
all_results <- list()

# Loop over scenarios
for (scenario in ext_scenarios) {
  
  # Loop over replicates
  for (i in seq_along(sensitivity_list_pen)) {
    
    df <- sensitivity_list_pen[[i]]
    
    # Filter by extinction scenario
    df_sub <- df %>%
      filter(ext_fraction == scenario) 
    
    # Fit GLM, wrapped in try() to handle errors
    model <- try(glmmTMB(cbind(threatened_species,
                               richness.x - threatened_species) ~ mean_age + rates,
                         family = betabinomial(),
                         data = df_sub),
                 silent = TRUE)
    
    # Skip model if error
    if (inherits(model, "try-error")) next
    
    # Extract summary
    smry <- summary(model)
    
    # Extract coefficient table and convert to data.frame
    coef_df <- as.data.frame(smry$coefficients$cond)
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add replicate and scenario metadata
    coef_df$replicate <- i
    coef_df$scenario <- scenario
    
    # Append to results list
    all_results[[length(all_results) + 1]] <- coef_df
  }
}

####GLM Summaries for each extinction scenarios
glm_summary_pen_df <- bind_rows(all_results)

##factors
glm_summary_pen_df$scenario <- factor(glm_summary_pen_df$scenario,
                                  levels = c("low_ex",
                                             "int_ex",
                                             "high_ex"),
                                  labels = c("Low",
                                             "Intermediate",
                                             "High"),
                                  ordered = TRUE)

glm_summary_pen_df$term <- factor(glm_summary_pen_df$term,
                              levels = c("(Intercept)",
                                         "mean_age",
                                         "rates"),
                              labels = c("Intercept",
                                         "Corrected age",
                                         "Diversification rates"),
                              ordered = TRUE)

##GLM stats
glm_summary_pen_stats <- glm_summary_pen_df %>%
  group_by(term, scenario) %>%
  summarize(
    mean_estimate = mean(Estimate),
    sd_estimate = sd(Estimate),
    median_estimate = median(Estimate),
    mean_p = mean(`Pr(>|z|)`, na.rm = TRUE),
    prop_significant = mean(`Pr(>|z|)` < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

##saving
write_csv(glm_summary_pen_stats, file = "Results/glm_sensitivity_Peninsula.csv")

write_xlsx(glm_summary_pen_stats, path = "Results/glm_sensitivity_Peninsula_excel.xlsx")


####plotting
png("Figures/Supplementary/Sensitivity/glm_random_estimates_Peninsula.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

##plots
 glm_summary_pen_df %>%
  filter(term != "Intercept") %>%
  ggplot(aes(x = scenario, y = Estimate, fill = scenario)) +
  geom_boxplot() +
  facet_wrap(~ term, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Random assignment of DD & NE sp. [Peninsular Spain]",
       y = "Estimate of GLM Coefficients", x = "Extinction Scenario") +
  scale_fill_brewer(palette = "Set2")+
  mynamestheme+
  theme(legend.position = "none")

dev.off()
  


##proportion of significant terms
png("Figures/Supplementary/Sensitivity/glm_random_pvalue_Peninsula.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

glm_summary_pen_df %>%
  filter(term != "Intercept") %>%
  mutate(significant = `Pr(>|z|)` < 0.05) %>%
  group_by(term, scenario) %>%
  summarize(prop_sig = mean(significant), .groups = "drop") %>%
  ggplot(aes(x = scenario, y = prop_sig, fill = scenario)) +
  geom_col() +
  facet_wrap(~ term) +
  ylim(0, 1)+
  theme_minimal(base_size = 14) +
  labs(y = "Proportion of Significant Results (p < 0.05)",
       x = "Extinction scenario",
       title = "Random assignment of DD & NE sp. [Peninsular Spain]") +
  scale_fill_brewer(palette = "Set2") +
  mynamestheme+
  theme(legend.position = "none")

dev.off()


# Andalusia ---------------------------------------------------------------

andalucia <- read_excel(path,
                        sheet = "FLORANDOR")

##species age and speciation rates data
andalucia_ages_rates <- read_csv("Data/Processed/andalucia_merged_iucn_clads.csv")

#extracting variables
andalucia_ages_rates <- andalucia_ages_rates %>% select(Genus, ext_fraction, mean_age,
                                            richness, rates)

### Ordering the conservation status for sensitivity analyses
not_threatened_ss <- c("LC", "NT")
threatened_ss <- c("DD", "NE", "VU", "EN", "CR", "EX", "RE", "EW")


##organizing threatening status
andalucia_ss <- andalucia %>%
  mutate(
    Status = factor(
      Status,
      ordered = TRUE
    ),
    threat_category = case_when(
      Status %in% not_threatened_ss ~ "not_threatened",
      Status %in% threatened_ss ~ "threatened"
    )
  )

### Genus arrangements
andalucia_ss <- andalucia_ss  %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(
    col = Species,
    into = c("Genus", "Epithet"),
    sep = " ",
    extra = "merge"
  )

#factors
andalucia_ss$threat_category <- factor(andalucia_ss$threat_category)

andalucia_ss$Status <- factor(andalucia_ss$Status)


# Filter out NE and DD to calculate proportions
andalucia_assessed <- andalucia_ss %>%
  filter(!Status %in% c("DD", "NE"))

##calculating proportions
prop_threatened <- mean(andalucia_assessed$threat_category == "threatened")
prop_not_threatened <- 1 - prop_threatened


##randomizing the DD and NE species
n_replicates <- 100
sensitivity_list_andalucia <- vector("list", n_replicates)

set.seed(42)

for (i in 1:n_replicates) {
  sensitivity_list_andalucia[[i]] <- andalucia_ss %>%
    mutate(threat_category_adjusted = threat_category) %>%
    mutate(threat_category_adjusted = ifelse(
      Status %in% c("DD", "NE"),
      sample(c("threatened", "not_threatened"),
             size = sum(Status %in% c("DD", "NE")), 
             replace = TRUE, prob = c(prop_threatened, prop_not_threatened)),
      as.character(threat_category_adjusted)
    )) %>% group_by(Genus) %>%
    summarize(
      richness = n(),
      threatened_species = sum(threat_category_adjusted == "threatened"),
      proportion_threatened = threatened_species / richness
    ) %>% 
    left_join(andalucia_ages_rates, #merging ages and speciation rates
              by = "Genus") %>% 
    drop_na(ext_fraction)
}

##saving
save(sensitivity_list_andalucia, file = "Data/Processed/Sensitivity/sensitivity_list_andalucia.RData")
## loading
load("Data/Processed/Sensitivity/sensitivity_list_andalucia.RData")

##phylogenetic signal

##calling genus tree
# Load the genus_phylo obtained from 02_phylo_manipulation 
load("Data/Processed/plant_genus_phylo.RData")

# Create a vector of tip labels once
tip_labels <- genus_phylo$tip.label

# Initialize result lists
K_and_list <- vector("list", length = 100)
lambda_and_list <- vector("list", length = 100)

# Start loop
for (i in seq_len(100)) {
  message("Processing replicate ", i)
  
  df <- sensitivity_list_andalucia[[i]] %>%
    dplyr::filter(ext_fraction == "low_ex") 
  
  ##filtering data
  df_signal <- df[df$Genus %in% genus_phylo$tip.label, ]
  
  # reorder data
  df_signal <- df_signal[match(genus_phylo$tip.label,
                               df_signal$Genus), ]
  
  # extracting the proportion of threatened species
  threatened_vector <- df_signal$proportion_threatened
  
  names(threatened_vector) <- df_signal$Genus
  
  
  # Try-catch in case phylosig fails
  K_and_list[[i]] <- tryCatch(
    phylosig(genus_phylo, 
             threatened_vector, nsim = 100, test = TRUE, method = "K"),
    error = function(e) NA
  )
  
  lambda_and_list[[i]] <- tryCatch(
    phylosig(genus_phylo, 
             threatened_vector, nsim = 100, test = TRUE, method = "lambda"),
    error = function(e) NA
  )
}

##converting to df

#Blomberg
K_and_df <- purrr::map_dfr(K_and_list, function(x) {
  if (is.list(x)) {
    tibble(K = x$K, P = x$P)
  } else {
    tibble(K = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(K_and_df,
          file = "Data/Processed/Sensitivity/Randomizations/K_and_df.csv")

#lambda
lambda_and_df <- purrr::map_dfr(lambda_and_list, function(x) {
  if (is.list(x)) {
    tibble(lambda = x$lambda, logL = x$logL, P = x$P)
  } else {
    tibble(lambda = NA, logL = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(lambda_and_df,
          file = "Data/Processed/Sensitivity/Randomizations/lambda_and_df.csv")


###merging and plotting########
phylo_signal_and <- dplyr::left_join(K_and_df, lambda_and_df,
                                     by = "replicate") %>% 
  rename(p_value_K = P.x,
         p_value_lambda = P.y)

#saving 
write_csv(phylo_signal_and, "Results/phylo_signal_andalucia_random.csv")


##summarising
phylo_signal_and_summary <- phylo_signal_and %>%
  summarise(
    lambda_mean = mean(lambda),
    lambda_sd = sd(lambda),
    lambda_p_mean = mean(p_value_lambda),
    K_mean = mean(K),
    K_sd = sd(K),
    K_p_mean = mean(p_value_K)
  ) %>%
  tibble::tibble(
    parameter = c("lambda", "K"),
    mean = c(.$lambda_mean, .$K_mean),
    sd = c(.$lambda_sd, .$K_sd),
    p_value_mean = c(.$lambda_p_mean, .$K_p_mean)
  )

#removing extra columns
phylo_signal_and_summary <- phylo_signal_and_summary[,-c(1:6)]

#assigning regions
phylo_signal_and_summary$region <- "Andalusia"

##merging Phylo signal regions
phylo_signal_random_summary <- rbind(phylo_signal_pen_summary,
                                     phylo_signal_and_summary)

write_xlsx(phylo_signal_random_summary,
           path = "Results/phylo_signal_random_summary.xlsx")

#pivoting for plotting

long_and_signal <- phylo_signal_and %>%
  select(K, lambda) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         K = "Blomberg's K",
                         lambda = "Pagel's λ"))

# Define reference values from main analysis (DD and NE as Not Threatened)
ref_main_and <- data.frame(
  Metric = c("Blomberg's K", "Pagel's λ"),
  ref = c(0.043, 0.403)
)


####plotting
png("Figures/Supplementary/Sensitivity/Phylo_signal_Andalucia.png", 
    width = 20, height = 12,
    units = "cm", pointsize = 8, res = 300)

ggplot(long_and_signal, aes(x = Value)) +
  geom_histogram(fill = "skyblue", alpha = 0.6, bins = 30) +
  facet_wrap(~Metric, scales = "free") +
  
  
  # Main analysis vertical line (red, dashed)
  geom_vline(data = ref_main_and, aes(xintercept = ref), 
             linetype = "dashed", size = 1.2, color = "red") +
  
  
  labs(x = NULL, y = NULL, 
       title = "Phylogenetic Signal in\n Eastern Andalusia") +
  theme_minimal(base_size = 14) +
  mynamestheme+
  theme(
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 15)
  ) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.0001))


dev.off()

###### fitting the glm 

# List of extinction scenarios
ext_scenarios <- c("low_ex", "int_ex", "high_ex")

# Initialize empty list to collect results
all_results_andalucia <- list()

# Loop over scenarios
for (scenario in ext_scenarios) {
  
  # Loop over replicates
  for (i in seq_along(sensitivity_list_andalucia)) {
    
    df <- sensitivity_list_andalucia[[i]]
    
    # Filter by extinction scenario
    df_sub <- df %>%
      filter(ext_fraction == scenario) 
    
    model <- try(glmmTMB(cbind(threatened_species,
                      richness.x - threatened_species) ~ mean_age + rates,
                family = betabinomial(),
                data = df_sub),
                silent = TRUE)
    
    # Skip model if error
    if (inherits(model, "try-error")) next
    
    # Extract summary
    smry <- summary(model)
    
    # Extract coefficient table and convert to data.frame
    coef_df <- as.data.frame(smry$coefficients$cond)
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add replicate and scenario metadata
    coef_df$replicate <- i
    coef_df$scenario <- scenario
    
    # Append to results list
    all_results_andalucia[[length(all_results_andalucia) + 1]] <- coef_df
  }
}

####GLM Summaries for each extinction scenarios
glm_summary_df_andalucia <- bind_rows(all_results_andalucia)

##factors
glm_summary_df_andalucia$scenario <- factor(glm_summary_df_andalucia$scenario,
                                  levels = c("low_ex",
                                             "int_ex",
                                             "high_ex"),
                                  labels = c("Low",
                                             "Intermediate",
                                             "High"),
                                  ordered = TRUE)

glm_summary_df_andalucia$term <- factor(glm_summary_df_andalucia$term,
                              levels = c("(Intercept)",
                                         
                                         "mean_age",
                                         "rates"),
                              labels = c("Intercept",
                                         
                                         "Corrected age",
                                         "Diversification rates"),
                              ordered = TRUE)

##GLM stats
glm_summary_stats_andalucia <- glm_summary_df_andalucia %>%
  group_by(term, scenario) %>%
  summarize(
    mean_estimate = mean(Estimate),
    sd_estimate = sd(Estimate),
    median_estimate = median(Estimate),
    mean_p = mean(`Pr(>|z|)`, na.rm = TRUE),
    prop_significant = mean(`Pr(>|z|)` < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

##saving
write_csv(glm_summary_stats_andalucia, 
          file = "Results/glm_sensitivity_andalucia.csv")

write_xlsx(glm_summary_stats_andalucia, 
           path = "Results/glm_sensitivity_andalucia_excel.xlsx")


####plotting
png("Figures/Supplementary/Sensitivity/glm_random_estimates_andalucia.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

##plots
glm_summary_df_andalucia %>%
filter(term != "Intercept") %>%
  ggplot(aes(x = scenario, y = Estimate, fill = scenario)) +
  geom_boxplot() +
  facet_wrap(~ term, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Random assignment of DD & NE sp. [Eastern Andalusia]",
       y = "Estimate of GLM Coefficients", x = "Extinction Scenario") +
  scale_fill_brewer(palette = "Set2")+
  mynamestheme+
  theme(legend.position = "none")

dev.off()



##proportion of significant terms
png("Figures/Supplementary/Sensitivity/glm_random_significant_Andalusia.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

glm_summary_df_andalucia %>%
  filter(term != "Intercept") %>%
  mutate(significant = `Pr(>|z|)` < 0.05) %>%
  group_by(term, scenario) %>%
  summarize(prop_sig = mean(significant), .groups = "drop") %>% 
  ggplot(aes(x = scenario, y = prop_sig, fill = scenario)) +
  geom_col() +
  facet_wrap(~ term) +
  theme_minimal(base_size = 14) +
  labs(y = "Proportion of Significant Results (p < 0.05)",
       x = "Extinction scenario",
       title = "Random assignment of DD & NE sp. [Eastern Andalusia]") +
  scale_fill_brewer(palette = "Set2") +
  mynamestheme+
  theme(legend.position = "none")

dev.off()


# Probability of extinction -----------------------------------------------

### Reading mega plant phylogeny
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

############### Peninsula #################################

#reading peninsula iucn data
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

## determining proportions as probabilities for the randomization

#removing DD and NE for obtaining proportions
peninsula_iucn_ss <- peninsula_iucn %>% filter(Status != c("DD", "NE"))

prop_threatened <- mean(peninsula_iucn_ss$threat_category == "threatened")
prop_not_threatened <- 1 - prop_threatened


##randomizing the DD and NE species
n_replicates <- 100
sensitivity_list_pen_EDGE <- vector("list", n_replicates)

set.seed(42)

##Implementing EDGE across the randomizations
for (i in 1:n_replicates) {
  pen_ss <- peninsula_iucn %>%
    mutate(Status_adjusted = Status)
  
  dd_ne <- pen_ss$Status %in% c("DD", "NE")
  
  pen_ss$Status_adjusted[dd_ne] <- sample(
    c("VU", "LC"),
    size = sum(dd_ne),
    replace = TRUE,
    prob = c(prop_threatened, prop_not_threatened))
  
  pen_ss <- pen_ss %>%
    select(species, Status_adjusted) %>% 
      rename(RL.cat = Status_adjusted) %>% 
      filter(species %in% common_species_pen)
    
    ## Regionally extinct as Extinct
    pen_ss$RL.cat[pen_ss$RL.cat == "RE"] <- "EX"
    
    pen_ss_clean <- data.frame(
      species = as.character(pen_ss$species),
      RL.cat  = as.character(pen_ss$RL.cat))
    
    #### applying EDGE2
    sensitivity_list_pen_EDGE[[i]] <- calculate_EDGE2(tree = peninsula_phylo,
                                       table = pen_ss_clean,
                                       verbose = FALSE)
    
}

#saving
save(sensitivity_list_pen_EDGE, file = "Data/Processed/Sensitivity/sensitivity_list_peninsula_EDGE.RData")

#loading
load("Data/Processed/Sensitivity/sensitivity_list_peninsula_EDGE.RData")

###########phylogenetic signal ##############

# Initialize result lists
K_pen_EDGE_list <- vector("list", length = 100)
lambda_pen_EDGE_list <- vector("list", length = 100)

# Start loop for measuring phylo signal
for (i in seq_len(100)) {
  message("Processing replicate ", i)
  
  df <- sensitivity_list_pen_EDGE[[i]]
  
  ##filtering data
  df_signal <- df[df$species %in% peninsula_phylo$tip.label, ]
  
  # reorder data
  df_signal <- df_signal[match(peninsula_phylo$tip.label,
                               df_signal$species), ]
  
  # extracting the probability of extinction
  pext_vector <- df_signal$pext
  
  names(pext_vector) <- df_signal$species
  
  
  # Try-catch in case phylosig fails
  K_pen_EDGE_list[[i]] <- tryCatch(
    phylosig(peninsula_phylo, 
             pext_vector, nsim = 100, test = TRUE, method = "K"),
    error = function(e) NA
  )
  
  lambda_pen_EDGE_list[[i]] <- tryCatch(
    phylosig(peninsula_phylo, 
             pext_vector, nsim = 100, test = TRUE, method = "lambda"),
    error = function(e) NA
  )
}


##converting to df

###Blomberg
K_pen_EDGE_df <- purrr::map_dfr(K_pen_EDGE_list, function(x) {
  if (is.list(x)) {
    tibble(K = x$K, P = x$P)
  } else {
    tibble(K = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(K_pen_EDGE_df,
          file = "Data/Processed/Sensitivity/Randomizations/K_pen_EDGE_df.csv")

#reading
K_pen_EDGE_df <- read_csv("Data/Processed/Sensitivity/Randomizations/K_pen_EDGE_df.csv")


###lambda
lambda_pen_EDGE_df <- purrr::map_dfr(lambda_pen_EDGE_list, function(x) {
  if (is.list(x)) {
    tibble(lambda = x$lambda, logL = x$logL, P = x$P)
  } else {
    tibble(lambda = NA, logL = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(lambda_pen_EDGE_df,
          file = "Data/Processed/Sensitivity/Randomizations/lambda_pen_EDGE_df.csv")

#reading
lambda_pen_EDGE_df <- read_csv("Data/Processed/Sensitivity/Randomizations/lambda_pen_EDGE_df.csv")



############### Andalusia ########################

## reading Andalusia IUCN data
andalusia_iucn <- read_csv("Data/Processed/andalucia_iucn.csv")

#unifying species name
andalusia_iucn$species <- str_replace_all(andalusia_iucn$species, " ", "_")

##matching names
common_species_andalusia <- intersect(plant_phylo$tip.label, andalusia_iucn$species)

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

## determining proportions as probabilities for the randomization

#removing DD and NE for obtaining proportions
andalusia_iucn_ss <- andalusia_iucn %>% filter(Status != c("DD", "NE"))

prop_threatened <- mean(andalusia_iucn_ss$threat_category == "threatened")
prop_not_threatened <- 1 - prop_threatened


##randomizing the DD and NE species
n_replicates <- 100
sensitivity_list_andalusia_EDGE <- vector("list", n_replicates)

set.seed(42)

##Implementing EDGE across the randomizations
for (i in 1:n_replicates) {
  andalusia_ss <- andalusia_iucn %>%
    mutate(Status_adjusted = Status)
  
  dd_ne <- andalusia_ss$Status %in% c("DD", "NE")
  
  andalusia_ss$Status_adjusted[dd_ne] <- sample(
    c("VU", "LC"),
    size = sum(dd_ne),
    replace = TRUE,
    prob = c(prop_threatened, prop_not_threatened))
  
  andalusia_ss <- andalusia_ss %>%
    select(species, Status_adjusted) %>% 
    rename(RL.cat = Status_adjusted) %>% 
    filter(species %in% common_species_andalusia)
  
  ## Regionally extinct as Extinct
  andalusia_ss$RL.cat[andalusia_ss$RL.cat == "RE"] <- "EX"
  
  andalusia_ss_clean <- data.frame(
    species = as.character(andalusia_ss$species),
    RL.cat  = as.character(andalusia_ss$RL.cat))
  
  ##there are some duplicated species
  andalusia_ss_clean <- andalusia_ss_clean[
    !duplicated(andalusia_ss_clean$species), ]
  
  #### applying EDGE2
  sensitivity_list_andalusia_EDGE[[i]] <- calculate_EDGE2(tree = andalusia_phylo,
                                                    table = andalusia_ss_clean,
                                                    verbose = FALSE)
  
}

#saving
save(sensitivity_list_andalusia_EDGE, file = "Data/Processed/Sensitivity/sensitivity_list_andalusia_EDGE.RData")

#loading
load("Data/Processed/Sensitivity/sensitivity_list_andalusia_EDGE.RData")

###########phylogenetic signal ##############

# Initialize result lists
K_andalusia_EDGE_list <- vector("list", length = 100)
lambda_andalusia_EDGE_list <- vector("list", length = 100)

# Start loop for measuring phylo signal
for (i in seq_len(100)) {
  message("Processing replicate ", i)
  
  df <- sensitivity_list_andalusia_EDGE[[i]]
  
  ##filtering data
  df_signal <- df[df$species %in% andalusia_phylo$tip.label, ]
  
  # reorder data
  df_signal <- df_signal[match(andalusia_phylo$tip.label,
                               df_signal$species), ]
  
  # extracting the probability of extinction
  pext_vector <- df_signal$pext
  
  names(pext_vector) <- df_signal$species
  
  
  # Try-catch in case phylosig fails
  K_andalusia_EDGE_list[[i]] <- tryCatch(
    phylosig(andalusia_phylo, 
             pext_vector, nsim = 100, test = TRUE, method = "K"),
    error = function(e) NA
  )
  
  lambda_andalusia_EDGE_list[[i]] <- tryCatch(
    phylosig(andalusia_phylo, 
             pext_vector, nsim = 100, test = TRUE, method = "lambda"),
    error = function(e) NA
  )
}


##converting to df

###Blomberg
K_andalusia_EDGE_df <- purrr::map_dfr(K_andalusia_EDGE_list, function(x) {
  if (is.list(x)) {
    tibble(K = x$K, P = x$P)
  } else {
    tibble(K = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(K_andalusia_EDGE_df,
          file = "Data/Processed/Sensitivity/Randomizations/K_andalusia_EDGE_df.csv")

#reading
K_andalusia_EDGE_df <- read_csv("Data/Processed/Sensitivity/Randomizations/K_andalusia_EDGE_df.csv")


###lambda
lambda_andalusia_EDGE_df <- purrr::map_dfr(lambda_andalusia_EDGE_list, function(x) {
  if (is.list(x)) {
    tibble(lambda = x$lambda, logL = x$logL, P = x$P)
  } else {
    tibble(lambda = NA, logL = NA, P = NA)
  }
}, .id = "replicate")

#saving
write_csv(lambda_andalusia_EDGE_df,
          file = "Data/Processed/Sensitivity/Randomizations/lambda_andalusia_EDGE_df.csv")

#reading
lambda_andalusia_EDGE_df <- read_csv("Data/Processed/Sensitivity/Randomizations/lambda_andalusia_EDGE_df.csv")

