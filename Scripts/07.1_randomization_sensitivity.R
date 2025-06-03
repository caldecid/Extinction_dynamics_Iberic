
# Randomization of DD and NE according to proportions --------------------------


## Libraries
library(tidyverse)
library(readxl)
library(phytools)
library(remotes)
library(SpeciesAge) #install_github("thauffe/SpeciesAge")
library(picante)
library(ggcorrplot)
library(DescTools)
library(gridExtra)
library(purrr)
library(broom)



# Peninsula  ----------------------------


## Iberic Peninsula data
path <- "Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls"
excel_sheets(path)

peninsula <- read_excel(path,
                        sheet = "Lista_Peninsula_Final.txt" )

##species age and speciation rates data
pen_ages_rates <- read_csv("Data/Processed/Sensitivity/peninsula_merged_iucn_clads_ss.csv")

pen_ages_rates <- pen_ages_rates %>% select(Genus, ext_fraction, mean_age,
                                            richness, mean_lambda)

### Ordering the conservation status for sensitivity analyses
not_threatened_ss <- c("LC", "NT")
threatened_ss <- c("DD", "NE", "VU", "EN", "CR", "EX", "RE", "EW")

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

#factos
peninsula_ss$threat_category <- factor(peninsula_ss$threat_category)

peninsula_ss$Status <- factor(peninsula_ss$Status)

### Genus arrangements
peninsula_ss <- peninsula_ss %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")


# Filter out NE and DD to calculate proportions
peninsula_assessed <- peninsula_ss %>%
  filter(!Status %in% c("DD", "NE"))

##calculating proportions
prop_threatened <- mean(peninsula_assessed$threat_category == "threatened")
prop_not_threatened <- 1 - prop_threatened


##randomizing the DD and NE species
n_replicates <- 100
sensitivity_list <- vector("list", n_replicates)

set.seed(42)

for (i in 1:n_replicates) {
  sensitivity_list[[i]] <- peninsula_ss %>%
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
    left_join(peninsula_ages_rates, #merging ages and speciation rates
              by = "Genus") %>% 
      drop_na(ext_fraction)
}

####### fitting the glm

# List of extinction scenarios
ext_scenarios <- c("low_ex", "int_ex", "high_ex")

# Initialize empty list to collect results
all_results <- list()

# Loop over scenarios
for (scenario in ext_scenarios) {
  
  # Loop over replicates
  for (i in seq_along(sensitivity_list)) {
    
    df <- sensitivity_list[[i]]
    
    # Filter by extinction scenario
    df_sub <- df %>%
      filter(ext_fraction == scenario) 
    
    # Fit GLM, wrapped in try() to handle errors
    model <- try(glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                     data = df_sub,
                     family = quasibinomial()),
                 silent = TRUE)
    
    # Skip model if error
    if (inherits(model, "try-error")) next
    
    # Extract summary
    smry <- summary(model)
    
    # Extract coefficient table and convert to data.frame
    coef_df <- as.data.frame(smry$coefficients)
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
glm_summary_df <- bind_rows(all_results)

##factors
glm_summary_df$scenario <- factor(glm_summary_df$scenario,
                                  levels = c("low_ex",
                                             "int_ex",
                                             "high_ex"),
                                  labels = c("Low",
                                             "Intermediate",
                                             "High"),
                                  ordered = TRUE)

glm_summary_df$term <- factor(glm_summary_df$term,
                              levels = c("(Intercept)",
                                         "richness",
                                         "mean_age",
                                         "mean_lambda"),
                              labels = c("Intercept",
                                         "Richness",
                                         "Corrected age",
                                         "lambda"),
                              ordered = TRUE)

##GLM stats
glm_summary_stats <- glm_summary_df %>%
  group_by(term, scenario) %>%
  summarize(
    mean_estimate = mean(Estimate),
    sd_estimate = sd(Estimate),
    median_estimate = median(Estimate),
    mean_p = mean(`Pr(>|t|)`, na.rm = TRUE),
    prop_significant = mean(`Pr(>|t|)` < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

##saving
write_csv(glm_summary_stats, file = "Results/glm_sensitivity_Peninsula.csv")


####plotting
png("Figures/Supplementary/Sensitivity/Random_estimates_Peninsula.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

##plotsglm_summary_df %>%
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
png("Figures/Supplementary/Sensitivity/Random_prop_Peninsula.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

glm_summary_df %>%
  filter(term != "Intercept") %>%
  mutate(significant = `Pr(>|t|)` < 0.05) %>%
  group_by(term, scenario) %>%
  summarize(prop_sig = mean(significant), .groups = "drop") %>%
  ggplot(aes(x = scenario, y = prop_sig, fill = scenario)) +
  geom_col() +
  facet_wrap(~ term) +
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
andalucia_ages_rates <- read_csv("Data/Processed/Sensitivity/andalucia_merged_iucn_clads_ss.csv")

andalucia_ages_rates <- andalucia_ages_rates %>% select(Genus, ext_fraction, mean_age,
                                             mean_lambda)

##organizing threatening status
andalucia_ss <- andalucia %>%
  mutate(
    Status = factor(
      Status,
      labels = c(not_threatened_ss, threatened_ss[-c(2,7,8)]), #removing absent categories
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
sensitivity_list <- vector("list", n_replicates)

set.seed(42)

for (i in 1:n_replicates) {
  sensitivity_list[[i]] <- andalucia_ss %>%
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

####### fitting the glm

# List of extinction scenarios
ext_scenarios <- c("low_ex", "int_ex", "high_ex")

# Initialize empty list to collect results
all_results <- list()

# Loop over scenarios
for (scenario in ext_scenarios) {
  
  # Loop over replicates
  for (i in seq_along(sensitivity_list)) {
    
    df <- sensitivity_list[[i]]
    
    # Filter by extinction scenario
    df_sub <- df %>%
      filter(ext_fraction == scenario) 
    
    # Fit GLM, wrapped in try() to handle errors
    model <- try(glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                     data = df_sub,
                     family = quasibinomial()),
                 silent = TRUE)
    
    # Skip model if error
    if (inherits(model, "try-error")) next
    
    # Extract summary
    smry <- summary(model)
    
    # Extract coefficient table and convert to data.frame
    coef_df <- as.data.frame(smry$coefficients)
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
glm_summary_df_andalucia <- bind_rows(all_results)

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
                                         "richness",
                                         "mean_age",
                                         "mean_lambda"),
                              labels = c("Intercept",
                                         "Richness",
                                         "Corrected age",
                                         "lambda"),
                              ordered = TRUE)

##GLM stats
glm_summary_stats_andalucia <- glm_summary_df_andalucia %>%
  group_by(term, scenario) %>%
  summarize(
    mean_estimate = mean(Estimate),
    sd_estimate = sd(Estimate),
    median_estimate = median(Estimate),
    mean_p = mean(`Pr(>|t|)`, na.rm = TRUE),
    prop_significant = mean(`Pr(>|t|)` < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

##saving
write_csv(glm_summary_stats_andalucia, 
          file = "Results/glm_sensitivity_andalucia.csv")


####plotting
png("Figures/Supplementary/Sensitivity/Random_estimates_andalucia.png", 
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
png("Figures/Supplementary/Sensitivity/Random_prop_Andalusia.png", 
    width = 20, height = 15,
    units = "cm", pointsize = 8, res = 300)

glm_summary_df_andalucia %>%
  filter(term != "Intercept") %>%
  mutate(significant = `Pr(>|t|)` < 0.05) %>%
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
