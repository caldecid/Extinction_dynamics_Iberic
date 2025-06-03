
# Sensitivity analyses (DD and NE as threatened) --------------------------


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



# IUCN data  ----------------------------


## Iberic Peninsula data
path <- "Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls"
excel_sheets(path)

peninsula <- read_excel(path,
                        sheet = "Lista_Peninsula_Final.txt" )


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

peninsula_ss$threat_category <- factor(peninsula_ss$threat_category)

peninsula_ss$Status <- factor(peninsula_ss$Status)

### Genus arrangements
peninsula_ss <- peninsula_ss %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")


### Save
write_csv(peninsula_ss, file = "Data/Processed/Sensitivity/peninsula_iucn_ss.csv")


## Andalucia data

andalucia <- read_excel(path,
                        sheet = "FLORANDOR")


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

##saving
write_csv(andalucia_ss, file = "Data/Processed/Sensitivity/andalucia_iucn_ss.csv")



# Phylogenetic manipulation -----------------------------------------------


## PENINSULA

### String from the peninsula species name
peninsula_sp_names <- str_replace(peninsula_ss$species, " ", "_")

### Eliminate subspecies
peninsula_without_sub <- unique(str_replace(peninsula_sp_names, " .*", ""))

### Reading mega plant phylogeny#############
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

### Eliminating species not present in phylogeny
n_pen <- peninsula_without_sub[peninsula_without_sub %in% plant_phylo$tip.label]

### Pruning tree preserving peninsula sp
peninsula_phy <- keep.tip(plant_phylo, n_pen)


### Determining the proportion of threatened and non-threatened species
peninsula_genus_ss <- peninsula_ss %>% group_by(Genus) %>%
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
peninsula_ages_total_ss <- left_join(peninsula_ages_rates, peninsula_genus_ss,
                                  by = "Genus") %>%
  left_join(ev.distinc.peninsula,
            by = "Genus") 


##saving
write_csv(peninsula_ages_total_ss, 
          file = "Data/Processed/Sensitivity/peninsula_ages_total_ss.csv")


##calculating Bloomberg's K and Pagel's lambda
##filtering data
data <- peninsula_genus_ss[peninsula_genus_ss$Genus %in% genus_phylo$tip.label, ]

# reorder data
data <- data[match(genus_phylo$tip.label, data$Genus), ]

# extracting the proportion of threatened species
threatened_vector <- data$proportion_threatened

names(threatened_vector) <- data$Genus

# Calculate Blomberg's K or Pagel's lambda
K_result_ss <- phylosig(genus_phylo, threatened_vector,
                     nsim = 100, test = TRUE, method = "K")

lambda_result_ss <- phylosig(genus_phylo,
                          threatened_vector, nsim = 100,
                          test = TRUE, method = "lambda")


#########################ANDALUCIA##############################################

##string from the andalucia species name
andalucia_sp_names <- str_replace(andalucia_ss$species, " ", "_")

## eliminate subspecies
andalucia_without_sub <- unique(str_replace(andalucia_sp_names, " .*", ""))

##reading mega plant phylogeny#############
plant_phylo <- read.tree("Data/Raw/PhytoPhylo.tre")

##eliminating species not present in phylogeny
n_and <- andalucia_without_sub[andalucia_without_sub %in% plant_phylo$tip.label]

##pruning tree preserving andalucia sp
andalucia_phy <- keep.tip(plant_phylo, n_and)

##saving
save(andalucia_phy,
     file = "Data/Processed/Sensitivity/andalucia_phy.RData")


##determining the proportion of threatened and non-threatened species
andalucia_genus_ss <- andalucia_ss %>% group_by(Genus) %>%
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
andalucia_ages_total_ss <- left_join(andalucia_ages_rates, andalucia_genus_ss,
                                  by = "Genus") %>%
  left_join(ev.distinc.andalucia,
            by = "Genus")


##saving
write_csv(andalucia_ages_total_ss, 
          file = "Data/Processed/Sensitivity/andalucia_ages_total_ss.csv")


##calculating Bloomberg's K and Pagel's lambda
##filtering data
data <- andalucia_genus_ss[andalucia_genus_ss$Genus %in% genus_phylo$tip.label, ]

# reorder data
data <- data[match(genus_phylo$tip.label, data$Genus), ]

# extracting the proportion of threatened species
threatened_vector <- data$proportion_threatened

names(threatened_vector) <- data$Genus

# Calculate Blomberg's K and Pagel's lambda
K_result_and_ss <- phylosig(genus_phylo, threatened_vector,
                         nsim = 100, test = TRUE, method="K")

lambda_result_and_ss <- phylosig(genus_phylo,
                              threatened_vector, nsim = 100,
                              test = TRUE, method="lambda")


# GLM analyses ------------------------------------------------------------


###############################Peninsula#######################################

###reading ClaDs results
pen_rates <- read_csv(file = "Data/Processed/peninsula_ClaDS_rates.csv")

##averaging to genus rates
pen_genus_rates <- pen_rates %>%
  separate(Species, into = c("Genus", "epithet"), sep = "_") %>%
  group_by(Genus) %>%
  summarise(
    mean_lambda = mean(rates, na.rm = TRUE),
    sd_lambda = sd(rates, na.rm = TRUE),
    n_species = n()
  )


pen_merged_ss <- peninsula_ages_total_ss %>% left_join(pen_genus_rates,
                                                       by = "Genus") %>% 
  select(-c( n_species))


##saving 
write_csv(pen_merged_ss, 
          file = "Data/Processed/Sensitivity/peninsula_merged_iucn_clads_ss.csv")

##excluding corrected age and diversification rate
peninsula_traits_non <- pen_merged_ss %>% filter(ext_fraction == "low_ex") %>% 
  select(c("age", "richness", "ev.dist",
           "mean_lambda",
           "proportion_threatened"))

##evaluating correlation in the predictors

# Selecting the predictor columns
predictor_data <- peninsula_traits_non[, c("richness", "ev.dist", "mean_lambda",
                                           "age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_pen_non.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ##excluding evo distinct due to correlation with age

dev.off()
##fitting GLM
model_pen_non <- glm(proportion_threatened ~ richness + age + mean_lambda,
                     data = peninsula_traits_non,
                     family = quasibinomial())

summary(model_pen_non)

##pseudo R2
pseudo_R2_pen_non <- PseudoR2(model_pen_non, which = "all")


# low extinction ----------------------------------------------------------

##only the predictors (corrected age)
peninsula_traits_low <- pen_merged_ss %>% filter(ext_fraction == "low_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))


# Selecting the predictor columns
predictor_data <-peninsula_traits_low[, c("richness", "ev.dist", "mean_lambda",
                                          "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_pen_low.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ##excluding evo distinct due to correlation with age


dev.off()

##fitting GLM
model_pen_low <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                     data = peninsula_traits_low,
                     family = quasibinomial())

summary(model_pen_low)

##pseudo R2
pseudo_R2_pen_low <- PseudoR2(model_pen_low, which = "all")



# intermediate extinction -------------------------------------------------

##only the predictors (corrected age)
peninsula_traits_int <- pen_merged_ss %>% filter(ext_fraction == "int_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <-peninsula_traits_int[, c("richness", "ev.dist", "mean_lambda",
                                          "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_pen_int.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ##excluding evo distinct due to correlation with age

dev.off()

##fitting GLM
model_pen_int <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                     data = peninsula_traits_int,
                     family = quasibinomial())

summary(model_pen_int)

##pseudo R2
pseudo_R2_pen_int <- PseudoR2(model_pen_int, which = "all")



# high extinction ---------------------------------------------------------


##only the predictors (corrected age)
peninsula_traits_high <- pen_merged_ss %>% filter(ext_fraction == "high_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <-peninsula_traits_high[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_pen_high.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ##not excluding evo.dist due to the rule of thumb
#######the question here is whether I include the evo.dist or not. 

dev.off()

##fitting GLM
model_pen_high <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                      data = peninsula_traits_high,
                      family = quasibinomial())

summary(model_pen_high)

##pseudo R2
pseudo_R2_pen_high <- PseudoR2(model_pen_high, which = "all")



######################ANDALUCIA##################################################

#reading clads output for Andalucia
and_rates <- read_csv("Data/Processed/andalucia_ClaDS_rates.csv")

#averaging rates per genus
and_genus_rates <- and_rates %>%
  separate(Species, into = c("Genus", "epithet"), sep = "_") %>%
  group_by(Genus) %>%
  summarise(
    mean_lambda = mean(rates, na.rm = TRUE),
    sd_lambda = sd(rates, na.rm = TRUE),
    n_species = n()
  )

# merging with andalucia iucn
and_merged_ss <- andalucia_ages_total_ss%>% left_join(and_genus_rates, by = "Genus") 

#saving
write_csv(and_merged_ss,
          file = "Data/Processed/Sensitivity/andalucia_merged_iucn_clads_ss.csv")


##only the predictors (not corrected age)
andalucia_traits_non <- and_merged_ss %>% filter(ext_fraction == "low_ex") %>% 
  select(c("age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_non[, c("richness", "ev.dist", "mean_lambda",
                                           "age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_andalucia_non.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)


ggcorrplot(cor_matrix, lab = TRUE) ## excluding ev.dist due to correlation with age




dev.off()

##fitting GLM
model_anda_non <- glm(proportion_threatened ~ richness + age + mean_lambda,
                      data = andalucia_traits_non,
                      family = quasibinomial())

summary(model_anda_non) ##no significant terms

##pseudo R2
pseudo_R2_anda_non <- PseudoR2(model_anda_non, which = "all")



# low extinction ----------------------------------------------------------

##only the predictors (corrected age)
andalucia_traits_low <- and_merged_ss %>% filter(ext_fraction == "low_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_low[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap

png("Figures/Supplementary/Sensitivity/Figure_cor_andalucia_low.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ## excluding ev.dist due to correlation with age

dev.off()
##fitting GLM
model_anda_low <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                      data = andalucia_traits_low,
                      family = quasibinomial())

summary(model_anda_low) ##no significant terms

##pseudo R2
pseudo_R2_anda_low <- PseudoR2(model_anda_low, which = "all")



# intermediate extinction -------------------------------------------------

#filtering and selecting
andalucia_traits_int <- and_merged_ss %>% filter(ext_fraction == "int_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_int[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_andalucia_int.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ## excluding ev.dist due to correlation with age

dev.off()
##fitting GLM
model_anda_int <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                      data = andalucia_traits_int,
                      family = quasibinomial())

summary(model_anda_int) ##no significant terms

##pseudo R2
pseudo_R2_anda_int <- PseudoR2(model_anda_int, which = "all")


# high extinction ---------------------------------------------------------

##filtering and selecting
andalucia_traits_high <- and_merged_ss %>% filter(ext_fraction == "high_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))
# Selecting the predictor columns
predictor_data <- andalucia_traits_high[, c("richness", "ev.dist", "mean_lambda",
                                            "mean_age")]


# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Sensitivity/Figure_cor_andalucia_high.png", 
    width = 15, height = 15,
    units = "cm", pointsize = 8, res = 300)

ggcorrplot(cor_matrix, lab = TRUE) ## excluding ev.dist due to correlation with age

dev.off()
##fitting GLM
model_anda_high <- glm(proportion_threatened ~ richness + mean_age + mean_lambda,
                       data = andalucia_traits_high,
                       family = quasibinomial())

summary(model_anda_high) ##no significant terms

##pseudo R2
pseudo_R2_anda_high <- PseudoR2(model_anda_high, which = "all")


# Figures -----------------------------------------------------------------


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


# Peninsula ---------------------------------------------------------------

#############################GLM coeficients############################

##not corrected

summary_model_pen_non <- summary(model_pen_non)

coefs_pen_non <- data.frame(
  Term = rownames(summary_model_pen_non$coefficients),
  Estimate = summary_model_pen_non$coefficients[, "Estimate"],
  StdError = summary_model_pen_non$coefficients[, "Std. Error"],
  LowerCI = summary_model_pen_non$coefficients[,
                                               "Estimate"] - 1.96 * summary_model_pen_non$coefficients[, "Std. Error"],
  UpperCI = summary_model_pen_non$coefficients[,
                                               "Estimate"] + 1.96 * summary_model_pen_non$coefficients[, "Std. Error"],
  pValue = summary_model_pen_non$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_pen_non <- coefs_pen_non[-1,]

coefs_pen_non$Term <- factor(coefs_pen_non$Term, levels = c("richness",
                                                            "age", "mean_lambda"),
                             labels = c("richness", "mean_age",
                                        "mean_lambda"),
                             ordered = TRUE)

##corrected
coefs_pen_non$corrected <- "not_corrected"


##########low extinction

summary_model_pen_low <- summary(model_pen_low)

coefs_pen_low <- data.frame(
  Term = rownames(summary_model_pen_low$coefficients),
  Estimate = summary_model_pen_low$coefficients[, "Estimate"],
  StdError = summary_model_pen_low$coefficients[, "Std. Error"],
  LowerCI = summary_model_pen_low$coefficients[,
                                               "Estimate"] - 1.96 * summary_model_pen_low$coefficients[, "Std. Error"],
  UpperCI = summary_model_pen_low$coefficients[,
                                               "Estimate"] + 1.96 * summary_model_pen_low$coefficients[, "Std. Error"],
  pValue = summary_model_pen_low$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_pen_low <- coefs_pen_low[-1,]

##ordering
coefs_pen_low$Term <- factor(coefs_pen_low$Term, levels = c("richness",
                                                            "mean_age",
                                                            "mean_lambda"),
                             ordered = TRUE)

##corrected
coefs_pen_low$corrected <- "low_ext"


########intermediate extinction

summary_model_pen_int <- summary(model_pen_int)

coefs_pen_int <- data.frame(
  Term = rownames(summary_model_pen_int$coefficients),
  Estimate = summary_model_pen_int$coefficients[, "Estimate"],
  StdError = summary_model_pen_int$coefficients[, "Std. Error"],
  LowerCI = summary_model_pen_int$coefficients[,
                                               "Estimate"] - 1.96 * summary_model_pen_int$coefficients[, "Std. Error"],
  UpperCI = summary_model_pen_int$coefficients[,
                                               "Estimate"] + 1.96 * summary_model_pen_int$coefficients[, "Std. Error"],
  pValue = summary_model_pen_int$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_pen_int <- coefs_pen_int[-1,]

##ordering
coefs_pen_int$Term <- factor(coefs_pen_int$Term, levels = c("richness",
                                                            "mean_age",
                                                            "mean_lambda"),
                             ordered = TRUE)

##corrected
coefs_pen_int$corrected <- "int_ext"


########high extinction

summary_model_pen_high <- summary(model_pen_high)

coefs_pen_high <- data.frame(
  Term = rownames(summary_model_pen_high$coefficients),
  Estimate = summary_model_pen_high$coefficients[, "Estimate"],
  StdError = summary_model_pen_high$coefficients[, "Std. Error"],
  LowerCI = summary_model_pen_high$coefficients[,
                                                "Estimate"] - 1.96 * summary_model_pen_high$coefficients[, "Std. Error"],
  UpperCI = summary_model_pen_high$coefficients[,
                                                "Estimate"] + 1.96 * summary_model_pen_high$coefficients[, "Std. Error"],
  pValue = summary_model_pen_high$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_pen_high <- coefs_pen_high[-1,]

##ordering
coefs_pen_high$Term <- factor(coefs_pen_high$Term, levels = c("richness",
                                                              "mean_age",
                                                              "mean_lambda"),
                              ordered = TRUE)
##corrected
coefs_pen_high$corrected <- "high_ext"


###binding dataframes
coefs_pen_total <- rbind(coefs_pen_non, coefs_pen_low, 
                         coefs_pen_int, coefs_pen_high)

##organizing
coefs_pen_total$Term <- factor(coefs_pen_total$Term, levels = c("richness",
                                                                "mean_age",
                                                                "mean_lambda"),
                               labels = c("Richness",
                                          "Genus age",
                                          "Speciation rates"),
                               ordered = TRUE)

coefs_pen_total$corrected <- factor(coefs_pen_total$corrected,
                                    levels = c("not_corrected",
                                               "low_ext",
                                               "int_ext",
                                               "high_ext"),
                                    labels = c("Not corrected",
                                               "Low extinction",
                                               "Intermediate extinction",
                                               "High extinction"),
                                    ordered = TRUE)



########plot coeficients###############

png("Figures/Supplementary/Sensitivity/Figure_pen_coef.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)


pen_coefs <- ggplot(coefs_pen_total, aes(x = Term, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Phylogenetic traits", y = "Coefficient Estimate")+
  theme_bw() +
  #ylim(-0.2, 0.2)+
  scale_x_discrete(labels = c("Richness","Genus age", expression(~lambda)))+
  facet_wrap(~ corrected, nrow = 1) +
  ggtitle("Iberian Peninsula [GLM predictors] Sensitivity")+
  mynamestheme

pen_coefs
dev.off()

#############################GLM AGE############################

###not corrected
# Prepare peninsula_traits_non for prediction in order to plot
pred_peninsula_non_age <- data.frame(age = seq(min(peninsula_traits_non$age),
                                               max(peninsula_traits_non$age),
                                               length.out = 100),
                                     richness = mean(peninsula_traits_non$richness,
                                                     na.rm = TRUE),
                                     mean_lambda = mean(peninsula_traits_non$mean_lambda, na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_non_age <- cbind(pred_peninsula_non_age,
                                predict(model_pen_non,
                                        newdata = pred_peninsula_non_age,
                                        type = "link", se.fit = TRUE))
pred_peninsula_non_age <- transform(pred_peninsula_non_age,
                                    predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                    conf.low = plogis(fit - 1.96 * se.fit),
                                    conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_non_age$corrected <- "not_corrected"

##changing names for binding
pred_peninsula_non_age <- pred_peninsula_non_age %>% rename(mean_age = age)


############low extinction
# Prepare for prediction in order to plot
pred_peninsula_low_age <- data.frame(mean_age = seq(min(peninsula_traits_low$mean_age),
                                                    max(peninsula_traits_low$mean_age),
                                                    length.out = 100),
                                     richness = mean(peninsula_traits_low$richness,
                                                     na.rm = TRUE),
                                     mean_lambda = mean(peninsula_traits_low$mean_lambda,
                                                        na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_low_age <- cbind(pred_peninsula_low_age,
                                predict(model_pen_low,
                                        newdata = pred_peninsula_low_age,
                                        type = "link", se.fit = TRUE))
pred_peninsula_low_age <- transform(pred_peninsula_low_age,
                                    predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                    conf.low = plogis(fit - 1.96 * se.fit),
                                    conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_low_age$corrected <- "low_ext"


########intermediate extinction

pred_peninsula_int_age <- data.frame(mean_age = seq(min(peninsula_traits_int$mean_age),
                                                    max(peninsula_traits_int$mean_age),
                                                    length.out = 100),
                                     richness = mean(peninsula_traits_int$richness,
                                                     na.rm = TRUE),
                                     mean_lambda = mean(peninsula_traits_int$mean_lambda,
                                                        na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_int_age <- cbind(pred_peninsula_int_age,
                                predict(model_pen_int,
                                        newdata = pred_peninsula_int_age,
                                        type = "link", se.fit = TRUE))

pred_peninsula_int_age <- transform(pred_peninsula_int_age,
                                    predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                    conf.low = plogis(fit - 1.96 * se.fit),
                                    conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_int_age$corrected <- "int_ext"


########high extinction

pred_peninsula_high_age <- data.frame(mean_age = seq(min(peninsula_traits_high$mean_age),
                                                     max(peninsula_traits_high$mean_age),
                                                     length.out = 100),
                                      richness = mean(peninsula_traits_high$richness,
                                                      na.rm = TRUE),
                                      mean_lambda = mean(peninsula_traits_high$mean_lambda,
                                                         na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_high_age <- cbind(pred_peninsula_high_age,
                                 predict(model_pen_high,
                                         newdata = pred_peninsula_high_age,
                                         type = "link", se.fit = TRUE))
pred_peninsula_high_age <- transform(pred_peninsula_high_age,
                                     predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                     conf.low = plogis(fit - 1.96 * se.fit),
                                     conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_high_age$corrected <- "high_ext"


#######binding data
pred_pen_age <- rbind(pred_peninsula_non_age,
                      pred_peninsula_low_age,
                      pred_peninsula_int_age,
                      pred_peninsula_high_age)


##order
pred_pen_age$corrected <- factor(pred_pen_age$corrected,
                                 levels = c("not_corrected",
                                            "low_ext",
                                            "int_ext",
                                            "high_ext"),
                                 labels = c("Not corrected",
                                            "Low extinction",
                                            "Intermediate extinction",
                                            "High extinction"),
                                 ordered = TRUE)

##save
write_csv(pred_pen_age,
          file = "Data/Processed/peninsula_predict_ages.csv")



# Predicted response vs predictor plot
########plot coeficients###############

png("Figures/Supplementary/Sensitivity/Figure_pen_age.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)

pen_ages <- ggplot(pred_pen_age, aes(x = mean_age, y = predicted)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Genus age", y = "Prop. of threatened species") +
  theme_bw()+
  facet_wrap(~corrected, nrow = 1, scales = "free_x")+
  #ggtitle("Peninsula")+
  mynamestheme
pen_ages
dev.off()


###########################GLM Richness########################

###not corrected
# Prepare peninsula_traits_non for prediction in order to plot
pred_peninsula_non_richness <- data.frame(richness = seq(min(peninsula_traits_non$richness),
                                                         max(peninsula_traits_non$richness),
                                                         length.out = 100),
                                          age = mean(peninsula_traits_non$age,
                                                     na.rm = TRUE),
                                          mean_lambda = mean(peninsula_traits_non$mean_lambda, na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_non_richness <- cbind(pred_peninsula_non_richness,
                                     predict(model_pen_non,
                                             newdata = pred_peninsula_non_richness,
                                             type = "link", se.fit = TRUE))

pred_peninsula_non_richness <- transform(pred_peninsula_non_richness,
                                         predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                         conf.low = plogis(fit - 1.96 * se.fit),
                                         conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_non_richness$corrected <- "not_corrected"

##changing names for binding
pred_peninsula_non_richness <- pred_peninsula_non_richness %>% rename(mean_age = age)


############low extinction

# Prepare for prediction in order to plot
pred_peninsula_low_richness <- data.frame(richness = seq(min(peninsula_traits_low$richness),
                                                         max(peninsula_traits_low$richness),
                                                         length.out = 100),
                                          mean_age = mean(peninsula_traits_low$mean_age,
                                                          na.rm = TRUE),
                                          mean_lambda = mean(peninsula_traits_low$mean_lambda,
                                                             na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_low_richness <- cbind(pred_peninsula_low_richness,
                                     predict(model_pen_low,
                                             newdata = pred_peninsula_low_richness,
                                             type = "link", se.fit = TRUE))

pred_peninsula_low_richness <- transform(pred_peninsula_low_richness,
                                         predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                         conf.low = plogis(fit - 1.96 * se.fit),
                                         conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_low_richness$corrected <- "low_ext"


########intermediate extinction

pred_peninsula_int_richness <- data.frame(richness = seq(min(peninsula_traits_int$richness),
                                                         max(peninsula_traits_int$richness),
                                                         length.out = 100),
                                          mean_age = mean(peninsula_traits_int$mean_age,
                                                          na.rm = TRUE),
                                          mean_lambda = mean(peninsula_traits_int$mean_lambda,
                                                             na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_int_richness <- cbind(pred_peninsula_int_richness,
                                     predict(model_pen_int,
                                             newdata = pred_peninsula_int_richness,
                                             type = "link", se.fit = TRUE))

pred_peninsula_int_richness <- transform(pred_peninsula_int_richness,
                                         predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                         conf.low = plogis(fit - 1.96 * se.fit),
                                         conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_int_richness$corrected <- "int_ext"


########high extinction

pred_peninsula_high_richness <- data.frame(richness = seq(min(peninsula_traits_high$richness),
                                                          max(peninsula_traits_high$richness),
                                                          length.out = 100),
                                           mean_age = mean(peninsula_traits_high$mean_age,
                                                           na.rm = TRUE),
                                           mean_lambda = mean(peninsula_traits_high$mean_lambda,
                                                              na.rm = TRUE))

# Get predicted values with confidence intervals
pred_peninsula_high_richness <- cbind(pred_peninsula_high_richness,
                                      predict(model_pen_high,
                                              newdata = pred_peninsula_high_richness,
                                              type = "link", se.fit = TRUE))

pred_peninsula_high_richness <- transform(pred_peninsula_high_richness,
                                          predicted = plogis(fit),  # Transform from log-odds to probability for binomial model
                                          conf.low = plogis(fit - 1.96 * se.fit),
                                          conf.high = plogis(fit + 1.96 * se.fit))

##factor
pred_peninsula_high_richness$corrected <- "high_ext"


#######binding data
pred_pen_richness <- rbind(pred_peninsula_non_richness,
                           pred_peninsula_low_richness,
                           pred_peninsula_int_richness,
                           pred_peninsula_high_richness)


##order
pred_pen_richness$corrected <- factor(pred_pen_richness$corrected,
                                      levels = c("not_corrected",
                                                 "low_ext",
                                                 "int_ext",
                                                 "high_ext"),
                                      labels = c("Not corrected",
                                                 "Low extinction",
                                                 "Intermediate extinction",
                                                 "High extinction"),
                                      ordered = TRUE)

##save
write_csv(pred_pen_richness,
          file = "Data/Processed/peninsula_predict_richness.csv")




# Predicted response vs predictor plot
########plot coeficients###############

png("Figures/Supplementary/Sensitivity/Figure_pen_richness.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)

pen_richness <- ggplot(pred_pen_richness, aes(x = richness, y = predicted)) +
  geom_line(color = "deepskyblue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Richness", y = "Prop. of threatened species") +
  theme_bw()+
  facet_wrap(~corrected, nrow = 1)+
  
  mynamestheme

pen_richness
dev.off()


## merging the three figures of Iberian Peninsula

##Plotting
png("Figures/Supplementary/Sensitivity/Figure_Pen_sensitivity.png", 
    width = 30, height = 30,
    units = "cm", pointsize = 8, res = 300)

grid.arrange(pen_coefs,
             pen_richness,
             pen_ages,
             nrow = 3)

dev.off()


# Andalucia ---------------------------------------------------------------

###not corrected########

summary_model_anda_non <- summary(model_anda_non)

coefs_anda_non <- data.frame(
  Term = rownames(summary_model_anda_non$coefficients),
  Estimate = summary_model_anda_non$coefficients[, "Estimate"],
  StdError = summary_model_anda_non$coefficients[, "Std. Error"],
  LowerCI = summary_model_anda_non$coefficients[,
                                                "Estimate"] - 1.96 * summary_model_anda_non$coefficients[, "Std. Error"],
  UpperCI = summary_model_anda_non$coefficients[,
                                                "Estimate"] + 1.96 * summary_model_anda_non$coefficients[, "Std. Error"],
  pValue = summary_model_anda_non$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_anda_non <- coefs_anda_non[-1,]

##ordering
coefs_anda_non$Term <- factor(coefs_anda_non$Term, levels = c("richness",
                                                              "age",
                                                              "mean_lambda"),
                              labels = c("richness",
                                         "mean_age",
                                         "mean_lambda"),
                              ordered = TRUE)
##factor
coefs_anda_non$corrected <- "Not corrected"


###########low extinction

summary_model_anda_low <- summary(model_anda_low)

coefs_anda_low <- data.frame(
  Term = rownames(summary_model_anda_low$coefficients),
  Estimate = summary_model_anda_low$coefficients[, "Estimate"],
  StdError = summary_model_anda_low$coefficients[, "Std. Error"],
  LowerCI = summary_model_anda_low$coefficients[,
                                                "Estimate"] - 1.96 * summary_model_anda_low$coefficients[, "Std. Error"],
  UpperCI = summary_model_anda_low$coefficients[,
                                                "Estimate"] + 1.96 * summary_model_anda_low$coefficients[, "Std. Error"],
  pValue = summary_model_anda_low$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_anda_low <- coefs_anda_low[-1,]

##ordering
coefs_anda_low$Term <- factor(coefs_anda_low$Term, levels = c("richness",
                                                              "mean_age",
                                                              "mean_lambda"),
                              ordered = TRUE)


##factor
coefs_anda_low$corrected <- "Low extinction"



######intermediate extinction
summary_model_anda_int <- summary(model_anda_int)

coefs_anda_int <- data.frame(
  Term = rownames(summary_model_anda_int$coefficients),
  Estimate = summary_model_anda_int$coefficients[, "Estimate"],
  StdError = summary_model_anda_int$coefficients[, "Std. Error"],
  LowerCI = summary_model_anda_int$coefficients[,
                                                "Estimate"] - 1.96 * summary_model_anda_int$coefficients[, "Std. Error"],
  UpperCI = summary_model_anda_int$coefficients[,
                                                "Estimate"] + 1.96 * summary_model_anda_int$coefficients[, "Std. Error"],
  pValue = summary_model_anda_int$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_anda_int <- coefs_anda_int[-1,]

##ordering
coefs_anda_int$Term <- factor(coefs_anda_int$Term, levels = c("richness",
                                                              "mean_age",
                                                              "mean_lambda"),
                              ordered = TRUE)

#factor
coefs_anda_int$corrected <- "Intermediate extinction"

####high extinction

summary_model_anda_high <- summary(model_anda_high)

coefs_anda_high <- data.frame(
  Term = rownames(summary_model_anda_high$coefficients),
  Estimate = summary_model_anda_high$coefficients[, "Estimate"],
  StdError = summary_model_anda_high$coefficients[, "Std. Error"],
  LowerCI = summary_model_anda_high$coefficients[,
                                                 "Estimate"] - 1.96 * summary_model_anda_high$coefficients[, "Std. Error"],
  UpperCI = summary_model_anda_high$coefficients[,
                                                 "Estimate"] + 1.96 * summary_model_anda_high$coefficients[, "Std. Error"],
  pValue = summary_model_anda_high$coefficients[, "Pr(>|t|)"]
)

##eliminating the intercept
coefs_anda_high <- coefs_anda_high[-1,]

##ordering
coefs_anda_high$Term <- factor(coefs_anda_high$Term, levels = c("richness",
                                                                "mean_age",
                                                                "mean_lambda"),
                               ordered = TRUE)

##factor
coefs_anda_high$corrected <- "High extinction"


###binding dataframes
coefs_anda_total <- rbind(coefs_anda_non, coefs_anda_low,
                          coefs_anda_int, coefs_anda_high)


##organizing
coefs_anda_total$corrected <- factor(coefs_anda_total$corrected, 
                                     levels = c("Not corrected",
                                                "Low extinction",
                                                "Intermediate extinction",
                                                "High extinction"),
                                     ordered = TRUE)


########plot coeficients############### 
png("Figures/Supplementary/Sensitivity/Figure_and_coef.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)

ggplot(coefs_anda_total, aes(x = Term, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Phylogenetic traits", y = "Coefficient Estimate", 
       title = "Eastern Andalusia") +
  theme_bw() +
  
  #ylim(-0.3, 0.3)+
  scale_x_discrete(labels = c("Richness","Genus age", expression(~lambda)))+
  facet_wrap(~ corrected, nrow = 1)+
  mynamestheme

dev.off()




