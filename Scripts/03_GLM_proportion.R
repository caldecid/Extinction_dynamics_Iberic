# Threatened species vs predictors (GLM) -------------------------------------------------------------

##Libraries
library(tidyverse)
library(ggcorrplot)
library(DescTools)



###############################Peninsula#######################################

##calling IUCN data (not selecting lambdas calculated by the anterior method)
peninsula_iucn <- read_csv("Data/Processed/peninsula_ages_total.csv") %>% 
                                                  select(-c(lambda_a, lambda_mean))

##loading ClaDS 
load("Data/Processed/claDS_results_peninsula_phy.tre.RData") ##it is loaded as CladsOutput

##ClaDS output
pen_rates <- data.frame(Species = CladsOutput[["tree"]][["tip.label"]],
                        rates = CladsOutput[["lambdatip_map"]]) 
pen_genus_rates <- pen_rates %>%
  separate(Species, into = c("Genus", "epithet"), sep = "_") %>%
  group_by(Genus) %>%
  summarise(
    mean_lambda = mean(rates, na.rm = TRUE),
    sd_lambda = sd(rates, na.rm = TRUE),
    n_species = n()
  )

pen_merged <- peninsula_iucn %>% left_join(pen_genus_rates, by = "Genus")

##excluding corrected age and diversification rate
peninsula_traits_non <- pen_merged %>% filter(ext_fraction == "low_ex") %>% 
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
png("Figures/Supplementary/Figure_cor_pen_non.png", 
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
peninsula_traits_low <- pen_merged %>% filter(ext_fraction == "low_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))


# Selecting the predictor columns
predictor_data <-peninsula_traits_low[, c("richness", "ev.dist", "mean_lambda",
                                          "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_pen_low.png", 
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
peninsula_traits_int <- pen_merged %>% filter(ext_fraction == "int_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <-peninsula_traits_int[, c("richness", "ev.dist", "mean_lambda",
                                          "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_pen_int.png", 
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
peninsula_traits_high <- pen_merged %>% filter(ext_fraction == "high_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <-peninsula_traits_high[, c("richness", "ev.dist", "mean_lambda",
                                          "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_pen_high.png", 
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

##calling data
andalucia_iucn <- read_csv("Data/Processed/andalucia_ages_total.csv")

##loading ClaDS 
load("Data/Processed/claDS_results_andalucia_phy.tre.RData") ##it is loaded as CladsOutput

##ClaDS output
and_rates <- data.frame(Species = CladsOutput[["tree"]][["tip.label"]],
                        rates = CladsOutput[["lambdatip_map"]]) 
and_genus_rates <- and_rates %>%
  separate(Species, into = c("Genus", "epithet"), sep = "_") %>%
  group_by(Genus) %>%
  summarise(
    mean_lambda = mean(rates, na.rm = TRUE),
    sd_lambda = sd(rates, na.rm = TRUE),
    n_species = n()
  )

and_merged <- andalucia_iucn %>% left_join(and_genus_rates, by = "Genus") %>% 
            select(-c(lambda_a, lambda_mean))


##only the predictors (not corrected age)
andalucia_traits_non <- and_merged %>% filter(ext_fraction == "low_ex") %>% 
  select(c("age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_non[, c("richness", "ev.dist", "mean_lambda",
                                           "age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_andalucia_non.png", 
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
andalucia_traits_low <- and_merged %>% filter(ext_fraction == "low_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_low[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap

png("Figures/Supplementary/Figure_cor_andalucia_low.png", 
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
andalucia_traits_int <- and_merged %>% filter(ext_fraction == "int_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))

# Selecting the predictor columns
predictor_data <- andalucia_traits_int[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_andalucia_int.png", 
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
andalucia_traits_high <- and_merged %>% filter(ext_fraction == "high_ex") %>% 
  select(c("mean_age", "richness", "ev.dist",
           "mean_lambda", "proportion_threatened"))
# Selecting the predictor columns
predictor_data <- andalucia_traits_high[, c("richness", "ev.dist", "mean_lambda",
                                           "mean_age")]


# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")
print(cor_matrix)

# Visualize the correlation matrix with a heatmap
png("Figures/Supplementary/Figure_cor_andalucia_high.png", 
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


save(model_pen_non, model_pen_low, model_pen_int, model_pen_high,
     model_anda_non, model_anda_low, model_anda_int, model_anda_high,
     peninsula_traits_non, peninsula_traits_low, peninsula_traits_int,
     peninsula_traits_high,
     andalucia_traits_non, andalucia_traits_low, andalucia_traits_int,
     andalucia_traits_high, 
     file = "Data/Processed/data1_4.RData")

