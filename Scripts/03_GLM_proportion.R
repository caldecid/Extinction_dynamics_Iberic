# Threatened species vs predictors (GLM) -------------------------------------------------------------

##Libraries
library(tidyverse)
library(DescTools)
library(glmmTMB) #for using betabinomial model for overdispersion
library(DHARMa)

##loading ClaDS from my desktop (due to memory it was not uploaded to github)
load("Data/Processed/claDS_plant_phylo.RData") 

##ClaDS output
plant_genus_rates <- data.frame(Species = CladsOutput[["tree"]][["tip.label"]],
                        rates = CladsOutput[["lambdatip_map"]]) %>% 
                      rename(Genus = Species)

#saving
write_csv(plant_genus_rates, file = "Data/Processed/plant_genus_rates.csv")

###############################Peninsula#######################################

##calling peninsula data
peninsula_ages_proportion <- read_csv("Data/Processed/peninsula_ages_proportion.csv")

##merging ages and diversification rates
peninsula_total <- peninsula_ages_proportion %>%
                          left_join(plant_genus_rates, by = "Genus") 

##saving 
write_csv(peninsula_total, 
                       file = "Data/Processed/peninsula_merged_iucn_clads.csv")

#reading
peninsula_total <- read_csv("Data/Processed/peninsula_merged_iucn_clads.csv")



# Low extinction ----------------------------------------------------------

peninsula_low <- peninsula_total %>% filter(ext_fraction == "low_ex")

##null model 
model_pen_low_null <- glmmTMB(cbind(threatened_species,
                                    richness - threatened_species) ~ 1,
                              family = betabinomial(),
                              data = peninsula_low)

summary(model_pen_low_null)


#### complete model
model_pen_low_1 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age + rates,
                           family = betabinomial(),
                           data = peninsula_low)
summary(model_pen_low_1)

#r squared
# McFadden R2
R2_pen_low_1 <- 1 - (logLik(model_pen_low_1) / logLik(model_pen_low_null))
R2_pen_low_1 ## 0.14% no variation explained

res_pen_low_1 <- simulateResiduals(model_pen_low_1)
plot(res_pen_low_1) #only high dispersion of the residuals detected
testDispersion(res_pen_low_1) ## high dispersion of the residuals
testZeroInflation(res_pen_low_1) ## no problems detected with the zero


### model with rates only 
model_pen_low_2 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ rates,
                           family = betabinomial(),
                           data = peninsula_low)
summary(model_pen_low_2)

#r squared
R2_pen_low_2 <- 1 - (logLik(model_pen_low_2) / logLik(model_pen_low_null))
R2_pen_low_2 

##Diagnostics
res_pen_low_2 <- simulateResiduals(model_pen_low_2)
plot(res_pen_low_2)#only high dispersion of the residuals detected
testDispersion(res_pen_low_2) ## high dispersion of the residuals
testZeroInflation(res_pen_low_2) # no problems detected with the zero

### model with age only 
model_pen_low_3 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age,
                           family = betabinomial(),
                           data = peninsula_low)
summary(model_pen_low_3)

#r squared
R2_pen_low_3 <- 1 - (logLik(model_pen_low_3) / logLik(model_pen_low_null))
R2_pen_low_3 

## Diagnostics
res_pen_low_3 <- simulateResiduals(model_pen_low_3)
plot(res_pen_low_3)#only high dispersion of the residuals detected
testDispersion(res_pen_low_3) ## high dispersion of the residuals
testZeroInflation(res_pen_low_3)

## alternative model using richess as predictor (artifact)
model_pen_low_alt_1 <- glmmTMB(cbind(threatened_species,
              richness - threatened_species) ~ richness + mean_age + rates,
        family = betabinomial(),
        data = peninsula_low)

summary(model_pen_low_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response

#r squared
R2_pen_low_alt <- 1 - (logLik(model_pen_low_alt_1) / logLik(model_pen_low_null))
R2_pen_low_alt 


### alternative model using richness as offset to remove its effect
model_pen_low_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = peninsula_low
)

summary(model_pen_low_alt_2)


# intermediate extinction ----------------------------------------------------

peninsula_int <- peninsula_total %>% filter(ext_fraction == "int_ex")

##null model 
model_pen_int_null <- glmmTMB(cbind(threatened_species,
                                    richness - threatened_species) ~ 1,
                              family = betabinomial(),
                              data = peninsula_int)
summary(model_pen_int_null) 


#### complete model
model_pen_int_1 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age + rates,
                           family = betabinomial(),
                           data = peninsula_int)
summary(model_pen_int_1)

#r squared
R2_pen_int_1 <- 1 - (logLik(model_pen_int_1) / logLik(model_pen_int_null))
R2_pen_int_1 

## Diagnostics
res_pen_int_1 <- simulateResiduals(model_pen_int_1)
plot(res_pen_int_1)#only high dispersion of the residuals detected
testDispersion(res_pen_int_1) ## high dispersion of the residuals
testZeroInflation(res_pen_int_1)


### model with rates only 
model_pen_int_2 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ rates,
                           family = betabinomial(),
                           data = peninsula_int)
summary(model_pen_int_2)

#r squared
R2_pen_int_2 <- 1 - (logLik(model_pen_int_2) / logLik(model_pen_int_null))
R2_pen_int_2 

## Diagnostics
res_pen_int_2 <- simulateResiduals(model_pen_int_2)
plot(res_pen_int_2) #only high dispersion of the residuals detected
testDispersion(res_pen_int_2) ## high dispersion of the residuals
testZeroInflation(res_pen_int_2)

### model with age only 
model_pen_int_3 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age,
                           family = betabinomial(),
                           data = peninsula_int)
summary(model_pen_int_3)

#r squared
R2_pen_int_3 <- 1 - (logLik(model_pen_int_3) / logLik(model_pen_int_null))
R2_pen_int_3 


## Diagnostics
res_pen_int_3 <- simulateResiduals(model_pen_int_3)
plot(res_pen_int_3) #only high dispersion of the residuals detected
testDispersion(res_pen_int_3) ## high dispersion of the residuals
testZeroInflation(res_pen_int_3)


## alternative model using richess as predictor (artifact)
model_pen_int_alt_1 <- glmmTMB(cbind(threatened_species,
                       richness - threatened_species) ~ richness + mean_age + rates,
                               family = betabinomial(),
                               data = peninsula_int)

summary(model_pen_int_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response
#r squared
R2_pen_int_alt <- 1 - (logLik(model_pen_int_alt_1) / logLik(model_pen_int_null))
R2_pen_int_alt 

### alternative model using richness as offset to remove its effect
model_pen_int_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = peninsula_int
)

summary(model_pen_int_alt_2)


# high extinction ---------------------------------------------------------

peninsula_high <- peninsula_total %>% filter(ext_fraction == "high_ex")

##null model 
model_pen_high_null <- glmmTMB(cbind(threatened_species,
                                     richness - threatened_species) ~ 1,
                               family = betabinomial(),
                               data = peninsula_high)
summary(model_pen_high_null) 


#### complete model
model_pen_high_1 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age + rates,
                           family = betabinomial(),
                           data = peninsula_high)
summary(model_pen_high_1)

#r squared
R2_pen_high_1 <- 1 - (logLik(model_pen_high_1) / logLik(model_pen_high_null))
R2_pen_high_1

## Diagnostics
res_pen_high_1 <- simulateResiduals(model_pen_high_1)
plot(res_pen_high_1)#only high dispersion of the residuals detected
testDispersion(res_pen_hig_1) ## high dispersion of the residuals
testZeroInflation(res_pen_high_1) ## zero no problem

### model with rates only 
model_pen_high_2 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ rates,
                           family = betabinomial(),
                           data = peninsula_high)
summary(model_pen_high_2)

R2_pen_high_2 <- 1 - (logLik(model_pen_high_2) / logLik(model_pen_high_null))
R2_pen_high_2

## Diagnostics
res_pen_high_2 <- simulateResiduals(model_pen_high_2)
plot(res_pen_high_2)#only high dispersion of the residuals detected
testDispersion(res_pen_high_2) ## high dispersion of the residuals
testZeroInflation(res_pen_high_2) ## zero no problem

### model with age only 
model_pen_high_3 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age,
                           family = betabinomial(),
                           data = peninsula_high)
summary(model_pen_high_3) 

## Diagnostics
res_pen_high_3 <- simulateResiduals(model_pen_high_3)
plot(res_pen_high_3)#only high dispersion of the residuals detected
testDispersion(res_pen_high_3) ## high dispersion of the residuals
testZeroInflation(res_pen_high_3) ## zero no problem

#r squared
R2_pen_high_3 <- 1 - (logLik(model_pen_high_3) / logLik(model_pen_high_null))
R2_pen_high_3

## alternative model using richess as predictor (artifact)
model_pen_high_alt_1 <- glmmTMB(cbind(threatened_species,
                                     richness - threatened_species) ~ richness + mean_age + rates,
                               family = betabinomial(),
                               data = peninsula_high)

summary(model_pen_high_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response

#r squared
R2_pen_high_alt <- 1 - (logLik(model_pen_high_alt_1) / logLik(model_pen_high_null))
R2_pen_high_alt 


### alternative model using richness as offset to remove its effect
model_pen_high_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = peninsula_high
)

summary(model_pen_high_alt_2)


###############################Andalusia#######################################

##calling andalucia data
andalucia_ages_proportion <- read_csv("Data/Processed/andalucia_ages_proportion.csv")

##merging ages and diversification rates
andalucia_total <- andalucia_ages_proportion %>%
                left_join(plant_genus_rates, by = "Genus") 

##saving 
write_csv(andalucia_total, 
          file = "Data/Processed/andalucia_merged_iucn_clads.csv")

##reading
andalucia_total <- read_csv("Data/Processed/andalucia_merged_iucn_clads.csv")



# Low extinction ----------------------------------------------------------
andalucia_low <- andalucia_total %>% filter(ext_fraction == "low_ex")

##null model 
model_andalucia_low_null <- glmmTMB(cbind(threatened_species,
                                          richness - threatened_species) ~ 1,
                                    family = betabinomial(),
                                    data = andalucia_low)
summary(model_andalucia_low_null) 


#### complete model
model_andalucia_low_1 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age + rates,
                           family = betabinomial(),
                           data = andalucia_low)
summary(model_andalucia_low_1)

#r squared
R2_andalucia_low_1 <- 1 - (logLik(model_andalucia_low_1) / logLik(model_andalucia_low_null))
R2_andalucia_low_1

## Diagnostics
res_andalucia_low_1 <- simulateResiduals(model_andalucia_low_1)
plot(res_andalucia_low_1)#no problems detected
testDispersion(res_andalucia_low_1) ## no overdispersion of the residuals
testZeroInflation(res_andalucia_low_1) ## zero no problem


### model with rates only 
model_andalucia_low_2 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ rates,
                           family = betabinomial(),
                           data = andalucia_low)

summary(model_andalucia_low_2)

#r squared
R2_andalucia_low_2 <- 1 - (logLik(model_andalucia_low_2) / logLik(model_andalucia_low_null))
R2_andalucia_low_2

## Diagnostics
res_andalucia_low_2 <- simulateResiduals(model_andalucia_low_2)
plot(res_andalucia_low_2)#no problems detected
testDispersion(res_andalucia_low_2) ## no problem with the dispersion of the residuals
testZeroInflation(res_andalucia_low_2) ## zero no problem

### model with age only 
model_andalucia_low_3 <- glmmTMB(cbind(threatened_species,
                                 richness - threatened_species) ~ mean_age,
                           family = betabinomial(),
                           data = andalucia_low)

summary(model_andalucia_low_3) 

#r squared
R2_andalucia_low_3 <- 1 - (logLik(model_andalucia_low_3) / logLik(model_andalucia_low_null))
R2_andalucia_low_3

## Diagnostics
res_andalucia_low_3 <- simulateResiduals(model_andalucia_low_3)
plot(res_andalucia_low_3)#no problems detected
testDispersion(res_andalucia_low_3) ## no problem with the dispersion of the residuals
testZeroInflation(res_andalucia_low_3) ## zero no problem


## alternative model using richess as predictor (artifact)
model_andalucia_low_alt_1 <- glmmTMB(cbind(threatened_species,
                                     richness - threatened_species) ~ richness + mean_age + rates,
                               family = betabinomial(),
                               data = andalucia_low)

summary(model_andalucia_low_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response
#r squared
R2_andalucia_low_alt_1 <- 1 - (logLik(model_andalucia_low_alt_1) / logLik(model_andalucia_low_null))
R2_andalucia_low_alt_1


### alternative model using richness as offset to remove its effect
model_andalucia_low_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = andalucia_low
)

summary(model_andalucia_low_alt_2) ### in any case mean age is significant


# Intermediate extinction -------------------------------------------------

andalucia_int <- andalucia_total %>% filter(ext_fraction == "int_ex")

##null model 
model_andalucia_int_null <- glmmTMB(cbind(threatened_species,
                                          richness - threatened_species) ~ 1,
                                    family = betabinomial(),
                                    data = andalucia_int)
summary(model_andalucia_int_null)

#### complete model
model_andalucia_int_1 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ mean_age + rates,
                                 family = betabinomial(),
                                 data = andalucia_int)
summary(model_andalucia_int_1)

#r squared
R2_andalucia_int_1 <- 1 - (logLik(model_andalucia_int_1) / logLik(model_andalucia_int_null))
R2_andalucia_int_1

## Diagnostics
res_andalucia_int_1 <- simulateResiduals(model_andalucia_int_1)
plot(res_andalucia_int_1)#no problems detected
testDispersion(res_andalucia_int_1) ## no problem with the dispersion of the residuals
testZeroInflation(res_andalucia_int_1) ## zero no problem


### model with rates only 
model_andalucia_int_2 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ rates,
                                 family = betabinomial(),
                                 data = andalucia_int)

summary(model_andalucia_int_2)

#r squared
R2_andalucia_int_2 <- 1 - (logLik(model_andalucia_int_2) / logLik(model_andalucia_int_null))
R2_andalucia_int_2

## Diagnostics
res_andalucia_int_2 <- simulateResiduals(model_andalucia_int_2)
plot(res_andalucia_int_2)#no problems detected
testDispersion(res_andalucia_int_2) ## no problem with the dispersion of the residuals
testZeroInflation(res_andalucia_int_2) ## zero no problem

### model with age only 
model_andalucia_int_3 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ mean_age,
                                 family = betabinomial(),
                                 data = andalucia_int)

summary(model_andalucia_int_3) 

## Diagnostics
res_andalucia_int_3 <- simulateResiduals(model_andalucia_int_3)
plot(res_andalucia_int_3)#no problems detected
testDispersion(res_andalucia_int_3) ## no problem with the dispersion of the residuals
testZeroInflation(res_andalucia_int_3) ## zero no problem

#r squared
R2_andalucia_int_3 <- 1 - (logLik(model_andalucia_int_3) / logLik(model_andalucia_int_null))
R2_andalucia_int_3

## alternative model using richess as predictor (artifact)
model_andalucia_int_alt_1 <- glmmTMB(cbind(threatened_species,
                                           richness - threatened_species) ~ richness + mean_age + rates,
                                     family = betabinomial(),
                                     data = andalucia_int)

summary(model_andalucia_int_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response
#r squared
R2_andalucia_int_alt_1 <- 1 - (logLik(model_andalucia_int_alt_1) / logLik(model_andalucia_int_null))
R2_andalucia_int_alt_1

### alternative model using richness as offset to remove its effect
model_andalucia_int_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = andalucia_int
)

summary(model_andalucia_int_alt_2) ### in any case mean age is significant


# High extinction ---------------------------------------------------------
andalucia_high <- andalucia_total %>% filter(ext_fraction == "high_ex")

##null model 
model_andalucia_high_null <- glmmTMB(cbind(threatened_species,
                                           richness - threatened_species) ~ 1,
                                     family = betabinomial(),
                                     data = andalucia_high)
summary(model_andalucia_high_null)

#### complete model
model_andalucia_high_1 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ mean_age + rates,
                                 family = betabinomial(),
                                 data = andalucia_high)
summary(model_andalucia_high_1)

# McFadden R2
R2_andalucia_high_1 <- 1 - (logLik(model_andalucia_high_1) / logLik(model_andalucia_high_null))
R2_andalucia_high_1

### model with rates only 
model_andalucia_high_2 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ rates,
                                 family = betabinomial(),
                                 data = andalucia_high)

summary(model_andalucia_high_2)

# McFadden R2
R2_andalucia_high_2 <- 1 - (logLik(model_andalucia_high_2) / logLik(model_andalucia_high_null))
R2_andalucia_high_2

### model with age only 
model_andalucia_high_3 <- glmmTMB(cbind(threatened_species,
                                       richness - threatened_species) ~ mean_age,
                                 family = betabinomial(),
                                 data = andalucia_high)

summary(model_andalucia_high_3) 

# McFadden R2
R2_andalucia_high_3 <- 1 - (logLik(model_andalucia_high_3) / logLik(model_andalucia_high_null))
R2_andalucia_high_3

## alternative model using richess as predictor (artifact)
model_andalucia_high_alt_1 <- glmmTMB(cbind(threatened_species,
                                           richness - threatened_species) ~ richness + mean_age + rates,
                                     family = betabinomial(),
                                     data = andalucia_high)

summary(model_andalucia_high_alt_1) 
##richness is a strong predictor, however we suspect that it could be due
#a statistical artifact, given that richness is part of the response
#r squared
R2_andalucia_high_alt_1 <- 1 - (logLik(model_andalucia_high_alt_1) / logLik(model_andalucia_high_null))
R2_andalucia_high_alt_1

### alternative model using richness as offset to remove its effect
model_andalucia_high_alt_2 <- glmmTMB(
  threatened_species ~ mean_age + rates + offset(log(richness)),
  family = nbinom2,
  data = andalucia_high
)

summary(model_andalucia_high_alt_2) ### in any case mean age is significant


########## saving ############
save(model_pen_low_1, model_pen_low_2, model_pen_low_3, model_pen_low_null, 
     model_pen_low_alt_1, model_pen_low_alt_2,
     model_pen_int_1, model_pen_int_2, model_pen_int_3, model_pen_int_null, 
     model_pen_int_alt_1, model_pen_int_alt_2,
     model_pen_high_1, model_pen_high_2, model_pen_high_3, model_pen_high_null, 
     model_pen_high_alt_1, model_pen_high_alt_2,
     model_andalucia_low_1, model_andalucia_low_2, model_andalucia_low_3,
     model_andalucia_low_null, 
     model_andalucia_low_alt_1, model_andalucia_low_alt_2,
     model_andalucia_int_1, model_andalucia_int_2, model_andalucia_int_3,
     model_andalucia_int_null, 
     model_andalucia_int_alt_1, model_andalucia_int_alt_2,
     model_andalucia_high_1, model_andalucia_high_2, model_andalucia_high_3,
     model_andalucia_high_null, 
     model_andalucia_high_alt_1, model_andalucia_high_alt_2,
     file = "Data/Processed/data1_4.RData")

