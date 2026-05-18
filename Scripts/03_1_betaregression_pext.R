# GLM using as response the probability of extinction ---------------------

##Libraries
library(tidyverse)
library(DescTools)
library(glmmTMB)
library(DHARMa)
library(phytools)
library(purrr)
library(patchwork)
library(betareg)
library(performance)
library(ggeffects)



# Peninsula ---------------------------------------------------------------

#reading
peninsula_total <- read_csv("Data/Processed/peninsula_merged_iucn_clads.csv")

##calling EDGE 
EDGE2_Peninsula <- read_csv(file = "Data/Processed/EDGE2_Peninsula.csv")

##Adding genus
EDGE2_Peninsula <- EDGE2_Peninsula %>% mutate(genus = str_extract(species, "^[^_]+"))

##grouping genus and obtainig the mean probability of extinction
EDGE2_peninsula_genera <- EDGE2_Peninsula %>% group_by(genus) %>% 
                                        summarise(mean_pext = mean(pext))

##left join
peninsula_total_pext <- peninsula_total %>% left_join(EDGE2_peninsula_genera,
                                                      by = c("Genus" = "genus")) %>% 
                        drop_na()

####### Betaregression 

############# low extinction #############
peninsula_total_pext_low <- peninsula_total_pext %>% filter(ext_fraction == "low_ex")

n_pen_low <- nrow(peninsula_total_pext_low)

peninsula_total_pext_low$mean_pext_adj <-
  (peninsula_total_pext_low$mean_pext * (n_pen_low - 1) + 0.5) / n_pen_low


## null
pen_pext_low_null <- betareg(mean_pext_adj ~ 1, data = peninsula_total_pext_low)

summary(pen_pext_low_null)

r2(pen_pext_low_null)


## full 
pen_pext_low_full <- betareg(mean_pext_adj ~ rates + mean_age,
                                        data = peninsula_total_pext_low)

summary(pen_pext_low_full)

r2(pen_pext_low_full)

## only rates
pen_pext_low_rates <- betareg(mean_pext_adj ~ rates,
                             data = peninsula_total_pext_low)

summary(pen_pext_low_rates)

r2(pen_pext_low_rates)

##only ages
pen_pext_low_ages <- betareg(mean_pext_adj ~ mean_age,
                              data = peninsula_total_pext_low)

summary(pen_pext_low_ages)

r2(pen_pext_low_ages)


## measuring AIC
AIC(pen_pext_low_full, pen_pext_low_rates, pen_pext_low_ages, pen_pext_low_null)
## best supported model was the pen_pext_low_ages

############# intermediate extinction #############

peninsula_total_pext_int <- peninsula_total_pext %>% filter(ext_fraction == "int_ex")

n_pen_int <- nrow(peninsula_total_pext_int)

peninsula_total_pext_int$mean_pext_adj <-
  (peninsula_total_pext_int$mean_pext * (n_pen_int - 1) + 0.5) / n_pen_int


## null
pen_pext_int_null <- betareg(mean_pext_adj ~ 1, data = peninsula_total_pext_int)

summary(pen_pext_int_null)

r2(pen_pext_int_null)


## full 
pen_pext_int_full <- betareg(mean_pext_adj ~ rates + mean_age,
                             data = peninsula_total_pext_int)

summary(pen_pext_int_full)

r2(pen_pext_int_full)

## only rates
pen_pext_int_rates <- betareg(mean_pext_adj ~ rates,
                              data = peninsula_total_pext_int)

summary(pen_pext_int_rates)

r2(pen_pext_int_rates)

##only ages
pen_pext_int_ages <- betareg(mean_pext_adj ~ mean_age,
                             data = peninsula_total_pext_int)

summary(pen_pext_int_ages)

r2(pen_pext_int_ages)


## measuring AIC
AIC(pen_pext_int_full, pen_pext_int_rates, pen_pext_int_ages, pen_pext_int_null)
## best supported model was the pen_pext_int_ages

############# high extinction #############
peninsula_total_pext_high <- peninsula_total_pext %>% filter(ext_fraction == "high_ex")

n_pen_high <- nrow(peninsula_total_pext_high)

peninsula_total_pext_high$mean_pext_adj <-
  (peninsula_total_pext_high$mean_pext * (n_pen_high - 1) + 0.5) / n_pen_high


## null
pen_pext_high_null <- betareg(mean_pext_adj ~ 1, data = peninsula_total_pext_high)

summary(pen_pext_high_null)

r2(pen_pext_high_null)


## full 
pen_pext_high_full <- betareg(mean_pext_adj ~ rates + mean_age,
                             data = peninsula_total_pext_high)

summary(pen_pext_high_full)

r2(pen_pext_high_full)

## only rates
pen_pext_high_rates <- betareg(mean_pext_adj ~ rates,
                              data = peninsula_total_pext_high)

summary(pen_pext_high_rates)

r2(pen_pext_high_rates)

##only ages
pen_pext_high_ages <- betareg(mean_pext_adj ~ mean_age,
                             data = peninsula_total_pext_high)

summary(pen_pext_high_ages)

r2(pen_pext_high_ages)


## measuring AIC
AIC(pen_pext_high_full, pen_pext_high_rates, pen_pext_high_ages, pen_pext_high_null)
## best supported model was the pen_pext_high_ages


# Andalusia ---------------------------------------------------------------

#reading
andalusia_total <- read_csv("Data/Processed/andalucia_merged_iucn_clads.csv")

##calling EDGE 
EDGE2_andalusia <- read_csv(file = "Data/Processed/EDGE2_Andalusia.csv")

##Adding genus
EDGE2_andalusia <- EDGE2_andalusia %>% mutate(genus = str_extract(species, "^[^_]+"))

##grouping genus and obtainig the mean probability of extinction
EDGE2_andalusia_genera <- EDGE2_andalusia %>% group_by(genus) %>% 
  summarise(mean_pext = mean(pext))

##left join
andalusia_total_pext <- andalusia_total %>% left_join(EDGE2_andalusia_genera,
                                                      by = c("Genus" = "genus")) %>% 
  drop_na()

####### Betaregression 


# low extinction ----------------------------------------------------------

andalusia_total_pext_low <- andalusia_total_pext %>% filter(ext_fraction == "low_ex")

n_andalusia_low <- nrow(andalusia_total_pext_low)

andalusia_total_pext_low$mean_pext_adj <-
  (andalusia_total_pext_low$mean_pext * (n_andalusia_low- 1) + 0.5) / n_andalusia_low


## null
andalusia_pext_low_null <- betareg(mean_pext_adj ~ 1, data = andalusia_total_pext_low)

summary(andalusia_pext_low_null)

r2(andalusia_pext_low_null)


## full 
andalusia_pext_low_full <- betareg(mean_pext_adj ~ rates + mean_age,
                             data = andalusia_total_pext_low)

summary(andalusia_pext_low_full)

r2(andalusia_pext_low_full)

## only rates
andalusia_pext_low_rates <- betareg(mean_pext_adj ~ rates,
                              data = andalusia_total_pext_low)

summary(andalusia_pext_low_rates)

r2(andalusia_pext_low_rates)

##only ages
andalusia_pext_low_ages <- betareg(mean_pext_adj ~ mean_age,
                             data = andalusia_total_pext_low)

summary(andalusia_pext_low_ages)

r2(andalusia_pext_low_ages)


## measuring AIC
AIC(andalusia_pext_low_full, andalusia_pext_low_rates, andalusia_pext_low_ages, andalusia_pext_low_null)
## best supported model was the andalusia_pext_low_ages

############# intermediate extinction #############

andalusia_total_pext_int <- andalusia_total_pext %>% filter(ext_fraction == "int_ex")

n_andalusia_int <- nrow(andalusia_total_pext_int)

andalusia_total_pext_int$mean_pext_adj <-
  (andalusia_total_pext_int$mean_pext * (n_andalusia_int - 1) + 0.5) / n_andalusia_int


## null
andalusia_pext_int_null <- betareg(mean_pext_adj ~ 1, data = andalusia_total_pext_int)

summary(andalusia_pext_int_null)

r2(andalusia_pext_int_null)


## full 
andalusia_pext_int_full <- betareg(mean_pext_adj ~ rates + mean_age,
                             data = andalusia_total_pext_int)

summary(andalusia_pext_int_full)

r2(andalusia_pext_int_full)

## only rates
andalusia_pext_int_rates <- betareg(mean_pext_adj ~ rates,
                              data = andalusia_total_pext_int)

summary(andalusia_pext_int_rates)

r2(andalusia_pext_int_rates)

##only ages
andalusia_pext_int_ages <- betareg(mean_pext_adj ~ mean_age,
                             data = andalusia_total_pext_int)

summary(andalusia_pext_int_ages)

r2(andalusia_pext_int_ages)


## measuring AIC
AIC(andalusia_pext_int_full, andalusia_pext_int_rates, andalusia_pext_int_ages, andalusia_pext_int_null)
## best supported model was the andalusia_pext_int_ages

############# high extinction #############
andalusia_total_pext_high <- andalusia_total_pext %>% filter(ext_fraction == "high_ex")

n_andalusia_high <- nrow(andalusia_total_pext_high)

andalusia_total_pext_high$mean_pext_adj <-
  (andalusia_total_pext_high$mean_pext * (n_andalusia_high - 1) + 0.5) / n_andalusia_high


## null
andalusia_pext_high_null <- betareg(mean_pext_adj ~ 1, data = andalusia_total_pext_high)

summary(andalusia_pext_high_null)

r2(andalusia_pext_high_null)


## full 
andalusia_pext_high_full <- betareg(mean_pext_adj ~ rates + mean_age,
                              data = andalusia_total_pext_high)

summary(andalusia_pext_high_full)

r2(andalusia_pext_high_full)

## only rates
andalusia_pext_high_rates <- betareg(mean_pext_adj ~ rates,
                               data = andalusia_total_pext_high)

summary(andalusia_pext_high_rates)

r2(andalusia_pext_high_rates)

##only ages
andalusia_pext_high_ages <- betareg(mean_pext_adj ~ mean_age,
                              data = andalusia_total_pext_high)

summary(andalusia_pext_high_ages)

r2(andalusia_pext_high_ages)


## measuring AIC
AIC(andalusia_pext_high_full, andalusia_pext_high_rates, andalusia_pext_high_ages, andalusia_pext_high_null)
## best supported model was the andalusia_pext_high_ages



# Plotting ----------------------------------------------------------------

### theme
my_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(family = "serif",
                             size = 10, color = "black"),
    axis.title = element_text(family = "serif", size = 11),
    strip.text = element_text(family = "serif", size = 11, face = "bold"),
    plot.title = element_text(family = "serif", size = 12, face = "bold"),
    legend.position = "none"
  )

# Predictions ---------------------------------------------------------

# Peninsula
pred_pen_low  <- ggpredict(pen_pext_low_ages,  terms = "mean_age")
pred_pen_int  <- ggpredict(pen_pext_int_ages,  terms = "mean_age")
pred_pen_high <- ggpredict(pen_pext_high_ages, terms = "mean_age")

# Andalucía
pred_and_low  <- ggpredict(andalusia_pext_low_ages,  terms = "mean_age")
pred_and_int  <- ggpredict(andalusia_pext_int_ages,  terms = "mean_age")
pred_and_high <- ggpredict(andalusia_pext_high_ages, terms = "mean_age")

# Plot function -------------------------------------------------------

plot_fun <- function(pred, rawdata, title_text, color){
  
  ggplot() +
    
    # Raw data
    geom_point(data = rawdata,
               aes(x = mean_age, y = mean_pext_adj),
               alpha = 0.15,
               size = 1,
               color = "gray40") +
    
    # Confidence ribbon
    geom_ribbon(data = pred,
                aes(x = x,
                    ymin = conf.low,
                    ymax = conf.high),
                fill = color,
                alpha = 0.15,
                color = NA) +
    
    # Prediction line
    geom_line(data = pred,
              aes(x = x, y = predicted),
              color = color,
              linewidth = 1.2) +
    
    labs(
      x = "Genus age",
      y = "Probability of extinction",
      title = title_text
    ) +
    
    coord_cartesian(ylim = c(0, 0.25)) +
    
    my_theme
}

# Peninsula plots -----------------------------------------------------

p1 <- plot_fun(pred_pen_low,
               peninsula_total_pext_low,
               "Peninsula – Low",
               color = "blue")

p2 <- plot_fun(pred_pen_int,
               peninsula_total_pext_int,
               "Peninsula – Intermediate",
               color = "blue")

p3 <- plot_fun(pred_pen_high,
               peninsula_total_pext_high,
               "Peninsula – High",
               color = "blue")

# Andalusia plots -----------------------------------------------------

p4 <- plot_fun(pred_and_low,
               andalusia_total_pext_low,
               "Andalusia – Low",
               color = "red")

p5 <- plot_fun(pred_and_int,
               andalusia_total_pext_int,
               "Andalusia – Intermediate",
               color = "red")

p6 <- plot_fun(pred_and_high,
               andalusia_total_pext_high,
               "Andalusia – High",
               color = "red")

# Combine plots -------------------------------------------------------

final_plot <- (p1 | p2 | p3) /
  (p4 | p5 | p6)

# Show
final_plot

# Save publication-quality SVG ---------------------------------------

#saving
svg("Figures/Figure_betaregression.svg",
    width = 24/2.54,   # convert cm → inches
    height = 16/2.54)

print(final_plot)

dev.off()
