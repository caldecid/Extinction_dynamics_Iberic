# Quantile Regressions

# packages
library(quantreg)
library(broom)
library(dplyr)
library(purrr)

# Data

# Peninsula
path_pen <- "Data/Processed/peninsula_ages_total.csv"
# Andalucia
path_and <- "Data/Processed/andalucia_ages_total.csv"

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
                      legend.position = "right")

# Function to avoid repeating the code for both areas
quant_analysis <- function(path, name_file) {
  # Data
  x <- read_csv(path)
  non_data <- x %>%  filter(ext_fraction == "low_ex") %>%
    mutate(mean_age = age, ext_fraction = "Not Corrected")
  x <- non_data %>%  bind_rows(x)
  # Plot for quantile regression
  lab_lev <- c("Not corrected",
               "Low extinction",
               "Intermediate extinction",
               "High extinction")
  data_quant_pen <- x %>%
    mutate(Corrected = factor(
      ext_fraction,
      levels = unique(ext_fraction),
      labels = lab_lev,
      ordered = TRUE
    ))
  
  # Quantiles
  q10 <- seq(.1, .9, by = .1)
  
  
  # Quantile regression
  n <- length(lab_lev)
  rq_pen <- list()
  for (i in 1:n) {
    data_quant_pen_temp <- filter(data_quant_pen, Corrected == lab_lev[i])
    rq_pen_temp <- quantreg::rq(proportion_threatened ~ mean_age,
                                data = data_quant_pen_temp,
                                tau = q10)
    rq_pen_temp2 <- summary(rq_pen_temp, se = "boot")
    rq_pen[[i]] <- map_df(rq_pen_temp2, tidy_rq_summary, .id = "tau") %>%
      mutate(tau = q10) %>%
      select(-term)
    
  }
  names(rq_pen) <- lab_lev
  q_pen <- bind_rows(rq_pen, .id = "Correction")
  write_csv(q_pen, paste0("Results/Quantile_reg_", name_file, ".csv"))
  
  # Figure

  
data_quant_pen %>%
  ggplot(aes(x = mean_age, y = proportion_threatened)) +
  geom_point(col = "Gray") +
  stat_quantile(
    aes(color = after_stat(quantile)),
    quantiles = q10,
    size = 1,
    alpha = .5,
    show.legend = TRUE) +
  scale_color_gradient(low = "#000000", high = "#d9d9d9" )  + 
  labs(x = "Genus age", y = "Prop. of threatened species") +
  scale_x_log10() +
  facet_wrap(. ~ Corrected, scale = "free") +
  theme_bw() +
  mynamestheme
  
  ggsave(filename = paste0("Figures/Figure_", name_file, "_quant.png"),
         width = 30,
         height = 20,
         units = "cm")
  
  return(list(q_pen))
}

# Results

res_and <- quant_analysis(path_and, "and")
res_pen <- quant_analysis(path_pen, "pen")
