############# GLM figures##########################################

###Libraries
library(tidyverse)
library(ggeffects)
library(gridExtra)
library(ggpubr)
library(patchwork)

# Data after running scripts from 1 to 4
load("Data/Processed/data1_4.RData")


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


##only plotting the supported model for each region and extinction scenario

# Peninsula ---------------------------------------------------------------

### low extinction (Div. rates)

##predicting
pred_pen_low <- ggpredict(model_pen_low_2, terms = "rates") 

f_pred_pen_low <- ggplot(pred_pen_low, aes(x = x, y = predicted)) +
  geom_line(color = "blue", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "blue", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Diversification rate", y = "Pred. proportion of threatened sp.")+
  ggtitle("Low extinction")+
  mynamestheme +
  theme(plot.title = element_text(family = "serif", size = (13),
                            face = "italic", hjust = 0.5))

## intermediate extinction (Div.rates)
pred_pen_int <- ggpredict(model_pen_int_2, terms = "rates")

f_pred_pen_int <- ggplot(pred_pen_int, aes(x = x, y = predicted)) +
  geom_line(color = "blue", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "blue", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Diversification rate", y = NULL
       #  "Pred. proportion of threatened sp."
       )+
  ggtitle("Intermediate extinction")+
  mynamestheme +
  theme(plot.title = element_text(family = "serif", size = (13),
                                  face = "italic", hjust = 0.5))

## high extinction (Genus ages)
pred_pen_high <- ggpredict(model_pen_high_3, terms = "mean_age")

f_pred_pen_high <- ggplot(pred_pen_high, aes(x = x, y = predicted)) +
  geom_line(color = "blue", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "blue", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Genus Age", y = NULL
         #"Pred. proportion of threatened sp." 
       )+
  ggtitle("High extinction")+
  mynamestheme +
  theme(plot.title = element_text(family = "serif", size = (13),
                                  face = "italic", hjust = 0.5))

##combining plots
peninsula_figure <- f_pred_pen_low + f_pred_pen_int + f_pred_pen_high

y <- peninsula_figure +
  plot_annotation(
    title = "Iberian Peninsula",
    theme =  theme(plot.title = element_text(family = "serif", size = (15),
                                       face = "bold", hjust = 0.5)))

##saving 
ggsave("Figures/Supplementary/fig_Peninsula_new.png",
       plot = y,
       width = 250,
       height = 110,
       units = "mm",
       dpi = 600)


# Andalusia ---------------------------------------------------------------

## low extinction (genus ages)
pred_and_low <- ggpredict(model_andalucia_low_3, terms = "mean_age")

f_pred_and_low <- ggplot(pred_and_low, aes(x = x, y = predicted)) +
  geom_line(color = "red", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Genus Age", y = "Pred. proportion of threatened sp.")+
  #ggtitle("Low extinction")+
  mynamestheme

## Intermediate extinction (genus ages)
pred_and_int <- ggpredict(model_andalucia_int_3, terms = "mean_age")

f_pred_and_int <- ggplot(pred_and_int, aes(x = x, y = predicted)) +
  geom_line(color = "red", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Genus Age", y = NULL
         #"Pred. proportion of threatened sp."
         )+
#  ggtitle("Intermediate extinction")+
  mynamestheme

## High extinction (genus ages)
pred_and_high <- ggpredict(model_andalucia_high_3, terms = "mean_age")

f_pred_and_high <- ggplot(pred_and_high, aes(x = x, y = predicted)) +
  geom_line(color = "red", linewidth= 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.15, color = NA)+
  theme_classic()+
  labs(x = "Genus Age", y = NULL
         #"Pred. proportion of threatened sp."
       )+
#  ggtitle("High extinction")+
  mynamestheme


##combininig Eastern Andalusica plots

andalucia_figure <- f_pred_and_low + f_pred_and_int + f_pred_and_high

########### plotting coefficients ################################

###########low extinction

summary_model_anda_low <- summary(model_andalucia_low_3)

coefs_anda_low <- data.frame(
  Term = rownames(summary_model_anda_low$coefficients$cond),
  Estimate = summary_model_anda_low$coefficients$cond[, "Estimate"],
  StdError = summary_model_anda_low$coefficients$cond[, "Std. Error"],
  LowerCI = summary_model_anda_low$coefficients$cond[,
                                                "Estimate"] - 1.96 * summary_model_anda_low$coefficients$cond[, "Std. Error"],
  UpperCI = summary_model_anda_low$coefficients$cond[,
                                                "Estimate"] + 1.96 * summary_model_anda_low$coefficients$cond[, "Std. Error"],
  pValue = summary_model_anda_low$coefficients$cond[, "Pr(>|z|)"]
)

##eliminating the intercept
coefs_anda_low <- coefs_anda_low[-1,]

##ordering
coefs_anda_low$Term <- factor(coefs_anda_low$Term, levels = c("mean_age"),
                              labels = c("Genus Age"),
                              ordered = TRUE)


##factor
coefs_anda_low$corrected <- "Low extinction"



######intermediate extinction
summary_model_anda_int <- summary(model_andalucia_int_3)

coefs_anda_int <- data.frame(
  Term = rownames(summary_model_anda_int$coefficients$cond),
  Estimate = summary_model_anda_int$coefficients$cond[, "Estimate"],
  StdError = summary_model_anda_int$coefficients$cond[, "Std. Error"],
  LowerCI = summary_model_anda_int$coefficients$cond[,
                                                "Estimate"] - 1.96 * summary_model_anda_int$coefficients$cond[, "Std. Error"],
  UpperCI = summary_model_anda_int$coefficients$cond[,
                                                "Estimate"] + 1.96 * summary_model_anda_int$coefficients$cond[, "Std. Error"],
  pValue = summary_model_anda_int$coefficients$cond[, "Pr(>|z|)"]
)

##eliminating the intercept
coefs_anda_int <- coefs_anda_int[-1,]

##ordering
coefs_anda_int$Term <- factor(coefs_anda_int$Term, levels = c("mean_age"),
                              labels = c("Genus Age"),
                              ordered = TRUE)

#factor
coefs_anda_int$corrected <- "Intermediate extinction"

####high extinction

summary_model_anda_high <- summary(model_andalucia_high_3)

coefs_anda_high <- data.frame(
  Term = rownames(summary_model_anda_high$coefficients$cond),
  Estimate = summary_model_anda_high$coefficients$cond[, "Estimate"],
  StdError = summary_model_anda_high$coefficients$cond[, "Std. Error"],
  LowerCI = summary_model_anda_high$coefficients$cond[,
                                                 "Estimate"] - 1.96 * summary_model_anda_high$coefficients$cond[, "Std. Error"],
  UpperCI = summary_model_anda_high$coefficients$cond[,
                                                 "Estimate"] + 1.96 * summary_model_anda_high$coefficients$cond[, "Std. Error"],
  pValue = summary_model_anda_high$coefficients$cond[, "Pr(>|z|)"]
)

##eliminating the intercept
coefs_anda_high <- coefs_anda_high[-1,]

##ordering
coefs_anda_high$Term <- factor(coefs_anda_high$Term, levels = c("mean_age"),
                               labels = c("Genus Age"),
                               ordered = TRUE)

##factor
coefs_anda_high$corrected <- "High extinction"


###binding dataframes
coefs_anda_total <- rbind(coefs_anda_low,
                          coefs_anda_int, coefs_anda_high)


##organizing
coefs_anda_total$corrected <- factor(coefs_anda_total$corrected, 
                                     levels = c("Low extinction",
                                                "Intermediate extinction",
                                                "High extinction"),
                                     ordered = TRUE)


########plot coeficients############### 
fig_coef_andalucia <- ggplot(coefs_anda_total, aes(x = Term, y = Estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = NULL, y = "Coefficient Estimate", 
       title = "Eastern Andalusia") +
  theme_bw() +
  facet_wrap(~ corrected, nrow = 1)+
  mynamestheme+
  theme(axis.text.x = element_text(family = "serif", size = (13),
                                                face = "bold", hjust = 0.5))


##combining plots

x <- fig_coef_andalucia/andalucia_figure

##saving 
ggsave("Figures/figure_andalucia_new.png",
       plot = x,
       width = 200,
       height = 170,
       units = "mm",
       dpi = 600)
