###################Figures##########################################

###Libraries
library(tidyverse)
library(gridExtra)
library(ggpubr)

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
                                                            "age", "lambda_a"),
                             labels = c("richness", "mean_age",
                                        "lambda_mean"),
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
                                                            "lambda_mean"),
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
                                                            "lambda_mean"),
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
                                                              "lambda_mean"),
                              ordered = TRUE)
##corrected
coefs_pen_high$corrected <- "high_ext"


###binding dataframes
coefs_pen_total <- rbind(coefs_pen_non, coefs_pen_low, 
                         coefs_pen_int, coefs_pen_high)

##organizing
coefs_pen_total$Term <- factor(coefs_pen_total$Term, levels = c("richness",
                                                              "mean_age",
                                                              "lambda_mean"),
                               labels = c("Richness",
                                          "Genus age",
                                          "lambda"),
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

png("text/Figures/Figure_pen_coef.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)


pen_coefs <- ggplot(coefs_pen_total, aes(x = Term, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Phylogenetic traits", y = "Coefficient Estimate")+
  theme_bw() +
  ylim(-0.2, 0.2)+
  scale_x_discrete(labels = c("Richness","Genus age", expression(~lambda)))+
  facet_wrap(~ corrected, nrow = 1) +
  ggtitle("Iberian Peninsula [GLM predictors]")+
  mynamestheme
  

dev.off()

#############################GLM AGE############################

###not corrected
# Prepare peninsula_traits_non for prediction in order to plot
pred_peninsula_non_age <- data.frame(age = seq(min(peninsula_traits_non$age),
                                        max(peninsula_traits_non$age),
                                        length.out = 100),
                                        richness = mean(peninsula_traits_non$richness,
                                                        na.rm = TRUE),
                                        lambda_a = mean(peninsula_traits_non$lambda_a, na.rm = TRUE))

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
pred_peninsula_non_age <- pred_peninsula_non_age %>% rename(mean_age = age,
                                                            lambda_mean = lambda_a)


############low extinction
# Prepare for prediction in order to plot
pred_peninsula_low_age <- data.frame(mean_age = seq(min(peninsula_traits_low$mean_age),
                                               max(peninsula_traits_low$mean_age),
                                               length.out = 100),
                                     richness = mean(peninsula_traits_low$richness,
                                                     na.rm = TRUE),
                                     lambda_mean = mean(peninsula_traits_low$lambda_mean,
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
                                     lambda_mean = mean(peninsula_traits_int$lambda_mean,
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
                                     lambda_mean = mean(peninsula_traits_high$lambda_mean,
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

png("text/Figures/Figure_pen_age.png", 
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

dev.off()



###########################GLM Richness########################

###not corrected
# Prepare peninsula_traits_non for prediction in order to plot
pred_peninsula_non_richness <- data.frame(richness = seq(min(peninsula_traits_non$richness),
                                               max(peninsula_traits_non$richness),
                                               length.out = 100),
                                     age = mean(peninsula_traits_non$age,
                                                     na.rm = TRUE),
                                     lambda_a = mean(peninsula_traits_non$lambda_a, na.rm = TRUE))

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
pred_peninsula_non_richness <- pred_peninsula_non_richness %>% rename(mean_age = age,
                                                            lambda_mean = lambda_a)


############low extinction

# Prepare for prediction in order to plot
pred_peninsula_low_richness <- data.frame(richness = seq(min(peninsula_traits_low$richness),
                                                    max(peninsula_traits_low$richness),
                                                    length.out = 100),
                                     mean_age = mean(peninsula_traits_low$mean_age,
                                                     na.rm = TRUE),
                                     lambda_mean = mean(peninsula_traits_low$lambda_mean,
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
                                     lambda_mean = mean(peninsula_traits_int$lambda_mean,
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
                                      lambda_mean = mean(peninsula_traits_high$lambda_mean,
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

png("text/Figures/Figure_pen_richness.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)

pen_richness <- ggplot(pred_pen_richness, aes(x = richness, y = predicted)) +
  geom_line(color = "deepskyblue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Richness", y = "Prop. of threatened species") +
  theme_bw()+
  facet_wrap(~corrected, nrow = 1)+
  
  mynamestheme

dev.off()


## merging the three figures of Iberian Peninsula

##Plotting
png("text/Figures/Figure_2.png", 
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
                                                              "lambda_a"),
                              labels = c("richness",
                                         "mean_age",
                                         "lambda_mean"),
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
                                                              "lambda_mean"),
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
                                                              "lambda_mean"),
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
                                                                "lambda_mean"),
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
png("text/Figures/Figure_3.png", 
    width = 30, height = 12,
    units = "cm", pointsize = 8, res = 300)

ggplot(coefs_anda_total, aes(x = Term, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Phylogenetic traits", y = "Coefficient Estimate", 
       title = "Eastern Andalusia") +
  theme_bw() +
  
  ylim(-0.3, 0.3)+
  scale_x_discrete(labels = c("Richness","Genus age", expression(~lambda)))+
  facet_wrap(~ corrected, nrow = 1)+
  mynamestheme

dev.off()
