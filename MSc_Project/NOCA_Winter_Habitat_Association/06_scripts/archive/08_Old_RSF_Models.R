# Script based on Kristen Lalla, Vanessa Poirier and Barbara Frei

#.libPaths("D:/R Packages")

library(lme4)
library(ggplot2)
library(car)
library(tidyverse)
library(corrplot)
library(MuMIn)
library(cowplot)
library(ResourceSelection)
library(corrplot)
library(dplyr)
library(performance) # for Nakagawa conditional/marginal R2
library(partR2)


# =============================================================================================== #
# 1. LOADING DATA & CHECKING FOR COLINEARITY FOR 2023-2024 DATA                                   #
# =============================================================================================== #
#### Load data and check for colinearity ####

# Anne-Marie's filepaths.
obs <- read.csv("03_processed_data/noca_obs_df_mea_sd_2324.csv")
random <- read.csv("03_processed_data/noca_ran_df_mea_sd_2324.csv")

random <- random %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)

# Combine used and available data
obs <- obs %>% 
  mutate(Used = 1)
random <- random %>% 
  mutate(Used = 0)

rsf_data <- bind_rows(obs, random)

# Prepare 4 models for each of the scale including the vegetation variables.
rsf_model_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                 family = binomial(link = "logit"), data = rsf_data)


rsf_model_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                    family = binomial(link = "logit"), data = rsf_data)

summary(rsf_model_25)
summary(rsf_model_50)
summary(rsf_model_100)
summary(rsf_model_200)
# AIC_25 = 627.04
# AIC_50 = 615.33 <---
# AIC_100 = 634.65
# AIC_200 = 639.19

vif(rsf_model_25) #Check for correlation between covariates
vif(rsf_model_50)
vif(rsf_model_100)
vif(rsf_model_200)
# VIF_25 = 2.14
# VIF_50 = 2.94
# VIF_100 = 3.27
# VIF_200 = 3.72
# All VIF < 10 - no multicollinearity in the variables used.

## prepare data by scaling
noca_scaled.df <- rsf_data

noca_scaled.df[,c(6:70)]  <-  scale(noca_scaled.df[,c(6:70)], center = TRUE, scale = TRUE)

noca_scaled.df$ID <- as.factor(noca_scaled.df$ID)

#### test for scale of effect by first running full model at each scale ####

noca_full_25 <- glmer(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m + (1|ID) + (1|Winter),  data = noca_scaled.df, family = binomial)

noca_full_50 <- glmer(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m + (1|ID) + (1|Winter), data = noca_scaled.df, family = binomial)

noca_full_100 <- glmer(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m + (1|ID) + (1|Winter), data = noca_scaled.df, family = binomial)

noca_full_200 <- glmer(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m + (1|ID)+ (1|Winter),  data = noca_scaled.df, family = binomial)

# For all of the full models I get an error: boundary (singluar) fit : see help ('isSingular').
# Yet everything still runs.
noca_SoE <- AIC(noca_full_25, noca_full_50, noca_full_100, noca_full_200)

# Comparing candidate set - the full model is the best one.

noca_data_AIC <- AIC(noca_data_full_25, noca_tree_height, noca_tree_div, noca_data_null_25)
noca_data_AIC
summary(noca_data_full_25)

library(AICcmodavg) # more info AICc

models_data_25 <- list(noca_data_full_25, noca_tree_height, noca_tree_div, noca_data_null_25) #define list of models
mod.names <- c('global', 'height', 'structure', 'null') #specify model names

aictab(cand.set = models_data_25, modnames = mod.names) 

# summary of each model
summary(noca_data_full_25)
summary(noca_tree_height)
summary(noca_tree_div)
summary(noca_data_null_25)

# compare marginal vs conditional
r.squaredGLMM(noca_data_full_25, nullfx =noca_data_null_25)

#H-S GOF test 
hoslem.test(rsf_data$Used, fitted(noca_data_full_25), g = 10)
hoslem.test(rsf_data$Used, fitted(noca_tree_height), g = 10)
hoslem.test(rsf_data$Used, fitted(noca_tree_div), g = 10)
# Seems that 50m is by far the best scale (> 10.0 AIC between the other scales) so will continue at 
# that buffer

# Exploring if we do the same process but by separating the sites.
rsf_data_MBO <- rsf_data %>% 
  filter(Site %in% c(1,2))

rsf_data_BDU <- rsf_data %>% 
  filter(Site %in% c(3,4))

rsf_data_CON <- rsf_data %>% 
  filter(Site %in% c(5,6))

# Prepare 4 models for each of the scale including the vegetation variables for MBO.
rsf_MBO_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                    family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                    family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                     family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                     family = binomial(link = "logit"), data = rsf_data_MBO)

summary(rsf_MBO_25)
summary(rsf_MBO_50)
summary(rsf_MBO_100)
summary(rsf_MBO_200)
# AIC_25 = 429.11
# AIC_50 = 420.83 <---
# AIC_100 = 451.19
# AIC_200 = 454.19

vif(rsf_model_25) #Check for correlation between covariates
vif(rsf_model_50)
vif(rsf_model_100)
vif(rsf_model_200)
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_BDU_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)

summary(rsf_BDU_25)
summary(rsf_BDU_50)
summary(rsf_BDU_100)
summary(rsf_BDU_200)
# AIC_25 = 79.331
# AIC_50 = 77.184
# AIC_100 = 76.316 <---
# AIC_200 = 80.698

vif(rsf_BDU_25) #Check for correlation between covariates
vif(rsf_BDU_50)
vif(rsf_BDU_100)
vif(rsf_BDU_200)
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_CON_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)

summary(rsf_CON_25)
summary(rsf_CON_50)
summary(rsf_CON_100)
summary(rsf_CON_200)
# AIC_25 = 98.071
# AIC_50 = 99.976
# AIC_100 = 97.635
# AIC_200 = 92.986 <---

vif(rsf_CON_25) #Check for correlation between covariates
vif(rsf_CON_50)
vif(rsf_CON_100)
vif(rsf_CON_200)
# All VIF < 10 but 100 m and 200 m have strong and moderately high correlation 


# =============================================================================================== #
# 1. LOADING DATA & CHECKING FOR COLINEARITY FOR 2024-2025 DATA                                   #
# =============================================================================================== #
#### Load data and check for colinearity ####

# Anne-Marie's filepaths.
obs <- read.csv("03_processed_data/noca_obs_df_mea_sd_2425.csv")
random <- read.csv("03_processed_data/noca_ran_df_mea_sd_2425.csv")

random <- random %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)

# Combine used and available data
obs <- obs %>% 
  mutate(Used = 1)
random <- random %>% 
  mutate(Used = 0)

rsf_data <- bind_rows(obs, random)

# Prepare 4 models for each of the scale including the vegetation variables.
rsf_model_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                     family = binomial(link = "logit"), data = rsf_data)


rsf_model_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                     family = binomial(link = "logit"), data = rsf_data)

summary(rsf_model_25)
summary(rsf_model_50)
summary(rsf_model_100)
summary(rsf_model_200)
# AIC_25 = 950.56 <-
# AIC_50 = 967.69
# AIC_100 = 966.13
# AIC_200 = 967.44
# The model that explains the variation the best amongst those 4 is the one for the 25m buffers.
# The difference between 25m and all the others is consistently > 10, meaning it is way better.

vif(rsf_model_25) #Check for correlation between covariates
vif(rsf_model_50)
vif(rsf_model_100)
vif(rsf_model_200)
# VIF_25 = 2.643861
# VIF_50 = 2.385655
# VIF_100 = 2.694027
# VIF_200 = 2.777495
# All VIF < 10 - no multicollinearity in the variables used.

# Exploring if we do the same process but by separating the sites.
rsf_data_MBO <- rsf_data %>% 
  filter(Site %in% c(1,2))

rsf_data_BDU <- rsf_data %>% 
  filter(Site %in% c(3,4))

rsf_data_CON <- rsf_data %>% 
  filter(Site %in% c(5,6))

# Prepare 4 models for each of the scale including the vegetation variables for MBO.
rsf_MBO_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_MBO)

summary(rsf_MBO_25)
summary(rsf_MBO_50)
summary(rsf_MBO_100)
summary(rsf_MBO_200)
# AIC_25 = 71.671
# AIC_50 = 68.981 <- 
# AIC_100 = 73.565
# AIC_200 = 72.403
# Not a big difference between any of the models, 50m seems a little better.

vif(rsf_MBO_25) #Check for correlation between covariates
vif(rsf_MBO_50)
vif(rsf_MBO_100)
vif(rsf_MBO_200)
# VIF_25 = 4.757315
# VIF_50 = 3.825303
# VIF_100 = 3.626452
# VIF_200 = 7.447256
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_BDU_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)

summary(rsf_BDU_25)
summary(rsf_BDU_50)
summary(rsf_BDU_100)
summary(rsf_BDU_200)
# AIC_25 = 533.62 <-
# AIC_50 = 547.03
# AIC_100 = 544.13
# AIC_200 = 546.74
# The model that best explains the variation is the 25m one

vif(rsf_BDU_25) #Check for correlation between covariates
vif(rsf_BDU_50)
vif(rsf_BDU_100)
vif(rsf_BDU_200)
# VIF_25 = 2.925958
# VIF_50 = 2.786459
# VIF_100 = 3.244581
# VIF_200 = 5.394469
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_CON_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)

summary(rsf_CON_25)
summary(rsf_CON_50)
summary(rsf_CON_100)
summary(rsf_CON_200)
# AIC_25 = 349.21
# AIC_50 = 344.56
# AIC_100 = 341.37 <-
# AIC_200 = 345.24
# Not a big difference but the 100m seems to be the best followed by 50m.


vif(rsf_CON_25) #Check for correlation between covariates
vif(rsf_CON_50)
vif(rsf_CON_100)
vif(rsf_CON_200)
# VIF_25 = 2.509554
# VIF_50 = 2.995298
# VIF_100 = 2.359233
# VIF_200 = 1.481259
# All VIF < 10. 


# RENDUE ICI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## candidate models:

noca_full_50 <- glmer(Used ~ X50m_Agricole + X50m_Anthropique + X50m_Forestier + 
                        X50m_Humide + forest_edge_50 + (1|habitat) + (1|ID),  data = swth_scaled.df, family = binomial)

#swth_full <- glm(Type ~ X50m_Agricole + X50m_Anthropique + X50m_Forestier + 
#X50m_Humide + forest_edge_50,  data = swth_scaled.df, family = binomial)

swth_edge <- glmer(Type ~  forest_edge_50 + 
                     (1|habitat) + (1|ID),  data = swth_scaled.df, family = binomial)

swth_composition <- glmer(Type ~ X50m_Agricole + X50m_Anthropique + X50m_Forestier + 
                            X50m_Humide + (1|habitat) + (1|ID),  data = swth_scaled.df, family = binomial)

swth_null_50 <- glmer(Type ~ (1|habitat) + (1|ID),  data = swth_scaled.df, family = binomial)


# comparing candidate set - full model is 'best' model
swth_AIC <- AIC(swth_full_50, swth_edge, swth_composition, swth_null_50)
swth_AIC
summary(swth_full_50)


library(AICcmodavg) # more info AICc

models <- list(swth_full_50, swth_edge, swth_composition, swth_null_50) #define list of models
mod.names <- c('global', 'edge', 'composition', 'null') #specify model names

aictab(cand.set = models, modnames = mod.names) 


# summary of each model
summary(swth_full_50)
summary(swth_composition)
summary(swth_edge)
summary(swth_null_50)

# compare marginal vs conditional
r.squaredGLMM(swth_full_50, nullfx =swth_null_50) 

#H-S GOF test 

hoslem.test(swth.df$Type, fitted(swth_full_50), g = 10)
hoslem.test(swth.df$Type, fitted(swth_composition), g = 10)
hoslem.test(swth.df$Type, fitted(swth_edge), g = 10)




### creating prediction plots ###

#creating databases to predict from for each variable
pred.df <- data.frame(     
  X50m_Forestier = rep(0, 500), 
  X50m_Humide = 0,
  X50m_Agricole = 0,
  X50m_Anthropique = 0,
  forest_edge_50 = 0,
  habitat = sample(swth.df$habitat,500), # factor("ISB", levels(swth.df$habitat))
  ID = sample(swth.df$ID,500))

#use the 99% quantile instead of the full range to deal with the outliers
#data is scaled and centered so anything under or above 2.33 will be ignored  
#For lof scaled variabel UPredictions is exp() to put back on original scale 

cut_min <- function(x){ifelse(min(x) < -2.33, -2.33, min(x))}
cut_max <- function(x){ifelse(max(x) > 2.33, 2.33, max(x))}


# the uPredictions is what is used on the x axis for plotting (non-logged, non-scaled values)


Forest_pred.df <- pred.df %>% mutate(
  X50m_Forestier = seq(cut_min(swth_scaled.df$X50m_Forestier), cut_max(swth_scaled.df$X50m_Forestier), length.out=500),
  uPredictions = (X50m_Forestier * sd(swth.df$X50m_Forestier) + mean(swth.df$X50m_Forestier)))


Agro_pred.df <- pred.df %>% mutate(
  X50m_Agricole = seq(cut_min(swth_scaled.df$X50m_Agricole), cut_max(swth_scaled.df$X50m_Agricole), length.out=500),
  uPredictions = (X50m_Agricole * sd(swth.df$X50m_Agricole) + mean(swth.df$X50m_Agricole)))


Anthro_pred.df <- pred.df %>% mutate(
  X50m_Anthropique = seq(cut_min(swth_scaled.df$X50m_Anthropique), cut_max(swth_scaled.df$X50m_Anthropique), length.out=500),
  uPredictions = (X50m_Anthropique * sd(swth.df$X50m_Anthropique) + mean(swth.df$X50m_Anthropique)))


#max(swth.df$forest_edge_50) generate values to the max observed (~700)
ForestE_pred.df <- pred.df %>% mutate(
  forest_edge_50 = seq(cut_min(swth_scaled.df$forest_edge_50), cut_max(swth_scaled.df$forest_edge_50), length.out=500),
  uPredictions = (forest_edge_50 * sd(swth.df$forest_edge_50) + mean(swth.df$forest_edge_50)))


# Creating the prediction for Type (1 = Used, 0 = Random)

Forest_pred.df$Type <-  predict(swth_full_50, newdata = Forest_pred.df, type="response", re.form = NA)

ForestE_pred.df$Type <-  predict(swth_full_50, newdata = ForestE_pred.df, type="response", re.form = NA)

Agro_pred.df$Type <-  predict(swth_full_50, newdata = Agro_pred.df, type="response", re.form = NA)

Anthro_pred.df$Type <-  predict(swth_full_50, newdata = Anthro_pred.df, type="response", re.form = NA)


# create column called presence just so that will be the label for graphing later
ForestE_pred.df$Presence <- ForestE_pred.df$Type
Agro_pred.df$Presence <- Agro_pred.df$Type
Anthro_pred.df$Presence <- Anthro_pred.df$Type



# plotting parameter effects of the full model

summod<- summary(swth_full_50)

signifVariables <- rownames(summod$coefficients)[summod$coefficients[,4]<0.05]


#dataList <- list(Forest_pred.df, ForestE_pred.df, Wetland_pred.df,  Agro_pred.df, Anthro_pred.df) # all variables
dataList <- list(ForestE_pred.df,  Agro_pred.df, Anthro_pred.df) # only signif variables


#Variables <- c("Forest", "Forest_edge", "Wetland", "Agriculture", "Anthropogenic")
Variables <- c("Forest_edge", "Agriculture", "Anthropogenic")


#axisLabels <- c("% Forest", "Forest edge (m)", "% Wetland", "% Agriculture", "% Anthropogenic")
axisLabels <- c("Forest edge (m)", "% Agriculture", "% Anthropogenic")


plotList <- list()

for(i in 1:length(dataList)){
  significant <- any(grepl(Variables[i], signifVariables))
  dat <- dataList[[i]]
  plotList[[i]] <- ggplot(data=dat, aes(x=uPredictions,y=Presence)) +
    geom_line(colour = "dodgerblue", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
    scale_x_continuous(axisLabels[i], n.breaks=4) +
    scale_y_continuous(n.breaks = 2, limits=c(0,1), labels = scales::number_format(accuracy = 1))+
    theme(axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

plot_grid(plotlist=plotList, ncol=3, labels=c("A)", "B)", "C)"))


### with observed points
library(ggpubr)

foreste_plot <- ggplot(data= ForestE_pred.df, aes(x=uPredictions,y=Presence)) +
  geom_line(colour = "dodgerblue", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[1], n.breaks=4, limits = c(0,500)) +
  scale_y_continuous(n.breaks = 2, limits=c(0,1), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = forest_edge_50, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = forest_edge_plot_id, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

agro_plot <- ggplot(data= Agro_pred.df, aes(x=uPredictions,y=Presence)) +
  geom_line(colour = "dodgerblue", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[2], n.breaks=4, limits = c(0,1)) +
  scale_y_continuous(n.breaks = 2, limits=c(0,1), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = X50m_Agricole, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = X50m_Agricole, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

anthro_plot <- ggplot(data= Anthro_pred.df, aes(x=uPredictions,y=Presence)) +
  geom_line(colour = "dodgerblue", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[3], n.breaks=4, limits = c(0,0.5)) +
  scale_y_continuous(n.breaks = 2, limits=c(0,1), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = X50m_Anthropique, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = X50m_Anthropique, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggarrange(foreste_plot, agro_plot, anthro_plot, ncol=3, labels=c("A)", "B)", "C)"))


hoslem.test(swth.df$Type, fitted(swth_full), g = 10)


library(ggeffects)


ggpredict(swth_full_50, "X50m_Forestier[all]") %>% plot()





swth_full_50 <- glmer(Type ~ X50m_Agricole + X50m_Anthropique + X50m_Forestier + 
                        X50m_Humide + forest_edge_50 + (1|habitat) + (1|ID),  data = swth_scaled.df, family = binomial)


# =============================================================================================== #
# 2. LOADING DATA & CHECKING FOR COLINEARITY FOR FIRST WINTER                                     #
# =============================================================================================== #
#### Load data and check for colinearity ####

# Anne-Marie's filepaths.
obs <- read.csv("03_processed_data/noca_obs_df_mea_sd_2223.csv")
random <- read.csv("03_processed_data/noca_ran_df_mea_sd_2223.csv")

random <- random %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)

# Combine used and available data
obs <- obs %>% 
  mutate(Used = 1)
random <- random %>% 
  mutate(Used = 0)

rsf_data <- bind_rows(obs, random)

# Prepare 4 models for each of the scale including the vegetation variables.
rsf_model_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                     family = binomial(link = "logit"), data = rsf_data)


rsf_model_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                     family = binomial(link = "logit"), data = rsf_data)

summary(rsf_model_25)
summary(rsf_model_50)
summary(rsf_model_100)
summary(rsf_model_200)
# AIC_25 = 942.69
# AIC_50 = 948.2
# AIC_100 = 976.19
# AIC_200 = 983.71
# The model that explains the variation the best amongst those 4 is the one for the 25m buffers.
# It's closely followed by the 50m buffer.

vif(rsf_model_25) #Check for correlation between covariates
vif(rsf_model_50)
vif(rsf_model_100)
vif(rsf_model_200)
# VIF_25 = 2.35
# VIF_50 = 2.61
# VIF_100 = 3.61
# VIF_200 = 5.98
# All VIF < 10 - no multicollinearity in the variables used, but the 200m model is higher than the rest.

# =============================================================================================== #
# 2. LOADING DATA & CHECKING FOR COLINEARITY FOR EVERYTHING                                       #
# =============================================================================================== #
#### Load data and check for colinearity ####

# Anne-Marie's filepaths.
obs.2223 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2223.csv")
random.2223 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2223.csv")

obs.2324 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2324.csv")
random.2324 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2324.csv")

obs.2425 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2425.csv")
random.2425 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2425.csv")


random.2223 <- random.2223 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)
random.2324 <- random.2324 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)
random.2425 <- random.2425 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)

# Combine used and available data
obs.2223 <- obs.2223 %>% 
  mutate(Used = 1)
random.2223 <- random.2223 %>% 
  mutate(Used = 0)
# Combine used and available data
obs.2324 <- obs.2324 %>% 
  mutate(Used = 1)
random.2324 <- random.2324 %>% 
  mutate(Used = 0)
# Combine used and available data
obs.2425 <- obs.2425 %>% 
  mutate(Used = 1)
random.2425 <- random.2425 %>% 
  mutate(Used = 0)

obs.2223$ID <- as.factor(obs.2223$ID)
obs.2324$ID <- as.factor(obs.2324$ID)
obs.2425$ID <- as.factor(obs.2425$ID)
random.2223$ID <- as.factor(random.2223$ID)
random.2324$ID <- as.factor(random.2324$ID)
random.2425$ID <- as.factor(random.2425$ID)


rsf_data <- bind_rows(obs.2223, obs.2324, obs.2425, random.2223, random.2324, random.2425)

# Prepare 4 models for each of the scale including the vegetation variables.
rsf_model_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                    family = binomial(link = "logit"), data = rsf_data)


rsf_model_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                     family = binomial(link = "logit"), data = rsf_data)


rsf_model_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                     family = binomial(link = "logit"), data = rsf_data)

summary(rsf_model_25)
summary(rsf_model_50)
summary(rsf_model_100)
summary(rsf_model_200)
# AIC_25 = 2511.7 <-
# AIC_50 = 2533.9 
# AIC_100 = 2580.9
# AIC_200 = 2593.9
# The model that explains the variation the best amongst those 4 is the one for the 25m buffers.


vif(rsf_model_25) #Check for correlation between covariates
vif(rsf_model_50)
vif(rsf_model_100)
vif(rsf_model_200)
# VIF_25 = 2.295925
# VIF_50 = 2.497666
# VIF_100 = 2.969782
# VIF_200 = 3.670893
# All VIF < 10 - no multicollinearity in the variables used, but the 200m model is higher than the rest.

## prepare data by scaling

noca_scaled.df <- rsf_data

noca_scaled.df[,c(8:71)]  <-  scale(noca_scaled.df[,c(8:71)], center = TRUE, scale = TRUE)

noca_scaled.df$ID <- as.factor(noca_scaled.df$ID)

#### test for scale of effect by first running full model at each scale ####

noca_full_25 <- glmer(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m + Winter + Site + (1|ID),  data = noca_scaled.df, family = binomial)

noca_full_50 <- glmer(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m + (1|ID) + (1|Winter) + (1|Site), data = noca_scaled.df, family = binomial)

noca_full_100 <- glmer(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m + (1|ID) + (1|Winter) + (1|Site), data = noca_scaled.df, family = binomial)

noca_full_200 <- glmer(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m + (1|ID) + (1|Winter) + (1|Site),  data = noca_scaled.df, family = binomial)

# For all of the full models I get an error: boundary (singluar) fit : see help ('isSingular').
# Yet everything still runs.

noca_SoE <- AIC(noca_full_25, noca_full_50, noca_full_100, noca_full_200)

summary(noca_full_25)
summary(noca_full_50)
summary(noca_full_100)
summary(noca_full_200)
# Seems that 50m is by far the best scale (> 10.0 AIC between the other scales) so will continue at 
# that buffer

# Checking 
r2_nakagawa(noca_full_25)














# Exploring if we do the same process but by separating the sites.
rsf_data_MBO <- rsf_data %>% 
  filter(Site %in% c(1,2))

rsf_data_BDU <- rsf_data %>% 
  filter(Site %in% c(3,4))

rsf_data_CON <- rsf_data %>% 
  filter(Site %in% c(5,6))

# Prepare 4 models for each of the scale including the vegetation variables for MBO.
rsf_MBO_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_MBO)


rsf_MBO_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_MBO)

summary(rsf_MBO_25)
summary(rsf_MBO_50)
summary(rsf_MBO_100)
summary(rsf_MBO_200)
# AIC_25 = 437.74
# AIC_50 = 438.62
# AIC_100 = 469.92
# AIC_200 = 467.97
# The model that explains the variation the best amongst those 4 seems to be either 25 or 50 m.
# The difference between 50m and 25m is < 2, but the difference between these 2 models and the rest
# is > 10.

vif(rsf_MBO_25) #Check for correlation between covariates
vif(rsf_MBO_50)
vif(rsf_MBO_100)
vif(rsf_MBO_200)
# VIF_25 = 2.14
# VIF_50 = 2.94
# VIF_100 = 3.27
# VIF_200 = 3.72
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_BDU_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)


rsf_BDU_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_BDU)

summary(rsf_BDU_25)
summary(rsf_BDU_50)
summary(rsf_BDU_100)
summary(rsf_BDU_200)
# AIC_25 = 76.68
# AIC_50 = 76.04
# AIC_100 = 77.46
# AIC_200 = 82.79
# The only models that seem to be less relevant is the 200 m one as it has a difference of 5 with
# the others.


vif(rsf_BDU_25) #Check for correlation between covariates
vif(rsf_BDU_50)
vif(rsf_BDU_100)
vif(rsf_BDU_200)
# VIF_25 = 2.54
# VIF_50 = 3.30
# VIF_100 = 3.42
# VIF_200 = 3.34
# All VIF < 10 - no multicollinearity in the variables used.

# Prepare 4 models for each of the scale including the vegetation variables for BDU.
rsf_CON_25 <- glm(Used ~ Mean_Veg_Height_25_m + Sd_Veg_Height_25_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_50 <- glm(Used ~ Mean_Veg_Height_50_m + Sd_Veg_Height_50_m,
                  family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_100 <- glm(Used ~ Mean_Veg_Height_100_m + Sd_Veg_Height_100_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)


rsf_CON_200 <- glm(Used ~ Mean_Veg_Height_200_m + Sd_Veg_Height_200_m,
                   family = binomial(link = "logit"), data = rsf_data_CON)

summary(rsf_CON_25)
summary(rsf_CON_50)
summary(rsf_CON_100)
summary(rsf_CON_200)
# AIC_25 = 94.86
# AIC_50 = 77.44
# AIC_100 = 63.46
# AIC_200 = 40.66
# The only models that seem to be less relevant is the 200 m one as it has a difference of 5 with
# the others.


vif(rsf_CON_25) #Check for correlation between covariates
vif(rsf_CON_50)
vif(rsf_CON_100)
vif(rsf_CON_200)
# VIF_25 = 4.54
# VIF_50 = 4.84
# VIF_100 = 9.81
# VIF_200 = 7.01
# All VIF < 10 but 100 m and 200 m have strong and moderately high correlation



