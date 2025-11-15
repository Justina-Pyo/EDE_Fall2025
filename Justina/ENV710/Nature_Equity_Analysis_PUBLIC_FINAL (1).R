

## Code to reproduce analyses and figures in 
# Spotswood et al. 2021 'Nature inequity and higher COVID-19 case rates in less green neighbourhoods in the United States' 
# DOI: https://doi.org/10.1038/s41893-021-00781-9
# Corresponding author: Erica Spotswood - San Francisco Estuary Institute
# email: ericas@sfei.org

# Nature- Equity Analysis for all urbanized areas in the US 

# Load packages ---------------------------------------------------------------

if (!require("rspatial")) devtools::install_github('rspatial/rspatial')
library(rspatial) 
library(raster)
library(sp)
library(rgdal)
library(spdep)
library(spatialreg)
library(dplyr)
library(car) 
library(ggplot2)
library(magrittr)
library(janitor)
library(showtext)

# Load data 

USdat<-readOGR("Nature_Equity_Data_FINAL.shp")

#====================
# data Dictionary 
#====================

# Data compiled from multiple sources. See manuscript main text for additional detail. 

# Data from US Census Bureau, park database, and other sources. 

# UA_ID: Urbanized Area ID - from US Census Bureau 
# STATE: State
# STATEFP: State code - from US Census Bureau 
# COUNTY: County 
# TRACTCE: Census tract code - from US census Bureau 
# BLKGRP: Block group code - from US census bureau 
# NAME: Block group name - from US census bureau 
# INTLAT: Latitude 
# INTLON: Longitude 
# OBJID: Unique object ID
# NDVI_US: Average NDVI value for block group  
# ParkAc: Acres of park within 1,000 m of the centroid of census blocks (population weighted and averaged to block group) 
# ParkPct: Percent of 1,000 m circle around centroid of block containing parks 
# BUI = Built up intensity value 
# Arid = aridity index 

# ACS 2018 5 year data: 

# ALAND: Land area of block group 
# ToPop_e: Total population estimate 
# ToPop_m: Total population margin of error 
# PCInc_e: Per capita income estimate 
# PCInc_m: Per capita income margine of error 
# White_e: Proportion of block group non-hispanic white - estimate 
# White_m: Proportion of block group non-hispanic white - margin of erroro  
# POC_e: Proportion of block group people of color - estimate 
# POC_m: Proportion of block group people of color - margin of error  

#===============
# DATA PREP 
#===============

## Race 

# calculate proportion white
USdat$propWhite<-USdat$White_e/USdat$ToPop_e
USdat$propWhite[is.na(USdat$propWhite)]<-0
# replace NA values with 0 = NAs produced when BG has 0population (divide by zero issue)

# Make BUI numeric 
USdat$BUI<-as.numeric(USdat$BUI)

# create categorical variable for race (majority white, or majority non-white)
USdat$raceCat<-ifelse(USdat$propWhite > 0.5, "majority white", "majority POC")

# create proportion POC (inverse of white)
USdat$propPOC<- 1 - USdat$propWhite

## Population density 

#Convert ALAND (in square meters) to Ha (we want population density in people/ha)
USdat$areaHA<-USdat$ALAND/10000

# calculate population density 
USdat$popDens<-USdat$ToPop_e/USdat$areaHA
USdat$popDens[is.na(USdat$popDens)]<-0

# Scale all variables 
USdat$BUI_scaled<-scale(USdat$BUI)
USdat$propWhite_scaled<-scale(USdat$propWhite)
USdat$propPOC_scaled<-scale(USdat$propPOC)
USdat$popDens_scaled<-scale(USdat$popDens)
USdat$medInc_scaled<-scale(USdat$PCInc_e)
USdat$Arid_scaled<-scale(USdat$Arid)

# make field converting park area from acres to ha
USdat$ParkHa<-USdat$ParkAc/2.471

# make quantiles for income unique to each urbanized area 

USdat@data <- USdat@data %>%
  group_by(UA_ID) %>%   mutate(IncomeQuant = findInterval(PCInc_e, 
  quantile(PCInc_e, prob=c(0, 0.25, 0.5, 0.75))))

USdat$IncomeQuant<-as.factor(USdat$IncomeQuant)


# ================================  
# NDVI Equity Analysis  
# =================================
  

# Checking Correlations

USdat_pairs<-USdat@data[,c(12:14,16, 18, 20,26,24:25, 29)]
round(cor(USdat_pairs, use = "pairwise.complete.obs"), digits = 2)

# Building Ordinary Least Squares model  

f1<-NDVI_US ~ propPOC_scaled + medInc_scaled + popDens_scaled + Arid_scaled
m <- lm(f1, data=USdat)
summary(m)

# Variance inflation factor 
vif(m) 

# get residuals from OLS model
USdat$residuals <- residuals(m)

# create neighborhood matrix with queen criteria which uses adjacent units to obtain weights matrix 
nb <- poly2nb(USdat, queen=TRUE)

# Add spatial weights to neighbors list. option "W" standardizes the matrix 
lw <- nb2listw(nb, style="W", zero.policy = T)

# Moran's I - alternative version - looks like it's actually not appropriate to 
# use the mc method for residuals (it's for raw variables instead)
# https://r-spatial.github.io/spdep/reference/lm.morantest.html
lm.morantest(m, lw, alternative="two.sided", zero.policy=TRUE)

# lagrange multiplier test -- to determine which SAR model to use 

LM<-lm.LMtests(m, lw, test="all", zero.policy=TRUE)

# significant for all -- conclusion is to run all models, and use fit to decide which one is best

lag <- lagsarlm(f1, data=USdat, lw, tol.solve=1.0e-30, method="LU", zero.policy=TRUE, na.action=na.exclude) # spatial lag model
error <- errorsarlm(f1, data=USdat, lw, method="LU", zero.policy=TRUE) # spatial error model

# compare models (OLS, lag, error)
AIC(m, lag, error)

# error model has best fit
summary(error)

#=====================================
# PARK ANALYSIS 
# =================================

# Start with OLS model  

f2<-ParkHa ~ propPOC_scaled + medInc_scaled + popDens_scaled 
m <- lm(f2, data=USdat)
summary(m)

# calculate variance inflation factor 
vif(m) 

# get residuals from OLS model
USdat$residuals2 <- residuals(m)

# build a neighbors list 
nb <- poly2nb(USdat, queen=TRUE)

# add spatial weights to neighbors list. 
# option "W" standardizes the matrix 
lw <- nb2listw(nb, style="W", zero.policy = T)

# Calculate Moran's I to evaluate presence of spatial autocorrelation 
lm.morantest(m, lw, alternative="two.sided", zero.policy=TRUE)

# lagrange multiplier test -- to determine which SAR model to use 

LM<-lm.LMtests(m, lw, test="all", zero.policy=TRUE)

# significant for all -- conclusion is to run all models, and use fit to decide which one is best

# lag and mixed models have convergence issues. Use error model only.  

error_park <- errorsarlm(f2, data=USdat, lw, method="LU", zero.policy=TRUE) # spatial error model

summary(error_park)

############################################# 
# FIGURES FOR PAPER 
############################################

# Fig 3b: NDVI vs. Race    
#-----------------

my_sum<-USdat@data %>% group_by(raceCat) %>% 
  summarise(n=n(), meanPk=mean(NDVI_US), sdPk=sd(NDVI_US)) %>% 
  mutate(se=sdPk/sqrt(n), CI=2*se) %>% 
  ungroup()

ggplot(my_sum, aes(x=raceCat, y=meanPk)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  geom_errorbar( aes(x=raceCat, ymin=meanPk-CI, ymax=meanPk+CI), width=0.3, 
                 colour="black", size=1, position=position_dodge(0.9)) +
  labs(x="", y="Greenness (NDVI)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.title=element_text(size=18))

# Fig 3c: Park Acres vs. Race    
#-----------------

my_sum<-USdat@data %>% group_by(raceCat) %>% 
  summarise(n=n(), meanPk=mean(ParkHa), sdPk=sd(ParkHa)) %>% 
  mutate(se=sdPk/sqrt(n), CI=2*se) %>% 
  ungroup()

ggplot(my_sum, aes(x=raceCat, y=meanPk)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  geom_errorbar( aes(x=raceCat, ymin=meanPk-CI, ymax=meanPk+CI), width=0.3, 
                 colour="black", size=1, position=position_dodge(0.9)) +
  labs(x="", y="Park access (ha)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.title=element_text(size=18))


# Fig 3d: NDVI vs. Income    
#-----------------

my_sum<-USdat@data %>% group_by(IncomeQuant) %>% 
  summarise(n=n(), meanPk=mean(NDVI_US), sdPk=sd(NDVI_US)) %>% 
  mutate(se=sdPk/sqrt(n), CI=2*se) %>% 
  ungroup()

ggplot(my_sum, aes(x=IncomeQuant, y=meanPk)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  geom_errorbar( aes(x=IncomeQuant, ymin=meanPk-CI, ymax=meanPk+CI), width=0.3, 
                 colour="black", size=1, position=position_dodge(0.9)) +
  labs(x="Income quantile", y="Greenness (NDVI)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.title=element_text(size=18))


# Fig 3e: Park access vs. Income    
#-----------------

my_sum<-USdat@data %>% group_by(IncomeQuant) %>% 
  summarise(n=n(), meanPk=mean(ParkHa), sdPk=sd(ParkHa)) %>% 
  mutate(se=sdPk/sqrt(n), CI=2*se) %>% 
  ungroup()

ggplot(my_sum, aes(x=IncomeQuant, y=meanPk)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  geom_errorbar( aes(x=IncomeQuant, ymin=meanPk-CI, ymax=meanPk+CI), width=0.3, 
                 colour="black", size=1, position=position_dodge(0.9)) +
  labs(x="Income quantile", y="Park access (ha)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.title=element_text(size=18))


# Extended Data Fig. 2: Coefficient values from SAR models  
#-----------------

# Extended Data Fig. 2a: NDVI equity 

# Extract coefficients

coefs <- broom::tidy(error) %>% 
  janitor::clean_names()
coefs$upper <- coefs$estimate + coefs$std_error * 1.96
coefs$lower <- coefs$estimate - coefs$std_error * 1.96

coefs %<>% 
  `[`(c(2, 3,4,5), )

coefs$coef_num <- 1:4
coefs$labs <- c("Proportion POC", "Median income", "Population density", "Aridity")
cols <- rep(c("grey", "white"),2)
sig_cols <- rep("#000000", 4)
sig_cols[] <- "#8D021F"

# Make figure  

ggplot(coefs) +
  geom_pointrange(aes(x = reorder(term, rev(coef_num)), y = estimate, ymin = lower, ymax = upper),
                  size = 1.25,  fatten = 1.25,
                  position = position_dodge(width = 1), show.legend = F, color = sig_cols) +
  geom_rect(aes(xmin = coef_num - 0.5, xmax = coef_num + 0.5, ymin = -0.05, ymax = 0.1, fill = reorder(term, rev(coef_num))), 
            alpha = 0.1, show.legend = F) +
  scale_y_continuous(limits = c(-0.05, 0.1), expand = c(0, 0)) +
  scale_x_discrete(labels = rev(coefs$labs)) +
  scale_fill_manual(values = cols)+
  ylab("Coefficient value") +
  xlab("") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  theme_classic() + 
  # Rescale text to whatever dpi is being used if exporting as a PNG, e.g., size = <font size>*<dpi>/72
  # Not applicable if exporting as PDF
  theme(text = element_text(family = "Myriad Pro", size = 12), 
        panel.border = element_rect(colour = "black"), rect = element_rect(fill = "transparent")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip(clip = "off")


# Extended Data Fig. 2b: Park access equity 

# Extract coefficients

coefs <- broom::tidy(error_park) %>% 
  janitor::clean_names()
coefs$upper <- coefs$estimate + coefs$std_error * 1.96
coefs$lower <- coefs$estimate - coefs$std_error * 1.96

coefs %<>% 
  `[`(c(2, 3,4), )

coefs$coef_num <- 1:3
coefs$labs <- c("Proportion POC", "Median income", "Population density")
cols <- rep(c("grey", "white"),2)[-4]
sig_cols <- rep("#000000", 3)
sig_cols[] <- "#8D021F"

# Make figure  

ggplot(coefs) +
  geom_pointrange(aes(x = reorder(term, rev(coef_num)), y = estimate, ymin = lower, ymax = upper),
          size = 1.25,  fatten = 1.25,
          position = position_dodge(width = 1), show.legend = F, color = sig_cols) +
  geom_rect(aes(xmin = coef_num - 0.5, xmax = coef_num + 0.5, ymin = -0.5, ymax = 1.1, 
                fill = reorder(term, rev(coef_num))), 
            alpha = 0.1, show.legend = F) +
  scale_y_continuous(limits = c(-0.5, 1.1), expand = c(0, 0)) +
  scale_x_discrete(labels = rev(coefs$labs)) +
  scale_fill_manual(values = cols)+
  ylab("Coefficient value") +
  xlab("") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  theme_classic() + 
  # Rescale text to whatever dpi is being used if exporting as a PNG, e.g., size = <font size>*<dpi>/72
  # Not applicable if exporting as PDF
  theme(text = element_text(family = "Myriad Pro", size = 12), 
        panel.border = element_rect(colour = "black"), rect = element_rect(fill = "transparent")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip(clip = "off")



