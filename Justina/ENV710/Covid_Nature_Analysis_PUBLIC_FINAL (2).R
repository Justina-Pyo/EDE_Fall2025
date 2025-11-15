

## Code to reproduce analyses and figures in 
# Spotswood et al. 2021 'Nature inequity and higher COVID-19 case rates in less green neighbourhoods in the United States' 
# DOI: https://doi.org/10.1038/s41893-021-00781-9
# Corresponding author: Erica Spotswood - San Francisco Estuary Institute
# email: ericas@sfei.org

## COVID-19 Data Analysis : 17 states by ZIP Code  

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
library(lubridate)
library(maptools)
library(MASS) # Note MASS makes dplyr 'select' inactive. Need to use dplyr::select everywhere 
library(lme4)
library(rgeos)
library(AER)
library(scales)
library(ciTools)

# Load Data ---------------------------------------------------------------

states<-readOGR("Covid_Nature_Data_FINAL.shp")

# Data Dictionary ---------------------------------------------------------

# County: County name
# STATE: State name 
# ZCTA5CE10: ZIP Code numerical code
# INTPTLAT10: Latitude
# INTPTLON10; Longitude 
# BUI: Built Up Intensity 
# PctPark: Percentage of ZIP Code occupied by publicly accessible open space 
  # Park area from original dataset created for this project combining publicly available datasets 
 # See Extended Data Table 4 for details. Area of ZIP from ACS 5 year 2018: ALAND10
 # values are population weighted 
# fips: County code 
# PopDens: Population density (from ACS 5 year 2018: AJWBE001/ALAND10) 
# propPOC: Proportion of population in ZIP Code people of color (all non-white and hispanic)
  # from ACS 5 year 2018: 1 -  (AJWVE003/AJWBE001) 
# POC: Is ZIP Code majority white or majority people of color? >50% white classified as majority white 
# Race: what is majority (>50%) race? If no race is >50%, classified as mixed
# medAge: median age (from ACS 5 year 2018: AJWCE001)
# totPop: total population size of ZIP Code (from ACS 5 year 2018: AJWBE001)
# Urban: Census designations "Urbanized Area" and "Urban Clusters" Collapsed to "Urban", "Non-Urban" classified as "Rural"
# Covid19: Covid-19 Cases per 100,000 people (downloaded and compiled from state health department websites)
# NDVI: Normalized Difference Vegetation Index 
# propWht: Proportion white population in ZIP Code (from ACS 5 year 2018: AJWVE003/AJWBE001)
# medInc: median Income (from ACS 5 year 2018: AJZAE001)
# Arid: Aridity Index
# Days: Number of days since first recorded case in county in which ZIP Code is located. Data taken from New York Times 

# Data Preparation ---------------------------

# Scale all variables 
states$NDVI_scaled <- scale(states$NDVI)
states$PctPark_scaled <- scale(states$PctPark)
states$popDens_scaled <- scale(states$popDens)
states$medAge_scaled <- scale(states$medAge)
states$Days_scaled <- scale(states$Days)
states$propWht_scaled <- scale(states$propWht)
states$propPOC_scaled <- scale(states$propPOC)
states$medIncome_scaled <- scale(states$medInc)
states$BUI_scaled <- scale(as.numeric(states$BUI))
states$AridIndex_scaled <- scale(states$Arid)

# Subset Urban only data 

# For spatial analysis (keeps data as spatial polygon layer)
urbDatSpatial <- subset(states, Urban!="Rural")

# for non-spatial analyses - makes data into dataframe object
urbDat<-states@data %>% filter(Urban!="Rural")

# [1] SAR models ---------------------------

# Summarizes correlations among variables, checks for spatial autocorrelation 
# and builds Spatial Autoregressive (SAR) models to account for spatial autocorrelation 

# Extended Data Table 5: First check correlation among variables
  round(cor(dplyr::select(urbDat, NDVI, PctPark,propPOC,medAge, medInc, Days, popDens,BUI,Arid), use = "pairwise.complete.obs"), digits = 2)

# Build Ordinary Least Squares (OLS) model  
  f1<-round(Covid19) ~ NDVI_scaled + PctPark_scaled +  
    propPOC_scaled + medIncome_scaled + medAge_scaled +
    + popDens_scaled + Days_scaled 

m <- lm(f1, data=urbDatSpatial)
summary(m)

# check variance inflation factor 
vif(m) 

# get residuals from OLS model
urbDatSpatial$residuals <- residuals(m)

# create neighborhood matrix. Queen criteria uses adjacent units to obtain weights matrix 
# Build a neighbors list 
nb <- poly2nb(urbDatSpatial, queen=TRUE)

# Add spatial weights to neighbors list. 
# option "W" standardizes the matrix 
lw <- nb2listw(nb, style="W", zero.policy = T)

# Test for spatial autocorrelation

# Moran's I
lm.morantest(m, lw, alternative="two.sided", zero.policy=TRUE)

# Geary test 
geary.test(urbDatSpatial$Covid19, lw, randomisation=TRUE, zero.policy=TRUE,
           alternative="greater", spChk=NULL, adjust.n=TRUE)
# Clear evidence of spatial autocorrelation using both Moran's I and Geary's c

# Build SAR models: lag, error, and mixed (Durban) models make 
# different assumptions about where spatial autocorrelation lies. Run all three and compare 
# Picking best model on basis of AIC 
lag <- lagsarlm(f1, data=urbDatSpatial, lw, tol.solve=1.0e-30, method="LU", zero.policy=TRUE, na.action=na.exclude) # spatial lag model
mixed <- lagsarlm(f1, data=urbDatSpatial, lw, method="LU", type="mixed", zero.policy=TRUE) # Durbin model - giving up on this one
error <- errorsarlm(f1, data=urbDatSpatial, lw, method="LU", zero.policy=TRUE) # spatial error model

AIC(lag, error, mixed) # AIC comparison - error model does best 

# Extended Data Table 7: Model summary - SAR error model 
summary(error) 


# [2] Mixed effects models ---------------------------------------------

# Build Generalized Linear Mixed Model (GLMM) with negative binomial error structure 
# and state as random effect 

# Negative binomial mixed effects model 

urban_mixed <- glmer.nb(round(Covid19) ~  NDVI_scaled +  
          PctPark_scaled + 
          propPOC_scaled +  medIncome_scaled +medAge_scaled  + 
         Days_scaled + popDens_scaled + (1|STATE), data=urbDat)

# Model has issues with identifiability. To fix this 
# Restart the fit at the previous optimum
  ss <- getME(urban_mixed,c("theta","fixef"))
  urban_mixed_2 <- update(urban_mixed,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

# Extended Data Table 1: Summary of negative binomial mixed effects  model 
  summary(urban_mixed_2)
  
# check variance inflation factor 
  vif(urban_mixed_2)
  
# Incident Rate Ratio calculation 
  
# calculate NDVI times 10 so 0.1 NDVI is the unit
urbDat$NDVI_x10 <- urbDat$NDVI*10
# run model with this instead of scaled NDVI
urban_mixed_3 <- glmer.nb(round(Covid19) ~  medIncome_scaled +
                           NDVI_x10 + medAge_scaled + PctPark_scaled +
                           propPOC_scaled + Days_scaled + popDens_scaled + (1|STATE), data=urbDat)

  ss <- getME(urban_mixed_3,c("theta","fixef"))
  urban_mixed_4 <- update(urban_mixed_3,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

fixed <- fixef(urban_mixed_4)
confnitfixed <- confint(urban_mixed_4, parm = "beta_", method = "Wald") 

# The exponentiated coefficients are also known as Incidence Rate Ratios (IRR)
IRR <- exp(cbind(fixed, confnitfixed ))

# [3] Instrument Variable Regression (2SLS) -----------------------------------

# Fitting IV model using two stage least squares 
# uses Built Up Intensity (BUI) as instrument 

sls <-ivreg(round(Covid19) ~ NDVI_scaled + PctPark_scaled +
              propPOC + medIncome_scaled + medAge_scaled + Days_scaled | 
              propPOC + medIncome_scaled + medAge_scaled + Days_scaled + 
              BUI_scaled + AridIndex_scaled, data=urbDat)

# Extended Data Table 8: summary of IV two stage least squares analysis 
summary(sls, diagnostics=TRUE)

# check variance inflation factor 
vif(sls)

# County vs ZIP code scale comparison -------------------------------------

# appears in Data section of Methods in manuscript 

# Calculate difference between Zip in county with lowest and highest values for NDVI, proportion people of color, and income  

  county<-urbDat %>% group_by(County, STATE) %>% 
    summarize(zipNum= n(),
    diffNDVI=max(NDVI)-min(NDVI), 
    diffPOC=max(propPOC)-min(propPOC), 
    diffInc=max(medInc)-min(medInc)) 

# Show mean of difference (highest - lowest) for NDVI, POC and Income 
  summary(county$diffNDVI, na.rm=TRUE)
  summary(county$diffPOC, na.rm=TRUE)
  summary(county$diffInc, na.rm=TRUE)

# Figures --------------------------------------------------------

## Fig 1d: Cases by  NDVI  

# Create NDVI quantiles  
  urbDat<- urbDat %>%
  mutate(NDVIquant = findInterval(NDVI, 
  quantile(NDVI, prob=c(0, 0.25, 0.5, 0.75)))) %>% 
  ungroup()

  urbDat$NDVIquant<-as.factor(urbDat$NDVIquant)

# summarize data by NDVI quantiles 
  my_sum<-urbDat %>% group_by(NDVIquant) %>% 
    summarise(n=n(), meanPk=mean(Covid19), sdPk=sd(Covid19)) %>% 
    mutate(se=sdPk/sqrt(n), CI=1.96*se)

# make plot 
  
  ggplot(my_sum, aes(x=NDVIquant, y=meanPk)) +
    geom_bar(stat="identity", position=position_dodge(), fill="plum3") +
    geom_errorbar( aes(x=NDVIquant, ymin=meanPk-CI, ymax=meanPk+CI), width=0.2, 
     size=1, position=position_dodge(0.9)) +
    labs(x="NDVI quantile", y="Cases per 100,000")+
    theme_minimal()+
      scale_y_continuous(labels = comma) +
      theme(axis.text=element_text(size=30), axis.title=element_text(size=36),
      legend.text=element_text(size=30),legend.title=element_text(size=36))

## Fig 1e: NDVI by  Race 

# summarize NDVI by race (POC majority vs. white majority zips) 
  my_sum<-urbDat %>% group_by(POC) %>% 
    summarise(n=n(), meanPk=mean(NDVI), sdPk=sd(NDVI)) %>% 
    mutate(se=sdPk/sqrt(n), CI=1.96*se)

# make plot 
  ggplot(my_sum, aes(x=POC, y=meanPk)) +
    geom_bar(stat="identity", position=position_dodge(), fill="#41AB5D") +
      scale_fill_brewer(palette="Set2")+
      geom_errorbar( aes(x=POC, ymin=meanPk-CI, ymax=meanPk+CI), width=0.2, 
      size=1, position=position_dodge(0.9)) +
    labs(x="", y="Greenness (NDVI)")+
    theme_minimal()+
      theme(axis.text=element_text(size=30), axis.title=element_text(size=36),
      legend.text=element_text(size=30),legend.title=element_text(size=36))

## Fig 1f: Cases by  Race 

# summarize Covid-19 cases by race (POC majority vs. white majority zips) 
  
  my_sum<-urbDat %>% group_by(POC) %>% 
    summarise(n=n(), meanPk=mean(Covid19), sdPk=sd(Covid19)) %>% 
    mutate(se=sdPk/sqrt(n), CI=1.96*se)

# make plot 
  
  ggplot(my_sum, aes(x=POC, y=meanPk)) +
    geom_bar(stat="identity", position=position_dodge(), fill="plum3") +
    geom_errorbar( aes(x=POC, ymin=meanPk-CI, ymax=meanPk+CI), width=0.2, 
      size=1, position=position_dodge(0.9)) +
    labs(x="", y="Cases per 100,000")+
    scale_y_continuous(labels = comma) +
    theme_minimal()+
    theme(axis.text=element_text(size=30), axis.title=element_text(size=36),
        legend.text=element_text(size=30),legend.title=element_text(size=36))

# Figure 2b: Coefficents for negative binomial mixed effects model 

# Extract coefficients 
v_cov <- vcov(urban_mixed_2, useScale = FALSE) # Variance-covariance matrix
betas <- fixef(urban_mixed_2)
se <- sqrt(diag(v_cov))
zval <- betas / se
pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
coefs <- (cbind(betas, se, zval, pval)[-1, ]) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "coef")
coefs$coef_num <- 1:7
coefs$labs <- c("NDVI", "Park access (ha)",  "Proportion POC", "Median income",
                "Median age", "Days since 1st Case", "Population Density")
cols <- rep(c("grey", "white"),4)[-8]

sig_cols <- rep("#000000", 7)
sig_cols[c(1, 3, 5)] <- "#8D021F"
names(sig_cols) <- coefs$coef

# make plot 
  ggplot(coefs) +
    geom_pointrange(aes(x = reorder(coef, rev(coef_num)), y = betas, ymin = betas - se, ymax = betas + se),
        size = 1.25,  fatten = 1.25,
        position = position_dodge(width = 1), show.legend = F, color = sig_cols) +
    geom_rect(aes(xmin = coef_num - 0.5, xmax = coef_num + 0.5, ymin = -0.16, ymax = 0.34, fill = reorder(coef, rev(coef_num))), 
            alpha = 0.1, show.legend = F) +
    scale_y_continuous(limits = c(-0.16, 0.34), expand = c(0, 0)) +
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

# Figure 2c and Extended Data Figure 1: Predicted Cases by NDVI, POC,and state 

# New dataset for predictions. 0 sets values to their average
  newdat_b <- expand.grid(
    NDVI_scaled=seq(-4.5, 1.2, by=0.1), # covering the full range of NDVI values from 0 to 0.8
    PctPark_scaled=0, 
    propPOC_scaled=c(-0.19, 1, 2.2), # low value (~20% POC), average (~50%), and high value (~80%)
    medIncome_scaled=0, 
    medAge_scaled=0, 
    Days_scaled=0, 
    popDens_scaled=0, 
    STATE=c("Arizona", "Arkansas", "Delaware", "Florida", 
          "Hawaii", "Illinois", "Indiana", "Maine", 
          "Maryland", "New Jersey", "New Mexico", "North Carolina", 
          "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", 
          "South Carolina"), 
  Covid19 = 0)


# Make prediction - includes random effect for state 
newdat_ci_b <- add_ci(newdat_b, urban_mixed, response = TRUE, nSims = 5) # nSims = 500 is ideal but slow

# Deal with scaling for visualization

# original variables for unscaling
cov_names <- fixef(urban_mixed) %>% 
  names() %>% 
  `[`(-1)
covs_scaled <- urbDat[ , cov_names]

# Backtransform the NDVI variable
cov_mean <- attributes(covs_scaled[["NDVI_scaled"]])$`scaled:center`
cov_sd <- attributes(covs_scaled[["NDVI_scaled"]])$`scaled:scale`

coef_values <- newdat_ci_b$NDVI_scaled
back_transformed <- coef_values * cov_sd + cov_mean
newdat_ci_b$NDVI_unscaled <- back_transformed

# sub in which states you want to see 
subber<-newdat_ci_b %>% group_by(STATE) %>%  filter(STATE==c("Arizona", "Florida"))

# Plot predictions vs NDVI for POC and White for two states 

ggplot(subber, aes(x = NDVI_unscaled, color=as.factor(propPOC_scaled))) +
  # adding the ribbon first so the line goes on top
  geom_ribbon(aes(ymax = LCB0.025, ymin = UCB0.975, fill=as.factor(propPOC_scaled)),
              alpha = 0.5) +
  geom_line(aes(y = pred)) +
  facet_wrap(~STATE) +
  labs(title = "", x = "NDVI", y = "Cases per 100k people") +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8000)) + 
  scale_color_discrete(name = "Proportion POC", breaks = c(-0.19, 1, 2.2), labels = c("Low", "Intermediate", "High")) +
  scale_fill_discrete(name = "Proportion POC", breaks = c(-0.19, 1, 2.2), labels = c("Low", "Intermediate", "High")) +
  theme(text = element_text(size = 12),
        legend.position = "top")



