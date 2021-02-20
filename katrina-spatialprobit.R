#
# Fit models using spatialprobit (Katrina dataset)
# Impacts are computed using the samples
#

#Load libraries
library(INLA)
library(INLABMA)
library(spdep)
library(parallel)
options(mc.cores = 3)


source("utils_slm.R")

load("katrina-slm.RData")


#Compute impacts using spatialprobit (used for comparison purposes)
library(spatialprobit)

#Add lagged covariates to Katrina for SDM model
colnames(mm2) <- c(colnames(mm), paste0("lag.", colnames(mm)[-1]))
mmlag <- data.frame(mm2[, 10:17])

Katrina <- cbind(Katrina, mmlag)


fit1 <- semprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
  low_status_customers +  high_status_customers + 
  owntype_sole_proprietor + owntype_national_chain, 
  W = W1, data = Katrina, ndraw = 25000, burn.in = 5000, showProgress = TRUE)
summary(fit1)

#impacts(fit1)

fit1lag <- semprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers + 
   owntype_sole_proprietor + owntype_national_chain +
   lag.flood_depth + lag.log_medinc + lag.small_size + lag.large_size +
   lag.low_status_customers +  lag.high_status_customers +
   lag.owntype_sole_proprietor + lag.owntype_national_chain, 
   W = W1, data = Katrina, ndraw = 25000, burn.in = 5000, showProgress = TRUE)
summary(fit1lag)

#impacts(fit1lag)



fit2 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers + 
   owntype_sole_proprietor + owntype_national_chain, 
   W = W1, data = Katrina, ndraw = 25000, burn.in = 5000, showProgress = TRUE)
summary(fit2)

impacts(fit2)



fit2lag <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers +
   owntype_sole_proprietor + owntype_national_chain +
   lag.flood_depth + lag.log_medinc + lag.small_size + lag.large_size +
   lag.low_status_customers +  lag.high_status_customers +
   lag.owntype_sole_proprietor + lag.owntype_national_chain,
   W = W1, data = Katrina, ndraw = 25000, burn.in = 5000, showProgress = TRUE)
summary(fit2lag)

impacts(fit2lag)




#save(file = "katrina-spatialprobit.RData", 
#  list = c("fit1", "fit1lag", "fit2", "fit2lag"))
