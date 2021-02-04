#
#Re-analysis of the Boston data set using the new 'slm' latent class in R-INLA
#

#Load libraries
library(INLA)
library(spdep)
library(parallel)

# Load some functions for the slm model
source("utils_slm.R")


#Load data
#data(boston)
library(rgdal)
boston.tr <- readOGR(system.file("shapes/boston_tracts.shp",
  package = "spData")[1])
boston_nb <- poly2nb(boston.tr)
censored <- boston.tr$CMEDV == 50
boston.c <- boston.tr[!censored,]
boston_nb_1 <- subset(boston_nb, !censored)
lw <- nb2listw(boston_nb_1, style = "W")


#Define some indices used in the models
n <- nrow(boston.c)
boston.c$idx<-1:n

#Adjcency matrix
#W<-nb2mat(boston.soi)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

#Model matrix for SLM models
f1 <- log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)+ I(RM^2) +  AGE + 
  log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)
mmatrix <- model.matrix(f1, boston.c)
mmatrix2 <- cbind(mmatrix, create_WX(mmatrix, lw, prefix = "lag"))


#Compute some Spatial Econometrics models using ML estimation
summary(m1 <- errorsarlm(f1, boston.c, lw))
summary(m2 <- lagsarlm(f1, boston.c, lw))
summary(m3 <- lagsarlm(f1, boston.c, lw, type = "mixed"))


#DEFINE PRIORS TO BE USED WITH R-INLA

#Zero-variance for Gaussian error term
zero.variance <- list(prec = list(initial = 15, fixed = TRUE))

#Compute eigenvalues for SLM model (as in Havard's code)
e <- eigenw(lw)
re.idx <- which(abs(Im(e)) < 1e-6)
rho.max <- 1 / max(Re(e[re.idx]))
rho.min <- 1 / min(Re(e[re.idx]))
rho <- mean(c(rho.min, rho.max))



#
#Variance-covarinace matrix for beta coeffients' prior
#
betaprec <- 0.0001
#Standard regression model
Q.beta <- Diagonal(n = ncol(mmatrix), x = 1)
Q.beta <- betaprec * Q.beta 
#Regression model with lagged covariates
Q.beta2 <- Diagonal(n = ncol(mmatrix2), x = 1)
Q.beta2 <- betaprec * Q.beta2 


#This defines the range for the spatial autocorrelation parameters
# the spatial adjacenty matrix W, the associated
#matrix of covariates X and the precision matrix Q.beta for the prior
#on the coefficients of the covariates 
#
args.slm <- list(
   rho.min = rho.min ,
   rho.max = rho.max,
   W = W,
   X = matrix(0, nrow(mmatrix), 0),
   Q.beta = matrix(1, 0, 0)
)

#Definition of priors for precision of the error effect in the slm latent
#effect and the spatial autocorrelation parameter (in the (0,1) interval).
#
hyper.slm <- list(
   prec = list(
      prior = "loggamma", param = c(0.01, 0.01)),
      rho = list(initial = 0, prior = "logitbeta", param = c(1, 1))
)

#SEM model
semm1 <- inla(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)+ I(RM^2) +  
   AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT) +
   f(idx, model = "slm", args.slm = args.slm, hyper = hyper.slm),
   data = as(boston.c, "data.frame"), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
)

# SLM model
slmm1 <- inla( log(CMEDV) ~ -1 + 
   f(idx, model = "slm", 
      args.slm = list(rho.min = rho.min, rho.max = rho.max, W = W, X = mmatrix, 
         Q.beta = Q.beta), 
      hyper = hyper.slm),
   data = as(boston.c, "data.frame"), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
)
#slmm1 <- inla.hyperpar(slmm1)

# SDM model
mmatrix2 <- cbind(mmatrix, create_WX(mmatrix, lw, prefix = "lag"))
# Regression model with lagged covariates
Q.beta2 <- Diagonal(n = ncol(mmatrix2), x = 1)
Q.beta2 <- betaprec * Q.beta2

sdmm1 <- inla( log(CMEDV) ~ -1 +
   f(idx, model = "slm",
      args.slm = list(rho.min = rho.min, rho.max = rho.max, W = W, X = mmatrix2,
         Q.beta = Q.beta2),
      hyper = hyper.slm),
   data = as(boston.c, "data.frame"), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
)
#sdmm1 <- inla.hyperpar(sdmm1)



# Use sampling to estimate fitted values
# SLM model
sampm1 <- inla.posterior.sample(2000, slmm1)

fitted_values <- inla.posterior.sample.eval(fun_fitted, sampm1, W = W)
 
summary(apply(fitted_values, 1, mean))
summary(apply(fitted_values, 1, sd))

# Plot fitted values
par(mfrow = c(1, 2))
plot(slmm1$summary.fitted.values[, "mean"], apply(fitted_values, 1, mean),
  xlab = "INLA", ylab = "SAMPLING")
abline(0,1)
plot(density(fitted_values[, 1]))
lines(slmm1$marginals.fitted.values[[1]], lty = 2)


# Compute impacts
# Load MCMC data
load("Boston_MCMC_Matlab/sar_g.RData")

# Eigenvalues of W
e.values <- eigen(W)$values

for(n.var in 1:13) {
#n.var <- 13
var.name <- colnames(mmatrix)[1 + n.var]

# Compute impacts
impact.SLM <- inla.posterior.sample.eval(compute_impacts_slm, sampm1,
 n.areas = 490, e.values = e.values, n.var = n.var, intercept = TRUE)

#Summary of impacts
apply(impact.SLM, 1, mean)
apply(impact.SLM, 1, sd)


# Compute impacts using ML
impacts.ML <- impacts(m2, listw = lw)


# Display INLA posterior marginals and ML estimates of IMPACTS
par(mfrow = c(1, 3))
plot(density(impact.SLM[1, ]), main = paste0("Avg. dir. imp., ", var.name))
lines(density(sar_g_direct[, n.var]), col = "red")
abline(v = impacts.ML$direct[n.var])
plot(density(impact.SLM[2, ]), main = paste0("Avg. indir. imp., ", var.name))
lines(density(-sar_g_direct[, n.var] + sar_g_total[, n.var]), col = "red")
abline(v = impacts.ML$indirect[n.var])
plot(density(impact.SLM[3, ]), main = paste0("Avg. total imp., ", var.name))
lines(density(sar_g_total[, n.var]), col = "red")
abline(v = impacts.ML$total[n.var])
}

# Display INLA, ML and MCMC estimates of model parameters

# Fixed effects
par(mfrow = c(4,4))
for(i in 1:(13 + 1)) {
  plot(slmm1$marginals.random$idx[[n + i]], main = colnames(mmatrix)[i],
    type = "l")
  lines(density(sar_g[, i]), col = "red")
}

# Spatial autocorrelation
marg.rho <- inla.tmarginal(function(x) { rho.min + x * (rho.max - rho.min)},
  slmm1$marginals.hyperpar[[2]])
plot(marg.rho, type = "l", main = "Spat. autoc.")
lines(density(sar_g[, 15]), col = "red")


# Precision
plot(slmm1$marginals.hyperpar[[1]], type = "l", main = "Precision")
lines(density(1 / sar_g[, 16]), col = "red")


# Impact for SDM model

# MCMC
load("Boston_MCMC_Matlab/sdm_g.RData")

# ML estimates
impacts.ML.SDM <- impacts(m3, listw = lw)
sampsdmm1 <- inla.posterior.sample(2000, sdmm1)

for(n.var in 1:13) {
#n.var <- 13
var.name <- colnames(mmatrix)[1 + n.var]


# Compute impacts
n.var <- 13
impact.SDM <- inla.posterior.sample.eval(compute_impacts_slm, sampsdmm1,
 n.areas = 490, e.values = e.values, n.var = n.var, intercept = TRUE,
  lag = TRUE)

#Summary of impacts
apply(impact.SDM, 1, mean)
apply(impact.SDM, 1, sd)
impacts.ML.SDM


# Display INLA posterior marginals and ML estimates of IMPACTS
par(mfrow = c(1, 3))
plot(density(impact.SDM[1, ]), main = paste0("Avg. dir. imp., ", var.name))
lines(density(sdm_g_direct[, n.var]), col = "red")
abline(v = impacts.ML.SDM$direct[n.var])
plot(density(impact.SDM[2, ]), main = paste0("Avg. indir. imp., ", var.name))
lines(density(-sdm_g_direct[, n.var] + sdm_g_total[, n.var]), col = "red")
abline(v = impacts.ML.SDM$indirect[n.var])
plot(density(impact.SDM[3, ]), main = paste0("Avg. total imp., ", var.name))
lines(density(sdm_g_total[, n.var]), col = "red")
abline(v = impacts.ML.SDM$total[n.var])
}

# Display INLA, ML and MCMC estimates of model parameters

# Fixed effects
par(mfrow = c(4,4))
for(i in 1:(13 + 1)) {
  plot(sdmm1$marginals.random$idx[[n + i]], main = colnames(mmatrix)[i],
    type = "l")
  lines(density(sdm_g[, i]), col = "red")
}

# Spatial autocorrelation
marg.rho <- inla.tmarginal(function(x) { rho.min + x * (rho.max - rho.min)},
  sdmm1$marginals.hyperpar[[2]])
plot(marg.rho, type = "l", main = "Spat. autoc.")
lines(density(sdm_g[, 28]), col = "red")


# Precision
plot(sdmm1$marginals.hyperpar[[1]], type = "l", main = "Precision")
lines(density(1 / sdm_g[, 29]), col = "red")

