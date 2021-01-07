#
#Re-analysis of the Boston data set using the new 'slm' latent class in R-INLA
#

#Load libraries
library(INLA)
library(spdep)
library(parallel)


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

#SLM model
slmm1 <- inla( log(CMEDV) ~ -1 + 
   f(idx, model = "slm", 
      args.slm = list(rho.min = rho.min, rho.max = rho.max, W = W, X = mmatrix, 
         Q.beta = Q.beta), 
      hyper = hyper.slm),
   data = as(boston.c, "data.frame"), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
)


# Use sampling to estimate fitted values
# SLM model
sampm1 <- inla.posterior.sample(1000, slmm1)

#n.areas: Number of areas
#W: adjacency amtrix (as Matrix)
fun_fitted <- function(..., n.areas, W) {
   coeffs <- idx[-c(1:n.areas)]

   rho <- rho.min + theta[2] * (rho.max - rho.min)

   # Add white noise
   XBe <- (mmatrix %*% matrix(coeffs, ncol = 1)) +
     matrix(rnorm(n.areas, 0, exp(- theta[1] / 2)), ncol = 1)
   print(dim(XBe))
   print(summary(XBe[, 1]))
   res <- as.vector(solve(Diagonal(n.areas, x = 1) - rho * W, XBe))

   return(res)
}

fitted_values <- inla.posterior.sample.eval(fun_fitted, sampm1, n.areas = 490,
  W = W)
summary(apply(fitted_values, 1, mean))
summary(apply(fitted_values, 1, sd))

# Plot fitted values
par(mfrow = c(1, 2))
plot(slmm1$summary.fitted.values[, "mean"], apply(fitted_values, 1, mean),
  xlab = "INLA", ylab = "SAMPLING")
abline(0,1)
plot(density(fitted_values[, 1]))
lines(slmm1$marginals.fitted.values[[1]], lty = 2)


# Compute impacts using sampling
#n.areas: Number of areas
#e.values: Eigenvalues of the adjacency amtrix (as Matrix) to compute impacts
#n.var: Variable to compute the impats (0 = Intercept, etc.)
#intercept: Has the model an intercept?
compute_impacts_slm <- function(..., n.areas, e.values,
  n.var, intercept = TRUE) {
  # Number of covariates
  n.cov <- length(idx) - n.areas - intercept #-1 is the intercept

  # Coefficient
  coeff <- idx[n.areas + intercept + n.var]
  
 # Spat. autocorr. coefficient
  rho <- rho.min + theta[2] * (rho.max - rho.min)

  total.impact <- coeff / (1 - rho)
  dir.impact <- sum(1 / (1 - rho * e.values)) * coeff / n.areas
  indir.impact <- total.impact - dir.impact

  return(c(dir.impact, indir.impact, total.impact))
}

# Compute impacts
n.var <- 5
impact.SLM <- inla.posterior.sample.eval(compute_impacts_slm, sampm1,
 n.areas = 490, e.values = eigen(W)$values, n.var = n.var, intercept = TRUE)

#Summary of impacts
apply(impact.SLM, 1, mean)
apply(impact.SLM, 1, sd)


# Compute impacts using ML
example(boston.c, run.dontrun=TRUE)
impacts.ML <- impacts(gp2, listw = nb2listw(boston.soi))

# Display INLA posterior marginals and ML estimates
par(mfrow = c(1, 3))
plot(density(impact.SLM[1, ]), main = "Avg. dir. imp.")
abline(v = impacts.ML$direct[n.var])
plot(density(impact.SLM[2, ]), main = "Avg. indir. imp.")
abline(v = impacts.ML$indirect[n.var])
plot(density(impact.SLM[3, ]), main = "Avg. total imp.")
abline(v = impacts.ML$total[n.var])


