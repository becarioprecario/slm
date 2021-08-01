# Simulations using a Gaussian model
# The SLM is uses as the aim is to assess the 'slm' latent effect in INLA


# Simulations based on the structure of the Boston housing dataset

#Load libraries
library(INLA)
library(spatialreg)
library(parallel)


#Load data
library(sf)
boston.tr <- st_read(system.file("shapes/boston_tracts.shp", package="spData")[1])
boston_nb <- spdep::poly2nb(boston.tr)
lw <- spdep::nb2listw(boston_nb, style="W")

#Define some indices used in the models
n <- nrow(boston.tr)
boston.tr$idx <- 1:n

#Adjcency matrix
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Simulate covariates
set.seed(1)
boston.tr$COV1 <- runif(n, -3, 3)

mmatrix <- model.matrix( ~ COV1, boston.tr)

# Simulate response
#coef_cov <- matrix(c(0, 0.5, -0.5), ncol = 1)
coef_cov <- matrix(c(-1, 1), ncol = 1)


lin_pred <- as.vector(mmatrix %*% coef_cov)

set.seed(1)
n_samples <- 100
sd_err <- 1
lin_pred_err <- sapply(1:n_samples, function(X) {
   lin_pred + rnorm(n, sd = sd_err)
})

# Response
rho <- 0.5
Y <- solve(Diagonal(x = 1, n) - rho * W, lin_pred_err)




#Model matrix for SLM models
f1 <- Y ~ COV1


#Compute some Spatial Econometrics models using ML estimation
#summary(m1 <- errorsarlm(f1, boston.c, lw))
#summary(m2 <- lagsarlm(f1, boston.c, lw))
#summary(m3 <- lagsarlm(f1, boston.c, lw, type = "mixed"))



#DEFINE PRIORS TO BE USED WITH R-INLA

#Zero-variance for Gaussian error term
zero.variance <- list(prec = list(initial = 15, fixed = TRUE))

#Compute eigenvalues for SLM model (as in Havard's code)
#e <- eigenw(lw)
#re.idx <- which(abs(Im(e)) < 1e-6)
#rho.max <- 1 / max(Re(e[re.idx]))
#rho.min <- 1 / min(Re(e[re.idx]))
#rho <- mean(c(rho.min, rho.max))

rho.min <- 0
rho.max <- 1


#
#Variance-covarinace matrix for beta coeffients' prior
#
betaprec <- 0.0001
#Standard regression model
Q.beta <- Diagonal(n = ncol(mmatrix), x = 1)
Q.beta <- betaprec * Q.beta

#This defines the range for the spatial autocorrelation parameters
# the spatial adjacenty matrix W, the associated
#matrix of covariates X and the precision matrix Q.beta for the prior
#on the coefficients of the covariates 
#
args.slm <- list(
   rho.min = rho.min ,
   rho.max = rho.max,
   W = W,
   X = mmatrix,
   Q.beta = Q.beta
)



#Definition of priors for precision of the error effect in the slm latent
#effect and the spatial autocorrelation parameter (in the (0,1) interval).
#
hyper.slm <- list(
   prec = list(
      prior = "loggamma", param = c(0.01, 0.01)),
      rho = list(initial = 0, prior = "logitbeta", param = c(1, 1))
)


# Run simulations in parallel
library(parallel)
#options(mc.cores = detectCores() - 1)
options(mc.cores = 4)

res_sims <- mclapply(1:n_samples, function(X) {
  print(X)

  # Add simulated response
  boston.tr$Y <- Y[, X]

  #SLM model
  slmm1 <- try(inla( Y ~ -1 +
    f(idx, model="slm", args.slm = args.slm, hyper = hyper.slm),
     num.threads = 2,
     data = as.data.frame(boston.tr), family = "gaussian",
     control.family = list(hyper = zero.variance) #,
     #control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
  ))

  return(slmm1)
})


# Check results for error
table(unlist(lapply(res_sims, class)))

res_sims <- res_sims [unlist(lapply(res_sims, class)) == "inla"]

# Compute coverage
params <- c(coef_cov, 
  1 / sd_err^2, # precision
  rho # rho
)

# Check coverage of a vector of values
# values: Vector of true values
# CIS: 2column matrix with credible intervals
check_coverage <- function(values, CIs) {
  sapply(1:length(values), function(i) {
    v <- values[i]
    CI <- CIs[i, ]

    return ( (v > CI[1] ) & (v < CI[2]))
 })
}

# Differences 
res_diff <- lapply(res_sims, function(obj) {
  estimates <- c(obj$summary.random$idx[-c(1:n), "mean"],
    obj$summary.hyperpar[, "mean"])

  params - estimates
})
res_diff <- do.call(rbind, res_diff)

# Average absolute differences
apply(res_diff, 2, function(X) {
  mean(abs(X))
})

# Average relative differences
sapply(1:ncol(res_diff), function(X) {
  mean(abs(res_diff[, X] / params[X]))
})


# % of coverage
res_coverage <- lapply(res_sims, function(obj) {
  CIs <- rbind(obj$summary.random$idx[-c(1:n), c(4, 6)],
    obj$summary.hyperpar[, c(3,5)])

  check_coverage(params, CIs)

})

res_coverage <- do.call (rbind, res_coverage) 

apply(res_coverage, 2, mean)


# Binomial simulations using probit link

# Compute probabilities
# IMPORTANT: We tale lin_pred_err + 1 so that the intercpet is zero;i.e.,
#  withoout covariates we have ~50% of 1's and ~50% of ~2's. 
Yprobs <- apply(Y, 2, pnorm)
# Simulate response
set.seed(1)
Ybin <- apply(Yprobs, 2, function(X) {
  rbinom(n, 1, Yprobs)
})

# Number of 1's
table(apply(Ybin, 2, sum))


# Set prec to a fixed value
hyper.slm_probit <- list(
  prec = list(initial = log(1), fixed = TRUE),#prior = "loggamma", param = c(0.01, 0.01)),
  theta2 = list(prior = "logitbeta", param = c(1, 1))#list(initial=0, prior = "logitbeta", param = c(1,1))
)

# Rmeove intercept in the model
args.slm_probit <- args.slm


res_bin_sims <- mclapply(1:n_samples, function(X) {
  print(X)

  # Add simulated response
  boston.tr$Ybin <- Ybin[, X]

  #SLM model
  slmm1 <- try(inla( Ybin ~ -1 +
    f(idx, model = "slm", args.slm = args.slm_probit, hyper = hyper.slm_probit),
     num.threads = 1,
     data = as.data.frame(boston.tr), family = "binomial",
     control.family = list(link = "probit") #,
     #control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE)
  ))

  return(slmm1)
})


# Check results for error
table(unlist(lapply(res_bin_sims, class)))

# Update params; intercept is zero
# Remove precision of slm term (is fixed to 1).
#params <- c(0, coef_cov[2], rho)
params <- c(coef_cov, rho)


# Differences 
res_bin_diff <- lapply(res_bin_sims, function(obj) {
  estimates <- c(obj$summary.random$idx[-c(1:n), "mean"],
    obj$summary.hyperpar[, "mean"])

  params - estimates
})
res_bin_diff <- do.call(rbind, res_bin_diff)

# Average absolute differences
apply(res_bin_diff, 2, function(X) {
  mean(abs(X))
})

# Average relative differences
sapply(1:ncol(res_bin_diff), function(X) {
  mean(abs(res_bin_diff[, X] / params[X]))
})



res_bin_coverage <- lapply(res_bin_sims, function(obj) {
  CIs <- rbind(obj$summary.random$idx[-c(1:n), c(4, 6)],
    obj$summary.hyperpar[, c(3,5)])

  check_coverage(params, CIs)

})

res_bin_coverage <- do.call (rbind, res_bin_coverage)

apply(res_bin_coverage, 2, mean)

