# Some functions to compute fitted values and impacts from inla models
#   fitted using the slm latent effect

# Function to re-estimate the fitted values using sampling
#
#W: adjacency amtrix (as Matrix)
fun_fitted <- function(..., W) {

   n.areas <- nrow(W)

   coeffs <- idx[-c(1:n.areas)]

   rho <- rho.min + theta[2] * (rho.max - rho.min)

   # Add white noise
   XBe <- (mmatrix %*% matrix(coeffs, ncol = 1)) +
     matrix(rnorm(n.areas, 0, exp(- theta[1] / 2)), ncol = 1)
   
   
   res <- as.vector(solve(Diagonal(n.areas, x = 1) - rho * W, XBe))

   return(res)
}

# Compute impacts using sampling
#
#n.areas: Number of areas
#e.values: Eigenvalues of the adjacency matrix (as Matrix) to compute impacts
#n.var: Variable to compute the impats (0 = Intercept, etc.)
#intercept: Has the model an intercept?
#lag: Whether lagged variables are included; in the linear predictor: covariates, lagged covariates (in the same order)
compute_impacts_slm <- function(..., n.areas, e.values,
  n.var, intercept = TRUE, lag = FALSE) {
  # Number of covariates + lagged covariates; used for lagged covariates
  n.cov <- length(idx) - n.areas - intercept #-1 is the intercept

  # Coefficient
  coeff <- idx[n.areas + intercept + n.var]
  
 # Spat. autocorr. coefficient
  rho <- rho.min + theta[2] * (rho.max - rho.min)

  # Impacts
  total.impact <- coeff / (1 - rho)
  dir.impact <- sum(1 / (1 - rho * e.values)) * coeff / n.areas

  # If lagged covariates
  if(lag) {
    # VAlue of the coefficient of the lagged covariate
    lag.coeff <- idx[n.areas + intercept + n.cov / 2  + n.var]

    # Compute trace of (I - rho * W)^{-1} W using a Taylor expansion of order 10 
    rho.e.values <- e.values
    trIrhoWW <- sum(e.values)
    for(i in 1:10) {
      rho.e.values <- rho * e.values * rho.e.values
      trIrhoWW  <- trIrhoWW  + sum(rho.e.values)
    }

    # Update impacts
    total.impact <- total.impact + lag.coeff / (1 - rho)
    dir.impact <- dir.impact + trIrhoWW * lag.coeff / n.areas
  }

  # Indirect impacts
  indir.impact <- total.impact - dir.impact

  return(c(dir.impact, indir.impact, total.impact))
}



# Compute impacts for first covariate
compute_impacts_sem_probit <- function(..., n.areas,
  n.var, intercept = TRUE, lag = FALSE, mmatrix) {

  # Number of covariates + lagged covariates; used for lagged covariates
  #n.cov <- length(idx) - n.areas - intercept #-1 is the intercept
  n.cov <- ncol(mmatrix) - intercept

  # Coefficient
  coeff <- eval(parse(text = colnames(mmatrix)[intercept + n.var]))
  # Covariate
  #X <- mmatrix[, intercept + n.var]

  eta <- Predictor[1:n.areas]

  # Spat. autocorr. coefficient (NOT USED for SEM)
  #rho <- rho.min + theta[2] * (rho.max - rho.min)

  if(!lag) {
    # Impacts
    #eta <- X * coeff
    total.impact <- mean(dnorm(eta)) * coeff
    dir.impact <- total.impact
    #dir.impact <- sum(1 / (1 - rho * e.values)) * coeff / n.areas
  } else { # If lagged covariates
    # Value of the coefficient of the lagged covariate
    coeff_lag <- eval(parse(text = colnames(mmatrix)[intercept + n.cov / 2 + n.var]))
    #X_lag <- mmatrix[, intercept + n.cov / 2 + n.var]

    # 'linear pedictor' eta
    #eta <- X * coeff + X_lag * coeff_lag
    # Impacts
    dir.impact <- mean(dnorm(eta)) * coeff
    total.impact <- mean(dnorm(eta)) * (coeff + coeff_lag)
  }

  # Indirect impacts
  indir.impact <- total.impact - dir.impact

  return(c(dir.impact, indir.impact, total.impact))
}


# Compute all impacts
# inla_samples: Samples from inla.posterior.sample
# n.var: NOW the number of covariates to compute the impacts
compute_impacts_sem_probit_all <- function(inla_samples, n.areas,
  n.var, intercept = TRUE, lag = FALSE, mmatrix) {

  res <- data.frame(sapply(1:n.var, function(X) {
    impacts <- inla.posterior.sample.eval(compute_impacts_sem_probit,
      inla_samples, n.areas = n.areas, n.var = X, lag = lag,
      mmatrix = mmatrix)

    impacts <- c(apply(impacts, 1, mean), apply(impacts, 1, sd))
    return(impacts)
  }))
  

  colnames(res) <- colnames(mmatrix)[intercept + 1:n.var]
  rownames(res) <- c("direct (mean)", "indirect (mean)", "total (mean)",
    "direct (s.d.)", "indirect (s.d.)", "total (s.d.)")

  return(res)
}

# Compute impacts for first covariate
compute_impacts_slm_probit <- function(..., n.areas, W,
  n.var, intercept = TRUE, lag = FALSE, mmatrix) {

  # Number of covariates + lagged covariates; used for lagged covariates
  #n.cov <- length(idx) - n.areas - intercept #-1 is the intercept
  n.cov <- ncol(mmatrix) - intercept

  # Coefficient
  coeff <- idx[n.areas + intercept + n.var]

  # Covariate
  #X <- mmatrix[, intercept + n.var]

  eta <- Predictor[1:n.areas]
  deta <- dnorm(eta)


  # Spat. autocorr. coefficient (NOT USED for SEM)
  rho <- rho.min + theta[1] * (rho.max - rho.min)


  # Traces of d(f(eta)) %*% W^k, k=0, 1, ..., 11
  # Note that traces are already divided by n.areas
  DETA <- Diagonal(n.areas, x = deta)
  TRACES <- rep(NA, 11 + lag)
  TRACES[1] <- mean(diag(DETA))
  TRACES[2] <- 0

  AUX <- DETA %*% W
  for(i in 3:(11 + lag)) {
    AUX <- AUX %*% W
    TRACES[i] <- mean(diag(AUX))
  }

  if(!lag) {
    # Impacts
    #eta <- X * coeff
    total.impact <- mean(deta) * coeff / (1 - rho)

    # Direct impacts
    dir.impact <- coeff * sum((rho^(0:10)) * TRACES[1:11])

  } else { # If lagged covariates
    # Value of the coefficient of the lagged covariate
    coeff_lag <- idx[n.areas + intercept + n.cov / 2  + n.var]

    #X_lag <- mmatrix[, intercept + n.cov / 2 + n.var]

    # 'linear pedictor' eta
    #eta <- X * coeff + X_lag * coeff_lag
    # Impacts
    total.impact <- mean(deta) * (coeff + coeff_lag) / (1 - rho)

    # Direct impacts
    dir.impact <- coeff * sum((rho^(0:10)) * TRACES[1:11])
    dir.impact <- dir.impact +
      coeff_lag * sum((rho^(0:10) * TRACES[1 + 1:11]))
  }

  # Indirect impacts
  indir.impact <- total.impact - dir.impact

  return(c(dir.impact, indir.impact, total.impact))
}


# Compute all impacts
# inla_samples: Samples from inla.posterior.sample
# n.var: NOW the number of covariates to compute the impacts
compute_impacts_slm_probit_all <- function(inla_samples, n.areas,
  W, n.var, intercept = TRUE, lag = FALSE, mmatrix) {

  res <- data.frame(sapply(1:n.var, function(X) {
    impacts <- inla.posterior.sample.eval(compute_impacts_slm_probit,
      inla_samples, n.areas = n.areas, W = W, n.var = X, lag = lag,
      mmatrix = mmatrix)

    impacts <- c(apply(impacts, 1, mean), apply(impacts, 1, sd))
    return(impacts)
  }))
  

  colnames(res) <- colnames(mmatrix)[intercept + 1:n.var]
  rownames(res) <- c("direct (mean)", "indirect (mean)", "total (mean)",
    "direct (s.d.)", "indirect (s.d.)", "total (s.d.)")

  return(res)
}

