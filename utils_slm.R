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
      rho.e.values <- rho * e.values * rho.e.values / i
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

