#
# Compute impacts using the models fitted to the Katrina dataset
# Impacts are computed using sampling
#

#Load libraries
library(INLA)
library(INLABMA)
library(spdep)
library(parallel)
options(mc.cores = 3)


source("utils_slm.R")

load("katrina-slm.RData")


# For testing, the impacts are computed for the first covariate (flood_depth)

# Obtain samples
samp_semm1 <- inla.posterior.sample(2000, semm1)
samp_sdemm1 <- inla.posterior.sample(2000, sdemm1)
samp_slxm1 <- inla.posterior.sample(2000, slxm1)


# Impacts SEM
imp_semm1 <- inla.posterior.sample.eval(compute_impacts_sem_probit, samp_semm1,
  n.areas = nrow(Katrina), n.var = 1, mmatrix = mm)
# Impacts SDEM
colnames(mm2) <- rownames(sdemm1$summary.fixed)
imp_sdemm1 <- inla.posterior.sample.eval(compute_impacts_sem_probit,
  samp_sdemm1, n.areas = nrow(Katrina), n.var = 1,
  lag = TRUE, mmatrix = mm2)
#Impacts SLX
colnames(mm2) <- rownames(slxm1$summary.fixed)
imp_slxm1 <- inla.posterior.sample.eval(compute_impacts_sem_probit,
  samp_slxm1, n.areas = nrow(Katrina), n.var = 1,
  lag = TRUE, mmatrix = mm2)

# Check
apply(imp_semm1, 1, mean)
apply(imp_sdemm1, 1, mean)
apply(imp_slxm1, 1, mean)

# Impacts SEM
imp_semm1_all <- compute_impacts_sem_probit_all(samp_semm1,
  n.areas = nrow(Katrina), n.var = 8, mmatrix = mm)
# Impacts SDEM
colnames(mm2) <- rownames(sdemm1$summary.fixed)
imp_sdemm1_all <- compute_impacts_sem_probit_all(samp_sdemm1,
  n.areas = nrow(Katrina), n.var = 8, mmatrix = mm2, lag = TRUE)
# Impacts SLX
colnames(mm2) <- rownames(slxm1$summary.fixed)
imp_slxm1_all <- compute_impacts_sem_probit_all(samp_slxm1,
  n.areas = nrow(Katrina), n.var = 8, mmatrix = mm2, lag = TRUE)

# Impacts for SLM and SDM spatial probit models

# Sampling
samp_slmm1 <- inla.posterior.sample(2000, slmm1)
samp_sdmm1 <- inla.posterior.sample(2000, sdmm1)


# Impacts SLM
imp_slmm1 <- inla.posterior.sample.eval(compute_impacts_slm_probit, samp_slmm1,
  W = W1, n.areas = nrow(Katrina), n.var = 1, mmatrix = mm)
apply(imp_slmm1, 1, mean)

# Impacts SDM
colnames(mm2) <- c(colnames(mm), paste0("lag.", colnames(mm)[-1]))
imp_sdmm1 <- inla.posterior.sample.eval(compute_impacts_slm_probit,
  samp_sdmm1, W = W1, n.areas = nrow(Katrina), e.values = NULL, n.var = 1,
  lag = TRUE, mmatrix = mm2)

# Check
apply(imp_slmm1, 1, mean)
apply(imp_sdmm1, 1, mean)



# Impacts SLM
imp_slmm1_all <- compute_impacts_slm_probit_all(samp_slmm1,
  n.areas = nrow(Katrina), W = W1, n.var = 8, mmatrix = mm)
# Impacts SDM
#colnames(mm2) <- rownames(sdmm1$summary.fixed)
colnames(mm2) <- c(colnames(mm), paste0("lag.", colnames(mm)[-1]))
imp_sdmm1_all <- compute_impacts_slm_probit_all(samp_sdmm1,
  n.areas = nrow(Katrina), W = W1, n.var = 8, mmatrix = mm2, lag = TRUE)

#Create table with summar results

#TOTAL IMPACTS
tabti <- data.frame(cbind(
  INLASEM = unlist(imp_semm1_all[3, ]),
  #MLSLM=rep(NA, nvars), 
  INLASLM = unlist(imp_slmm1_all[3, ]),
  #MLSDM=NA,
  INLASDM = unlist(imp_sdmm1_all[3, ]),
  INLASDEM = unlist(imp_sdemm1_all[3, ]),
  INLASLX = unlist(imp_slxm1_all[3, ])
))
rownames(tabti) <- colnames(imp_semm1_all)


#DIRECT IMPACTS
tabdi <- data.frame(cbind(
  INLASEM = unlist(imp_semm1_all[1, ]),
  #MLSLM=rep(NA, nvars), 
  INLASLM = unlist(imp_slmm1_all[1, ]),
  #MLSDM=NA,
  INLASDM = unlist(imp_sdmm1_all[1, ]),
  INLASDEM = unlist(imp_sdemm1_all[1, ]),
  INLASLX = unlist(imp_slxm1_all[1, ])
))
rownames(tabdi) <- colnames(imp_semm1_all)


#INDIRECT IMPACTS
tabii <- data.frame(cbind(
  INLASEM = unlist(imp_semm1_all[2, ]),
  #MLSLM=rep(NA, nvars), 
  INLASLM = unlist(imp_slmm1_all[2, ]),
  #MLSDM=NA,
  INLASDM = unlist(imp_sdmm1_all[2, ]),
  INLASDEM = unlist(imp_sdemm1_all[2, ]),
  INLASLX = unlist(imp_slxm1_all[2, ])
))
rownames(tabii) <- colnames(imp_semm1_all)



#Produce tables for the paper
library(xtable)

print(xtable(tabdi, digits = 3))
print(xtable(tabti, digits = 3))
print(xtable(tabii, digits = 3))

#Compute impacts using spatialprobit (used for comparison purposes)
library(spatialprobit)

#Add lagged covariates to Katrina for SDM model
#mmlag<- create_WX(mm, listw)
mmlag <- data.frame(mm2[, 10:17])
#names(mmlag) <- paste("lag.", names(mmlag), sep = "")

Katrina <- cbind(Katrina, mmlag)



if(!file.exists("katrina-spatialprobit.RData")){
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




  save(file = "katrina-spatialprobit.RData", 
  list = c("fit1", "fit1lag", "fit2", "fit2lag"))
} else {
  load("katrina-spatialprobit.RData")
}


#Compute direct impacts (flood_depth) for several models
idx.cov <- 2 #Flood depth
n.cov <- 8 #Number of covariates

#SEM model using output from spatialprobit
# S_r(W) = D(f(eta)) * I_n*beta_r
eta <- fit1$X %*% t(fit1$bdraw) #Fitted values 673 x 25000
Dfeta<- dnorm(eta)

dirimp.sem <- apply(Dfeta, 2, sum) * fit1$bdraw[, idx.cov] / n
totimp.sem <- dirimp.sem
#TOtal impacts = direct impacts
#Indirect impacts = 0

#Compute direct impacts (flood_depth) for 
#SDEM model using output from spatialprobit
# S_r(W) = D(f(eta)) * I_n*(beta_r + W gamma_r)
eta.sdem <- fit1lag$X %*% t(fit1lag$bdraw) #Fitted values 673 x 25000
Dfeta.sdem<- dnorm(eta.sdem)
dirimp.sdem <- apply(Dfeta.sdem, 2, sum) * fit1lag$bdraw[, idx.cov] / n

indirimp.sdem <- apply(Dfeta.sdem, 2, function(X) {
   sum(Diagonal(n, X) %*% fit1lag$W)
}) * fit1lag$bdraw[, idx.cov + n.cov] /  n

totimp.sdem <- dirimp.sdem + indirimp.sdem

#Compute impacts for SDM
# eta = (I - rho * W)^{-1} * X* beta
# S_r(W) = D(f(eta)) * (I_n - rho * W)^{-1} * (beta_r + W gamma_r)
#Fitted values 673 x 25000
eta.sdm <- sapply(1:nrow(fit2lag$bdraw), function (X) {
  solve(Diagonal(n, 1) - fit2lag$pdraw[X] * fit2lag$W,
    fit2lag$X %*% matrix(fit2lag$bdraw[X, ], ncol = 1))[, 1]
})

Dfeta.sdm<- dnorm(eta.sdm)
totimp.sdm <- apply(Dfeta.sdm, 2, sum) * (1 / (1 - fit2lag$pdraw)) *
  (fit2lag$bdraw[, idx.cov] + fit2lag$bdraw[, idx.cov + n.cov])/ n

# Direct impacts
#evalues <- eigen(fit2lag$W)$values

#dirimp.sdm <- sapply(1:nrow(fit2lag$bdraw), function(X) {
dirimp.sdm <- unlist(mclapply(sample(1:nrow(fit2lag$bdraw), 5000), function(X) {

  DETA <- Diagonal(n, x = Dfeta.sdm[, X])

  coeff <- fit2lag$bdraw[X, idx.cov]
  coeff_lag <- fit2lag$bdraw[X, idx.cov + n.cov]

  rho <- fit2lag$pdraw[X]

  TRACES <- rep(NA, 11 + 1)
  TRACES[1] <- mean(diag(DETA))
  TRACES[2] <- 0

  AUX <- DETA %*% fit2lag$W
  for(i in 3:(11 + 1)) {
    AUX <- AUX %*% fit2lag$W
    TRACES[i] <- mean(diag(AUX))
  }

  dir.impact <- coeff * sum((rho^(0:10)) * TRACES[1:11])
  dir.impact <- dir.impact +
      coeff_lag * sum((rho^(0:10) * TRACES[1 + 1:11]))

  return(dir.impact)
}))

indimp.sdm <- totimp.sdm - dirimp.sdm 





#And compare to those obtained with MCMC
#RESULTS: They look very close


# MCMC
summary((fit2$direct[, 1]))#Flood depth
summary((fit2$indirect[, 1]))#Flood depth
summary((fit2$total[, 1]))#Flood depth

# INLA
c(tabdi$INLASLM[1], tabii$INLASLM[1], tabti$INLASLM[1])


# Display posterior densities of impacts

pdf("katrina-impacts-flood-depth.pdf")
par(mfrow = c(1, 3))
plot(density(imp_slmm1[1,]), type = "l", main = "direct (f.d.)")
lines(density(fit2$direct[, 1]), col = "red")
plot(density(imp_slmm1[2,]), type = "l", main = "indirect (f.d.)")
lines(density(fit2$indirect[, 1]), col = "red")
plot(density(imp_slmm1[3,]), type = "l", main = "total (f.d.)")
lines(density(fit2$total[, 1]), col = "red")
dev.off()


# Compara all impacts
# INLA impacts are approximated using a Gaussian distribution
#  but the actual samples could be used

# DIRECT impacts
pdf(file = "katrina-directimpacts-comparison.pdf")
par(mfrow = c(3, 3))

for(i in 1:8)
{
  plot(density(fit2$direct[, i]), col = "red", main = rownames(tabti)[i])
  abline(v = tabdi$INLASLM[i])
  curve(dnorm(x, tabdi$INLASLM[i], imp_slmm1_all[4, i]), add = TRUE)
}
dev.off()


# INDIRECT impacts
pdf(file = "katrina-indirectimpacts-comparison.pdf")
par(mfrow = c(3, 3))

for(i in 1:8)
{
  plot(density(fit2$indirect[, i]), col = "red", main = rownames(tabti)[i])
  abline(v = tabii$INLASLM[i])
  curve(dnorm(x, tabii$INLASLM[i], imp_slmm1_all[5, i]), add = TRUE)
}
dev.off()


# TOTAL impacts
pdf(file = "katrina-totalimpacts-comparison.pdf")
par(mfrow = c(3, 3))

for(i in 1:8)
{
  plot(density(fit2$total[, i]), col = "red", main = rownames(tabti)[i])
  abline(v = tabti$INLASLM[i])
  curve(dnorm(x, tabti$INLASLM[i], imp_slmm1_all[6, i]), add = TRUE)
}
dev.off()


#TOTAL Impacts per model

pdf(file = "flooddepth-impacts.pdf")
par(mfrow = c(2, 2))

#SEM
plot(density(imp_semm1[3, ]), type = "l", main = "SEM: total (f.d.)",
  xlim = c(-0.16, 0),
  ylim = c(0, 32), xlab = "total impact")
lines(density(totimp.sem), lty = 1, lwd = 3)
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1, 3), bty = "n")

#SDEM
plot(density(imp_sdemm1[3,]), type = "l", main = "SDEM: total (f.d.)",
  xlim = c(-0.16, 0),
  ylim = c(0, 32), xlab = "total impact")
lines(density(totimp.sdem), lty = 1, lwd = 3)
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1, 3), bty = "n")

#SLM
plot(density(imp_slmm1[3,]), type = "l", main = "SLM: total (f.d.)",
  xlim = c(-0.16, 0),
  ylim = c(0, 32), xlab = "total impact")
lines(density(fit2$total[, idx.cov - 1]), lty = 1, lwd = 3)#SLM
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1, 3), bty = "n")

#SDM
plot(density(imp_sdmm1[3,]), type = "l", main = "SDM: total (f.d.)",
  xlim = c(-0.16, 0),
  ylim = c(0, 32), xlab = "total impact")
lines(density(totimp.sdm), lty = 1, lwd = 3)#SDM
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1, 3), bty = "n")

dev.off()


