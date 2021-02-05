#
#Compute impacts using the models fitted to the Boston dataset
#

#NOTE: Impacts rely on a series approximation to (I - rho * W) ^{-1}
# using the eigenvalues of W.
#
# FIXME: Use the implementation in LeSage and Pace (2005), Sect. 4.9, page 114-...
# 
# ALSO: Rspectre may be used to obtain the range of the spat. auto. param.
#

#Load libraries
library(INLA)
library(INLABMA)
library(spatialreg)
library(parallel)

# Load some function for the slm model
source("utils_slm.R")

load("boston-slm.RData")


#Compute some indices to be used later
#Position of beta's in the random effects vector
nvars <- ncol(mmatrix) - 1
idxbeta <- n + 1 + 1:nvars
idxgamma <-  n + 1 + nvars + 1:nvars


#SEM model (no need to re-compute anything)
#Avg. Direct impacts = \beta_r
#Avg. total impacts = \beta_r
#Avg. Indir. impacts = 0 


#SDEM and SLX models 
#Avg. Direct impacts = \beta_r
#Avg. total impacts = \beta_r+\gamma_r (re-compute using lin. comb.)
#Avg. Indir. impacts = \gamma_r

lcbase <- inla.make.lincomb(a = 1, b = 1)

nms <- rownames(sdemm1$summary.fixed)[-1]

lcsdem <- lapply(1:nvars, function(X) {
  lcaux <- lcbase
  names(lcaux$lc[[1]]) <- nms[X]
  names(lcaux$lc[[2]]) <- nms[X + nvars]

  lcaux
})
lcsdem <- do.call(c, lcsdem)
names(lcsdem) <- nms[1:nvars]

#Refit models
sdemm1 <- inla(fsdem,
  data = boston.c2, family = "gaussian",
  lincomb = lcsdem,
  control.family = list(hyper = zero.variance),
  control.compute = list(dic = TRUE, cpo = TRUE)
)

slxm1 <- inla(fslx,
   data = boston.c2, family = "gaussian",
   lincomb = lcsdem,
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)


# SLM and SDM models; for SLM model  \gamma_r=0
#   These are computed using sampling

#
#Avg. Direct impacts = n^{-1} tr[(I_n-\rholag W)^{-1}]\beta_r + 
#				n^{-1}tr[(I_n-\rholag W)^{-1}W]\gamma_r
#Avg. Total impacts = (1/(1-rho))* (\beta_r+\gamma_r) (re-compute using lin. comb.)
#Avg. Indir. impacts =  Avg. Total -  Avg. Direct

# SLM model
# Sampling
samp_slmm1 <- inla.posterior.sample(2000, slmm1)
# Compute eigenvalues of W
e.values <- eigen(W)$values


# Compute impacts and summary of impacts
# Returns impacts dir, ind, tot by column: means + s.d.'s
imp_slmm1 <- sapply(1:13, function(n.var) {
  res <- inla.posterior.sample.eval(compute_impacts_slm, samp_slmm1,
    n.areas = n, e.values = e.values, n.var = n.var, intercept = TRUE)

   return(c(apply(res, 1, mean), apply(res, 1, sd)))
})

# SDM model
samp_sdmm1 <- inla.posterior.sample(2000, sdmm1)

# Impacts
imp_sdmm1 <- sapply(1:13, function(n.var) {
  res <- inla.posterior.sample.eval(compute_impacts_slm, samp_sdmm1,
    n.areas = n, e.values = e.values, n.var = n.var, intercept = TRUE,
    lag = TRUE)

   return(c(apply(res, 1, mean), apply(res, 1, sd)))
})



#Compute impacts for ML models
trs <- trW(as(lw, "CsparseMatrix"))
impm2<-impacts(m2, tr=trs)#SLM
impm3<-impacts(m3, tr=trs)#SDM


#Create table with summary results

#TOTAL IMPACTS
tabti <- data.frame(MLSLM = impm2$total, INLASLM = imp_slmm1[3, ],
  MLSDM = impm3$total, INLASDM = imp_sdmm1[3, ])
rownames(tabti) <- attr(impm2, "bnames")

#DIRECT IMPACTS
tabdi <- data.frame(MLSLM = impm2$direct, INLASLM = imp_slmm1[1, ],
  MLSDM = impm3$direct, INLASDM = imp_sdmm1[1, ])
rownames(tabdi) <- attr(impm2, "bnames")

#DIRECT IMPACTS
tabin <- data.frame(MLSLM = impm2$indirect, INLASLM = imp_slmm1[2, ],
  MLSDM = impm3$indirect, INLASDM = imp_sdmm1[2, ])
rownames(tabin) <- attr(impm2, "bnames")

#ADD SDEM
tabti$INLASDEM <- sdemm1$summary.lincomb.derived[, "mean"]
tabdi$INLASDEM <-  sdemm1$summary.fixed[1 + 1:nvars, "mean"]
tabin$INLASDEM <- sdemm1$summary.fixed[1 + nvars + 1:nvars, "mean"]

#ADD SLX
tabti$INLASLX <- slxm1$summary.lincomb.derived[, "mean"]
tabdi$INLASLX <-  slxm1$summary.fixed[1 + 1:nvars, "mean"]
tabin$INLASLX <- slxm1$summary.fixed[1 + nvars + 1:nvars, "mean"]

#Produce tables for the paper
library(xtable)

print(xtable(tabdi, digits = 3))
print(xtable(tabti, digits = 3))
print(xtable(tabin, digits = 3))


#Create tables of the fixed effects so that thay are equal
#to the impacts reported for some models

#Create table to sumarise results
tabfixed <- data.frame(
  SEM = semm1$summary.fixed[, 1],
  SLM = slmm1$summary.random$idx[n + 1:(1 + 13), 2],
  SDM = sdmm1$summary.random$idx[n + 1:(1 + 13), 2],
  SDMlag = c(NA, sdmm1$summary.random$idx[n + 14 + 1:13, 2]),
  SDEM = sdemm1$summary.fixed[1:(1 + 13), 1],
  SDEMlag = c(NA, sdemm1$summary.fixed[14 + 1:13, 1]),
  SLX = slxm1$summary.fixed[1:(1 + 13), 1],
  SLXlag = c(NA, slxm1$summary.fixed[14 + 1:13, 1])
)

#Create table in LaTeX format
library(xtable)
xtable(round(tabfixed, 5), digits = 3)

#Transform Spatial autocorrelation parameters to be in  (rho.min, rho.max)
#
ff <- function(z) {
  z * (rho.max - rho.min) + rho.min
}
semmarg <- inla.tmarginal(ff, semm1$marginals.hyperpar[[2]])
slmmarg <- inla.tmarginal(ff, slmm1$marginals.hyperpar[[2]])
sdmmarg <- inla.tmarginal(ff, sdmm1$marginals.hyperpar[[2]])
sdemmarg <- inla.tmarginal(ff, sdemm1$marginals.hyperpar[[2]])

#Create table with summary values of spatial autocorrelations  in the
#R-INLA MODEL SCALE
#
tabrhotrans <- cbind(
   SEM = unlist(inla.zmarginal(semmarg)[c(1, 2, 3, 7)], FALSE),
   SLM = unlist(inla.zmarginal(slmmarg)[c(1, 2, 3, 7)], FALSE),
   SDM = unlist(inla.zmarginal(sdmmarg)[c(1, 2, 3, 7)], FALSE),
   SDEM = unlist(inla.zmarginal(sdemmarg)[c(1, 2, 3, 7)], FALSE)
#   SLX = as.numeric(NA)#slxm1$summary.hyper[2,1]
)
#xtable(round(tabrho, 4), digits=3)


xtable(round(tabrhotrans, 4), digits = 3)

#Plot of rho: R-INLA INTERNAL SCALE marginals
pdf(file = "Boston-rho-orig.pdf")

plot(semm1$marginals.hyperpar[[2]], type = "l", xlab = "rho", ylab = "density",
   xlim = c(0.68, 0.92), ylim = c(0, 30))
lines(slmm1$marginals.hyperpar[[2]], lty = 2)
lines(sdmm1$marginals.hyperpar[[2]], lty = 3)
lines(sdemm1$marginals.hyperpar[[2]], lty = 4)

legend(.80, 30, c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
#system("convert Boston-rho-orig.pdf Boston-rho-orig.eps")



#Plot of rho: TRANSFORMED marginals
pdf(file = "Boston-rho-trans.pdf")

plot(semmarg, type = "l", xlab = "rho", ylab = "density",
   xlim = c(0.4, 0.83), ylim = c(0,15))
lines(slmmarg, lty = 2)
lines(sdmmarg, lty = 3)
lines(sdemmarg, lty = 4)

legend(.60, 15, c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
#system("convert Boston-rho-trans.pdf Boston-rho-trans.eps")


# Plot some of the impacts for checking






#Use Gaussian approximation with mean and variance of the product of
#two random variables
#E[X * Y] =E[X] * E[Y]
#VAR[X * Y]=(E[X]^2) * VAR[Y]+(E[Y]^2) * VAR[X] + VAR[X] * VAR[Y]
#          =E[X^2]*E[Y^2]-(E[X]*E[Y])^2 
#X=1/(1-\rholag))
#Y=\beta_r

appimpacts<-function(meanX, meanY, sdX, sdY)
{
	meanXY=meanX*meanY
	varXY=(meanX*sdY)^2+(meanY*sdX)^2+(sdX*sdY)^2

	return(list(meanXY=meanXY, varXY=varXY))
}

#Approximation for SLM
#appslmnox<-appimpacts(
#   inla.zmarginal(invrhoslm,TRUE)$mean, slmm1$summary.random$idx[490+1+5,]$mean,
#   inla.zmarginal(invrhoslm,TRUE)$sd, slmm1$summary.random$idx[490+1+5,]$sd)
#Approximation for SDM
#appsdmnox<-appimpacts(
#   inla.zmarginal(invrhosdm,TRUE)$mean, sdmm1$summary.lincomb.derived[5,]$mean,
#   inla.zmarginal(invrhosdm,TRUE)$sd, sdmm1$summary.lincomb.derived[5,]$sd)

# Compute the impacts for NOX-squared
n.var <- 5
NOX2_imp_slm <- inla.posterior.sample.eval(compute_impacts_slm, samp_slmm1,
  n.areas = n, e.values = e.values, n.var = n.var, intercept = TRUE)
NOX2_imp_sdm <- inla.posterior.sample.eval(compute_impacts_slm, samp_sdmm1,
  n.areas = n, e.values = e.values, n.var = n.var, intercept = TRUE,
  lag = TRUE)

#Display total impacts for NOX-squared
pdf(file="NOXtotimpacts.pdf")
#SEM
plot(semm1$marginals.fixed$`I(NOX^2)`, type = "l", xlim = c(-2, .5),
  ylim = c(0, 3.1))
#SLM
lines(density(NOX2_imp_slm[, 2]), lty = 2)
#SDM
lines(density(NOX2_imp_sdm[, 2]), lty = 3)
#SDEM
lines(sdemm1$marginals.lincomb.derived$var5, lty = 4)
#SLX
lines(slxm1$marginals.lincomb.derived$var5, lty = 5)

legend(-2, 3, c("SEM", "SLM", "SDM", "SDEM", "SLX"), lty = 1:5, bty = "n")
dev.off()



#Check that the approximations for the SLM and SDM models are accurate
#Here we load and display the total impacts using Roger's run of the 
#SET Matlab code and compare to the approximations for all variables

#SLM
#load("sar_g.RData")
sar_spBreg_lag <- spBreg_lag(f1, data=boston.c, listw=lw, Durbin=FALSE, control=list(ndraw=25000L, nomit=5000L))
trs <- trW(as(lw, "CsparseMatrix"))
sar_imps <- impacts(sar_spBreg_lag, tr=trs)

pdf(file = "totimp-slm.pdf")
par(mfrow = c(4, 4))
for(i in 1:13) {
  imp_slm <- inla.posterior.sample.eval(compute_impacts_slm, samp_slmm1,
    n.areas = n, e.values = e.values, n.var = i, intercept = TRUE)
  plot(density(sar_imps$sres$total[, i]), main = colnames(mmatrix)[1 + i])
  lines(density(imp_slm[3, ]), lty = 2, col = "red")
}
dev.off()


#SDM
#load("sdm_g.RData")
sdm_spBreg_lag <- spBreg_lag(f1, data=boston.c, listw=lw, Durbin=TRUE, control=list(ndraw=25000L, nomit=5000L))
sdm_imps <- impacts(sdm_spBreg_lag, tr=trs)

pdf(file = "totimp-sdm.pdf")
par(mfrow = c(4, 4))
for(i in 1:13) {
  imp_sdm <- inla.posterior.sample.eval(compute_impacts_slm, samp_sdmm1,
    n.areas = n, e.values = e.values, n.var = i, intercept = TRUE,
    lag = TRUE)
  plot(density(sdm_imps$sres$total[,i]), main = colnames(mmatrix)[1 + i])
  lines(density(imp_sdm[3, ]), lty = 2, col = "red")
}
dev.off()



#Display total impacts for NOX-squared INCLUDING MCMC impacts for SLM and SDM
pdf(file="NOXtotimpacts-2.pdf")
#SEM
plot(semm1$marginals.fixed$`I(NOX^2)`, lty = 2, type="l", xlim = c(-2, .5), 
  ylim = c(0, 3.1))
#SLM
lines(density(NOX2_imp_slm[3, ]), lty = 2)
#SDM
lines(density(NOX2_imp_sdm[3, ]), lty = 3)
#SDEM
lines(sdemm1$marginals.lincomb.derived$var5, lty = 4)
#SLX
lines(slxm1$marginals.lincomb.derived$var5, lty = 5)
#MCMC
lines(density(sar_imps$sres$total[, 5]), lty = 1, lwd = 3)
lines(density(sdm_imps$sres$total[, 5]), lty = 3, lwd = 3)

legend(-2, 3, c("SEM", "SLM", "SDM", "SDEM", "SLX"), lty = c(2, 1, 3, 4, 5),
  bty = "n")
legend(-2, 2, c("INLA", "MCMC"), lwd = c(1, 3), bty = "n")
dev.off()


