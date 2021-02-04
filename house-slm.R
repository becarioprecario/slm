#
#Re-analysis of the Lucas County house data set using the new 'slm' latent class in R-INLA
#

#Load libraries
library(INLA)
library(spatialreg)
library(parallel)


#Load data
library(spData)
data(house)
lw <- spdep::nb2listw(LO_nb, style="W")


#Define some indices used in the models
n<-nrow(house)
house$idx<-1:n

#Adjcency matrix
#W<-nb2mat(boston.soi)
W <- as(lw, "CsparseMatrix")

#Model matrix for SLM models
f1<-log(price) ~ rooms + lotsize + yrbuilt + stories + syear + garage
mmatrix <- model.matrix(f1, house)


#Compute some Spatial Econometrics models using ML estimation
(ml_time <- system.time(m2<-lagsarlm(f1, house, lw, method="Matrix")))

(ml_imps_time <- system.time(ml_imps <- impacts(m2, tr=trMat, R=20000)))

(mcmc_time <- system.time(m2B <- spBreg_lag(f1, house, lw, control=list(ndraw=25000L, nomit=5000L))))

(mcmc_imps_time <- system.time(impacts(m2B, tr=trMat)))

#DEFINE PRIORS TO BE USED WITH R-INLA

(inla_time <- system.time({
#Zero-variance for Gaussian error term
zero.variance = list(prec=list(initial = 25, fixed=TRUE))

#Compute eigenvalues for SLM model (as in Havard's code)
#e = eigenw(lw)
#re.idx = which(abs(Im(e)) < 1e-6)
#rho.max = 1/max(Re(e[re.idx]))
#rho.min = 1/min(Re(e[re.idx]))
B <- as(similar.listw(lw), "CsparseMatrix")
(rho.max <- RSpectra::eigs(B, k=1, which="LR")$values)
(rho.min <- RSpectra::eigs(B, k=1, which="SR")$values)
rho = mean(c(rho.min, rho.max))



#
#Variance-covarinace matrix for beta coeffients' prior
#
betaprec<-.0001
#Standard regression model
Q.beta = Diagonal(n=ncol(mmatrix), x=1)
Q.beta = betaprec*Q.beta 
#Regression model with lagged covariates
#Q.beta2 = Diagonal(n=ncol(mmatrix2), x=1)
#Q.beta2 = betaprec*Q.beta2 


#This defines the range for the spatial autocorrelation parameters
# the spatial adjacenty matrix W, the associated
#matrix of covariates X and the precision matrix Q.beta for the prior
#on the coefficients of the covariates 
#
args.slm = list(
   rho.min = rho.min ,
   rho.max = rho.max,
   W = W,
   X = matrix(0, nrow(mmatrix),0),
   Q.beta = matrix(1,0,0)
)

#Definition of priors for precision of the error effect in the slm latent
#effect and the spatial autocorrelation parameter (in the (0,1) interval).
#
hyper.slm = list(
   prec = list(
      prior = "loggamma", param = c(0.01, 0.01)),
      rho = list(initial=0, prior = "logitbeta", param = c(1,1))
)



#SLM model
slmm1<-inla( log(price) ~ -1 + 
   f(idx, model="slm", 
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W=W, X=mmatrix, 
         Q.beta=Q.beta), 
      hyper=hyper.slm),
   data=as(house, "data.frame"), family="gaussian",
   control.family = list(hyper=zero.variance),
   control.compute=list(dic = TRUE, cpo = TRUE, config = TRUE)
)
}))

(inla_imps_time <- system.time({
  nvars<-ncol(mmatrix)-1
  idxbeta<-n+1+1:nvars
  ff<-function(z){z*(rho.max-rho.min)+rho.min}
  slmmarg<-inla.tmarginal(ff, slmm1$marginals.hyperpar[[2]])
  trIrhoWbslm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W, rho, order=5)
  }, slmmarg, n=30)
  res <- inla.zmarginal(trIrhoWbslm, TRUE)$mean * 
  slmm1$summary.random$idx[idxbeta,"mean"]/n
}))

#Save all objects
save(file="house-slm.RData", list=ls())
