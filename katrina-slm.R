#
#RE-analysis of the Katrina dataset using R-INLA
#
library(INLA)
library(parallel)
library(spdep)

#Here I use the katrina dataset from spatialprobit

library(spatialprobit)
data(Katrina)


#And index for slm model
n <- nrow(Katrina)
Katrina$idx <- 1:n

# (a) 0-3 months time horizon
# LeSage et al. (2011) use k=11 nearest neighbors in this case
nb <- knn2nb(knearneigh(cbind(Katrina$long, Katrina$lat), k=11, longlat=TRUE))
listw <- nb2listw(nb, style="W")
W1 <- as(listw, "CsparseMatrix")


#Model matrix for SLM models
f1<- y1 ~ 1 + flood_depth + log_medinc + small_size + large_size +
  low_status_customers + high_status_customers + owntype_sole_proprietor +
  owntype_national_chain
mm <- model.matrix(f1, Katrina)
#With lagged covariates
mm2 <-  cbind(mm, as.matrix(W1) %*% mm[, -1])

#DEFINE PRIORS TO BE USED WITH R-INLA

#Compute eigenvalues for SLM model (as in Havard's code)
e <- eigen(W1)$values
re.idx <- which(abs(Im(e)) < 1e-6)
rho.max <- 1 / max(Re(e[re.idx]))
rho.min <- 1 / min(Re(e[re.idx]))
rho <- mean(c(rho.min, rho.max))



#Variance-covarinace matrix for beta coeffients' prior
#
betaprec1 <- 0.0001
#Standard regression model
Q.beta1 <- Diagonal(n = ncol(mm), x = 1)
Q.beta1 <- betaprec1 * Q.beta1
#Regression model with lagged covariates
Q.beta2 <- Diagonal(n = ncol(mm2), x = 1)
Q.beta2 <- betaprec1 * Q.beta2

#This defines the range for the spatial autocorrelation parameters
# the spatial adjacenty matrix W, the associated
#matrix of covariates X and the precision matrix Q.beta for the prior
#on the coefficients of the covariates 
#
args.slm <- list(
   rho.min = rho.min ,
   rho.max = rho.max,
   W = W1,
   X = matrix(0, nrow(mm), 0),
   Q.beta = matrix(1, 0, 0)
)

#Definition of priors for precision of the error effect in the slm latent
#effect and the spatial autocorrelation parameter (in the (0,1) interval).
#
hyper.slm <- list(
   prec = list(initial = log(1), fixed = TRUE),#prior = "loggamma", param = c(0.01, 0.01)),
      rho = list(prior = "logitbeta", param = c(1, 1))#list(initial=0, prior = "logitbeta", param = c(1,1))
)

#SEM model
semm1<-inla(
   update(f1, ~.+f(idx, model="slm", args.slm=args.slm, hyper=hyper.slm)),
   data=Katrina, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE)
)

#SLM model
slmm1<-inla( y1 ~ -1 +
   f(idx, model="slm",
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W=W1, X=mm,
         Q.beta=Q.beta1),
      hyper=hyper.slm),
   data=Katrina, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE)
)


#SDM model
sdmm1<-inla( y1 ~ -1 +
   f(idx, model="slm",
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W=W1, X=mm2,
         Q.beta=Q.beta2),
      hyper=hyper.slm),
   data=Katrina, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE)
)

#SDEM model (this requires a previous setup)
#First of all, we add the lagged covariates in the new data.frame
Katrina2 <- as.data.frame(mm2[, -1])
names(Katrina2) <- paste("var", 1:ncol(Katrina2), sep = "")

Katrina2$y1 <- Katrina$y1
Katrina2$idx <- Katrina$idx

#Construct formula  for SDEM model with lagged covariates
fsdem <- as.formula(paste("y1 ~ ", paste(names(Katrina2)[1:16], collapse=" + "),
  "+f(idx, model = \"slm\", args.slm = args.slm, hyper = hyper.slm)",
    collapse = ""))

#Fit SDEM model
sdemm1<-inla(fsdem,
   data=Katrina2, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE)
)

#SLX model
fslx <- as.formula(paste("y1 ~ ",
  paste(names(Katrina2)[1:16], collapse = " + "),
    "+f(idx, model = \"iid\")", collapse = ""))

slxm1<-inla(fslx,
   data=Katrina2, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE)
)


#Create table to sumarise results
tabfixed <- data.frame(
  SEM = semm1$summary.fixed[, 1],
  SLM = slmm1$summary.random$idx[n + 1:(1 + 8), 2],
  SDM = sdmm1$summary.random$idx[n + 1:(1 + 8), 2],
  SDMlag = c(NA, sdmm1$summary.random$idx[n + 9 + 1:8, 2]),
  SDEM = sdemm1$summary.fixed[1:(1 + 8), 1],
  SDEMlag = c(NA, sdemm1$summary.fixed[9 + 1:8, 1]),
  SLX = slxm1$summary.fixed[1:(1 + 8), 1],
  SLXlag = c(NA, slxm1$summary.fixed[9 + 1:8, 1])
)
rownames(tabfixed) <- colnames(mm)

#Table in LaTeX format
library(xtable)
xtable(round(tabfixed, 5), digits = 3)



#Transform Spatial autocorrelation parameters to be in  (rho.min, rho.max)
#
ff <- function(z) {
  z * (rho.max - rho.min) + rho.min
}
semmarg <- inla.tmarginal(ff, semm1$marginals.hyperpar[[1]])
slmmarg <- inla.tmarginal(ff, slmm1$marginals.hyperpar[[1]])
sdmmarg <- inla.tmarginal(ff, sdmm1$marginals.hyperpar[[1]])
sdemmarg <- inla.tmarginal(ff, sdemm1$marginals.hyperpar[[1]])

#Transformed data
#Given that the transformation is linear this is a shortcut to get
#summary statistics
tabrhotrans <- data.frame(cbind(
   SEM = unlist(inla.zmarginal(semmarg)[c(1, 2, 3, 7)], FALSE),
   SLM = unlist(inla.zmarginal(slmmarg)[c(1, 2, 3, 7)], FALSE),
   SDM = unlist(inla.zmarginal(sdmmarg)[c(1, 2, 3, 7)], FALSE),
   SDEM = unlist(inla.zmarginal(sdemmarg)[c(1, 2, 3, 7)], FALSE)
#   SLX = as.numeric(NA)#slxm1$summary.hyper[2,1]
))

xtable(round(tabrhotrans, 4), digits=3)


#Plot of rho: ORIGINAL marginals
pdf(file = "Katrina-rho-orig.pdf")

plot(semm1$marginals.hyperpar[[1]], type = "l",
   xlab = expression(rho), ylab = expression(paste(pi, "(", rho, "|y)")),
   xlim = c(0.5, 1), ylim = c(0, 30))
lines(slmm1$marginals.hyperpar[[1]], lty = 2)
lines(sdmm1$marginals.hyperpar[[1]], lty = 3)
lines(sdemm1$marginals.hyperpar[[1]], lty = 4)

legend(0.80, 30, c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
system("convert Katrina-rho-orig.pdf Katrina-rho-orig.eps")



#Plot of rho: TRANSFORMED marginals
pdf(file = "Katrina-rho-trans.pdf")

plot(semmarg, type = "l", 
  xlab = expression(rho), ylab = expression(paste(pi, "(", rho, "|y)")),
  xlim = c(-0.2, 0.8), ylim = c(0, 6))
abline(v = tabrhotrans$SEM[1], lty = 1)
lines(slmmarg, lty = 2)
abline(v = tabrhotrans$SLM[1], lty = 2)
lines(sdmmarg, lty = 3)
abline(v = tabrhotrans$SDM[1], lty = 3)
lines(sdemmarg, lty = 4)
abline(v = tabrhotrans$SDEM[1], lty = 4)

legend("top", legend = c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
system("convert Katrina-rho-trans.pdf Katrina-rho-trans.eps")



save(file = "katrina-slm.RData", list = ls())

