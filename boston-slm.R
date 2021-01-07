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

#SDM model
sdmm1 <- inla( log(CMEDV) ~ -1 + 
   f(idx, model = "slm", 
      args.slm = list(rho.min = rho.min, rho.max = rho.max, W = W, X = mmatrix2,
         Q.beta = Q.beta2), 
      hyper = hyper.slm),
   data = as(boston.c, "data.frame"), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)

#SDEM model (this requires a previous setup)
#First of all, we add the lagged covariates in the new data.frame
boston.c2 <- as.data.frame(mmatrix2[, -1])
names(boston.c2) <- paste("var", 1:ncol(boston.c2), sep = "")

boston.c2$CMEDV <- boston.c$CMEDV
boston.c2$idx <- boston.c$idx

#Construct formula  for SDEM model with lagged covariates
fsdem <- as.formula(paste("log(CMEDV) ~ ", paste(names(boston.c2)[1:26],
    collapse=" + "),
  "+f(idx,model=\"slm\", args.slm=args.slm, hyper=hyper.slm)", collapse=""))

#Fit SDEM model
sdemm1 <- inla(fsdem,
   data = boston.c2, family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)

#SLX model
fslx <- as.formula(paste("log(CMEDV) ~ ",
  paste(names(boston.c2)[1:26], collapse=" + "), "+f(idx,model=\"iid\")",
  collapse=""))

slxm1 <- inla(fslx, 
   data = boston.c2, family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)


#Create table to sumarise results
tabfixed<-data.frame(
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
xtable(round(tabfixed, 5), digits = 5)


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
#R-INLA INTERNAL SCALE
#
tabrho <- data.frame(
   SEM = semm1$summary.hyper[2, c(1, 2, 3, 5)],
   SLM = slmm1$summary.hyper[2, c(1, 2, 3, 5)],
   SDM = sdmm1$summary.hyper[2, c(1, 2, 3, 5)],
   SDEM = sdemm1$summary.hyper[2, c(1, 2, 3, 5)]#,
#   SLX=as.numeric(NA)#slxm1$summary.hyper[2,1]
)
#xtable(round(tabrho, 4), digits=3)


#Transformed data
#Given that the transformation is linear this is a shortcut to get
#summary statistics
tabrhotrans <- tabrho
tabrhotrans[1, ] <- tabrhotrans[1, ] * (rho.max - rho.min) + rho.min
tabrhotrans[2, ] <- tabrhotrans[2, ] * (rho.max - rho.min)#+rho.min
tabrhotrans[3, ] <- tabrhotrans[3, ] * (rho.max - rho.min) + rho.min
tabrhotrans[4, ] <- tabrhotrans[4, ] * (rho.max - rho.min) + rho.min

xtable(round(tabrhotrans, 4), digits = 3)



#Plot of rho: R-INLA INTERNAL SCALE marginals
pdf(file = "Boston-rho-orig.pdf")

plot(semm1$marginals.hyperpar[[2]], type = "l", xlab = "rho", ylab = "density",
   xlim = c(0.68, 0.92), ylim = c(0,30))
lines(slmm1$marginals.hyperpar[[2]], lty = 2)
lines(sdmm1$marginals.hyperpar[[2]], lty = 3)
lines(sdemm1$marginals.hyperpar[[2]], lty = 4)

legend(0.80, 30, c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
system("convert Boston-rho-orig.pdf Boston-rho-orig.eps")



#Plot of rho: TRANSFORMED marginals
pdf(file = "Boston-rho-trans.pdf")

plot(semmarg, type = "l", xlab = "rho", ylab = "density",
   xlim = c(0.4, 0.83), ylim = c(0,15))
abline(v = tabrhotrans$SEM[1], lty = 1)
lines(slmmarg, lty = 2)
abline(v = tabrhotrans$SLM[1], lty = 2)
lines(sdmmarg, lty = 3)
abline(v = tabrhotrans$SDM[1], lty = 3)
lines(sdemmarg, lty = 4)
abline(v = tabrhotrans$SDEM[1], lty = 4)

legend("topleft", legend = c("SEM", "SLM", "SDM", "SDEM"), lty = 1:4, bty = "n")

dev.off()
system("convert Boston-rho-trans.pdf Boston-rho-trans.eps")





#Check other form for NOX using a continuous random walk. Will
#not work with SLM and SDM as covariates are inside the 'slm' class

#SEM model
semm1rw2 <- inla(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + 
   f(NOX, model = "rw2", hyper = list(theta = list(param = c(1, 1))))+
   I(RM^2) +
   AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)+
   f(idx, model = "slm", args.slm = args.slm, hyper = hyper.slm),
   data = as.data.frame(boston.c), family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)
plot(semm1rw2, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
dev.off()

#SDEM model
fsdemrw2 <- update(fsdem, . ~ . - var5 - var18)
fsdemrw2 <- update(fsdemrw2, . ~ . +
  f(var5, model = "rw2", hyper = list(theta = list(param =c(1, 1)))) )

sdemm1rw2 <- inla(fsdemrw2,
   data = boston.c2, family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)
plot(sdemm1rw2, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)


#SLX model
fslxrw2 <- update(fslx, . ~ . - var5 - var18)
fslxrw2 <- update(fslxrw2, . ~ . + 
#  f(var5, model = "rw2", hyper = list(theta = list(param = c(1, 1)))) )
  f(inla.group(var5), model = "rw2",
    hyper = list(theta = list(param = c(1, 1)))) )

slxm1rw2 <- inla(fslxrw2,
   data = boston.c2, family = "gaussian",
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)
plot(slxm1rw2, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)

#Create plot from these models
pdf(file = "NOX.pdf", width = 5, height = 2.5)

par(mfrow = c(1, 3))
plot(semm1rw2$summary.random$NOX[, 1:2], type = "l", lwd = 2,
  main = "SEM model", xlab = "NOX", ylab = "f(NOX)")
lines(semm1rw2$summary.random$NOX[, c(1, 4)], lty = 2)
lines(semm1rw2$summary.random$NOX[, c(1, 6)], lty = 2)

plot(sdemm1rw2$summary.random$var5[, 1:2], type = "l", lwd = 2,
  main = "SDEM model", xlab = "NOX", ylab = "f(NOX)")
lines(sdemm1rw2$summary.random$var5[, c(1, 4)], lty = 2)
lines(sdemm1rw2$summary.random$var5[, c(1, 6)], lty = 2)

plot(slxm1rw2$summary.random$var5[, 1:2], type = "l", lwd = 2,
  main = "SLX model", xlab = "NOX", ylab = "f(NOX)")
lines(slxm1rw2$summary.random$var5[, c(1, 4)], lty = 2)
lines(slxm1rw2$summary.random$var5[, c(1, 6)], lty = 2)

dev.off()


#Add map of SEM and SDEM random effects
boston.c$SEM <- semm1$summary.random$idx$mean
boston.c$SDEM <- sdemm1$summary.random$idx$mean

ats <- seq(-0.8, 0.8, by = 0.1)
cls <- grey.colors(length(ats) + 1)
pdf(file = "slm-effects.pdf", width = 8, height = 4)
print(spplot(boston.c, c("SEM", "SDEM"), at = ats, col.regions = cls, 
   col = "transparent"))
dev.off()



#Save all objects
save(file = "boston-slm.RData", list = ls())
