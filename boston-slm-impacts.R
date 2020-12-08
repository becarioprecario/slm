#
#Compute impacts using the models fitted to the Boston dataset
#

#Load libraries
library(INLA)
library(INLABMA)
library(spdep)
library(parallel)

load("boston-slm.RData")


#Compute some indices to be used later
#Position of beta's in the random effects vector
nvars<-ncol(mmatrix)-1
idxbeta<-n+1+1:nvars
idxgamma<-n+1+nvars+1:nvars


#SEM model (no need to re-compute anything)
#Avg. Direct impacts = \beta_r
#Avg. total impacts = \beta_r
#Avg. Indir. impacts = 0 


#SDEM and SLX models 
#Avg. Direct impacts = \beta_r
#Avg. total impacts = \beta_r+\gamma_r (re-compute using lin. comb.)
#Avg. Indir. impacts = \gamma_r

lcbase<-inla.make.lincomb(a=1, b=1)

nms<-rownames(sdemm1$summary.fixed)[-1]

lcsdem<-lapply(1:nvars, function(X){
	lcaux<-lcbase
	names(lcaux$lc[[1]])<-nms[X]
	names(lcaux$lc[[2]])<-nms[X+nvars]

	lcaux
})
lcsdem<-do.call(c, lcsdem)
names(lcsdem)<-nms[1:nvars]

#Refit models
sdemm1<-inla(fsdem,
   data=boston.c2, family="gaussian",
   lincomb=lcsdem,
   control.family = list(hyper=zero.variance),
   control.compute=list(dic=TRUE, cpo=TRUE)
)

slxm1<-inla(fslx,
   data=boston.c2, family="gaussian",
   lincomb=lcsdem,
   control.family = list(hyper=zero.variance),
   control.compute=list(dic=TRUE, cpo=TRUE)
)


#SLM and SDM models; for SLM model  \gamma_r=0
#
#Avg. Direct impacts = n^{-1} tr[(I_n-\rholag W)^{-1}]\beta_r + 
#				n^{-1}tr[(I_n-\rholag W)^{-1}W]\gamma_r
#Avg. Total impacts = (1/(1-rho))* (\beta_r+\gamma_r) (re-compute using lin. comb.)
#Avg. Indir. impacts =  Avg. Total -  Avg. Direct

#Some linear combinations need to be recomputed for the SDM model
lcsdm<-lapply(1:nvars, function(X){
	auxw<-rep(0, nrow(sdemm1$summary.random$idx))
	auxw[c(idxbeta[X], idxgamma[X])]<-1
	lc<-inla.make.lincomb(idx=auxw)

})
lcsdm<-do.call(c, lcsdm)
names(lcsdm)<-nms[1:nvars]


#Refit SDM model
sdmm1<-inla( log(CMEDV) ~ -1 +
   f(idx, model="slm",
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W=W, X=mmatrix2,
         Q.beta=Q.beta2),
      hyper=hyper.slm),
   lincomb=lcsdm,
   data=as.data.frame(boston.c), family="gaussian",
   control.family = list(hyper=zero.variance),
   control.compute=list(dic=TRUE, cpo=TRUE)
)


#Compute function on rho
invrhoslm<-inla.tmarginal(function(rho){1/(1-rho)}, slmmarg)
invrhosdm<-inla.tmarginal(function(rho){1/(1-rho)}, sdmmarg)

#This takes a long time...
#  user  system elapsed 
#199.016 126.748 339.910 
trIrhoWbslm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W, rho, order=5)
}, slmmarg, n=30)

trIrhoWbsdm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W, rho, order=5)
}, sdmmarg, n=30)

trIrhoWgsdm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W, rho, order=5, offset=1)
}, sdmmarg, n=30)



#Compute impacts for ML models
impm2<-impacts(m2, listw=lw)#SLM
impm3<-impacts(m3, listw=lw)#SDM



#Create table with summar results

#TOTAL IMPACTS
tabti<-data.frame(MLSLM=impm2$total, INLASLM=NA, MLSDM=impm3$total, 
   INLASDM=NA)
rownames(tabti)<-attr(impm2, "bnames")

#Add SLM
tabti$INLASLM<-inla.zmarginal(invrhoslm, TRUE)$mean * slmm1$summary.random$idx[idxbeta,"mean"]
#Add SDM
tabti$INLASDM<-inla.zmarginal(invrhosdm, TRUE)$mean * sdmm1$summary.lincomb.derived[, "mean"]



#DIRECT IMPACTS
tabdi<-data.frame(MLSLM=impm2$direct, INLASLM=NA, MLSDM=impm3$direct, 
   INLASDM=NA)
rownames(tabdi)<-attr(impm2, "bnames")

#Add SLM
tabdi$INLASLM<-inla.zmarginal(trIrhoWbslm, TRUE)$mean * slmm1$summary.random$idx[idxbeta,"mean"]/n
tabdi$INLASDM<-(inla.zmarginal(trIrhoWbsdm, TRUE)$mean * sdmm1$summary.random$idx[idxbeta,"mean"]+inla.zmarginal(trIrhoWgsdm, TRUE)$mean * sdmm1$summary.random$idx[idxgamma,"mean"])/n


#ADD SDEM
tabti$INLASDEM<- sdemm1$summary.lincomb.derived[, "mean"]
tabdi$INLASDEM<-  sdemm1$summary.fixed[1+1:nvars, "mean"]

#ADD SLX
tabti$INLASLX<- slxm1$summary.lincomb.derived[, "mean"]
tabdi$INLASLX<-  slxm1$summary.fixed[1+1:nvars, "mean"]

#Produce tables for the paper
library(xtable)

print(xtable(tabdi, digits=5))
print(xtable(tabti, digits=5))
print(xtable(tabti-tabdi, digits=5))


#Create tables of the fixed effects so that thay are equal
#to the impacts reported for some models

#Create table to sumarise results
tabfixed<-data.frame(
   SEM=semm1$summary.fixed[,1],
   SLM=slmm1$summary.random$idx[n+1:(1+13),2],
   SDM=sdmm1$summary.random$idx[n+1:(1+13),2],
   SDMlag=c(NA, sdmm1$summary.random$idx[n+14+1:13,2]),
   SDEM=sdemm1$summary.fixed[1:(1+13),1],
   SDEMlag=c(NA, sdemm1$summary.fixed[14+1:13,1]),
   SLX=slxm1$summary.fixed[1:(1+13),1],
   SLXlag=c(NA, slxm1$summary.fixed[14+1:13,1])
)

#Create table in LaTeX format
library(xtable)
xtable(round(tabfixed, 5), digits=5)

#Transform Spatial autocorrelation parameters to be in  (rho.min, rho.max)
#
ff<-function(z){z*(rho.max-rho.min)+rho.min}
semmarg<-inla.tmarginal(ff, semm1$marginals.hyperpar[[2]])
slmmarg<-inla.tmarginal(ff, slmm1$marginals.hyperpar[[2]])
sdmmarg<-inla.tmarginal(ff, sdmm1$marginals.hyperpar[[2]])
sdemmarg<-inla.tmarginal(ff, sdemm1$marginals.hyperpar[[2]])

#Create table with summary values of spatial autocorrelations  in the
#R-INLA INTERNAL SCALE
#
tabrho<-data.frame(
   SEM=semm1$summary.hyper[2,c(1,2,3,5)],
   SLM=slmm1$summary.hyper[2,c(1,2,3,5)],
   SDM=sdmm1$summary.hyper[2,c(1,2,3,5)],
   SDEM=sdemm1$summary.hyper[2,c(1,2,3,5)]#,
#   SLX=as.numeric(NA)#slxm1$summary.hyper[2,1]
)
#xtable(round(tabrho, 4), digits=3)


#Transformed data
#Given that the transformation is linear this is a shortcut to get
#summary statistics
tabrhotrans<-tabrho
tabrhotrans[1,]<-tabrhotrans[1,]*(rho.max-rho.min)+rho.min
tabrhotrans[2,]<-tabrhotrans[2,]*(rho.max-rho.min)#+rho.min
tabrhotrans[3,]<-tabrhotrans[3,]*(rho.max-rho.min)+rho.min
tabrhotrans[4,]<-tabrhotrans[4,]*(rho.max-rho.min)+rho.min

xtable(round(tabrhotrans, 4), digits=3)

#Plot of rho: R-INLA INTERNAL SCALE marginals
pdf(file="Boston-rho-orig.pdf")

plot(semm1$marginals.hyperpar[[2]], type="l", xlab="rho", ylab="density",
   xlim=c(0.68, 0.92), ylim=c(0,30))
lines(slmm1$marginals.hyperpar[[2]], lty=2)
lines(sdmm1$marginals.hyperpar[[2]], lty=3)
lines(sdemm1$marginals.hyperpar[[2]], lty=4)

legend(.80, 30, c("SEM", "SLM", "SDM", "SDEM"), lty=1:4, bty="n")

dev.off()
system("convert Boston-rho-orig.pdf Boston-rho-orig.eps")



#Plot of rho: TRANSFORMED marginals
pdf(file="Boston-rho-trans.pdf")

plot(semmarg, type="l", xlab="rho", ylab="density",
   xlim=c(0.4, 0.83), ylim=c(0,15))
lines(slmmarg, lty=2)
lines(sdmmarg, lty=3)
lines(sdemmarg, lty=4)

legend(.60, 15, c("SEM", "SLM", "SDM", "SDEM"), lty=1:4, bty="n")

dev.off()
system("convert Boston-rho-trans.pdf Boston-rho-trans.eps")








#Use Gaussian approximation with mean and variance of the product of
#two random variables
#E[X * Y] =E[X] * E[Y]
#VAR[X * Y]=(E[X]^2) * VAR[Y]+(E[Y]^2) * VAR[X] + VAR[X]¯* VAR[Y]
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
appslmnox<-appimpacts(
   inla.zmarginal(invrhoslm,TRUE)$mean, slmm1$summary.random$idx[490+1+5,]$mean,
   inla.zmarginal(invrhoslm,TRUE)$sd, slmm1$summary.random$idx[490+1+5,]$sd)
#Approximation for SDM
appsdmnox<-appimpacts(
   inla.zmarginal(invrhosdm,TRUE)$mean, sdmm1$summary.lincomb.derived[5,]$mean,
   inla.zmarginal(invrhosdm,TRUE)$sd, sdmm1$summary.lincomb.derived[5,]$sd)


#Display total impacts for NOX-squared
pdf(file="NOXtotimpacts.pdf")
#SEM
plot(semm1$marginals.fixed$`I(NOX^2)`, type="l", xlim=c(-2, .5), ylim=c(0,3.1))
#SLM
curve(dnorm(x, mean=appslmnox$meanXY, sd=sqrt(appslmnox$varXY)), lty=2, add=TRUE)
#SDM
curve(dnorm(x, mean=appsdmnox$meanXY, sd=sqrt(appsdmnox$varXY)), lty=3, add=TRUE)
#SDEM
lines(sdemm1$marginals.lincomb.derived$var5, lty=4)
#SLX
lines(slxm1$marginals.lincomb.derived$var5, lty=5)

legend(-2, 3, c("SEM", "SLM", "SDM", "SDEM", "SLX"), lty=1:5, bty="n")
dev.off()



#Check that the approximations for the SLM and SDM models are accurate
#Here we load and display the total impacts using Roger's run of the 
#SET Matlab code and compare to the approximations for all variables

#SLM
load("Roger_files/Boston/sar_g.RData")
pdf(file="totimp-slm.pdf")
par(mfrow=c(4,4))
for(i in 1:13)
{
appslm<-appimpacts(
   inla.zmarginal(invrhoslm,TRUE)$mean, slmm1$summary.random$idx[490+1+i,]$mean,
   inla.zmarginal(invrhoslm,TRUE)$sd, slmm1$summary.random$idx[490+1+i,]$sd)
plot(density(sar_g_total[,i]), main=colnames(mmatrix)[1+i])
curve(dnorm(x, mean=appslm$meanXY, sd=sqrt(appslm$varXY)), lty=2, col="red",
   add=TRUE)
}
dev.off()


#SDM
load("Roger_files/Boston/sdm_g.RData")
pdf(file="totimp-sdm.pdf")
par(mfrow=c(4,4))
for(i in 1:13)
{
appsdm<-appimpacts(
   inla.zmarginal(invrhosdm,TRUE)$mean, sdmm1$summary.lincomb.derived[i,]$mean,
   inla.zmarginal(invrhosdm,TRUE)$sd, sdmm1$summary.lincomb.derived[i,]$sd)
plot(density(sdm_g_total[,i]), main=colnames(mmatrix)[1+i])
curve(dnorm(x, mean=appsdm$meanXY, sd=sqrt(appsdm$varXY)), lty=2, col="red",
   add=TRUE)
}
dev.off()



#Display total impacts for NOX-squared INCLUDING MCMC impacts for SLM and SDM
pdf(file="NOXtotimpacts-2.pdf")
#SEM
plot(semm1$marginals.fixed$`I(NOX^2)`, lty = 2, type="l", xlim=c(-2, .5), 
  ylim=c(0,3.1))
#SLM
curve(dnorm(x, mean=appslmnox$meanXY, sd=sqrt(appslmnox$varXY)), lty=1, add=TRUE)
lines(density(sar_g_total[,5]), lty=1, lwd=3)
#SDM
curve(dnorm(x, mean=appsdmnox$meanXY, sd=sqrt(appsdmnox$varXY)), lty=3, add=TRUE)
lines(density(sdm_g_total[,5]), lty=3, lwd=3)
#SDEM
lines(sdemm1$marginals.lincomb.derived$var5, lty=4)
#SLX
lines(slxm1$marginals.lincomb.derived$var5, lty=5)

legend(-2, 3, c("SEM", "SLM", "SDM", "SDEM", "SLX"), lty=c(2,1,3,4,5), bty="n")
legend(-2, 2, c("INLA", "MCMC"), lwd=c(1, 3),  bty="n")
dev.off()


