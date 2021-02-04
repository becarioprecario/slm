#
#Compute impacts using the models fitted to the Katrina dataset
#

#Load libraries
library(INLA)
library(INLABMA)
library(spdep)
library(parallel)

load("katrina-slm.RData")


#Compute some indices to be used later
#Position of beta's in the random effects vector
nvars<-ncol(mm)-1
idxbeta<-n+1+1:nvars
idxgamma<-n+1+nvars+1:nvars


#Compute approximation to D(f) for several models
#
#In our approximation, mean(D(f)) is constant and it is used
#as weights in the linear combinations

Dfsem<- mean(dnorm(semm1$summary.linear.predictor$mean))
Dfslm<- mean(dnorm(slmm1$summary.linear.predictor$mean))
Dfslmvec<- Diagonal(x=dnorm(slmm1$summary.linear.predictor$mean))
Dfsdm<- mean(dnorm(sdmm1$summary.linear.predictor$mean))
Dfsdmvec<- Diagonal(x=dnorm(sdmm1$summary.linear.predictor$mean))
Dfsdem<- mean(dnorm(sdemm1$summary.linear.predictor$mean))
Dfslx<- mean(dnorm(slxm1$summary.linear.predictor$mean))

#Alternative approximation using the mode
Dfsem2<- mean(dnorm(semm1$summary.linear.predictor$mode))
Dfslm2<- mean(dnorm(slmm1$summary.linear.predictor$mode))
Dfsdm2<- mean(dnorm(sdmm1$summary.linear.predictor$mode))
Dfsdem2<- mean(dnorm(sdemm1$summary.linear.predictor$mode))
Dfslx2<- mean(dnorm(slxm1$summary.linear.predictor$mode))

#Compare both approximations
plot(c(Dfsem, Dfslm, Dfsdm, Dfsdem, Dfslx),
   c(Dfsem2, Dfslm2, Dfsdm2, Dfsdem2, Dfslx2),
   xlab="f(Mean)", ylab="f(Mode)")
abline(0,1)


#Alternative approximation using a point estimate of f(\eta_i)
#NOTE: dnorm() is not a one-to-one function, so inla.tmarginal() cannot be used
#In fact, I am not sure how  the marginal of f(\eta_i) could be derived.
#
#Dfsem3<-mclapply(semm1$marginals.linear.predictor[1:10], function(X){inla.tmarginal(dnorm, X)})
#Dfsem3<-unlist(mclapply(Dfsem3, function(X){inla.zmarginal(X, TRUE)$mean}))



#SEM model (no need to re-compute anything)
#Avg. Direct impacts = mean(D(f))*\beta_r
#Avg. total impacts = mean(D(f))*\beta_r
#Avg. Indir. impacts = 0 


#SDEM and SLX models 
#Avg. Direct impacts = mean(D(f))*\beta_r
#Avg. total impacts = mean(D(f))*(\beta_r+\gamma_r) (re-compute using lin. comb.)
#Avg. Indir. impacts = mean(D(f))*\gamma_r


nms<-colnames(mm)[-1]

#Function to define the linear combinations used to compute the impacts
#Note that the term on D(f) is included as weights
#
impactlclag<-function(Df, nms, nvars)
{
	lcbase<-inla.make.lincomb(a=Df, b=Df)

	lcs<-lapply(1:nvars, function(X){
		lcaux<-lcbase
		names(lcaux$lc[[1]])<-nms[X]
		names(lcaux$lc[[2]])<-nms[X+nvars]

		lcaux
	})
lcs<-do.call(c, lcs)
names(lcs)<-nms[1:nvars]

	return(lcs)

}

#Refit models
#Fit SDEM model
lcsdem<-impactlclag(Dfsdem, rownames(sdemm1$summary.fixed)[-1], nvars)

sdemm1<-inla(fsdem,
   data=Katrina2, family="binomial",
   lincomb=lcsdem,
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE)
)


lcslx<-impactlclag(Dfslx, rownames(slxm1$summary.fixed)[-1], nvars)
slxm1<-inla(fslx,
   data=Katrina2, family="binomial",
   lincomb=lcslx,
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE)
)



#SLM and SDM models; for SLM model  \gamma_r=0
#
#Avg. Direct impacts = n^{-1} tr[D(f) (I_n-\rholag W)^{-1}]\beta_r + 
#				n^{-1}tr[D(f) (I_n-\rholag W)^{-1}W]\gamma_r
#Avg. Total impacts = mean(D(f))*[(1/(1-rho))* (\beta_r+\gamma_r)] (re-compute using lin. comb.)
#Avg. Indir. impacts =  Avg. Total -  Avg. Direct

#Some linear combinations need to be recomputed for the SDM model
lcsdm<-lapply(1:nvars, function(X){
	auxw<-rep(0, nrow(sdmm1$summary.random$idx))
	auxw[c(idxbeta[X], idxgamma[X])]<-Dfsdm
	lc<-inla.make.lincomb(idx=auxw)

})
lcsdm<-do.call(c, lcsdm)
names(lcsdm)<-nms[1:nvars]


#Refit SDM model
sdmm1<-inla( y1 ~ -1 +
   f(idx, model="slm",
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W=W1, X=mm2,
         Q.beta=Q.beta2),
      hyper=hyper.slm),
   lincomb=lcsdm,
   data=Katrina, family="binomial",
   control.family = list(link="probit"),
   control.compute=list(dic=TRUE, cpo=TRUE)
)



#Compute function on rho
invrhoslm<-inla.tmarginal(function(rho){1/(1-rho)}, slmmarg)
invrhosdm<-inla.tmarginal(function(rho){1/(1-rho)}, sdmmarg)

#This takes a long time...
#    user   system  elapsed 
# 869.168  244.716 1149.683 
trIrhoWbslm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W1, rho, order=10, Df=Dfslmvec)
}, slmmarg, n=30)

trIrhoWbsdm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W1, rho, order=10, Df=Dfslmvec)
}, sdmmarg, n=30)

trIrhoWgsdm<-inla.tmarginal(function(rho){
	INLABMA:::trIrhoWinv(W1, rho, order=10, offset=1, Df=Dfslmvec)
}, sdmmarg, n=30)



#Compute impacts for ML models
#
#No ML estimates in this case...






#Create table with summar results

#TOTAL IMPACTS
tabti<-data.frame(INLASEM=Dfsem*semm1$summary.fixed[,1][-1], 
   #MLSLM=rep(NA, nvars), 
   INLASLM=NA, 
   #MLSDM=NA,
   INLASDM=NA)
rownames(tabti)<-nms

#Add SLM
tabti$INLASLM<-Dfslm*inla.zmarginal(invrhoslm, TRUE)$mean * slmm1$summary.random$idx[idxbeta,"mean"]
#Add SDM
tabti$INLASDM<-inla.zmarginal(invrhosdm, TRUE)$mean * sdmm1$summary.lincomb.derived[, "mean"]



#DIRECT IMPACTS
tabdi<-data.frame(INLASEM=tabti$INLASEM, 
   #MLSLM=rep(NA, nvars), 
   INLASLM=NA, 
   #MLSDM=NA,
   INLASDM=NA)
rownames(tabdi)<-nms

#Add SEM
tabdi$INLASEM<-Dfsem*semm1$summary.fixed[,1][-1]

#Add SLM
tabdi$INLASLM<-inla.zmarginal(trIrhoWbslm, TRUE)$mean * slmm1$summary.random$idx[idxbeta,"mean"]/n
tabdi$INLASDM<-(inla.zmarginal(trIrhoWbsdm, TRUE)$mean * sdmm1$summary.random$idx[idxbeta,"mean"]+inla.zmarginal(trIrhoWgsdm, TRUE)$mean * sdmm1$summary.random$idx[idxgamma,"mean"])/n


#ADD SDEM
tabti$INLASDEM<- sdemm1$summary.lincomb.derived[, "mean"]
tabdi$INLASDEM<-  Dfsdem*sdemm1$summary.fixed[1+1:nvars, "mean"]

#ADD SLX
tabti$INLASLX<- slxm1$summary.lincomb.derived[, "mean"]
tabdi$INLASLX<-  slxm1$summary.fixed[1+1:nvars, "mean"]

#Produce tables for the paper
library(xtable)

print(xtable(tabdi, digits=5))
print(xtable(tabti, digits=5))
print(xtable(tabti-tabdi, digits=5))


#Compute impacts using spatialprobit (used for comparison purposes)
library(spatialprobit)

#Add lagged covariates to Katrina for SDM model
#mmlag<- create_WX(mm, listw)
mmlag <- data.frame(mm2[,10:17])
names(mmlag) <- paste("lag.", names(mmlag), sep ="")

Katrina <- cbind(Katrina, mmlag)



if(!file.exists("katrina-spatialprobit.RData")){
fit1 <- semprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers + 
   owntype_sole_proprietor + owntype_national_chain, 
   W=W1, data=Katrina, ndraw=25000, burn.in = 5000, showProgress=TRUE)
summary(fit1)

#impacts(fit1)

fit1lag <- semprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers + 
   owntype_sole_proprietor + owntype_national_chain +
   lag.flood_depth + lag.log_medinc + lag.small_size + lag.large_size +
   lag.low_status_customers +  lag.high_status_customers +
   lag.owntype_sole_proprietor + lag.owntype_national_chain, 
   W=W1, data=Katrina, ndraw=25000, burn.in = 5000, showProgress=TRUE)
summary(fit1lag)

#impacts(fit1lag)



fit2 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers + 
   owntype_sole_proprietor + owntype_national_chain, 
   W=W1, data=Katrina, ndraw=25000, burn.in = 5000, showProgress=TRUE)
summary(fit2)

impacts(fit2)



fit2lag <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
   low_status_customers +  high_status_customers +
   owntype_sole_proprietor + owntype_national_chain +
   lag.flood_depth + lag.log_medinc + lag.small_size + lag.large_size +
   lag.low_status_customers +  lag.high_status_customers +
   lag.owntype_sole_proprietor + lag.owntype_national_chain,
   W=W1, data=Katrina, ndraw=25000, burn.in = 5000, showProgress=TRUE)
summary(fit2lag)

impacts(fit2lag)




save(file="katrina-spatialprobit.RData", 
  list=c("fit1", "fit1lag", "fit2", "fit2lag"))
}else{
load("katrina-spatialprobit.RData")
}


#Compute direct impacts (flood_depth) for 
#SEM model using output from spatialprobit
# S_r(W) = D(f(eta)) * I_n*beta_r
eta <- fit1$X%*%t(fit1$bdraw) #Fitted values 673 x 25000
Dfeta<- dnorm(eta)

dirimp.sem <- apply(Dfeta, 2, sum) * fit1$bdraw[,2]/n
totimp.sem <- dirimp.sem
#TOtal impacts = direct impacts
#Indirect impacts = 0

#Compute direct impacts (flood_depth) for 
#SDEM model using output from spatialprobit
# S_r(W) = D(f(eta)) * I_n*(beta_r + W gamma_r)
eta.sdem <- fit1lag$X%*%t(fit1lag$bdraw) #Fitted values 673 x 25000
Dfeta.sdem<- dnorm(eta.sdem)
dirimp.sdem <- apply(Dfeta.sdem, 2, sum) * fit1lag$bdraw[,2]/n

indirimp.sdem <- apply(Dfeta.sdem, 2, function(X){
   sum(Diagonal(n, X) %*% fit1lag$W)
}) * fit1lag$bdraw[,2+8]/n

totimp.sdem <- dirimp.sdem + indirimp.sdem



#Compute marginal of total impacts for SLM
#And compare to those obtained with MCMC
#RESULTS: They look very close
muX<-inla.zmarginal(invrhoslm, TRUE)$mean
sdX<-inla.zmarginal(invrhoslm, TRUE)$sd

#Load data provided by Roger using Matlab
load("A174p4lesage/sarp.RData")
summary((sarp$results$total[,1]))#Flood depth

pdf(file="katrina-totalimpacts-comparison.pdf")
par(mfrow=c(3,3))

for(i in 1:8)
{
	muY<-slmm1$summary.random$idx[idxbeta[i],"mean"]
	sdY<-slmm1$summary.random$idx[idxbeta[i],"sd"]

	muDfXY<-Dfslm*muX*muY
	sdDfXY<-Dfslm*sqrt( (muX*sdY)^2+(muY*sdX)^2+(sdX*sdY)^2  )
	#WARNING!!!! Note the 'n *' in the next line
	plot(density(n*fit2$total[,i]), main=nms[i], col="red")
	lines(density(sarp$results$total[,i]), main=nms[i], col="blue")
	curve(dnorm(x, mean=muDfXY, sd=sdDfXY), add=TRUE, col="black")
	#abline(v=muDfXY, col="red")

	legend(x="topleft", legend=c("Gaussian", "spprobit", "SET"), lty=c(1,1,1), 
	   col=c("black", "red", "blue"), bty="n", cex=.85 )
	
}
dev.off()


#Create a plot with the total impacts for fllod_depth

fdsem<-inla.tmarginal(function(x){Dfsem*x}, semm1$marginals.fixed$flood_depth)
fdsdem<- sdemm1$marginals.lincomb.derived$var1
fdslx<- slxm1$marginals.lincomb.derived$var1

#SLM
muX<-inla.zmarginal(invrhoslm, TRUE)$mean
sdX<-inla.zmarginal(invrhoslm, TRUE)$sd
muY<-slmm1$summary.random$idx[idxbeta[1],"mean"]
sdY<-slmm1$summary.random$idx[idxbeta[1],"sd"]
muDfXY<-Dfslm*muX*muY
sdDfXY<-Dfslm*sqrt( (muX*sdY)^2+(muY*sdX)^2+(sdX*sdY)^2  )

#SDM
muX2<-inla.zmarginal(invrhosdm, TRUE)$mean
sdX2<-inla.zmarginal(invrhosdm, TRUE)$sd
muY2<-sdmm1$summary.lincomb.derived[1, "mean"]
sdY2<-sdmm1$summary.lincomb.derived[1, "sd"]
muDfXY2<-muX2*muY2
sdDfXY2<-sqrt( (muX2*sdY2)^2+(muY2*sdX2)^2+(sdX2*sdY2)^2  )




pdf(file="flooddepth-impacts.pdf")
plot(fdsem, xlim=c(-.20,0), ylim=c(0, 30), type="l")
curve(dnorm(x, mean=muDfXY, sd=sdDfXY), type="l", lty=2, add=TRUE)
curve(dnorm(x, mean=muDfXY2, sd=sdDfXY2), type="l", lty=3, add=TRUE)
lines(fdsdem, lty=4)
lines(fdslx, lty=5)

#MCMC: spatialprobit
#lines(density(sarp$results$total[,1]),  lty=2, lwd=3)#Matlab SET model
lines(density(totimp.sem), lty = 1, lwd = 3)#SEM
lines(density(totimp.sdem), lty = 4, lwd = 3)#SDEM
lines(density(fit2$total[,1]), lty=2, lwd = 3)#SLM
lines(density((fit2lag$total[,1]+fit2lag$total[,1+8])), lty=3, lwd = 3)#SDM

legend(-.2, 30, c("SEM", "SLM", "SDM", "SDEM", "SLX"), lty=1:5, bty="n")
legend(-.2, 20, c("INLA", "MCMC"), lwd=c(1,3), bty="n")
dev.off()

#TOTAL Impacts per model

pdf(file="flooddepth-impacts2.pdf")
par(mfrow=c(2,2))

#SEM
plot(fdsem, xlim=c(-.20,0), ylim=c(0, 30), type="l", main = "SEM Model")
lines(density(totimp.sem), lty = 1, lwd = 3)
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1,3), bty = "n")

#SDEM
plot(fdsdem, xlim=c(-.20,0), ylim=c(0, 30), type="l", main ="SDEM Model")
lines(density(totimp.sdem), lty = 1, lwd = 3)
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1,3), bty = "n")

#SLM
curve(dnorm(x, mean=muDfXY, sd=sdDfXY), from = -.2, to = 0, type="l", 
  main = "SLM Model", ylab = "y", ylim=c(0, 30))
lines(density(fit2$total[,1]), lty=1, lwd = 3)#SLM
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1,3), bty = "n")

#SDM
curve(dnorm(x, mean=muDfXY2, sd=sdDfXY2), from = -.2, to = 0, 
  type="l", main = "SDM Model", ylab = "y", ylim=c(0, 30))
lines(density((fit2lag$total[,1]+fit2lag$total[,1+8])), lty=1, lwd = 3)#SDM
legend("topleft", legend = c("INLA", "MCMC"), lwd = c(1,3), bty = "n")

dev.off()


#TOTAL Impacts per model WITH CORRECTED STANDADRD DEVIATIONS
#
#St. Dev. have been corrected to match those of the MCMC
#This probably means that the approximation of the mean is correct
#but that the st. dev. are biased. Probably because of autocorrelation
#between the spatial autocorrelation parameter and the covariate coefficient.


#Correct marginal by setting s.d. to new s.d.
sdmarg <- function(x, orig.mean, orig.sd, new.sd){

   z <- (x-orig.mean)/orig.sd
   z <- z*new.sd+orig.mean
}

#SEM
fdsem.zmarg <- inla.zmarginal(fdsem, TRUE)
fdsem.new <- inla.tmarginal(sdmarg, fdsem, orig.mean = fdsem.zmarg$mean, 
  orig.sd = fdsem.zmarg$sd, new.sd = sd(totimp.sem))

#SDEM
fdsdem.zmarg <- inla.zmarginal(fdsdem, TRUE)
fdsdem.new <- inla.tmarginal(sdmarg, fdsdem, orig.mean = fdsdem.zmarg$mean,
  orig.sd = fdsdem.zmarg$sd, new.sd = sd(totimp.sdem))


#Ratio of standard deviations MCMC/INLA
sd(totimp.sem)/fdsem.zmarg$sd#SLM
sd(totimp.sdem)/fdsdem.zmarg$sd#SDM
sdDfXY/sd(fit2$total[,1])#SEM
sdDfXY2/sd(fit2lag$total[,1])#SDEM

#Correlation bewtween spatial autocorrelation and beta_r
cor(fit1$pdraw, fit1$bdra[,2])#SEM
c(cor(fit1lag$pdraw, fit1lag$bdra[,2]),#SDEM beta_r
  cor(fit1lag$pdraw, fit1lag$bdra[,2+8])#SDEM gamma_r
)
cor(fit2$pdraw, fit2$bdra[,2])#SLM
c(cor(fit2lag$pdraw, fit2lag$bdra[,2]), #SDEM beta_r
  cor(fit2lag$pdraw, fit2lag$bdra[,2+8])#SDEM gamma_r
)


pdf(file="flooddepth-impacts-correctedsd.pdf")
par(mfrow=c(2,2))

#SEM
plot(fdsem, xlim=c(-.20,0), ylim=c(0, 30), type="l", main = "SEM Model")
lines(density(totimp.sem), lty = 1, lwd = 3)
lines(fdsem.new, lty = 3, lwd = 1)
legend("topleft", legend = c("INLA", "INLAsd", "MCMC"), 
   lty = c(1,3,1), lwd = c(1,1, 3), bty = "n")

#SDEM
plot(fdsdem, xlim=c(-.20,0), ylim=c(0, 30), type="l", main ="SDEM Model")
lines(density(totimp.sdem), lty = 1, lwd = 3)
lines(fdsem.new, lty = 3, lwd = 1)
legend("topleft", legend = c("INLA", "INLAsd", "MCMC"), 
   lty = c(1,3,1), lwd = c(1,1, 3), bty = "n")

#SLM
curve(dnorm(x, mean=muDfXY, sd=sdDfXY), from = -.2, to = 0, type="l", 
  main = "SLM Model", ylab = "y", ylim=c(0, 30))
lines(density(fit2$total[,1]), lty=1, lwd = 3)#SLM
curve(dnorm(x, mean=muDfXY, sd= sd(fit2$total[,1])), from = -.2, to = 0, 
  type="l",main = "SLM Model", ylab = "y", ylim=c(0, 30), add =TRUE,
  lty = 3)
legend("topleft", legend = c("INLA", "INLAsd", "MCMC"), 
   lty = c(1,3,1), lwd = c(1,1, 3), bty = "n")

#SDM
curve(dnorm(x, mean=muDfXY2, sd=sdDfXY2), from = -.2, to = 0, 
  type="l", main = "SDM Model", ylab = "y", ylim=c(0, 30))
lines(density((fit2lag$total[,1]+fit2lag$total[,1+8])), lty=1, lwd = 3)#SDM
curve(dnorm(x, mean=muDfXY2, sd= sd(fit2lag$total[,1]+fit2lag$total[,1+8])),
  type="l",add =TRUE,  lty = 3)
legend("topleft", legend = c("INLA", "INLAsd", "MCMC"),
   lty = c(1,3,1), lwd = c(1,1, 3), bty = "n")

dev.off()


