# ASSESSOR: An age-structured state-space stock-recruit model for Pacific salmon
[Mark D. Scheuerell](https://faculty.washington.edu/scheuerl/)  
`r Sys.Date()`  

## Overview
This document outlines how to fit an age-structured state-space stock-recruit model for semelparous species like Pacific salmon (_Oncorhynchus_ spp.). The general structure follows that of Fleischman _et al._ (2013), but it also allows for the inclusion of specific external drivers of productivity, such as climate variability. The model is composed of two primary pieces: a process model that governs the true population dynamics, and an observation model that relates the data in hand to the true process.

### Process component
We begin with our process model that describes the true, but unknown production of offspring from their parents. In any given year _t_, spawning adults produce some number of surviving offspring, which follows a general Ricker model, such that
	
$$\log(R_t) = \log(S_t) + a_t \ – bS_t + w_t.$$
	
Here $R_t$ is the total number of subsequent recruits (offspring) born in year _t_; $S_t$ is the true, but unobserved, number of spawning adults; $a_t$ is the annual density-independent productivity; $b$ is the strength of density dependence; and $w_t$ is a process error representing environmental stochasticity, which is autocorrelated over time according to $w_t \sim \text{N}(\phi w_{t-1}, q_a)$.

Previous applications of time-varying productivity (e.g., Dorner _et al._ 2008, Peterman _et al._ 2003) have used a Markov form where $a_t \sim \text{N}(a_{t-1}, \sigma_a)$, but we will model $(a_t)$ as a function of time-varying covariates. Specifically,

$$a_t = \bar{a} + \sum_{i=1}^{M} c_{i,t} \ X_{i,t+h} $$

Here $\bar{a}$ is the underlying mean productivity, and $c_{i,t}$ is the effect of covariate $i$ at time $t$, $X_{i,t+h}$. To allow for direct comparison of effect sizes, the covariates are typically standardized to have a zero mean and unit variance.

The estimated number of fish of age $a$ returning in year $t$ $(N_{a,t})$ is then product of the total number of brood-year recruits in year $t – a$ and the proportion of mature fish from that brood year that returned to spawn at age $a$ $(p_{a,t-a})$, such that

$$N_{a,t} = R_{t-a} \ p_{a,t-a}.$$

The vector of age-specific return rates for brood year $t$ $(\mathbf{p}_t)$ has length $A$, which equals the number of adult age classes. We modeled the vector of year specific maturity schedules $(\mathbf{p}_t)$ as random effects using a hierarchical form of the Dirichlet distribution, where

$$\mathbf{p}_t \sim \text{Dirichlet}(\mathbf{\mu},\pi).$$

In this formulation, the mean vector $\mathbf{\mu}$ is itself distributed as a Dirichlet, and therefore has a total of $A$ elements that are all greater than zero. The precision parameter $\pi$ affects each of the elements in $\mu$, such that large values of $\pi$ result in values of $\mathbf{p}_t$ very close to $\mathbf{\mu}$.

### Observation component

Estimates of the number of spawning adults necessarily contain some sampling or observation errors due to incomplete censuses, mis-identification, etc. Therefore, we will assume that the estimates of escapement $(E_t)$ are log-normally distributed about the true number of spawners $(S_t)$

$$\log(E_t) \sim \text{Normal}\left(\log(S_t), r_s\right).$$

We do not have the ability to estimate the observation variances for both the escapement and harvest without any additional prior information. Therefore, we will assume the harvest is recorded without error and calculate $S_t$ as the difference between the estimated total run size $(N_t)$ and harvest $(H_t)$

$$S_t = N_t - H_t.$$

and $N_t$ is the sum of $N_{a,t}$ over all age classes.

The age composition data include the number of fish in each age class $a$ in year $t$ $(O_{a,t})$. The age data are then modeled as a multinomial process with order $Y_t$ and proportion vector $\mathbf{d}_t$, such that

$$\mathbf{O}_t \sim \text{Multinomial}(Y_t, \mathbf{d}_t).$$

The order of the multinomial is simply the sum of the observed numbers of fish across all ages returning in year $t$:

$$Y_t = \sum_{a=min}^{max} O_{a,t}$$

The proportion vector $\mathbf{d}_t$ for the multinomial is based on the age-specific, model-derived estimates of adult returns in year $t$ $(N_{a,t})$, such that

$$d_{a,t} = \frac{N_{a,t}}{\sum_{a=min}^{max} N_{a,t}}.$$

## Requirements
Our analyses require the [R software](https://cran.r-project.org/) (v3.2.3) for data retrieval, data processing, and summarizing model results, and the [JAGS software](http://mcmc-jags.sourceforge.net/) (v4.2.0) for Markov chain Monte Carlo (MCMC) simulation in model fitting. Please note that some of the R code below may not work with older versions of JAGS due to some changes in the ways that arrays are handled.

We also need the `R2jags` and `gsl` packages, which are not included with base `R`, so we begin by installing them (if necessary) and then loading them.


```r
if(!require("R2jags")) {
  install.packages("R2jags")
  library("R2jags")
}
if(!require("gsl")) {
  install.packages("gsl")
  library("gsl")
}
```

## User inputs
We begin by specifying the names of four necessary data files that contain the following information:
 
 1. estimated total number of adult spawners (escapement) by year;
 2. estimated age composition of adult spawners by year;
 3. estimated total harvest by year;
 4. covariate(s) by year.

Let's also define the following parameters, which will be referenced throughout the analysis.

 * _TT_: number of calendar years of data
 * _AA_: number of age classes 
 * _MM_: number of covariates


```r
## file with escapement data
## [TT x 2] matrix of obs counts; 1st col is calendar yr
file.esc <- "SkagitSthdEsc.csv"

## file with age comp data
## [TT x (1+AA)]; 1st col is calendar yr
file.age <- "SkagitSthdAge.csv"
## min & max ages
age.min <- 3
age.max <- 8
## years, if any, of age-comp to skip; see below
age.skip <- 2

## file with catch data
## [TT x 2] matrix of obs catch; 1st col is calendar yr
file.harv <- "SkagitSthdCatch.csv"

## file with covariate data
## [TT x (1+MM)]; 1st col is calendar yr
file.cov <- "SkagitEnvCov.csv"

## number of years of forecasts
nFore <- 1

## file where to save JAGS model
file.JAGS <- "SkagitSthd_RR_JAGS.txt"

## upper threshold for Gelman & Rubin's (1992) potential scale reduction factor (Rhat).
Rhat_thresh <- 1.1
```

## Loading the data
Here we load in the four data files and do some simple calculations and manipulations. First the spawner data:


```r
## escapement
datEsc <- read.csv(file.esc)
## years of data
datYrs <- datEsc$year
## number of years of data
TT <- length(datYrs)
## get first & last years
yr.first <- min(datYrs)
yr.last <- max(datYrs)
## log of escapement
datLnEsc <- c(log(datEsc[,-1]),rep(NA,nFore))
```

Next the age composition data:


```r
## age comp data
datAge <- read.csv(file.age)
## drop year col & first age.min+age.skip rows
datAge <- datAge[-(1:(age.min+age.skip)),-1]
## num of age classes
AA <- age.max-age.min+1
## add row(s) of NA's for forecast years
datAge <- rbind(datAge,matrix(0,nFore,AA,dimnames=list(TT+seq(nFore),colnames(datAge))))
## total num of age obs by cal yr
datAge[,"sum"] <- apply(datAge,1,sum)
## row indices for any years with no obs age comp
iYNA <- which(datAge$sum<AA,TRUE)
## replace 0's in yrs w/o any obs with NA's
datAge[iYNA,(1:AA)] <- NA
## change total in yrs w/o any obs from 0 to AA to help dmulti()
datAge[iYNA,"sum"] <- AA
## convert class
datAge <- as.matrix(datAge)
```

Then the harvest data:


```r
## harvest
datHarv <- read.csv(file.harv)
## drop year col & first age.max rows
datHarv <- c(datHarv[,-1],rep(0,nFore))
```

And finally the covariates:


```r
## covariate(s)
datCov <- read.csv(file.cov)
## drop year col
datCov <- datCov[,-1]
## transform the covariates to z-scores
datCov <- scale(datCov)
```

## SSSR model in JAGS

Now we can specify the model in JAGS. Note that the code below is not written to be universally generic with respect to the number of covariates, but rather to emphasize how to incorporate the three in this specific application. The important point is that the number of covariate parameters in the `PRIORS` and `LIKELIHOOD` sections (i.e., there must be a unique `gamma` parameter for each of the _MM_ covariates).


```r
cat("

model {
	
	#--------
	# PRIORS
	#--------
	# alpha = intrinsic productivity
	muAlpha ~ dnorm(0,1e-3)I(-4,4);
	alpha <- exp(muAlpha);
	lnAlphaP <- muAlpha + var.Qr/(2 - 2*phi^2);

	# gamma = covariate effect
	gammaFlow ~ dnorm(0,0.001);
	gammaPDO ~ dnorm(0,0.001);
	gammaHrel ~ dnorm(0,0.001);

	# beta = strength of dens depend
	beta ~ dunif(0,0.1);

	# AR(1) param for proc errors
	phi ~ dunif(-0.999,0.999);
	
	# Qr = process variance for recruits model
	# diffuse gamma prior on precision
	# tau.Qr ~ dgamma(0.001,0.001);
	# var.Qr <- pow(tau.Qr,-1);
	# unif prior on SD
	tau.Qr <- pow(sd.Qr,-2);
	sd.Qr ~ dunif(0.001,20);
	var.Qr <- pow(sd.Qr,2)
	
	# innov in yr 1
	resLnRec0 ~ dnorm(0,tau.Qr*(1-phi*phi));
	
	# Rs = variance for Sp obs model
	# diffuse gamma prior on precision
	# tau.Rs ~ dgamma(0.001,0.001);
	# var.R <- pow(tau.Rs,-1);
	# unif prior on SD
	tau.Rs <- pow(sd.Rs,-2);
	sd.Rs ~ dunif(0.001,20);
	
	# unobservable early total run size
	tRunMu ~ dunif(1,5);
	tRunTau ~ dunif(1,20);
	
	# unprojectable early recruits;
	# hyper mean across all popns
	Rec.mu ~ dnorm(0,0.001);
	# hyper SD across all popns
	Rec.sig ~ dunif(0,100);
	# precision across all popns
	Rec.tau <- pow(Rec.sig,-2);

	# maturity schedule
	# unif vec for Dirch prior
	for(i in 1:AA) { theta[i] <- 1 }
	# hyper-mean for maturity
	muHD ~ ddirch(theta);
	# hyper-prec for maturity
	piHD ~ dunif(0.001,1e3);
	for(t in 1:(TT-age.min+nFore)) { pi[t,1:AA] ~ ddirch(muHD*piHD) }
	
	#------------
	# LIKELIHOOD
	#------------
	# 1st brood yr requires different innovation
	# predicted recruits in BY t
	lnAlpha[1] <- muAlpha + gammaFlow*datCov[1,1] + gammaPDO*datCov[1,2] + gammaHrel*datCov[1,3];
	muLnRec[1] <- lnSp[1] + lnAlpha[1] - beta*Sp[1];
	totLnRec[1] ~ dnorm(muLnRec[1] + phi*resLnRec0,tau.Qr);
	resLnRec[1] <- totLnRec[1] - muLnRec[1];
	# MEDIAN of total recruits
	totRec[1] <- exp(totLnRec[1]);
	# MEAN of total recuits
	# totRec[1] <- exp(totLnRec[1] + var.Qr/2);
	
	# R/S
	lnRS[1] <- totLnRec[1] - lnSp[1];
		
	# brood-yr recruits by age
	for(a in 1:AA) {
		Rec[1,a] <- max(1,totRec[1] * pi[1,a]);
		}
	
	# brood years 2:(TT-age.min)
	for(t in 2:(TT-age.min+nFore)) {
		# predicted recruits in BY t
		lnAlpha[t] <- muAlpha + gammaPDO*datCov[t,1] + gammaFlow*datCov[t,2] + gammaHrel*datCov[t,3];
		muLnRec[t] <- lnSp[t] + lnAlpha[t] - beta*Sp[t];
		totLnRec[t] ~ dnorm(muLnRec[t] + phi*resLnRec[t-1],tau.Qr);
		resLnRec[t] <- totLnRec[t] - muLnRec[t];
		# median of total recruits
		totRec[t] <- exp(totLnRec[t]);
		# R/S
		lnRS[t] <- totLnRec[t] - lnSp[t];
		# brood-yr recruits by age
		for(a in 1:AA) {
			Rec[t,a] <- max(1,totRec[t] * pi[t,a]);
			}
		} # end t loop over year

	# get total cal yr returns for first age.min yrs
	for(i in 1:(age.min+age.skip)) {
		lnTotRun[i] ~ dnorm(tRunMu*Rec.mu,Rec.tau/tRunTau);
		totRun[i] <- exp(lnTotRun[i]);
	}

	# get predicted calendar year returns by age
	# matrix Run has dim [(TT-age.min) x AA]
	# step 1: incomplete early broods
	# first cal yr of this grp is first brood yr + age.min + age.skip
	for(i in 1:(age.max-age.min-age.skip)) {
		# projected recruits
		for(a in 1:(i+age.skip)) {
			Run[i,a] <- Rec[(age.skip+i)-a+1,a];
			}
		# imputed recruits
		for(a in (i+1+age.skip):AA) {
			lnRec[i,a] ~ dnorm(Rec.mu,Rec.tau);
			Run[i,a] <- exp(lnRec[i,a]);
			}
		# total run size
		totRun[i+age.min+age.skip] <- sum(Run[i,1:AA]);
		# predicted age-prop vec for multinom
		for(a in 1:AA) {
			aVec[i,a] <- Run[i,a] / totRun[i+age.min];
			}
		# multinomial for age comp
		datAge[i,1:AA] ~ dmulti(aVec[i,1:AA],datAge[i,AA+1]);
		}
	
	# step 2: info from complete broods
	# first cal yr of this grp is first brood yr + age.max
	for(i in (AA-age.skip):(TT-age.min-age.skip+nFore)) {
		for(a in 1:AA) {
			Run[i,a] <- Rec[(age.skip+i)-a+1,a];
			}
		# total run size
		totRun[i+age.min+age.skip] <- sum(Run[i,1:AA]);
		# predicted age-prop vec for multinom
		for(a in 1:AA) {
			aVec[i,a] <- Run[i,a] / totRun[i+age.min];
			}
		# multinomial for age comp
		datAge[i,1:AA] ~ dmulti(aVec[i,1:AA],datAge[i,AA+1]);
		}
		
	# get predicted calendar year spawners
	# first cal yr is first brood yr
	for(t in 1:(TT+nFore)) {
		# obs model for spawners
		lnSp[t] <- log(max(1,totRun[t] - datHarv[t]));
		Sp[t] <- exp(lnSp[t]);
		datLnEsc[t] ~ dnorm(lnSp[t], tau.Rs);
		}
			
} # end model description

", file=file.JAGS)
```

***
## Fitting the model

The last thing we need to do before fitting the model in JAGS is to specify:

1. the data and indices that go into the model;
2. the model parameters and states that we want JAGS to return;
3. the MCMC control parameters.

Please note that the following code takes ~30 min to run on a quad-core machine.


```r
## data to pass to JAGS
data.JAGS <- c("datAge","datLnEsc","datHarv","datCov",
               "TT","AA","age.min","age.max","age.skip","nFore")

## 2. model params/states for JAGS to return
par.JAGS <- c("alpha","lnAlphaP","beta","Sp","Rec","totLnRec","lnRS",
              "gammaFlow","gammaPDO","gammaHrel",
              "var.Qr","pi","resLnRec")

## 3. MCMC control params
## function to create JAGS inits
init.vals <- function() {
	list(muAlpha=1, gammaFlow=0.1, gammaPDO=-0.1, gammaHrel=-0.2,
	     beta=1/exp(mean(datLnEsc, na.rm=TRUE)),
		piHD=1, muHD=rep(1,AA),
		 pi=matrix(c(0.01,0.3,0.48,0.15,0.05,0.01),TT-age.min+nFore,AA,byrow=TRUE),
		 Rec.mu=log(1000),
		 Rec.sig=0.1,
		 totLnRec=rep(log(1000),TT-age.min+nFore),
		 resLnRec0=0,
		 phi=0.5)
	}

mod.JAGS <- list(data=data.JAGS,
				 inits=init.vals,
				 parameters.to.save=par.JAGS,
				 model.file=file.JAGS,
				 n.chains=as.integer(4),
				 n.burnin=as.integer(5e5),
				 n.thin=as.integer(1000),
				 n.iter=as.integer(10e5),
				 DIC=TRUE)

## start timer
timer_start <- proc.time()

## fit the model in JAGS & store results
mod.fit <- do.call(jags.parallel, mod.JAGS)

## stop timer
(run_time_in_min <- round(((proc.time()-timer_start)/60)["elapsed"], 1))
```

## Model diagnostics


```r
## Rhat values for all parameters
rh <- mod.fit$BUGSoutput$summary[,"Rhat"]
## histogram of Rhat values for all parameters
par(mai=c(0.9,0.9,0.1,0.1))
hist(rh, breaks=seq(1,ceiling(max(rh)/0.01)*0.01,by=0.01),
     xlab=expression(italic(R[hat])))
## Rhat values > threshold
bad_Rhat <- rh[rh>Rhat_thresh]
## prop of params with Rhat > threshold
round(length(bad_Rhat)/length(rh),2)
## param names
par_names <- sub("\\[.*","",names(bad_Rhat))
## number of Rhat > threshold by param name
table(par_names)
## index values for offenders
idx <- as.integer(sub("(^.*\\[)([0-9]{1,3})(.*)","\\2",names(bad_Rhat)))
## data frame of offenders
(df <- data.frame(par=par_names, index=idx))
```

The convergence statistics indicate that some of the elements of the in the estimated proportions of the youngest and oldest age classes (i.e., 3 and 8, respectively) did not converge to our desired threshold. However, there is very little data to inform those parameters, so we should not be too concerned.

## Results

### Ricker _a_

Here is the mean instrinsic productivity (i.e., the slope of the curve near the origin).


```r
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
par(mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
RaDat <- sort(mod.fit$BUGSoutput$sims.list$alpha)
RaDat[RaDat>10] <- 10.2
alphaMean <- median(RaDat)
alphaCI <- RaDat[c(floor(mcmc.samp*0.025),ceiling(mcmc.samp*0.975))]
hist(RaDat, freq=FALSE, xlab="", main="", breaks=seq(0,max(RaDat),0.2), col=clr, border="blue3",
	 ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/10
arrows(alphaMean,par()$usr[3],alphaMean,par()$usr[3]-aHt, code=1, length=0.05, xpd=NA, col="blue3")
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt, code=1, length=0.05, xpd=NA, col="blue3")
mtext("Ricker a", 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)
```


### Ricker _b_

Here is the estimated strength of density dependence.


```r
par(mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
RbDat <- sort(mod.fit$BUGSoutput$sims.list$beta)
RbDat <- RbDat*10^abs(floor(log(max(RbDat),10)))
ylM <- max(RbDat)
brks <- seq(0,ceiling(ylM),0.1)
betaMean <- median(RbDat)
betaCI <- RbDat[c(floor(mcmc.samp*0.025),ceiling(mcmc.samp*0.975))]
hist(RbDat, freq=FALSE, breaks=brks, col=clr, border="blue3",
	 xlab="", xaxt="n", yaxt="n",
	 main="", ylab="Posterior density", cex.lab=1.2)
axis(1, at=seq(0,3))
aHt <- (par()$usr[4]-par()$usr[3])/10
arrows(betaMean,par()$usr[3]-0.005,betaMean,par()$usr[3]-aHt, code=1, length=0.05, xpd=NA, col="blue3")
arrows(betaCI,par()$usr[3]-0.005,betaCI,par()$usr[3]-aHt, code=1, length=0.05, xpd=NA, col="blue3")
mtext(expression(paste("Ricker b ",(10^{-4}),"")), 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)
```

### Covariate effects

Here are the effects of the three covariates on productivity.


```r
par(mfrow=c(length(covars),2), mai=c(0.4,0.2,0.1,0.1), omi=c(0.2,0.4,0,0))
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
offSet <- 0.07
datCov <- datCov[,c(2,3,1)]
covars <- grep("gamma",names(mod.fit$BUGSoutput$sims.list),value=TRUE)
ylN <- NULL
ylM <- NULL
for(i in 1:length(covars)) {
	ylN <- min(ylN,mod.fit$BUGSoutput$sims.list[[covars[i]]])
	ylM <- max(ylN,mod.fit$BUGSoutput$sims.list[[covars[i]]])
}
ylN <- floor(ylN*10)/10
ylM <- ceiling(ylM*10)/10
brks <- seq(ylN,ylM,length.out=41)
cov_names <- c("Flow","H releases","PDO")
tSeries <- seq(yr.first,length.out=TT-age.min)
for(i in 1:length(covars)) {
	# plot covar ts
	plot(tSeries, datCov[seq(length(tSeries)),i], xlab="", ylab="",
		 main="", cex.lab=1.3, pch=16, col="blue3", type="o")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
		 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i])
	mtext(side=2, cov_names[i], line=3)
	if(i==length(covars)) {
		mtext(side=1,"Brood year", line=3)
	}
	# plot covar effect
	covEff <- sort(mod.fit$BUGSoutput$sims.list[[covars[i]]])
	hist(covEff, freq=FALSE, breaks=brks, col=clr, border="blue3",
		 xlab="", yaxt="n",
		 main="", ylab="", cex.lab=1.2)
	abline(v=0, lty="dashed")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
		 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i+length(covars)])
	if(i==length(covars)) {
		mtext(side=1,"Effect size", line=3)
	}
}
```

### Spawners

Here is our estimate of the number of spawners over time, which includes a forecast for 2016. The black points are the data, the blue line is the median posterior estimate, and the shaded region is the 95% credible interval.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)[,(1:TT)]
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
ypMin <- min(pDat)
ypMax <- max(pDat)
tSeries <- seq(yr.first,length.out=TT)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", xaxt="n", yaxt="n",
	 xlab="Year", ylab="Spawners", main="", cex.lab=1.2)
polygon(c(tSeries,rev(tSeries)),c(pDat[3,1:TT],rev(pDat[1,1:TT])), col=clr, border=NA)
lines(tSeries, pDat[2,1:TT], col="blue3", lwd=2)
points(seq(yr.first,length.out=TT+nFore), exp(datLnEsc), pch=16, cex=1)
axis(1,at=seq(1980,2015,5))
axis(2,at=c(3000,6000,12000))
```

### Total run size

Here is our estimate of the total run size (i.e., catch + escapement) over time, which includes a forecast for 2016. The black points are the data, the blue line is the median posterior estimate, and the shaded region is the 95% credible interval.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
pDat <- pDat + matrix(datHarv,length(CI.vec),TT+nFore,byrow=TRUE)
tSeries <- seq(yr.first,length.out=TT+nFore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", xaxt="n", yaxt="n",
	 xlab="Year", ylab="Catch + escapement", main="", cex.lab=1.2)
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
points(tSeries, exp(datLnEsc)+datHarv, pch=16, cex=1)
axis(1,at=seq(1980,2015,5))
axis(2,at=c(4000,8000,16000))
```

### Recruits by age class

Here are the estimated number of recruits by brood year and age. Note that the uncertainty increases in more recent years as fewer complete age classes have been observed.


```r
CI.vec <- c(0.05,0.5,0.95)
par(mfrow=c(AA,1), mai=c(0.1,0.1,0.1,0.1), omi=c(0.5,0.5,0,0))
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
nRec <- TT-age.min
tSeries <- seq(yr.first,length.out=nRec+nFore)
pltTT <- seq(min(round(tSeries/5,0)*5),max(round(tSeries/5,0)*5),5)
for(i in rev(1:AA)) {
	pDat <- apply(mod.fit$BUGSoutput$sims.list$Rec[,,i],2,sort)
	pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })/100
	dd <- ifelse(max(pDat)<20,1,10)
	ypMax <- real2prec(max(pDat),prec=dd)
	while(ypMax %% 3 != 0) { ypMax <- ypMax + dd }
	plot(tSeries,pDat[3,], xlim=c(yr.first+1,yr.last-nFore-2), ylim=c(0.001,ypMax), type="n",
		 xaxt="n", yaxt="n", xlab="", ylab="", main="", las=1)
	polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
	lines(tSeries, pDat[2,], col="blue3", lwd=2)
	aHt <- (par()$usr[4]-par()$usr[3])/7
	ttl <- paste("Age-",i+age.min-1,sep="")
	text(tSeries[1]-0, par()$usr[4]-aHt, ttl, pos=4, cex=0.9)
	axis(2,seq(0,ypMax,length.out=4),las=1,cex=0.9)
	if(i!=1) {axis(1,at=pltTT,labels=FALSE)} else {axis(1,at=pltTT)}
}
mtext("Recruits (100s)", 2, line=2, outer=TRUE, cex=1.2)
mtext("Year", 1, line=2.5, outer=TRUE, cex=1.2)
```

### Total recruits

Here are the estimated total number of recruits by brood year.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Rec,c(1,2),sum)
pDat <- apply(apply(pDat,2,sort),2,function(x) { x[mcmc.samp*CI.vec] })
tSeries <- seq(yr.first,length.out=TT-age.min+nFore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", yaxt="n",
	 xlab="Brood year", ylab="Recruits", main="", cex.lab=1.2)
axis(2,at=c(3000,9000,27000))
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

### Recruits per spawner

Here is the time series of estimated recruits-per-spawner.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$lnRS,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
pDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$lnRS,2,median)
tSeries <- seq(yr.first,length.out=TT-age.min+nFore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="ln(R/S)", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

### Innovations

Here is the time series of the so-called "innovations", which are the residuals from the process model. They give some indications of population productivity after accounting for the effects of density dependence.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$resLnRec,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
pDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$resLnRec,2,median)
tSeries <- seq(yr.first,length.out=TT-age.min+nFore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="Innovations", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

### Age composition

Here are time series of the estimated proportions of each age class.


```r
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.1,0.2))
CI.vec <- c(0.025,0.5,0.975)
tSeries <- seq(yr.first,length.out=TT-age.min+nFore)
clr <- rgb(0, 0, 255, alpha = 40, maxColorValue = 255)
ageComp <- t(apply(apply(mod.fit$BUGSoutput$sims.list$pi,c(3,2),mean),2,cumsum))
plot(tSeries, rep(1,nRec+nFore), ylab="Proportion", xlab="Brood year", ylim=c(0,1), las=1,
		xaxs="i", yaxs="i", type="n", lty="solid", col="blue3", cex.lab=1.2)
for(i in c(1,2,3,4,6)) {
	polygon(c(tSeries,rev(tSeries)),c(ageComp[,i],rep(0,nRec+nFore)), col=clr, border=NA)
	}
lbl <- apply(cbind(c(0,ageComp[nRec+nFore,-AA]),ageComp[nRec+nFore,]),1,mean)
text(par()$usr[2],par()$usr[4]*1.05,"Age", xpd=NA, pos=4, offset=0.05, col="black", cex=0.8)
text(par()$usr[2],lbl[1:4],seq(3,6), xpd=NA, pos=4, col="black", cex=0.7)
text(par()$usr[2],lbl[5],"7&8", xpd=NA, pos=4, offset=0.15, col="black", cex=0.7)
```

### Spawner-recruit relationship

Here is the relationship between spawner and subsequent recruits.


```r
CI.vec <- c(0.025,0.5,0.975)
sDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)
sDat <- apply(sDat,2,function(x) { x[mcmc.samp*CI.vec] })
sDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,median)
sDat <- sDat[,1:(TT-age.min+nFore)]
rDat <- apply(mod.fit$BUGSoutput$sims.list$totLnRec,2,sort)
rDat <- apply(rDat,2,function(x) { x[mcmc.samp*CI.vec] })
rDat <- exp(rDat)
dd <- 3000
yM <- real2prec(max(rDat),"ceiling",dd)
yM <- 30000
xM <- real2prec(max(sDat),"ceiling",dd)
aa <- exp(matrix(mod.fit$BUGSoutput$sims.array[,,"lnAlpha"],ncol=1))
bb <- matrix(mod.fit$BUGSoutput$sims.array[,,"beta"], ncol=1)
idx <- sample(seq(1000),MC)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.1,0.2))
plot(sDat[2,],rDat[2,], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3", type="n",
	 xaxs="i", yaxs="i", ylab="Recruits (1000s)", xlab="Spawners (1000s)", cex.lab=1.2,
	 xaxt="n", yaxt="n")
axis(1, at=seq(0,xM,dd*2), labels=seq(0,xM,dd*2)/1000)
axis(2, at=seq(0,yM,dd*2), labels=seq(0,yM,dd*2)/1000)
for(i in 1:MC) { lines(simRS(aa[idx[i]],bb[idx[i]]), col="darkgray") }
abline(a=0,b=1,lty="dashed")
nCB <- TT-age.max
points(sDat[2,1:nCB],rDat[2,1:nCB], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3")
segments(sDat[2,1:nCB],rDat[1,1:nCB],sDat[2,1:nCB],rDat[3,1:nCB], col="blue3")
segments(sDat[1,1:nCB],rDat[2,1:nCB],sDat[3,1:nCB],rDat[2,1:nCB], col="blue3")
nTB <- dim(sDat)[2]
clr <- rgb(100, 0, 200, alpha = seq(200,100,length.out=age.max-age.min+nFore), maxColorValue = 255)
segments(sDat[2,(nCB+1):nTB],rDat[1,(nCB+1):nTB],sDat[2,(nCB+1):nTB],rDat[3,(nCB+1):nTB], col=clr)
segments(sDat[1,(nCB+1):nTB],rDat[2,(nCB+1):nTB],sDat[3,(nCB+1):nTB],rDat[2,(nCB+1):nTB], col=clr)
points(sDat[2,(nCB+1):nTB],rDat[2,(nCB+1):nTB],
		xlim=c(0,xM), ylim=c(0,yM), pch=16, col=clr)
```

### Management reference points

Here are a number of management reference points.


```r
# abbreviations for ref points
refNames <- c("Smsy","Umsy","MSY","Umax","MSR","Smsr","Seq")
# proportions of MSY to consider
yieldProps <- c(0.7,0.8,0.9)
propNames <- paste("OYP",yieldProps,sep="")
lnA <- matrix(mod.fit$BUGSoutput$sims.array[,,"lnAlphaP"],ncol=1)
bb <- matrix(mod.fit$BUGSoutput$sims.array[,,"beta"], ncol=1)
mcmc <- length(lnA)
# empty matrix for ref pts
ref.pts <- matrix(NA,mcmc,length(c(refNames,propNames)))
colnames(ref.pts) <- c(refNames,propNames)
# spawner series for optimal yield profile
SS <- seq(100,1e4,100)
# empty matrix for OLP
OYP <- matrix(0,length(SS),length(yieldProps))
for(i in 1:mcmc) {
	# spawners at MSY
	ref.pts[i,"Smsy"] <- (1 - lambert_W0(exp(1-lnA[i]))) / bb[i]
	# harvest rate at MSY
	ref.pts[i,"Umsy"] <- (1 - lambert_W0(exp(1-lnA[i])))
	# MSY
	ref.pts[i,"MSY"] <- ref.pts[i,"Smsy"]*((exp(lnA[i]-bb[i]*ref.pts[i,"Smsy"])) - 1)
	# maximum harvest rate
	ref.pts[i,"Umax"] <- 1 - 1/lnA[i]
	# maximum sustained production
	ref.pts[i,"MSR"] <- lnA[i]/(bb[i]*exp(1))
	# spawners at MSR
	ref.pts[i,"Smsr"] <- 1/bb[i]
	# equil spawning abundance
	ref.pts[i,"Seq"] <- lnA[i]/bb[i]
	# yield over varying S
	yield <- SS*(exp(lnA[i]-bb[i]*SS) - 1)
	for(j in 1:length(yieldProps)) {
		OYP[,j] <- OYP[,j] + 1*(yield > yieldProps[j]*ref.pts[i,"MSY"])
	}
}
OYP <- OYP/mcmc
```

### Optimal yield profiles

Here are the estimated optimal yield profiles based on the alternative states of nature.


```r
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0,0,0.2))
matplot(SS, OYP, type="l", lty="solid", las=1, col=c("slateblue","blue","darkblue"),
		xlab="Escapement", ylab="Probability", lwd=2)
points(x=c(4900,5700,6450), y=c(0.6,0.5,0.4), pch=21, cex=3.5, col="white", bg="white")
text(x=c(4900,5700,6450), y=c(0.6,0.5,0.4), c("90%","80%","70%"),
	 col=c("darkblue","blue","slateblue"), cex=0.7)
```

