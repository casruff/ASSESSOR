# ASSESSOR: An age-structured state-space stock-recruit model for Pacific salmon

***

### Contributors
[__Mark D. Scheuerell__](https://faculty.washington.edu/scheuerl/)  
_Fish Ecology Division, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, Seattle, WA USA, mark.scheuerell@noaa.gov_

__Casey P. Ruff__  
_Skagit River System Cooperative, Mt. Vernon, WA USA, cruff@skagitcoop.org_

__Joseph H. Anderson__  
_Washington Department of Fish and Wildlife, Olympia, WA USA, joseph.anderson@dfw.wa.gov_

__Eric M. Beamer__   
_Skagit River System Cooperative, Mt. Vernon, WA USA, ebeamer@skagitcoop.org_

### Acknowledgments

__Eric R. Buhle__  
_Ocean Associates, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, Seattle, WA USA, eric.buhle@noaa.gov_

[__James T. Thorson__](https://sites.google.com/site/thorsonresearch/)  
_Fishery Resource Analysis and Monitoring Division, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, Seattle, WA USA, james.thorson@noaa.gov_

### Version
This is version 0.16.05.02.

***

## Overview
ASSESSOR incorporates data on spawners (escapement), harvest, and age composition into a retrospective run reconstruction and probabilistic forecast under a Bayesian framework. The general structure follows that of Fleischman et al. [-@fleischman2013], but ASSESSOR also allows for the inclusion of specific external drivers of productivity, both natural (e.g., climate variability) and anthropogenic (e.g., flow alteration). The model is composed of two primary pieces: a process model that governs the true population dynamics, and an observation model that relates the data in hand to the true process.

### Process component
We begin with our process model that describes the true, but unknown production of offspring from their parents. In any given year _t_, spawning adults produce some number of surviving offspring, which follows a general Ricker model, such that
	
$$\log(R_t) = \log(S_t) + a_t \ – bS_t + w_t.$$
	
Here $R_t$ is the total number of subsequent recruits (offspring) born in year _t_; $S_t$ is the true, but unobserved, number of spawning adults; $a_t$ is the annual density-independent productivity; $b$ is the strength of density dependence; and $w_t$ is a process error representing environmental stochasticity, which is autocorrelated over time according to $w_t \sim \text{N}(\phi w_{t-1}, q_a)$.

Previous applications of time-varying productivity [e.g., @peterman2003; @dorner2008] have used a Markov form where $a_t \sim \text{N}(a_{t-1}, \sigma_a)$, but we will model $(a_t)$ as a function of time-varying covariates. Specifically,

$$a_t = \bar{a} + \sum_{i=1}^{M} c_{i,t} \ X_{i,t+h} $$

Here $\bar{a}$ is the underlying mean productivity, and $c_{i,t}$ is the effect of covariate $i$ at time $t$, $X_{i,t+h}$. To allow for direct comparison of effect sizes, the covariates are typically standardized to have a zero mean and unit variance.

The estimated number of fish of age $a$ returning in year $t$ $(N_{a,t})$ is then product of the total number of brood-year recruits in year $t – a$ and the proportion of mature fish from that brood year that returned to spawn at age $a$ $(p_{a,t-a})$, such that

$$N_{a,t} = R_{t-a} \ p_{a,t-a}.$$

The vector of age-specific return rates for brood year $t$ $(\mathbf{p}_t)$ has length $A$, which equals the number of adult age classes. That is, $\mathbf{p}_t$ is a combination of the probability of surviving to, and maturing in years $t + a_\min$ to $t + a_\max$. We modeled $(\mathbf{p}_t)$ as a random effect using a hierarchical form of the Dirichlet distribution, where

$$\mathbf{p}_t \sim \text{Dirichlet}(\boldsymbol{\mu},\pi).$$

In this formulation, the mean vector $\boldsymbol{\mu}$ is itself distributed as a Dirichlet, and therefore has a total of $A$ elements that are all greater than zero. The precision parameter $\pi$ affects each of the elements in $\boldsymbol{\mu}$, such that large values of $\pi$ results in values of $\mathbf{p}_t$ that are very close to $\boldsymbol{\mu}$.

### Observation component

Estimates of the number of spawning adults necessarily contain some sampling or observation errors due to incomplete censuses, mis-identification, etc. Therefore, we will assume that the estimates of escapement $(E_t)$ are log-normally distributed about the true number of spawners $(S_t)$

$$\log(E_t) \sim \text{Normal}\left(\log(S_t), r_s\right).$$

We do not have the ability to estimate the observation variances for both the escapement and harvest without any additional prior information. Therefore, we will assume the harvest is recorded without error and calculate $S_t$ as the difference between the estimated total run size $(N_t)$ and harvest $(H_t)$

$$S_t = N_t - H_t.$$

and $N_t$ is the sum of $N_{a,t}$ over all age classes.

The age composition data include the number of fish in each age class $a$ in year $t$ $(O_{a,t})$. The age data are then modeled as a multinomial process with order $Y_t$ and proportion vector $\mathbf{d}_t$, such that

$$\mathbf{O}_t \sim \text{Multinomial}(Y_t, \mathbf{d}_t).$$

The order of the multinomial is simply the sum of the observed numbers of fish across all ages returning in year $t$:

$$Y_t = \sum_{a=a_\min}^{a_\max} O_{a,t}$$

The proportion vector $\mathbf{d}_t$ for the multinomial is based on the age-specific, model-derived estimates of adult returns in year $t$ $(N_{a,t})$, such that

$$d_{a,t} = {N_{a,t} \over \displaystyle \sum_{a=a_\min}^{a_\max} N_{a,t}}.$$

## Requirements
ASSESSOR requires the [R software](https://cran.r-project.org/) (v3.2.3) for data retrieval, data processing, and summarizing model results, and the [JAGS software](http://mcmc-jags.sourceforge.net/) (v4.2.0) for Markov chain Monte Carlo (MCMC) simulation. Please note that some of the R code below may not work with older versions of JAGS due to some changes in the ways that arrays are handled.

We also need a few packages that are not included with base `R`, so we begin by installing them (if necessary) and then loading them.


```r
if(!require("R2jags")) {
  install.packages("R2jags")
  library("R2jags")
}
if(!require("RCurl")) {
  install.packages("RCurl")
  library("RCurl")
}
if(!require("gsl")) {
  install.packages("gsl")
  library("gsl")
}
```

The last thing we'll need is the following function for trimming real numbers to a desired precision when plotting output.


```r
real2prec <- function(x,map="round",prec=1) {
  if(prec==0) { stop("\"prec\" cannot be 0") }
  do.call(map,list(x/prec))*prec
}
```

## User inputs
We begin by specifying the names of four necessary data files that contain the following information:
 
 1. observed total number of adult spawners (escapement) by year;
 2. observed age composition of adult spawners by year;
 3. observed total harvest by year;
 4. measured covariate(s) by year.

Let's also define the following parameters, which will be referenced throughout the analysis.

 * `n_yrs`: number of calendar years of data
 * `A`: number of age classes 
 * `M`: number of covariates

The example data we use here are for steelhead trout (_Oncorhynchus mykiss_) from the Skagit River basin, which drains into northern Puget Sound.


```r
## 1. file with escapement data
## [n_yrs x 2] matrix of obs counts; 1st col is calendar yr
fn_esc <- "SkagitSthdEsc.csv"

## 2. file with age comp data
## [n_yrs x (1+A)]; 1st col is calendar yr
fn_age <- "SkagitSthdAge.csv"
## min & max ages
age_min <- 3
age_max <- 8
## years, if any, of age-comp to skip; see below
age_skip <- 2

## 3. file with harvest data
## [n_yrs x 2] matrix of obs catch; 1st col is calendar yr
fn_harv <- "SkagitSthdCatch.csv"

## 4. file with covariate data
## [n_yrs x (1+MM)]; 1st col is calendar yr
fn_cvrs <- "SkagitEnvCov.csv"

## number of years of forecasts
n_fore <- 1

## file where to save JAGS model
fn_jags <- "SkagitSthd_RR_JAGS.txt"

## upper threshold for Gelman & Rubin's (1992) potential scale reduction factor (Rhat).
Rhat_thresh <- 1.1

## URL for example data files
## set to NULL if using a local folder/directory
ex_url <- "https://raw.githubusercontent.com/mdscheuerell/ASSESSOR/master/"
```

## Loading the data
Here we load in the four data files and do some simple calculations and manipulations. First the spawner data:


```r
## escapement
#dat_esc <- read.table(paste0(ex_url,fn_esc), header=TRUE, sep=",")
dat_esc <- read.csv(text = getURL(paste0(ex_url,fn_esc)))
## years of data
dat_yrs <- dat_esc$year
## number of years of data
n_yrs <- length(dat_yrs)
## get first & last years
yr_frst <- min(dat_yrs)
yr_last <- max(dat_yrs)
## log of escapement
ln_dat_esc <- c(log(dat_esc[,-1]),rep(NA,n_fore))
```

Next the age composition data:


```r
## age comp data
dat_age <- read.csv(text = getURL(paste0(ex_url,fn_age)))
## drop year col & first age_min+age_skip rows
dat_age <- dat_age[-(1:(age_min+age_skip)),-1]
## num of age classes
A <- age_max-age_min+1
## add row(s) of NA's for forecast years
dat_age <- rbind(dat_age,matrix(0,n_fore,A,dimnames=list(n_yrs+seq(n_fore),colnames(dat_age))))
## total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age,1,sum)
## row indices for any years with no obs age comp
idx_NA_yrs <- which(dat_age$sum<A,TRUE)
## replace 0's in yrs w/o any obs with NA's
dat_age[idx_NA_yrs,(1:A)] <- NA
## change total in yrs w/o any obs from 0 to A to help dmulti()
dat_age[idx_NA_yrs,"sum"] <- A
## convert class
dat_age <- as.matrix(dat_age)
```

Then the harvest data:


```r
## harvest
dat_harv <- read.csv(text = getURL(paste0(ex_url,fn_harv)))
## drop year col & first age_max rows
dat_harv <- c(dat_harv[,-1],rep(0,n_fore))
```

And finally the covariates:


```r
## covariate(s)
dat_cvrs <- read.csv(text = getURL(paste0(ex_url,fn_cvrs)))
## drop year col
dat_cvrs <- dat_cvrs[,-1]
## transform the covariates to z-scores
dat_cvrs <- scale(dat_cvrs)
```

## SSSR model in JAGS

Now we can specify the model in JAGS. Note that the code below is not written to be universally generic with respect to the number of covariates, but rather to emphasize how to incorporate the three in this specific application. The important thing is the number of covariate parameters in the `PRIORS` and `LIKELIHOOD` sections (i.e., there must be a unique `c_` parameter for each of the _MM_ covariates).


```r
cat("

model {
	
	#--------
	# PRIORS
	#--------
	# alpha = intrinsic productivity
	Rkr_a ~ dnorm(0,1e-3)I(-4,4);
	alpha <- exp(Rkr_a);
	mu_Rkr_a <- Rkr_a + var_Qr/(2 - 2*phi^2);

	# gamma = covariate effect
	c_Flow ~ dnorm(0,0.001);
	c_PDO ~ dnorm(0,0.001);
	c_Hrel ~ dnorm(0,0.001);

	# beta = strength of dens depend
	Rkr_b ~ dunif(0,0.1);

	# AR(1) param for proc errors
	phi ~ dunif(-0.999,0.999);
	
	# Qr = process variance for recruits model
	sd_Qr ~ dunif(0.001,20);
	tau_Qr <- pow(sd_Qr,-2);
	var_Qr <- pow(sd_Qr,2)
	
	# innov in yr 1
	innov_1 ~ dnorm(0,tau_Qr*(1-phi*phi));
	
	# Rs = variance for Sp obs model
	sd_Rs ~ dunif(0.001,20);
	tau_Rs <- pow(sd_Rs,-2);
	var_Rs <- pow(sd_Rs,2)
	
	# unobservable early total run size
	ttl_run_mu ~ dunif(1,5);
	ttl_run_tau ~ dunif(1,20);
	
	# unprojectable early recruits;
	# hyper mean across all popns
	Rec_mu ~ dnorm(0,0.001);
	# hyper SD across all popns
	Rec_sig ~ dunif(0,100);
	# precision across all popns
	Rec_tau <- pow(Rec_sig,-2);

	# maturity schedule
	# unif vec for Dirch prior
	for(i in 1:A) { theta[i] <- 1 }
	# hyper-mean for maturity
	muHD ~ ddirch(theta);
	# hyper-prec for maturity
	piHD ~ dunif(0.001,1e3);
	for(t in 1:(n_yrs-age_min+n_fore)) { p_vec[t,1:A] ~ ddirch(muHD*piHD) }
	
	#------------
	# LIKELIHOOD
	#------------
	# 1st brood yr requires different innovation
	# predicted recruits in BY t
	ln_Rkr_a[1] <- Rkr_a + c_Flow*dat_cvrs[1,1] + c_PDO*dat_cvrs[1,2] + c_Hrel*dat_cvrs[1,3];
	E_ln_Rec[1] <- ln_Sp[1] + ln_Rkr_a[1] - Rkr_b*Sp[1];
	tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1] + phi*innov_1,tau_Qr);
	res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
	# MEDIAN of total recruits
	tot_Rec[1] <- exp(tot_ln_Rec[1]);

	# R/S
	ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
		
	# brood-yr recruits by age
	for(a in 1:A) {
		Rec[1,a] <- max(1,tot_Rec[1] * p_vec[1,a]);
		}
	
	# brood years 2:(n_yrs-age_min)
	for(t in 2:(n_yrs-age_min+n_fore)) {
		# predicted recruits in BY t
		ln_Rkr_a[t] <- Rkr_a + c_Flow*dat_cvrs[t,1] + c_PDO*dat_cvrs[t,2] + c_Hrel*dat_cvrs[t,3];
		E_ln_Rec[t] <- ln_Sp[t] + ln_Rkr_a[t] - Rkr_b*Sp[t];
		tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t] + phi*res_ln_Rec[t-1],tau_Qr);
		res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
		# median of total recruits
		tot_Rec[t] <- exp(tot_ln_Rec[t]);
		# R/S
		ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
		# brood-yr recruits by age
		for(a in 1:A) {
			Rec[t,a] <- max(1,tot_Rec[t] * p_vec[t,a]);
			}
		} # end t loop over year

	# get total cal yr returns for first age_min yrs
	for(i in 1:(age_min+age_skip)) {
		ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
		tot_Run[i] <- exp(ln_tot_Run[i]);
	}

	# get predicted calendar year returns by age
	# matrix Run has dim [(n_yrs-age_min) x A]
	# step 1: incomplete early broods
	# first cal yr of this grp is first brood yr + age_min + age_skip
	for(i in 1:(age_max-age_min-age_skip)) {
		# projected recruits
		for(a in 1:(i+age_skip)) {
			Run[i,a] <- Rec[(age_skip+i)-a+1,a];
			}
		# imputed recruits
		for(a in (i+1+age_skip):A) {
			lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
			Run[i,a] <- exp(lnRec[i,a]);
			}
		# total run size
		tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
		# predicted age-prop vec for multinom
		for(a in 1:A) {
			age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
			}
		# multinomial for age comp
		dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
		}
	
	# step 2: info from complete broods
	# first cal yr of this grp is first brood yr + age_max
	for(i in (A-age_skip):(n_yrs-age_min-age_skip+n_fore)) {
		for(a in 1:A) {
			Run[i,a] <- Rec[(age_skip+i)-a+1,a];
			}
		# total run size
		tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
		# predicted age-prop vec for multinom
		for(a in 1:A) {
			age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
			}
		# multinomial for age comp
		dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
		}
		
	# get predicted calendar year spawners
	# first cal yr is first brood yr
	for(t in 1:(n_yrs+n_fore)) {
		# obs model for spawners
		ln_Sp[t] <- log(max(1,tot_Run[t] - dat_harv[t]));
		Sp[t] <- exp(ln_Sp[t]);
		ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_Rs);
		}
			
} # end model description

", file=fn_jags)
```

***
## Fitting the model

The last thing we need to do before fitting the model in JAGS is to specify:

1. the data and indices that go into the model;
2. the model parameters and states that we want JAGS to return;
3. the MCMC control parameters.

Please note that the following code takes ~20 min to run on a quad-core machine with 3.5 GHz Intel Core i7 processors.




```r
## data to pass to JAGS
dat_jags <- c("dat_age","ln_dat_esc","dat_harv","dat_cvrs",
              "n_yrs","A","age_min","age_max","age_skip","n_fore")

## 2. model params/states for JAGS to return
par_jags <- c("alpha","mu_Rkr_a","Rkr_b","Sp","Rec","tot_ln_Rec","ln_RS",
              "c_Flow","c_PDO","c_Hrel",
              "var_Qr","var_Rs","p_vec","res_ln_Rec")

## 3. MCMC control params
# MCMC parameters
mcmc_chains <- 4
mcmc_length <- 10e5
mcmc_burn <- 5e5
mcmc_thin <- 1000
# total number of MCMC samples
mcmc_samp <- (mcmc_length-mcmc_burn)*mcmc_chains/mcmc_thin

## function to create JAGS inits
init_vals <- function() {
	list(Rkr_a=1, c_Flow=0, c_PDO=0, c_Hrel=0,
	     Rkr_b=1/exp(mean(ln_dat_esc, na.rm=TRUE)),
	     piHD=1, muHD=rep(1,A),
	     p_vec=matrix(c(0.01,0.3,0.48,0.15,0.05,0.01),n_yrs-age_min+n_fore,A,byrow=TRUE),
	     Rec_mu=log(1000),
	     Rec_sig=0.1,
	     tot_ln_Rec=rep(log(1000),n_yrs-age_min+n_fore),
	     innov_1=0,
	     phi=0.5)
	}

mod_jags <- list(data=dat_jags,
				 inits=init_vals,
				 parameters.to.save=par_jags,
				 model.file=fn_jags,
				 n.chains=as.integer(mcmc_chains),
				 n.iter=as.integer(mcmc_length),
				 n.burnin=as.integer(mcmc_burn),
				 n.thin=as.integer(mcmc_thin),
				 DIC=TRUE)

## fit the model in JAGS & store results
mod_fit <- do.call(jags.parallel, mod_jags)
```



## Model diagnostics

Here is a histogram of the Gelman & Rubin statistics $(R_{hat})$ for the estimated parameters. Recall that we set an upper threshold of 1.1, so values larger than that deserve some additional inspection.


```r
## Rhat values for all parameters
rh <- mod_fit$BUGSoutput$summary[,"Rhat"]
## histogram of Rhat values for all parameters
par(mai=c(0.9,0.9,0.3,0.1))
hist(rh, breaks=seq(1,ceiling(max(rh)/0.01)*0.01,by=0.01),main="",
     col=rgb(0,0,255,alpha=50,maxColorValue=255),border="blue3",xlab=expression(italic(R[hat])))
```

![](sssr_files/figure-html/model_diagnostics-1.png)<!-- -->

```r
## Rhat values > threshold
bad_Rhat <- rh[rh>Rhat_thresh]
## prop of params with Rhat > threshold
round(length(bad_Rhat)/length(rh),2)
```

```
## [1] 0.03
```

```r
## param names
par_names <- sub("\\[.*","",names(bad_Rhat))
## number of Rhat > threshold by param name
table(par_names)
```

```
## par_names
## p_vec 
##    17
```

```r
## index values for offenders
idx <- as.integer(sub("(^.*\\[)([0-9]{1,3})(.*)","\\2",names(bad_Rhat)))
## data frame of offenders
(df <- data.frame(par=par_names, index=idx))
```

```
##      par index
## 1  p_vec     2
## 2  p_vec     3
## 3  p_vec     6
## 4  p_vec     7
## 5  p_vec    12
## 6  p_vec    15
## 7  p_vec    17
## 8  p_vec    18
## 9  p_vec    19
## 10 p_vec    20
## 11 p_vec    21
## 12 p_vec    23
## 13 p_vec    24
## 14 p_vec    26
## 15 p_vec    29
## 16 p_vec    32
## 17 p_vec    33
```

The convergence statistics indicate that some of the elements in $p$ the estimated proportions of the youngest and oldest age classes (i.e., 3 and 8, respectively) did not converge to our desired threshold. However, there is very little data to inform those parameters, so we should not be too concerned.

## Results

Here is a table of summary statistics for some of the model parameters.


```r
print(mod_fit$BUGSoutput$summary[c("mu_Rkr_a","alpha","Rkr_b",
                                   "c_Flow","c_PDO","c_Hrel",
                                   "var_Qr","var_Rs"),
                                 c("mean","sd","2.5%","50%","97.5%")],
      digits=3,quote=FALSE,justify="right")
```

```
##               mean       sd      2.5%       50%     97.5%
## mu_Rkr_a  1.239865 2.93e-01  6.91e-01  1.244175  1.706699
## alpha     3.345830 9.80e-01  1.82e+00  3.278258  5.134761
## Rkr_b     0.000141 3.43e-05  6.99e-05  0.000143  0.000202
## c_Flow   -0.174193 7.04e-02 -3.15e-01 -0.174611 -0.035560
## c_PDO     0.106493 7.38e-02 -5.21e-02  0.109664  0.243990
## c_Hrel   -0.277020 1.08e-01 -4.60e-01 -0.289274 -0.028317
## var_Qr    0.094688 3.54e-02  4.33e-02  0.088710  0.180815
## var_Rs    0.021849 2.24e-02  4.70e-05  0.015965  0.077751
```

Here are a series of figures summarizing model output.

### Ricker _a_

Here are histograms of the posterior samples for the mean Ricker $a$ (left) and $\alpha = \exp{\alpha}$ (right) parameters. Note that for plotting purposes only, the density in the largest bin for each parameter contains counts for all values greater or equal to that. Vertical arrows under the x-axes indicate the 2.5^th^, 50^th^, and 97.5^th^ percentiles.


```r
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
par(mfrow=c(1,2), mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
## Ricker a
R_a_est <- mod_fit$BUGSoutput$sims.list$mu_Rkr_a
R_a_est[R_a_est>3] <- 3
alphaCI <- quantile(R_a_est,c(0.025,0.5,0.975))
brks <- seq(floor(min(R_a_est)/0.1)*0.1,ceiling(max(R_a_est)/0.1)*0.1,0.1)
hist(R_a_est,freq=FALSE,xlab="",main="",breaks=brks,
     col=clr, border="blue3", ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/20
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(paste("Ricker ",italic(a))), 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)
## Ricker alpha
#R_alpha_est <- mod_fit$BUGSoutput$sims.list$alpha
R_alpha_est <- exp(R_a_est)
R_alpha_est[R_alpha_est>9] <- 9
alphaCI <- quantile(R_alpha_est,c(0.025,0.5,0.975))
hist(R_alpha_est,freq=FALSE,xlab="",main="",breaks=seq(0,ceiling(max(R_alpha_est)/0.3)*0.3,0.3),
     col=clr, border="blue3", ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/20
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
#mtext("Ricker exp(a)", 1, line=3, cex=1.2)
mtext(expression(paste("Ricker ",alpha," ",(e^italic(a)))), 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)
```

![](sssr_files/figure-html/plot_Ricker_a-1.png)<!-- -->

### Ricker _b_

Here is a histogram of the estimated strength of density dependence. Vertical arrows under the x-axis indicate the 2.5^th^, 50^th^, and 97.5^th^ percentiles.


```r
par(mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
RbDat <- mod_fit$BUGSoutput$sims.list$Rkr_b
RbDat <- RbDat*10^abs(floor(log(max(RbDat),10)))
ylM <- max(RbDat)
brks <- seq(0,ceiling(ylM),0.1)
betaCI <- quantile(RbDat,c(0.025,0.5,0.975))
hist(RbDat, freq=FALSE, breaks=brks, col=clr, border="blue3",
	 xlab="", xaxt="n", yaxt="n",
	 main="", ylab="Posterior density", cex.lab=1.2)
axis(1, at=seq(0,3))
aHt <- (par()$usr[4]-par()$usr[3])/20
arrows(betaCI,par()$usr[3]-0.005,betaCI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(paste("Ricker ",italic(b)," ",(10^{-4}),"")), 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)
```

![](sssr_files/figure-html/plot_Ricker_b-1.png)<!-- -->

### Covariate effects

Here are plots of the covariate time series (_z_-scores, A-C) and histograms of their effects on productivity (D-F).


```r
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
offSet <- 0.07
covars <- mod_fit$BUGSoutput$sims.matrix[,grep("c_",names(mod_fit$BUGSoutput$sims.list),value=TRUE)]
par(mfrow=c(ncol(covars),2), mai=c(0.4,0.2,0.1,0.1), omi=c(0.2,0.4,0,0))
ylN <- floor(min(covars)*10)/10
ylM <- ceiling(max(covars)*10)/10
brks <- seq(ylN,ylM,length.out=diff(c(ylN,ylM))*40+1)
cov_names <- c("Flow","PDO","H releases")
tSeries <- seq(yr_frst,length.out=n_yrs-age_min)
for(i in 1:ncol(covars)) {
	## plot covar ts
	plot(tSeries, dat_cvrs[seq(length(tSeries)),i], xlab="", ylab="",
		 main="", cex.lab=1.3, pch=16, col="blue3", type="o")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i])
	mtext(side=2, cov_names[i], line=3)
	if(i==ncol(covars)) { mtext(side=1,"Brood year", line=3) }
	## plot covar effect
	hist(covars[,grep(colnames(dat_cvrs)[i],colnames(covars))],
	     freq=FALSE,breaks=brks,col=clr,border="blue3",
	     xlab="", yaxt="n", main="", ylab="", cex.lab=1.2)
	c_CI <- quantile(covars[,grep(colnames(dat_cvrs)[i],colnames(covars))],c(0.025,0.5,0.975))
	aHt <- (par()$usr[4]-par()$usr[3])/20
	arrows(c_CI,par()$usr[3]-0.005,c_CI,par()$usr[3]-aHt,
	       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
	abline(v=0, lty="dashed")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i+ncol(covars)])
	if(i==ncol(covars)) { mtext(side=1,"Effect size", line=3) }
}
```

![](sssr_files/figure-html/plot_cov_effects-1.png)<!-- -->

### Spawners

Here is the estimate of the number of spawners over time. The black points are the data, the blue line is the median posterior estimate, and the shaded region is the 95% credible interval.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,sort)[,(1:n_yrs)]
pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI.vec] })
ypMin <- min(pDat)
ypMax <- max(pDat)
tSeries <- seq(yr_frst,length.out=n_yrs)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", xaxt="n", yaxt="n",
	 xlab="Year", ylab="Spawners", main="", cex.lab=1.2)
polygon(c(tSeries,rev(tSeries)),c(pDat[3,1:n_yrs],rev(pDat[1,1:n_yrs])), col=clr, border=NA)
lines(tSeries, pDat[2,1:n_yrs], col="blue3", lwd=2)
points(seq(yr_frst,length.out=n_yrs+n_fore), exp(ln_dat_esc), pch=16, cex=1)
axis(1,at=seq(1980,2015,5))
axis(2,at=c(3000,6000,12000))
```

![](sssr_files/figure-html/plot_spawners-1.png)<!-- -->

### Total run size

Here is our estimate of the total run size (i.e., catch + escapement) over time, which includes a forecast for 2016. The black points are the data, the blue line is the median posterior estimate, and the shaded region is the 95% credible interval.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI.vec] })
pDat <- pDat + matrix(dat_harv,length(CI.vec),n_yrs+n_fore,byrow=TRUE)
tSeries <- seq(yr_frst,length.out=n_yrs+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", xaxt="n", yaxt="n",
	 xlab="Year", ylab="Catch + escapement", main="", cex.lab=1.2)
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
points(tSeries, exp(ln_dat_esc)+dat_harv, pch=16, cex=1)
axis(1,at=seq(1980,2015,5))
axis(2,at=c(4000,8000,16000))
```

![](sssr_files/figure-html/plot_run_size-1.png)<!-- -->

### Recruits by age class

Here are the estimated number of recruits by brood year and age. Note that the uncertainty increases in more recent years as fewer complete age classes have been observed.


```r
CI.vec <- c(0.05,0.5,0.95)
par(mfrow=c(A,1), mai=c(0.1,0.1,0.05,0.1), omi=c(0.5,0.5,0.1,0))
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
nRec <- n_yrs-age_min
tSeries <- seq(yr_frst,length.out=nRec+n_fore)
pltTT <- seq(min(round(tSeries/5,0)*5),max(round(tSeries/5,0)*5),5)
for(i in rev(1:A)) {
	pDat <- apply(mod_fit$BUGSoutput$sims.list$Rec[,,i],2,sort)
	pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI.vec] })/100
	dd <- ifelse(max(pDat)<20,1,10)
	ypMax <- real2prec(max(pDat),prec=dd)
	while(ypMax %% 3 != 0) { ypMax <- ypMax + dd }
	plot(tSeries,pDat[3,], xlim=c(yr_frst+1,yr_last-n_fore-2), ylim=c(0.001,ypMax),
	     type="n", xaxt="n", yaxt="n", xlab="", ylab="", main="", las=1)
	polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
	lines(tSeries, pDat[2,], col="blue3", lwd=2)
	aHt <- (par()$usr[4]-par()$usr[3])/7
	ttl <- paste("Age-",i+age_min-1,sep="")
	text(tSeries[1]-0, par()$usr[4]-aHt, ttl, pos=4, cex=0.9)
	axis(2,seq(0,ypMax,length.out=4),las=1,cex=0.9)
	if(i!=1) {axis(1,at=pltTT,labels=FALSE)} else {axis(1,at=pltTT)}
}
mtext("Recruits (100s)", 2, line=2, outer=TRUE, cex=1.2)
mtext("Year", 1, line=2.5, outer=TRUE, cex=1.2)
```

![](sssr_files/figure-html/plot_recruits_by_age-1.png)<!-- -->

### Total recruits

Here are the estimated total number of recruits by brood year.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod_fit$BUGSoutput$sims.list$Rec,c(1,2),sum)
pDat <- apply(apply(pDat,2,sort),2,function(x) { x[mcmc_samp*CI.vec] })
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", yaxt="n",
	 xlab="Brood year", ylab="Recruits", main="", cex.lab=1.2)
axis(2,at=c(3000,9000,27000))
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

![](sssr_files/figure-html/plot_total_recruits-1.png)<!-- -->

### Recruits per spawner

Here is the time series of estimated recruits-per-spawner.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod_fit$BUGSoutput$sims.list$ln_RS,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI.vec] })
pDat[2,] <- apply(mod_fit$BUGSoutput$sims.list$ln_RS,2,median)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="ln(R/S)", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

![](sssr_files/figure-html/plot_R_per_S-1.png)<!-- -->

### Innovations

Here is the time series of the so-called "innovations", which are the residuals from the process model. They give some indications of population productivity after accounting for the effects of density dependence.


```r
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod_fit$BUGSoutput$sims.list$res_ln_Rec,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI.vec] })
pDat[2,] <- apply(mod_fit$BUGSoutput$sims.list$res_ln_Rec,2,median)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="Innovations", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)
```

![](sssr_files/figure-html/plot_innovations-1.png)<!-- -->

### Age composition

Here are time series of the estimated proportions of each age class.


```r
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.2,0.2))
CI.vec <- c(0.025,0.5,0.975)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
clr <- rgb(0, 0, 255, alpha = 40, maxColorValue = 255)
ageComp <- t(apply(apply(mod_fit$BUGSoutput$sims.list$p_vec,c(3,2),mean),2,cumsum))
plot(tSeries, rep(1,nRec+n_fore), ylab="Proportion", xlab="Brood year", ylim=c(0,1), las=1,
     xaxs="i", yaxs="i", type="n", lty="solid", col="blue3", cex.lab=1.2)
for(i in c(1,2,3,4,6)) {
	polygon(c(tSeries,rev(tSeries)),c(ageComp[,i],rep(0,nRec+n_fore)), col=clr, border=NA)
	}
lbl <- apply(cbind(c(0,ageComp[nRec+n_fore,-A]),ageComp[nRec+n_fore,]),1,mean)
text(par()$usr[2],par()$usr[4]*1.05,"Age", xpd=NA, pos=4, offset=0.05, col="black", cex=0.8)
text(par()$usr[2],lbl[1:4],seq(3,6), xpd=NA, pos=4, col="black", cex=0.7)
text(par()$usr[2],lbl[5],"7&8", xpd=NA, pos=4, offset=0.15, col="black", cex=0.7)
```

![](sssr_files/figure-html/plot_age_comp-1.png)<!-- -->

### Spawner-recruit relationship

Here is the relationship between spawner and subsequent recruits. Gray lines show 100 possible relationships based on random draws from the posterior distributions for Ricker's _a_ and _b_ parameters.


```r
CI.vec <- c(0.025,0.5,0.975)
MC <- 100
idx <- sample(seq(mcmc_samp),MC)
sDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,sort)
sDat <- apply(sDat,2,function(x) { x[mcmc_samp*CI.vec] })
sDat[2,] <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,median)
sDat <- sDat[,1:(n_yrs-age_min+n_fore)]
rDat <- apply(mod_fit$BUGSoutput$sims.list$tot_ln_Rec,2,sort)
rDat <- apply(rDat,2,function(x) { x[mcmc_samp*CI.vec] })
rDat <- exp(rDat)
dd <- 3000
yM <- real2prec(max(rDat),"ceiling",dd)
yM <- 30000
xM <- real2prec(max(sDat),"ceiling",dd)
aa <- exp(matrix(mod_fit$BUGSoutput$sims.array[,,"mu_Rkr_a"],ncol=1))
bb <- matrix(mod_fit$BUGSoutput$sims.array[,,"Rkr_b"], ncol=1)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.1,0.2))
plot(sDat[2,],rDat[2,], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3", type="n",
	 xaxs="i", yaxs="i", ylab="Recruits (1000s)", xlab="Spawners (1000s)", cex.lab=1.2,
	 xaxt="n", yaxt="n")
axis(1, at=seq(0,xM,dd*2), labels=seq(0,xM,dd*2)/1000)
axis(2, at=seq(0,yM,dd*2), labels=seq(0,yM,dd*2)/1000)
for(i in 1:MC) { lines(aa[idx[i]]*seq(xM)*exp(-bb[idx[i]]*seq(xM)), col="darkgray") }
abline(a=0,b=1,lty="dashed")
nCB <- n_yrs-age_max
points(sDat[2,1:nCB],rDat[2,1:nCB], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3")
segments(sDat[2,1:nCB],rDat[1,1:nCB],sDat[2,1:nCB],rDat[3,1:nCB], col="blue3")
segments(sDat[1,1:nCB],rDat[2,1:nCB],sDat[3,1:nCB],rDat[2,1:nCB], col="blue3")
nTB <- dim(sDat)[2]
clr <- rgb(100, 0, 200, alpha = seq(200,100,length.out=age_max-age_min+n_fore), maxColorValue = 255)
segments(sDat[2,(nCB+1):nTB],rDat[1,(nCB+1):nTB],sDat[2,(nCB+1):nTB],rDat[3,(nCB+1):nTB], col=clr)
segments(sDat[1,(nCB+1):nTB],rDat[2,(nCB+1):nTB],sDat[3,(nCB+1):nTB],rDat[2,(nCB+1):nTB], col=clr)
points(sDat[2,(nCB+1):nTB],rDat[2,(nCB+1):nTB],
       xlim=c(0,xM), ylim=c(0,yM), pch=16, col=clr)
```

![](sssr_files/figure-html/plot_SR-1.png)<!-- -->

### Management reference points

Here are a number of management reference points.


```r
# abbreviations for ref points
refNames <- c("MSY","Smsy","Umsy")
# proportions of MSY to consider
yieldProps <- c(0.7,0.8,0.9)
propNames <- paste("OYP",yieldProps,sep="")
lnA <- matrix(mod_fit$BUGSoutput$sims.array[,,"mu_Rkr_a"],ncol=1)
bb <- matrix(mod_fit$BUGSoutput$sims.array[,,"Rkr_b"], ncol=1)
mcmc <- length(lnA)
# empty matrix for ref pts
ref.pts <- matrix(NA,mcmc,length(refNames))
colnames(ref.pts) <- refNames
# spawner series for optimal yield profile
SS <- seq(100,1e4,100)
# empty matrix for optimal yield profiles
OYP <- matrix(0,length(SS),length(yieldProps))
for(i in 1:mcmc) {
	# spawners at MSY
	ref.pts[i,"Smsy"] <- (1 - lambert_W0(exp(1-lnA[i]))) / bb[i]
	# MSY
	ref.pts[i,"MSY"] <- ref.pts[i,"Smsy"]*((exp(lnA[i]-bb[i]*ref.pts[i,"Smsy"])) - 1)
	# harvest rate at MSY
	ref.pts[i,"Umsy"] <- (1 - lambert_W0(exp(1-lnA[i])))
	# yield over varying S
	yield <- SS*(exp(lnA[i]-bb[i]*SS) - 1)
	for(j in 1:length(yieldProps)) {
		OYP[,j] <- OYP[,j] + 1*(yield > yieldProps[j]*ref.pts[i,"MSY"])
	}
}
OYP <- OYP/mcmc
```

Here are histograms of the estimated (A) maximum sustainable yield, (B) spawners at MSY, (C) harvest rate at MSY, and (D) the optimal yield profiles based on the alternative states of nature.


```r
ttl <- expression(bold(MSY), bold(S[MSY]), bold(U[MSY]))
x_lbl <- c("Escapement","Escapement","Rate")
plt_dat <- ref.pts
plt_dat[plt_dat>1e4] <- 1e4
offSet <- 0.07
par(mfcol=c(2,2), mai=c(0.8,0.8,0.4,0.1), omi=c(0,0,0,0.2))
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
for(i in 1:3) {
	hist(plt_dat[,refNames[i]],xlim=c(0,max(plt_dat[,refNames[i]])),
	     freq=FALSE,breaks=50,col=clr,border="blue3",
	     xlab=x_lbl[i], yaxt="n", main=ttl[i], ylab="", cex.lab=1.2)
	mtext("Posterior probability", 2, cex=1)
	c_CI <- quantile(ref.pts[,refNames[i]],c(0.025,0.5,0.975))
	aHt <- (par()$usr[4]-par()$usr[3])/20
	arrows(c_CI,par()$usr[3],c_CI,par()$usr[3]-aHt,
	       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i])
}
x_lp <- yieldProps
for(i in 1:length(x_lp)) {
	x_lp[i] <- SS[max(which(OYP[,i] == max(OYP[,i]) | abs(OYP[,i] - (yieldProps[i]-0.3)) <= 0.05))]
}
matplot(SS, OYP, type="l", lty="solid", las=1, col=c("slateblue","blue","darkblue"), lwd=2,
		xlab="Escapement", ylab="Prob. of meeting XX% of MSY", cex.lab=1.2, bty="n",
		main="Optimal yield profiles")
points(x=x_lp, y=yieldProps-0.3, pch=21, cex=3.5, col="white", bg="white")
text(x=x_lp, y=yieldProps-0.3, paste0(yieldProps*100,"%"),
	 col=c("slateblue","blue","darkblue"), cex=0.7)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"D")
```

![](sssr_files/figure-html/plot_OYP-1.png)<!-- -->

***

## References
