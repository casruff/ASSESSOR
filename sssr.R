## ----load_pkgs, message=FALSE, warning=FALSE-----------------------------
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

## ----real2prec-----------------------------------------------------------
real2prec <- function(x,map="round",prec=1) {
  if(prec==0) stop("\"prec\" cannot be 0")
    do.call(map,list(x/prec))*prec
}

## ----get_user_inputs-----------------------------------------------------
## file with escapement data
## [n_yrs x 2] matrix of obs counts; 1st col is calendar yr
fn_esc <- "SkagitSthdEsc.csv"

## file with age comp data
## [n_yrs x (1+A)]; 1st col is calendar yr
fn_age <- "SkagitSthdAge.csv"
## min & max ages
age_min <- 3
age_max <- 8
## years, if any, of age-comp to skip; see below
age_skip <- 2

## file with catch data
## [n_yrs x 2] matrix of obs catch; 1st col is calendar yr
fn_harv <- "SkagitSthdCatch.csv"

## file with covariate data
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

## ----get_escapement_data-------------------------------------------------
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

## ----get_age_data--------------------------------------------------------
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

## ----get_harvest---------------------------------------------------------
## harvest
dat_harv <- read.csv(text = getURL(paste0(ex_url,fn_harv)))
## drop year col & first age_max rows
dat_harv <- c(dat_harv[,-1],rep(0,n_fore))

## ----get_covariates------------------------------------------------------
## covariate(s)
dat_cvrs <- read.csv(text = getURL(paste0(ex_url,fn_cvrs)))
## drop year col
dat_cvrs <- dat_cvrs[,-1]
## transform the covariates to z-scores
dat_cvrs <- scale(dat_cvrs)

## ----SSSR_in_JAGS--------------------------------------------------------
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
	tau_Qr <- pow(sd_Qr,-2);
	sd_Qr ~ dunif(0.001,20);
	var_Qr <- pow(sd_Qr,2)
	
	# innov in yr 1
	innov_1 ~ dnorm(0,tau_Qr*(1-phi*phi));
	
	# Rs = variance for Sp obs model
	# diffuse gamma prior on precision
	# tau_Rs ~ dgamma(0.001,0.001);
	# var.R <- pow(tau_Rs,-1);
	# unif prior on SD
	tau_Rs <- pow(sd_Rs,-2);
	sd_Rs ~ dunif(0.001,20);
	
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
	for(t in 1:(n_yrs-age_min+n_fore)) { pi[t,1:A] ~ ddirch(muHD*piHD) }
	
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
		Rec[1,a] <- max(1,tot_Rec[1] * pi[1,a]);
		}
	
	# brood years 2:(n_yrs-age_min)
	for(t in 2:(n_yrs-age_min+n_fore)) {
		# predicted recruits in BY t
		ln_Rkr_a[t] <- Rkr_a + c_PDO*dat_cvrs[t,1] + c_Flow*dat_cvrs[t,2] + c_Hrel*dat_cvrs[t,3];
		E_ln_Rec[t] <- ln_Sp[t] + ln_Rkr_a[t] - Rkr_b*Sp[t];
		tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t] + phi*res_ln_Rec[t-1],tau_Qr);
		res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
		# median of total recruits
		tot_Rec[t] <- exp(tot_ln_Rec[t]);
		# R/S
		ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
		# brood-yr recruits by age
		for(a in 1:A) {
			Rec[t,a] <- max(1,tot_Rec[t] * pi[t,a]);
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

## ----JAGS_IO, message=FALSE, warning=FALSE, cache=TRUE-------------------
## data to pass to JAGS
data.JAGS <- c("dat_age","ln_dat_esc","dat_harv","dat_cvrs",
               "n_yrs","A","age_min","age_max","age_skip","n_fore")

## 2. model params/states for JAGS to return
par.JAGS <- c("alpha","mu_Rkr_a","Rkr_b","Sp","Rec","tot_ln_Rec","ln_RS",
              "c_Flow","c_PDO","c_Hrel",
              "var_Qr","pi","res_ln_Rec")

## 3. MCMC control params
# MCMC parameters
mcmc.chains <- 4
mcmc.length <- 10e5
mcmc.burn <- 5e5
mcmc.thin <- 1000
# total number of MCMC samples
mcmc.samp <- (mcmc.length-mcmc.burn)*mcmc.chains/mcmc.thin

## function to create JAGS inits
init.vals <- function() {
	list(Rkr_a=1, c_Flow=0.1, c_PDO=-0.1, c_Hrel=-0.2,
	     Rkr_b=1/exp(mean(ln_dat_esc, na.rm=TRUE)),
	     piHD=1, muHD=rep(1,A),
	     pi=matrix(c(0.01,0.3,0.48,0.15,0.05,0.01),n_yrs-age_min+n_fore,A,byrow=TRUE),
	     Rec_mu=log(1000),
	     Rec_sig=0.1,
	     tot_ln_Rec=rep(log(1000),n_yrs-age_min+n_fore),
	     innov_1=0,
	     phi=0.5)
	}

mod.JAGS <- list(data=data.JAGS,
				 inits=init.vals,
				 parameters.to.save=par.JAGS,
				 model.file=fn_jags,
				 n.chains=as.integer(mcmc.chains),
				 n.iter=as.integer(mcmc.length),
				 n.burnin=as.integer(mcmc.burn),
				 n.thin=as.integer(mcmc.thin),
				 DIC=TRUE)

## start timer
timer_start <- proc.time()

## fit the model in JAGS & store results
mod.fit <- do.call(jags.parallel, mod.JAGS)

## stop timer
(run_time_in_min <- round(((proc.time()-timer_start)/60)["elapsed"], 1))

## ----model_diagnostics, cache=TRUE, eval=TRUE----------------------------
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

## ----plot_Ricker_a, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
par(mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
RaDat <- mod.fit$BUGSoutput$sims.list$alpha
RaDat[RaDat>10] <- 10.2
alphaCI <- quantile(RaDat,c(0.025,0.5,0.975))
hist(RaDat,freq=FALSE,xlab="",main="",breaks=seq(0,max(RaDat),0.2),col=clr,border="blue3",
	 ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/10
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt,code=1,length=0.05,xpd=NA,col="blue3")
mtext("Ricker a", 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)

## ----plot_Ricker_b, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
par(mai=c(0.8,0.4,0.3,0.1), omi=c(0,0,0,0.2))
RbDat <- mod.fit$BUGSoutput$sims.list$Rkr_b
RbDat <- RbDat*10^abs(floor(log(max(RbDat),10)))
ylM <- max(RbDat)
brks <- seq(0,ceiling(ylM),0.1)
betaCI <- quantile(RbDat,c(0.025,0.5,0.975))
hist(RbDat, freq=FALSE, breaks=brks, col=clr, border="blue3",
	 xlab="", xaxt="n", yaxt="n",
	 main="", ylab="Posterior density", cex.lab=1.2)
axis(1, at=seq(0,3))
aHt <- (par()$usr[4]-par()$usr[3])/10
arrows(betaCI,par()$usr[3]-0.005,betaCI,par()$usr[3]-aHt,code=1,length=0.05,xpd=NA,col="blue3")
mtext(expression(paste("Ricker b ",(10^{-4}),"")), 1, line=3, cex=1.2)
mtext("Posterior probability", 2, cex=1.2)

## ----plot_cov_effects, fig.width=6.5, fig.height=8.5, fig.pos="placeHere", warnings=FALSE, messages=FALSE, eval=TRUE----
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
offSet <- 0.07
covars <- mod.fit$BUGSoutput$sims.matrix[,grep("c_",names(mod.fit$BUGSoutput$sims.list),
                                               value=TRUE)]
par(mfrow=c(ncol(covars),2), mai=c(0.4,0.2,0.1,0.1), omi=c(0.2,0.4,0,0))
ylN <- floor(min(covars)*10)/10
ylM <- ceiling(max(covars)*10)/10
brks <- seq(ylN,ylM,length.out=diff(c(ylN,ylM))*40+1)
cov_names <- c("Flow","H releases","PDO")
tSeries <- seq(yr_frst,length.out=n_yrs-age_min)
for(i in 1:ncol(covars)) {
	# plot covar ts
	plot(tSeries, dat_cvrs[seq(length(tSeries)),i], xlab="", ylab="",
		 main="", cex.lab=1.3, pch=16, col="blue3", type="o")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
		 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i])
	mtext(side=2, cov_names[i], line=3)
	if(i==ncol(covars)) {
		mtext(side=1,"Brood year", line=3)
	}
	# plot covar effect
	hist(covars[,grep(colnames(dat_cvrs)[i],colnames(covars))],
	     freq=FALSE,breaks=brks,col=clr,border="blue3",
		 xlab="", yaxt="n",
		 main="", ylab="", cex.lab=1.2)
	abline(v=0, lty="dashed")
	text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
		 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),LETTERS[i+ncol(covars)])
	if(i==ncol(covars)) {
		mtext(side=1,"Effect size", line=3)
	}
}

## ----plot_spawners, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)[,(1:n_yrs)]
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
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

## ----plot_run_size, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
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

## ----plot_recruits_by_age, fig.width=6, fig.height=6, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.05,0.5,0.95)
par(mfrow=c(A,1), mai=c(0.1,0.1,0.0,0.1), omi=c(0.5,0.5,0.1,0))
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
nRec <- n_yrs-age_min
tSeries <- seq(yr_frst,length.out=nRec+n_fore)
pltTT <- seq(min(round(tSeries/5,0)*5),max(round(tSeries/5,0)*5),5)
for(i in rev(1:A)) {
	pDat <- apply(mod.fit$BUGSoutput$sims.list$Rec[,,i],2,sort)
	pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })/100
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

## ----plot_total_recruits, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$Rec,c(1,2),sum)
pDat <- apply(apply(pDat,2,sort),2,function(x) { x[mcmc.samp*CI.vec] })
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", yaxt="n",
	 xlab="Brood year", ylab="Recruits", main="", cex.lab=1.2)
axis(2,at=c(3000,9000,27000))
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)

## ----plot_R_per_S, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$ln_RS,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
pDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$ln_RS,2,median)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="ln(R/S)", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)

## ----plot_innovations, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
pDat <- apply(mod.fit$BUGSoutput$sims.list$res_ln_Rec,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc.samp*CI.vec] })
pDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$res_ln_Rec,2,median)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(tSeries,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="Innovations", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(tSeries,rev(tSeries)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(tSeries, pDat[2,], col="blue3", lwd=2)

## ----plot_age_comp, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.2,0.2))
CI.vec <- c(0.025,0.5,0.975)
tSeries <- seq(yr_frst,length.out=n_yrs-age_min+n_fore)
clr <- rgb(0, 0, 255, alpha = 40, maxColorValue = 255)
ageComp <- t(apply(apply(mod.fit$BUGSoutput$sims.list$pi,c(3,2),mean),2,cumsum))
plot(tSeries, rep(1,nRec+n_fore), ylab="Proportion", xlab="Brood year", ylim=c(0,1), las=1,
		xaxs="i", yaxs="i", type="n", lty="solid", col="blue3", cex.lab=1.2)
for(i in c(1,2,3,4,6)) {
	polygon(c(tSeries,rev(tSeries)),c(ageComp[,i],rep(0,nRec+n_fore)), col=clr, border=NA)
	}
lbl <- apply(cbind(c(0,ageComp[nRec+n_fore,-A]),ageComp[nRec+n_fore,]),1,mean)
text(par()$usr[2],par()$usr[4]*1.05,"Age", xpd=NA, pos=4, offset=0.05, col="black", cex=0.8)
text(par()$usr[2],lbl[1:4],seq(3,6), xpd=NA, pos=4, col="black", cex=0.7)
text(par()$usr[2],lbl[5],"7&8", xpd=NA, pos=4, offset=0.15, col="black", cex=0.7)

## ----plot_SR, fig.width=6, fig.height=5, fig.pos="placeHere", eval=TRUE----
CI.vec <- c(0.025,0.5,0.975)
MC <- 100
idx <- sample(seq(mcmc.samp),MC)
sDat <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,sort)
sDat <- apply(sDat,2,function(x) { x[mcmc.samp*CI.vec] })
sDat[2,] <- apply(mod.fit$BUGSoutput$sims.list$Sp,2,median)
sDat <- sDat[,1:(n_yrs-age_min+n_fore)]
rDat <- apply(mod.fit$BUGSoutput$sims.list$tot_ln_Rec,2,sort)
rDat <- apply(rDat,2,function(x) { x[mcmc.samp*CI.vec] })
rDat <- exp(rDat)
dd <- 3000
yM <- real2prec(max(rDat),"ceiling",dd)
yM <- 30000
xM <- real2prec(max(sDat),"ceiling",dd)
aa <- matrix(mod.fit$BUGSoutput$sims.array[,,"alpha"],ncol=1)
bb <- matrix(mod.fit$BUGSoutput$sims.array[,,"Rkr_b"], ncol=1)
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

## ----ref_pts, eval=TRUE--------------------------------------------------
# abbreviations for ref points
refNames <- c("Smsy","Umsy","MSY","Umax","MSR","Smsr","Seq")
# proportions of MSY to consider
yieldProps <- c(0.7,0.8,0.9)
propNames <- paste("OYP",yieldProps,sep="")
lnA <- matrix(mod.fit$BUGSoutput$sims.array[,,"mu_Rkr_a"],ncol=1)
bb <- matrix(mod.fit$BUGSoutput$sims.array[,,"Rkr_b"], ncol=1)
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

## ----plot_OYP, fig.width=6, fig.height=4, fig.pos="placeHere", eval=TRUE----
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0,0,0.2))
matplot(SS, OYP, type="l", lty="solid", las=1, col=c("slateblue","blue","darkblue"),
		xlab="Escapement", ylab="Probability", lwd=2)
points(x=c(4900,5700,6450), y=c(0.6,0.5,0.4), pch=21, cex=3.5, col="white", bg="white")
text(x=c(4900,5700,6450), y=c(0.6,0.5,0.4), c("90%","80%","70%"),
	 col=c("darkblue","blue","slateblue"), cex=0.7)

