
# ====================== Simulation testing of Depletion based MPs ======================
# 
# Harford and Carruthers 2017
#
# Jan 20, 2017
#
# A modified version of runMSE is used. The following changes were made relative to the runMSE base from DLMtool v3.2
# 1) DLM_data@Dt is depletion of vulnerable biomass
# 2) DLM_data@Dep becomes 2D array ncol=2, nrow=nsim, (does not recursively update through time): [,1] depletion of vulnerable biomass in nyear, [,2] BMSY_B0 in vulnerable biomass
# 3) OM object @OM$Depletion is the depletion (SSB/SSB0) in the final historical year, "true" value. 
#
# Depletion used in DCAC, SP, and DBSRA is based on vulnerable biomass in the MPs coded below.
#
# ============================================================================================================================

rm(list=ls())
graphics.off()
library(DLMtool)

for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]])
sfInit(parallel=T,cpus=detectCores())


runMSE <- function(OM="1", MPs=NA, nsim=48, proyears=28, interval=4, pstar=0.5,
                   maxF=0.8, timelimit=1, reps=1, custompars=0, CheckMPs=TRUE){ 
  print("Loading operating model")
  
  flush.console()
  if(class(OM)!="OM")stop("You must specify an operating model")
  #if(!sfIsRunning())stop("You must initialize snowfall functions sfInit() see ??DLMtool")
  
  nyears <- OM@nyears  # number of  historical years
  maxage <- OM@maxage  # maximum age (no plus group)
  
  calcMax <- -log(0.01)/(min(OM@M)) # Age at which 1% of cohort survives
  maxage <- max(maxage, calcMax) # If maximum age is lower, increase it to calcMax
  OM@maxage <- maxage
  dep <- runif(nsim,OM@D[1],OM@D[2])  # sample from the range of user-specified depletion (Bcurrent/B0)
  Esd <- runif(nsim,OM@Fsd[1],OM@Fsd[2]) # interannual variability in fishing effort (log normal sd)
  
  Deriv <- getEffhist(Esd, nyears, EffYears=OM@EffYears, EffLower=OM@EffLower, EffUpper=OM@EffUpper) # Historical fishing effort
  Find <- Deriv[[1]] # Calculate fishing effort rate
  dFfinal <- Deriv[[2]] # Final gradient in fishing effort yr-1 
  
  # dep[order(dFfinal)]<-dep[order(dep,decreasing=T)] # robustifies 
  
  # matplot(t(Find), type="l")
  # plot(dep, dFfinal) 
  
  # Sample operating model parameters ===========================================================
  procsd<-runif(nsim,OM@Perr[1],OM@Perr[2])         # Process error standard deviation
  AC<-runif(nsim,OM@AC[1],OM@AC[2])    # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
  
  M<-runif(nsim,OM@M[1],OM@M[2]) # natural mortality rate 	
  Msd<-runif(nsim,OM@Msd[1],OM@Msd[2]) # sample inter annual variability in M from specified range
  Mgrad<-runif(nsim,OM@Mgrad[1],OM@Mgrad[2]) # sample gradient in M (M y-1)
  hs<-runif(nsim,OM@h[1],OM@h[2]) # sample of recruitment compensation (steepness - fraction of unfished recruitment at 20% of unfished biomass)
  Linf<-runif(nsim,OM@Linf[1],OM@Linf[2]) # sample of asymptotic length
  Linfsd<-runif(nsim,OM@Linfsd[1],OM@Linfsd[2]) # sample of interannual variability in Linf
  Linfgrad<-runif(nsim,OM@Linfgrad[1],OM@Linfgrad[2]) # sample of gradient in Linf (Linf y-1)
  recgrad<-runif(nsim,OM@recgrad[1],OM@recgrad[2]) # gradient in recent recruitment
  K<-runif(nsim,OM@K[1],OM@K[2])  # now predicted by a log-linear model
  Ksd<-runif(nsim,OM@Ksd[1],OM@Ksd[2])#runif(nsim,OM@Ksd[1],OM@Ksd[2])# sd is already added in the linear model prediction
  Kgrad<-runif(nsim,OM@Kgrad[1],OM@Kgrad[2]) # gradient in Von-B K parameter (K y-1)
  t0<-runif(nsim,OM@t0[1],OM@t0[2]) # a sample of theoretical age at length zero
  lenM <- array(runif(nsim*50,OM@L50[1],OM@L50[2]),c(nsim,50)) # length at 50% maturity
  lenM[lenM/Linf>0.8]<-NA
  lenM<-apply(lenM,1,function(x)x[!is.na(x)][1])
  len95 <- array(lenM + runif(nsim*50,OM@L50_95[1],OM@L50_95[2]),c(nsim,50)) # length at 95% maturity
  len95[len95/Linf>0.9]<-NA
  len95<-apply(len95,1,function(x)x[!is.na(x)][1])
  
  Spat_targ<-runif(nsim,OM@Spat_targ[1],OM@Spat_targ[2]) # spatial targetting Ba^targetting param 
  Frac_area_1<-runif(nsim,OM@Frac_area_1[1],OM@Frac_area_1[2]) # sampled fraction of unfished biomass in area 1 (its a two area model by default)
  Prob_staying<-runif(nsim,OM@Prob_staying[1],OM@Prob_staying[2]) # sampled probability of individuals staying in area 1 among years
  Size_area_1<-runif(nsim,OM@Size_area_1[1],OM@Size_area_1[2]) # currently redundant parameter for the habitat area size of area 1
  
  # Sample observation error model parameters ===============================================================
  Csd<-runif(nsim,OM@Cobs[1],OM@Cobs[2])                                    # Sampled catch observation error (lognormal sd)
  Cbias<-rlnorm(nsim,mconv(1,OM@Cbiascv),sdconv(1,OM@Cbiascv))              # Sampled catch bias (log normal sd)
  CAA_nsamp<-ceiling(runif(nsim,OM@CAA_nsamp[1],OM@CAA_nsamp[2]))                     # Number of catch-at-age observations
  CAA_ESS<-ceiling(runif(nsim,OM@CAA_ESS[1],OM@CAA_ESS[2]))                                          # Effective sample size
  CAL_nsamp<-runif(nsim,OM@CAL_nsamp[1],OM@CAL_nsamp[2])                              # Observation error standard deviation for single catch at age by area
  CAL_ESS<-ceiling(runif(nsim,OM@CAL_ESS[1],OM@CAL_ESS[2]))                                          # Effective sample size
  CALcv<-runif(nsim,OM@CALcv[1],OM@CALcv[2])                                    # Observation error standard deviation for single catch at age by area
  betas<-exp(runif(nsim,log(OM@beta[1]),log(OM@beta[2]))) # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  Isd<-runif(nsim,OM@Iobs[1],OM@Iobs[2])                    # Abundance index observation error (log normal sd)
  Derr<-runif(nsim,OM@Dcv[1],OM@Dcv[2])
  Dbias<-rlnorm(nsim,mconv(1,OM@Dbiascv),sdconv(1,OM@Dbiascv)) # sample of depletion bias
  Mbias<-rlnorm(nsim,mconv(1,OM@Mcv),sdconv(1,OM@Mcv))         # sample of M bias
  FMSY_Mbias<-rlnorm(nsim,mconv(1,OM@FMSY_Mcv),sdconv(1,OM@FMSY_Mcv)) # sample of FMSY/M bias
  ntest<-20                               # number of trials
  
  lenMbias<-rlnorm(nsim,mconv(1,OM@LenMcv),sdconv(1,OM@LenMcv))      # sample of length at maturity bias - assume same error as age based maturity
  LFCbias<-rlnorm(nsim,mconv(1,OM@LFCcv),sdconv(1,OM@LFCcv))        # sample of length at first capture bias
  LFSbias<-rlnorm(nsim,mconv(1,OM@LFScv),sdconv(1,OM@LFScv))        # sample of length at full selection bias
  Aerr<-runif(nsim,OM@Btcv[1],OM@Btcv[2])
  Abias<-exp(runif(nsim,log(OM@Btbias[1]),log(OM@Btbias[2])))#rlnorm(nsim,mconv(1,OM@Btbiascv),sdconv(1,OM@Btbiascv))    # smaple of current abundance bias
  Kbias<-rlnorm(nsim,mconv(1,OM@Kcv),sdconv(1,OM@Kcv))              # sample of von B. K parameter bias
  t0bias<-rlnorm(nsim,mconv(1,OM@t0cv),sdconv(1,OM@t0cv))           # sample of von B. t0 parameter bias
  Linfbias<-rlnorm(nsim,mconv(1,OM@Linfcv),sdconv(1,OM@Linfcv))     # sample of von B. maximum length bias
  Irefbias<-rlnorm(nsim,mconv(1,OM@Irefcv),sdconv(1,OM@Irefcv))     # sample of bias in reference (target) abundance index
  Crefbias<-rlnorm(nsim,mconv(1,OM@Crefcv),sdconv(1,OM@Crefcv))     # sample of bias in reference (target) catch index
  Brefbias<-rlnorm(nsim,mconv(1,OM@Brefcv),sdconv(1,OM@Brefcv))     # sample of bias in reference (target) biomass index
  Recsd<-runif(nsim,OM@Reccv[1],OM@Reccv[2])                        # Recruitment deviation 
  
  # Sample fishing efficiency parameters =======================================================
  qinc<-runif(nsim,OM@qinc[1],OM@qinc[2])  
  qcv<-runif(nsim,OM@qcv[1],OM@qcv[2])                 # interannual variability in catchability
  
  # dat<-as.data.frame(cbind(procsd,AC,M,Msd,Mgrad,hs,Linf,Linfsd,Linfgrad,recgrad,K,Ksd,Kgrad,t0,lenM,len95,L5,LFS,
  #           Vmaxlen,Spat_targ,Frac_area_1,Prob_staying,Size_area_1,Csd,Cbias,CAA_nsamp,CAA_ESS,CALcv,betas,
  #           Isd,Derr,Dbias,Mbias,FMSY_Mbias,lenMbias,LFCbias,LFSbias,Aerr,Abias,Kbias,
  #           t0bias,Linfbias,Irefbias,Crefbias,Brefbias,Recsd,qinc,qcv))
  
  #save(dat,file="F:/DLM/Operating models/Other/custompars")
  
  # Sample custom parameters ===================================================================
  if(sum(custompars)!=0){
    if(nrow(custompars)<nsim){
      ind<-sample(nrow(custompars),nsim,replace=T)
    }else{
      ind<-sample(nrow(custompars),nsim,replace=F)
    }
    for(i in 1:ncol(custompars))assign(names(custompars)[i],custompars[ind,i])
  }
  
  # Recruitment Deviations 
  procmu <- -0.5*(procsd)^2 # adjusted log normal mean
  Perr<-array(rnorm((nyears+proyears)*nsim,
                    rep(procmu,nyears+proyears),
                    rep(procsd,nyears+proyears)),
              c(nsim,nyears+proyears))
  for(y in 2:(nyears+proyears))Perr[,y]<-AC*Perr[,y-1]+Perr[,y]*(1-AC*AC)^0.5#2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
  Perr <-exp(Perr) # normal space (mean 1 on average)
  
  # Add cycle (phase shift) to recruitment deviations - if specified 
  if (is.finite(OM@Period[1]) & is.finite(OM@Amplitude[1])) {
    Shape <- "sin" # default sine wave - alternative - 'shift' for step changes
    recMulti <- sapply(1:nsim, SetRecruitCycle, Period=OM@Period, Amplitude=OM@Amplitude, 
                       TotYears=nyears+proyears, Shape=Shape)
    Perr <- Perr * t(recMulti) # Add cyclic pattern to recruitment
    print("Adding cyclic recruitment pattern")
    flush.console()
  }
  
  cumlRecDev <- apply(Perr[,1:nyears], 1, prod)
  dep[order(cumlRecDev)] <- dep[order(dep,decreasing=F)] # robustifies 
  
  R0 <- OM@R0  # Initial recruitment
  
  Marray<-gettempvar(M,Msd,Mgrad,nyears+proyears,nsim) # M by sim and year according to gradient and inter annual variability
  SRrel<-rep(OM@SRrel,nsim) # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  Linfarray<-gettempvar(Linf,Linfsd,Linfgrad,nyears+proyears,nsim) # Linf array
  Karray<-gettempvar(K,Ksd,Kgrad,nyears+proyears,nsim) # the K array
  
  Agearray<-array(rep(1:maxage,each=nsim),dim=c(nsim,maxage))   # Age array
  Len_age<-array(NA,dim=c(nsim,maxage,nyears+proyears)) # Length at age array
  ind<-as.matrix(expand.grid(1:nsim,1:maxage,1:(nyears+proyears))) # an index for calculating Length at age
  Len_age[ind]<-Linfarray[ind[,c(1,3)]]*(1-exp(-Karray[ind[,c(1,3)]]*(Agearray[ind[,1:2]]-t0[ind[,1]])))
  Wt_age<-array(NA,dim=c(nsim,maxage,nyears+proyears)) # Weight at age array
  Wt_age[ind] <- OM@a*Len_age[ind]^OM@b                  # Calculation of weight array
  
  ageM <- -((log(1-lenM/Linf))/K) + t0 # calculate ageM from L50 and growth parameters (non-time-varying)
  ageM[ageM < 1] <- 1 # age at maturity must be at least 1 
  age95 <- -((log(1-len95/Linf))/K) + t0 	
  age95[age95 < 1] <- 1.5 # must be greater than 0 and ageM
  
  ageMsd <- sapply(1:nsim,getroot,ageM,age95)
  ageMarray <- array(ageM,dim=c(nsim,maxage)) # Age at maturity array
  Mat_age <- 1/(1+exp((ageMarray-(Agearray))/(ageMarray*ageMsd)))  # Maturity at age array
  
  # Selectivity at Length ------------------------------------------------------
  # if (max(OM@L5) > 1.5) {
  # message("L5 set too high (maximum value of 1.5). \nDefaulting to L5 = 1.5")
  # OM@L5[OM@L5 > 1.5] <- 1.5
  # }
  
  Selnyears <- length(OM@SelYears)
  # are selectivity parameters relative to size at maturity?
  chk <- class(OM@isRel)
  if (length(OM@isRel) < 1) OM@isRel <- "true" 
  if (chk == "character") {
    chkRel <- tolower(OM@isRel)
    if (chkRel == "true" | OM@isRel == "1") multi <- lenM
    if (chkRel == "false" | OM@isRel == "0") multi <- 1
  }
  if (chk == "numeric") {
    if (OM@isRel == 1) multi <- lenM
    if (OM@isRel == 0) multi <- 1
  }
  if (Selnyears <= 1) { 
    tL5 <- runif(nsim, OM@L5[1], OM@L5[2]) * multi  # length at 0.05% selectivity ascending
    tLFS <- runif(nsim, OM@LFS[1], OM@LFS[2]) * multi   # first length at 100% selection
    tVmaxlen <- runif(nsim,OM@Vmaxlen[1],OM@Vmaxlen[2])   # selectivity at maximum length 
    L5 <- matrix(tL5, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
    LFS <- matrix(tLFS, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
    Vmaxlen <- matrix(tVmaxlen, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
  }
  
  if (Selnyears > 1) { # More than one break point in historical selection pattern
    L5 <- matrix(0, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
    LFS <- matrix(0, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
    Vmaxlen <- matrix(0, nrow=nyears+proyears, ncol=nsim, byrow=TRUE)
    SelYears <- OM@SelYears
    # length at 0.05% selectivity ascending
    L5_bk <- mapply(runif, n=nsim, min=OM@L5Lower, max=OM@L5Upper) * multi 
    # first length at 100% selection
    LFS_bk <- mapply(runif, n=nsim, min=OM@LFSLower, max=OM@LFSUpper) * multi
    # selectivity at maximum length
    Vmaxlen_bk <- mapply(runif, n=nsim, min=OM@VmaxLower, max=OM@VmaxUpper)
    
    for (X in 1:(Selnyears-1)) { 
      bkyears <- SelYears[X]:SelYears[X+1]
      L5[bkyears,] <- matrix(rep((L5_bk[,X]), length(bkyears)), 
                             ncol=nsim, byrow=TRUE)
      LFS[bkyears,] <- matrix(rep((LFS_bk[,X]), length(bkyears)), 
                              ncol=nsim, byrow=TRUE)
      Vmaxlen[bkyears,] <- matrix(rep((Vmaxlen_bk[,X]), length(bkyears)), 
                                  ncol=nsim, byrow=TRUE)
    }
    restYears <- max(SelYears):(nyears+proyears)
    L5[restYears,] <- matrix(rep((L5_bk[,Selnyears]), 
                                 length(restYears)), ncol=nsim, byrow=TRUE)
    LFS[restYears,] <- matrix(rep((LFS_bk[,Selnyears]), 
                                  length(restYears)), ncol=nsim, byrow=TRUE)
    Vmaxlen[restYears,] <- matrix(rep((Vmaxlen_bk[,Selnyears]), 
                                      length(restYears)), ncol=nsim, byrow=TRUE)  
  }
  
  ind <- which(LFS/Linf>1, arr.ind=T)
  if (length(ind) > 0) {
    message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
    LFS[ind] <- Linf[ind[,2]] * 0.9 
  }
  
  # LFS[LFS/Linf>1]<-NA
  # LFS<-apply(LFS,1,function(x)x[!is.na(x)][1])
  
  mod <- -((log(1-LFS[nyears,]/Linf))/K) + t0 # the age at modal (or youngest max) selectivity
  # age05 <- -((log(1-L5/Linf))/K) + t0# the highest age at %5 selectivity
  
  V <- array(NA, dim=c(nsim, maxage, nyears+proyears))
  for (Yr in 1:(nyears+proyears)) { # selectivity pattern for all years (possibly updated in projection)
    V[,,Yr] <- t(sapply(1:nsim, SelectFun, L5[Yr,], LFS[Yr,], Vmaxlen[Yr,], Linfs=Linfarray[,Yr], Lens=Len_age[,,Yr]))
  }
  
  Asize<-cbind(Size_area_1,1-Size_area_1)
  
  print("Optimizing for user-specified movement")  # Print a progress update
  flush.console()                                  # refresh the console
  
  if(sfIsRunning()){ # if the cluster is initiated 
    sfExport(list=c("Frac_area_1","Prob_staying")) # export some of the new arrays and ...
    mov<-array(t(sfSapply(1:nsim,getmov,Frac_area_1=Frac_area_1,Prob_staying=Prob_staying)),dim=c(nsim,2,2)) # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  }else{ # no cluster initiated
    mov<-array(t(sapply(1:nsim,getmov,Frac_area_1=Frac_area_1,Prob_staying=Prob_staying)),dim=c(nsim,2,2)) # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  }
  
  nareas<-2  # default is a two area model
  N<-array(NA,dim=c(nsim,maxage,nyears,nareas))        # stock numbers array
  Biomass<-array(NA,dim=c(nsim,maxage,nyears,nareas))  # stock biomass array
  VBiomass<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # vulnerable biomass array
  
  SSN<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # spawning stock numbers array
  SSB<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # spawning stock biomass array
  FM<-array(NA,dim=c(nsim,maxage,nyears,nareas))  # fishing mortality rate array
  Z<-array(NA,dim=c(nsim,maxage,nyears,nareas))   # total mortality rate array
  
  Agearray<-array(rep(1:maxage,each=nsim),dim=c(nsim,maxage))   # Age array
  surv<-exp(-Marray[,1])^(Agearray-1)                           # Survival array
  Nfrac<-surv*Mat_age                                           # predicted Numbers of mature ages
  initdist<-as.matrix(cbind(Frac_area_1,1-Frac_area_1))         # Get the initial spatial distribution of each simulated population
  
  R0a<-R0*initdist                                              # Unfished recruitment by area
  
  SAYR<-as.matrix(expand.grid(1:nareas,1,1:maxage,1:nsim)[4:1]) # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
  SAY<-SAYR[,1:3]
  SA<-SAYR[,1:2]
  SR<-SAYR[,c(1,4)]
  S<-SAYR[,1]
  SY<-SAYR[,c(1,3)]
  
  SSN[SAYR]<-Nfrac[SA]*R0*initdist[SR]                           # Calculate initial spawning stock numbers
  N[SAYR]<-R0*surv[SA]*initdist[SR]                              # Calculate initial stock numbers
  
  Biomass[SAYR]<-N[SAYR]*Wt_age[SAY]                             # Calculate initial stock biomass
  SSB[SAYR]<-SSN[SAYR]*Wt_age[SAY]                               # Calculate spawning stock biomass
  VBiomass[SAYR]<-Biomass[SAYR]*V[SAY]                            # Calculate vunerable biomass
  SSN0<-apply(SSN[,,1,],c(1,3),sum)                              # Calculate unfished spawning stock numbers
  SSB0<-apply(SSB[,,1,],1,sum)                                   # Calculate unfished spawning stock numbers
  SSBpR<-SSB0/R0                                                 # Spawning stock biomass per recruit
  SSB0a<-apply(SSB[,,1,],c(1,3),sum)                             # Calculate unfished spawning stock numbers
  bR<-log(5*hs)/(0.8*SSB0a)                                      # Ricker SR params
  aR<-exp(bR*SSB0a)/SSBpR                                        # Ricker SR params
  
  print("Optimizing for user-specified depletion")               # Print a progress update
  flush.console()                                                # update console
  
  if(sfIsRunning()){
    sfExport(list=c("dep","Find","Perr","Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
    qs<-sfSapply(1:nsim,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
  }else{
    qs <- sapply(1:nsim,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
  }
  
  # Check that depletion target is reached
  HighQ <- which(qs> 13 | qs < 0.005)
  if (length(HighQ) > 0) { # If q has hit bound, re-sample depletion and try again. Tries 30 times 
    # and then alerts user
    Err <- TRUE
    Nsec <- 10
    Nprob <- length(HighQ)
    message(Nprob," simulations have final biomass that is not close to specified depletion. \n qs have hit bounds.\n ")
    message("Re-sampling depletion and trying again")
    
    # Attempt again with different depletions
    count <- 0
    while (Err & count < 30) {
      count <- count + 1
      # message("Attempt ", count)
      cat(".")
      flush.console()
      Nprob <- length(HighQ)
      dep[HighQ] <- runif(Nprob,OM@D[1],OM@D[2])
      if(sfIsRunning()){
        sfExport(list=c("dep","Find","Perr","Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
        qs[HighQ]<-sfSapply(HighQ,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
      }else{
        qs[HighQ] <- sapply(HighQ,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
      }
      HighQ <- which(qs > 13 | qs < 0.005)
      if (length(HighQ) == 0) Err <- FALSE
    }
    if (!Err) {
      cat ("Success \n")
      cat("q range = ", range(qs),"\n")
      flush.console()
    }
    if (Err) { # still a problem
      # print(qs[HighQ])
      stop("Can't get to specified level of depletion. Try again for a complete new sample, or modify the input parameters")
      # print(qs[HighQ])
      # cat ("qs still very high \n")
      # cat ("Press ESC to quit now, or will continue in", Nsec, "seconds \n")
      # flush.console()
      # for (xx in 1:Nsec) {
      # Sys.sleep(1)
      # message(Nsec+1-xx)
      # flush.console()
      # }
    }  
  }
  
  print("Calculating historical stock and fishing dynamics")     # Print a progress update
  flush.console()                                                # update console
  
  fishdist<-(apply(VBiomass[,,1,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass[,,1,],c(1,3),sum)^Spat_targ,1,mean)  # spatial preference according to spatial biomass
  FM[SAYR]<-qs[S]*Find[SY]*V[SAY]*fishdist[SR]                    # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
  Z[SAYR]<-FM[SAYR]+Marray[SY]                                   # Total mortality rate                 
  
  for(y in 1:(nyears-1)){
    # set up some indices for indexed calculation
    SAYR<-as.matrix(expand.grid(1:nareas,y,1:maxage,1:nsim)[4:1]) # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
    SAY1R<-as.matrix(expand.grid(1:nareas,y+1,1:maxage,1:nsim)[4:1])
    SAY<-SAYR[,1:3]
    SA<-SAYR[,1:2]
    SR<-SAYR[,c(1,4)]
    S<-SAYR[,1]
    SY<-SAYR[,c(1,3)]
    SY1<-SAY1R[,c(1,3)]
    indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y+1,1:maxage,1:nsim)[5:1]) # Movement master index
    indMov2<-indMov[,c(1,2,3,4)]                                               # Movement from index
    indMov3<-indMov[,c(1,4,5)]                                                 # Movement to index
    
    if(SRrel[1]==1){
      N[,1,y+1,]<-Perr[,y]*(0.8*R0a*hs*apply(SSB[,,y,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB[,,y,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
    }else{ # most transparent form of the Ricker uses alpha and beta params
      N[,1,y+1,]<-Perr[,y+nyears]*aR*apply(SSB[,,y,],c(1,3),sum)*exp(-bR*apply(SSB[,,y,],c(1,3),sum))
    }
    
    fishdist<-(apply(VBiomass[,,y,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass[,,y,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
    FM[SAY1R]<-qs[S]*Find[SY1]*V[SAY]*fishdist[SR]                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Z[SAY1R]<-FM[SAY1R]+Marray[SY]                                          # Total mortality rate
    N[,2:maxage,y+1,]<-N[,1:(maxage-1),y,]*exp(-Z[,1:(maxage-1),y,])        # Total mortality
    temp<-array(N[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))   # Move individuals
    N[,,y+1,]<-apply(temp,c(4,3,1),sum)
    Biomass[SAY1R]<-N[SAY1R]*Wt_age[SAY]                                    # Calculate biomass
    VBiomass[SAY1R]<-Biomass[SAY1R]*V[SAY]                                   # Calculate vulnerable biomass
    SSN[SAY1R]<-N[SAY1R]*Mat_age[SA]                                        # Calculate spawning stock numbers
    SSB[SAY1R]<-SSN[SAY1R]*Wt_age[SAY]                                      # Calculate spawning stock biomass
    
  } # end of year
  
  CN<-apply(N*(1-exp(-Z))*(FM/Z),c(1,3,2),sum)  # Catch in numbers
  CN[is.na(CN)]<-0
  CB<-Biomass*(1-exp(-Z))*(FM/Z)                                            # Catch in biomass
  
  Cbiasa<-array(Cbias,c(nsim,nyears+proyears))                              # Bias array
  Cerr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Csd,(nyears+proyears))),sdconv(1,rep(Csd,nyears+proyears))),c(nsim,nyears+proyears)) # composite of bias and observation error
  Cobs<-Cbiasa[,1:nyears]*Cerr[,1:nyears]*apply(CB,c(1,3),sum)              # Simulated observed catch (biomass)
  
  CAA<-array(NA,dim=c(nsim,nyears,maxage))                                  # Catch  at age array
  cond<-apply(CN,1:2,sum,na.rm=T)<1                                         # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  fixind<-as.matrix(cbind(expand.grid(1:nsim,1:nyears),rep(floor(maxage/3),nyears))) # more fix
  CN[fixind[cond,]]<-1                                                      # puts a catch in the most vulnerable age class
  # for(i in 1:nsim)for(j in 1:nyears)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CN[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
  for(i in 1:nsim)for(j in 1:nyears)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_ESS[i],CN[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
  
  # Temporary fix til effdist simulator is fixed
  #for(i in 1:nsim)CAA[i,1,]<-CAA[i,2,]
  
  LatASD <- Len_age * 0.1 # This is currently fixed to cv of 10%
  MaxBin <- ceiling(max(Linfarray) + 2 * max(LatASD))
  binWidth <- ceiling(0.03 * MaxBin)
  CAL_bins <- seq(from=0, to=MaxBin+binWidth, by=binWidth) 
  CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(CAL_bins)-1)
  nCALbins <- length(CAL_binsmid)
  
  CAL <- array(NA,dim=c(nsim,nyears,nCALbins))                                # the catch at length array
  LFC <- rep(NA,nsim) # length at first capture
  
  for(i in 1:nsim){
    for(j in 1:nyears){
      tempCN <- ceiling(-0.5+rmultinom(1, size=CAL_ESS[i], prob=CN[i,j,])*CAL_nsamp[i]/CAL_ESS[i])
      #ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5)          # sample expected age
      lens <- unlist(sapply(1:maxage, function (X) rnorm(tempCN[X],  Len_age[i,X,j], LatASD[i,X,j])))
      lens[lens > (max(Linfarray) + 2 * max(LatASD))|lens>max(CAL_bins)] <- max(Linfarray) + 2 * max(LatASD) # truncate at 2 sd 
      CAL[i,j,] <- hist(lens,CAL_bins,plot=F)$counts                       # assign to bins
      LFC[i] <- min(c(lens,LFC[i]),na.rm=T)                                # get the smallest CAL observation
      #CAL[i,j,] <- ceiling(rmultinom(1, size=ESS[i], prob=tempCAL)*CAL_nsamp[i]*CAL_ESS[i]-0.5) # could replace with Dirichlet distribution
    }
  }
  
  Ierr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Isd,nyears+proyears)),sdconv(1,rep(Isd,nyears+proyears))),c(nsim,nyears+proyears))
  II<-(apply(Biomass,c(1,3),sum)*Ierr[,1:nyears])^betas     # apply hyperstability / hyperdepletion
  II<-II/apply(II,1,mean)                                   # normalize
  
  print("Calculating MSY reference points")                 # Print a progress update
  flush.console()                                           # update the console
  if(sfIsRunning()){
    sfExport(list=c("Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage")) # export some newly made arrays to the cluster
    MSYrefs<-sfSapply(1:nsim,getFMSY,Marray,hs,Mat_age,Wt_age,R0,V=V[,,nyears],maxage,nyears,proyears=200,Spat_targ,mov,SRrel,aR,bR) # optimize for MSY reference points	
  }else{
    MSYrefs<-sapply(1:nsim,getFMSY,Marray,hs,Mat_age,Wt_age,R0,V=V[,,nyears],maxage,nyears,proyears=200,Spat_targ,mov,SRrel,aR,bR) # optimize for MSY reference points
  }
  
  MSY <- MSYrefs[1,]  # record the MSY results (Vulnerable)
  FMSY <- MSYrefs[2,] # instantaneous apical FMSY  (Vulnerable)
  VBMSY <-(MSY/(1-exp(-FMSY))) # Biomass at MSY (Vulnerable)
  UMSY <- MSY/VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
  SSBMSY <- MSYrefs[3,] # Spawing Stock Biomass at MSY
  BMSY_B0 <- MSYrefs[4,] # SSBMSY relative to unfished (SSB)
  
  BMSY_B0bias<-array(rlnorm(nsim*ntest,mconv(1,OM@BMSY_B0cv),sdconv(1,OM@BMSY_B0cv)),dim=c(nsim,ntest)) # trial samples of BMSY relative to unfished
  
  print("Calculating reference yield - best fixed F strategy") # Print a progress update
  flush.console()                                              # update the console
  if(sfIsRunning()){ # Numerically optimize for F that provides highest long term yield
    RefY<-sfSapply(1:nsim,getFref,Marray=Marray,Wt_age=Wt_age,Mat_age=Mat_age,
                   Perr=Perr,N_s=N[,,nyears,],SSN_s=SSN[,,nyears,], Biomass_s=Biomass[,,nyears,],
                   VBiomass_s=VBiomass[,,nyears,],SSB_s=SSB[,,nyears,],Vn=V[,,(nyears+1):(nyears+proyears)],
                   hs=hs,R0a=R0a,nyears=nyears,proyears=proyears,nareas=nareas,maxage=maxage,mov=mov,SSBpR=SSBpR,
                   aR=aR,bR=bR,SRrel=SRrel, Spat_targ=Spat_targ)
  }else{
    RefY<-sapply(1:nsim,getFref,Marray=Marray,Wt_age=Wt_age,Mat_age=Mat_age,Perr=Perr,N_s=N[,,nyears,],SSN_s=SSN[,,nyears,],
                 Biomass_s=Biomass[,,nyears,],VBiomass_s=VBiomass[,,nyears,],SSB_s=SSB[,,nyears,],
                 Vn=V[,,(nyears+1):(nyears+proyears)],hs=hs,R0a=R0a,nyears=nyears,proyears=proyears,nareas=nareas,maxage=maxage,mov=mov,SSBpR=SSBpR,
                 aR=aR,bR=bR,SRrel=SRrel, Spat_targ=Spat_targ) 
  }
  

  Depletion<-(apply(VBiomass[,,nyears,],1,sum)/apply(VBiomass[,,1,],1,sum))#^betas   # apply hyperstability / hyperdepletion
  #Depletion<-(apply(SSB[,,nyears,],1,sum)/apply(SSB[,,1,],1,sum))#^betas   # apply hyperstability / hyperdepletion
  
  #cbind(dep,Depletion)
  FMSY_M<-FMSY/M                      # ratio of true FMSY to natural mortality rate M
  # LFS<-Linf*(1-exp(-K*(mod-t0)))      # Length at full selection
  A<-apply(VBiomass[,,nyears,],1,sum) # Abundance
  OFLreal<-A*FMSY                     # the true simulated Over Fishing Limit
  Recerr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Recsd,(nyears+proyears))),sdconv(1,rep(Recsd,nyears+proyears))),c(nsim,nyears+proyears))
  
  test<-array(BMSY_B0*BMSY_B0bias,dim=c(nsim,ntest)) # the simulated observed BMSY_B0 
  indy<-array(rep(1:ntest,each=nsim),c(nsim,ntest))  # index
  indy[test>0.9]<-NA                                 # interval censor
  BMSY_B0bias<-BMSY_B0bias[cbind(1:nsim,apply(indy,1,min,na.rm=T))] # sample such that BMSY_B0<90%
  
  I3<-apply(Biomass,c(1,3),sum)^betas     # apply hyperstability / hyperdepletion
  I3<-I3/apply(I3,1,mean)                 # normalize index to mean 1
  Iref<-apply(I3[,1:5],1,mean)*BMSY_B0    # return the real target abundance index corresponding to BMSY
  
  hsim<-rep(NA,nsim)                      # simulate values in steepness 
  cond<-hs>0.6
  hsim[cond]<-0.2+rbeta(sum(hs>0.6),alphaconv((hs[cond]-0.2)/0.8,(1-(hs[cond]-0.2)/0.8)*OM@hcv),betaconv((hs[cond]-0.2)/0.8,(1-(hs[cond]-0.2)/0.8)*OM@hcv))*0.8
  hsim[!cond]<-0.2+rbeta(sum(hs<0.6),alphaconv((hs[!cond]-0.2)/0.8,(hs[!cond]-0.2)/0.8*OM@hcv),betaconv((hs[!cond]-0.2)/0.8,(hs[!cond]-0.2)/0.8*OM@hcv))*0.8
  hbias<-hsim/hs                          # back calculate the simulated bias
  
  qmu<--0.5*qcv^2                                      # Mean
  qvar<-array(exp(rnorm(proyears*nsim,rep(qmu,proyears),rep(qcv,proyears))),c(nsim,proyears)) # Variations in interannual variation
  colnames(qvar) <- paste0("qvar", 1:proyears)
  FinF <-Find[,nyears] # Effort in final historical year
  
  DLM_data<-new('DLM_data',stock="MSE")             # create a blank DLM data object
  if(reps==1)DLM_data<-OneRep(DLM_data)             # make stochastic variables certain for only one rep
  DLM_data<-replic8(DLM_data,nsim)                  # make nsim sized slots in the DLM data object
  DLM_data@Name<-OM@Name
  DLM_data@Year<-1:nyears
  DLM_data@Cat<-Cobs
  DLM_data@Ind<-II
  DLM_data@Rec<-apply(N[,1,,],c(1,2),sum)*Recerr[,1:nyears] 
  DLM_data@t<-rep(nyears,nsim)
  DLM_data@AvC<-apply(Cobs,1,mean)
  DLM_data@Dt<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
  DLM_data@Mort<-M*Mbias
  DLM_data@FMSY_M<-FMSY_M*FMSY_Mbias
  DLM_data@BMSY_B0<-BMSY_B0*BMSY_B0bias
  DLM_data@Cref<-MSY*Crefbias
  DLM_data@Bref<-VBMSY*Brefbias
  DLM_data@Iref<-Iref*Irefbias
  DLM_data@LFC<-LFC*LFCbias
  DLM_data@LFS<-LFS*LFSbias
  DLM_data@CAA<-CAA
  #DLM_data@Dep<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
  DLM_data@Dep<-matrix(c(Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr)), VBMSY/apply(VBiomass[,,1,],1,sum)*BMSY_B0bias), nrow=nsim, ncol=2)
  DLM_data@Abun<-A*Abias*rlnorm(nsim,mconv(1,Aerr),sdconv(1,Aerr))
  DLM_data@vbK<-K*Kbias
  DLM_data@vbt0<-t0*t0bias
  DLM_data@vbLinf<-Linf*Linfbias
  DLM_data@L50 <- lenM * lenMbias
  DLM_data@L95 <- len95 * lenMbias
  DLM_data@L95[DLM_data@L95>0.9*DLM_data@vbLinf]<-0.9*DLM_data@vbLinf[DLM_data@L95>0.9*DLM_data@vbLinf] # Set a hard limit on ratio of L95 to Linf
  DLM_data@L50[DLM_data@L50>0.9*DLM_data@L95]<-0.9*DLM_data@L95[DLM_data@L50>0.9*DLM_data@L95] # Set a hard limit on ratio of L95 to Linf
  DLM_data@steep<-hs*hbias
  DLM_data@CAL_bins<-CAL_bins
  DLM_data@CAL<-CAL
  MLbin<-(CAL_bins[1:(length(CAL_bins)-1)]+CAL_bins[2:length(CAL_bins)])/2
  temp<-CAL*rep(MLbin,each=nsim*nyears)
  DLM_data@ML<-apply(temp,1:2,sum)/apply(CAL,1:2,sum) 
  DLM_data@Lc<-array(MLbin[apply(CAL,1:2,which.max)],dim=c(nsim,nyears))
  nuCAL<-CAL
  for(i in 1:nsim)for(j in 1:nyears)nuCAL[i,j,1:match(max(1,DLM_data@Lc[i,j]),MLbin)]<0
  temp<-nuCAL*rep(MLbin,each=nsim*nyears)
  DLM_data@Lbar<-apply(temp,1:2,sum)/apply(nuCAL,1:2,sum)
  DLM_data@MaxAge<-maxage
  DLM_data@Units<-"unitless"
  DLM_data@Ref<-OFLreal
  DLM_data@Ref_type<-'Simulated OFL'
  DLM_data@wla<-rep(OM@a,nsim)
  DLM_data@wlb<-rep(OM@b,nsim)
  DLM_data@OM<-as.data.frame(cbind(RefY,M,Depletion=apply(SSB[,,nyears,],1,sum)/apply(SSB[,,1,],1,sum),A,BMSY_B0,FMSY_M,Mgrad,Msd,
                                   procsd,Esd,dFfinal,MSY,qinc,qcv, FMSY,Linf,K,t0,hs,Linfgrad,Kgrad,Linfsd,
                                   recgrad,Ksd,ageM, L5[nyears,],LFS[nyears,],Vmaxlen[nyears,],LFC,OFLreal,
                                   Spat_targ,Frac_area_1,Prob_staying,AC, lenM, len95, qvar)) # put all the operating model parameters in one table
  
  names(DLM_data@OM)[26:28]<-c("L5","LFS","Vmaxlen") # These are missing labels in the line above
  
  DLM_data@Obs<-as.data.frame(cbind(Cbias,Csd,CAA_nsamp,CAA_ESS,CAL_nsamp,CAL_ESS,Isd,Dbias,Derr,Mbias,FMSY_Mbias,BMSY_B0bias,
                                    lenMbias,LFCbias,LFSbias,Abias,Aerr,Kbias,t0bias,Linfbias,hbias,Irefbias,Crefbias,Brefbias,betas))  # put all the observation error model parameters in one table
  
  DLM_data@LHYear<-OM@nyears # Last historical year is nyears (for fixed MPs)
  DLM_data@MPrec<-Cobs[,nyears]
  DLM_data@MPeff <- rep(1, nsim)
  DLM_data@Misc  <- vector("list", nsim)
  # assign("DLM_data",DLM_data,envir=.GlobalEnv) # for debugging fun
  
  # Run projections ===========================================================================
  
  if(is.na(MPs[1])) CheckMPs <- TRUE
  if (CheckMPs) {
    print("Determining available methods")  # print an progress report
    flush.console()                         # update the console
    PosMPs <- Can(DLM_data, timelimit=timelimit)  # list all the methods that could be applied to a DLM data object
    
    if(is.na(MPs[1])) {
      MPs<-PosMPs      # if the user does not supply an argument MPs run the MSE or all available methods
      print("No MPs specified: running all available")
    }	
    
    if(!is.na(MPs[1]))MPs<-MPs[MPs%in%PosMPs] # otherwise run the MSE for all methods that are deemed possible
    if(length(MPs)==0) {
      print(Cant(DLM_data, timelimit=timelimit))
      stop('MSE stopped: no viable methods \n\n') # if none of the user specied methods are possible stop the run
    }
  }	
  
  nMP <- length(MPs)                    # the total number of methods used
  
  MSElist<-list(DLM_data)[rep(1,nMP)]        # create a data object for each method (they have identical historical data and branch in projected years)
  
  B_BMSYa<-array(NA,dim=c(nsim,nMP,proyears))  # store the projected B_BMSY
  F_FMSYa<-array(NA,dim=c(nsim,nMP,proyears))  # store the projected F_FMSY
  Ba<-array(NA,dim=c(nsim,nMP,proyears))       # store the projected Biomass
  FMa<-array(NA,dim=c(nsim,nMP,proyears))      # store the projected fishing mortality rate
  Ca<-array(NA,dim=c(nsim,nMP,proyears))       # store the projected catch
  TACa<-array(NA,dim=c(nsim,nMP,proyears))     # store the projected TAC recommendation
  
  for(mm in 1:nMP){    # MSE Loop over methods
    
    print(paste(mm,"/",nMP," Running MSE for ",MPs[mm],sep=""))  # print a progress report
    flush.console()                                                  # update the console
    
    # projection arrays
    N_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    Biomass_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    VBiomass_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    SSN_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    SSB_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    FM_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    FM_nospace<-array(NA,dim=c(nsim,maxage,proyears,nareas)) # stores prospective F before reallocation to new areas
    FML<-array(NA,dim=c(nsim,nareas)) # last apical F
    Z_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    CB_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    
    # indexes
    SAYRL<-as.matrix(expand.grid(1:nsim,1:maxage,nyears,1:nareas))   # Final historical year
    SAYRt<-as.matrix(expand.grid(1:nsim,1:maxage,1+nyears,1:nareas)) # Trajectory year
    SAYR<-as.matrix(expand.grid(1:nsim,1:maxage,1,1:nareas))
    SYt<-SAYRt[,c(1,3)]
    SAYt<-SAYRt[,1:3]
    SR<-SAYR[,c(1,4)]
    SA1<-SAYR[,1:2]
    S1<-SAYR[,1]
    SY1<-SAYR[,c(1,3)]
    SAY1<-SAYR[,1:3] 
    SYA<-as.matrix(expand.grid(1:nsim,1,1:maxage))         # Projection year
    SY<-SYA[,1:2]
    SA<-SYA[,c(1,3)]
    SAY<-SYA[,c(1,3,2)]
    S<-SYA[,1]
    
    V_P <- V # Reset vulnerability array for MP 
    
    if(SRrel[1]==1){
      N_P[,1,1,]<-Perr[,nyears]*(0.8*R0a*hs*apply(SSB[,,nyears,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB[,,nyears,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
    }else{ # most transparent form of the Ricker uses alpha and beta params
      N_P[,1,1,]<-Perr[,nyears]*aR*apply(SSB[,,nyears,],c(1,3),sum)*exp(-bR*apply(SSB[,,nyears,],c(1,3),sum))
    }
    indMov<-as.matrix(expand.grid(1:nareas,1:nareas,1,1:maxage,1:nsim)[5:1])
    indMov2<-indMov[,c(1,2,3,4)]
    indMov3<-indMov[,c(1,4,5)]
    
    N_P[,2:maxage,1,]<-N[,1:(maxage-1),nyears,]*exp(-Z[,1:(maxage-1),nyears,])        # Total mortality
    temp<-array(N_P[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))   # Move individuals
    N_P[,,1,]<-apply(temp,c(4,3,1),sum)
    Biomass_P[SAYR]<-N_P[SAYR]*Wt_age[SAY1]                                    # Calculate biomass
    VBiomass_P[SAYR]<-Biomass_P[SAYR]*V_P[SAYt]                                  # Calculate vulnerable biomass
    SSN_P[SAYR]<-N_P[SAYR]*Mat_age[SA1]                                        # Calculate spawning stock numbers
    SSB_P[SAYR]<-SSN_P[SAYR]*Wt_age[SAY1]
    FML<-apply(FM[,,nyears,],c(1,3),max)
    
    if (class(match.fun(MPs[mm]))=="DLM_output") {
      DLM_data <- Sam(MSElist[[mm]],MPs=MPs[mm],perc=pstar,reps=reps)
      TACused <- apply(DLM_data@TAC,3,quantile,p=pstar,na.rm=T) 
      TACa[,mm,1] <- TACused
      
      fishdist<-(apply(VBiomass_P[,,1,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass_P[,,1,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
      
      CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-V_P[SAYt]*fishdist[SR]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
      
      temp<-CB_P[,,1,]/apply(CB_P[,,1,],1,sum)   # how catches are going to be distributed
      CB_P[,,1,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
      
      temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
      temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
      FM_P[SAYR]<--log(1-temp)
      
    }else{ # input control 
      runIn <- runInMP(MSElist[[mm]],MPs=MPs[mm], reps=reps) # Apply input control MP
      inc <- runIn[[1]]
      DLM_data <- runIn[[2]]
      
      Ai <- inc[1,,1]
      Ei <- inc[2,,1]
      Si <- t(inc[3:4,,1])
      newSel <- inc[5:6,,1]
      
      newUppLim <- inc[7,,1]
      newVmax <- inc[8,,1]
      y <- 1
      
      chngSel <- which(colSums(apply(newSel, 2, is.na))==0) # selectivity pattern changed 
      if (length(chngSel) > 0) {
        L5[y+nyears,chngSel] <- newSel[1,chngSel]
        LFS[y+nyears,chngSel] <- newSel[2,chngSel]
        if (any(!is.na(inc[7,,1]))) {
          ind <- which(!is.na(inc[7,,1]))
          Vmaxlen[y+nyears,ind] <- inc[7,ind,1]
        }
      }	  
      Vi <- t(sapply(1:nsim, SelectFun, L5[y+nyears,], LFS[y+nyears,], 
                     Vmaxlen[y+nyears,], Linfs=Linfarray[,y+nyears], 
                     Lens=Len_age[,,y+nyears]))
      
      # Maximum Size Limit 
      if (!all(is.na(newUppLim))) { # a upper size limit has been set 
        Vi[Len_age[,,y+nyears] >= newUppLim] <- 0   
      }
      # Vuln flag 
      Vchange <- any(!is.na(inc[5:8]))
      
      if(sum(Si!=1)==0){ # if there is no spatial closure
        if(!Vchange){ # if no vulnerability schedule is specified
          newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
          fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          FM_P[SAYR]<-FinF[S1]*Ei[S1]*V_P[SAYt]*fishdist[SR]*qvar[SY1]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        }else{
          V_P[,,(nyears+1):(proyears+nyears)] <- Vi # Update vulnerability schedule for all future years	  
          newVB<-apply(VBiomass_P[,,y,]*Vi[SA1],c(1,3),sum) # vulnerability modified
          fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA1]*fishdist[SR]*qvar[SY1]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        }
      }else{  # A spatial closure
        if(!Vchange){ # if no vulnerability schedule is specified
          newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
          fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
          FM_P[SAYR]<-FinF[S1]*Ei[S1]*V_P[SAYt]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY1]*qs[S1]^(1+qinc[S1]/100)^y 
        }else{
          V_P[,,(nyears+1):(proyears+nyears)] <- Vi # Update vulnerability schedule for all future years
          newVB<-apply(VBiomass_P[,,y,]*Vi[SA1],c(1,3),sum) # vulnerability modified
          fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
          FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA1]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY1]*qs[S1]^(1+qinc[S1]/100)^y 
        } # vulnerability specified
      }  # spatial closure specified  
    }   # input control  
    
    Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt]
    # CB_P[SAYR] <- Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
    CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] *  Biomass_P[SAYR]*(1-exp(-Z_P[SAYR])) 
    
    TACa[ ,mm,1] <- apply(CB_P[,,1,], 1, sum) # Adjust TAC to actual catch in the year 
    # To account for years where TAC is higher than catch
    
    upyrs<-1+(0:(floor(proyears/interval)-1))*interval  # the years in which there are updates (every three years)
    cat(".")
    flush.console()
    
    for(y in 2:proyears){
      cat(".")
      flush.console()
      if(class(match.fun(MPs[mm]))=="DLM_output") TACa[,mm,y]<-TACused
      SAYRt<-as.matrix(expand.grid(1:nsim,1:maxage,y+nyears,1:nareas)) # Trajectory year
      SAYt<-SAYRt[,1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt<-SAYRt[,c(1,3)]
      SAY1R<-as.matrix(expand.grid(1:nsim,1:maxage,y-1,1:nareas))
      SAYR<-as.matrix(expand.grid(1:nsim,1:maxage,y,1:nareas))
      SY<-SAYR[,c(1,3)]
      SA<-SAYR[,1:2]
      S1<-SAYR[,1]
      
      SAY<-SAYR[,1:3]
      S<-SAYR[,1]
      SR<-SAYR[,c(1,4)]
      SA2YR<-as.matrix(expand.grid(1:nsim,2:maxage,y,1:nareas))
      SA1YR<-as.matrix(expand.grid(1:nsim,1:(maxage-1),y-1,1:nareas))
      indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y,1:maxage,1:nsim)[5:1])
      indMov2<-indMov[,c(1,2,3,4)]
      indMov3<-indMov[,c(1,4,5)]
      
      N_P[SA2YR]<-N_P[SA1YR]*exp(-Z_P[SA1YR])         # Total mortality
      if(SRrel[1]==1){
        N_P[,1,y,]<-Perr[,y+nyears]*(0.8*R0a*hs*apply(SSB_P[,,y-1,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB_P[,,y-1,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
      }else{ # most transparent form of the Ricker uses alpha and beta params
        N_P[,1,y,]<-Perr[,y+nyears]*aR*apply(SSB_P[,,y-1,],c(1,3),sum)*exp(-bR*apply(SSB_P[,,y-1,],c(1,3),sum))
      }
      
      temp<-array(N_P[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))  # Move individuals
      N_P[,,y,]<-apply(temp,c(4,3,1),sum)
      
      Biomass_P[SAYR]<-N_P[SAYR]*Wt_age[SAYt]                                    # Calculate biomass
      VBiomass_P[SAYR]<-Biomass_P[SAYR]*V_P[SAYt]                       # Calculate vulnerable biomass
      SSN_P[SAYR]<-N_P[SAYR]*Mat_age[SA]                                       # Calculate spawning stock numbers
      SSB_P[SAYR]<-SSN_P[SAYR]*Wt_age[SAYt]                                      # Calculate spawning stock biomass
      
      if(y%in%upyrs){  # rewrite the DLM object and run the TAC function
        yind<-upyrs[match(y,upyrs)-1]:(upyrs[match(y,upyrs)]-1)
        CNtemp<-array(N_P[,,yind,]*exp(Z_P[,,yind,])*(1-exp(-Z_P[,,yind,]))*(FM_P[,,yind,]/Z_P[,,yind,]),c(nsim,maxage,interval,nareas))
        CBtemp<-array(Biomass_P[,,yind,]*exp(Z_P[,,yind,])*(1-exp(-Z_P[,,yind,]))*(FM_P[,,yind,]/Z_P[,,yind,]),c(nsim,maxage,interval,nareas))
        CNtemp[is.na(CNtemp)]<-0
        CBtemp[is.na(CBtemp)]<-0
        CNtemp[!is.finite(CNtemp)]<-0
        CBtemp[!is.finite(CBtemp)]<-0
        CNtemp<-apply(CNtemp,c(1,3,2),sum,na.rm=T)
        
        Cobs<-Cbiasa[,nyears+yind]*Cerr[,nyears+yind]*apply(CBtemp,c(1,3),sum,na.rm=T)
        Cobs[is.na(Cobs)]<-tiny
        Recobs<-Recerr[,nyears+yind]*apply(array(N_P[,1,yind,],c(nsim,interval,nareas)),c(1,2),sum)
        
        cond<-apply(CNtemp,1:2,sum,na.rm=T)<1                                         # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
        fixind<-as.matrix(cbind(expand.grid(1:nsim,1:interval),rep(floor(maxage/3),interval))) # more fix
        
        # assign("fixind",fixind,envir=.GlobalEnv) # for debugging fun
        # assign("CNtemp",CNtemp,envir=.GlobalEnv) # for debugging fun
        
        CNtemp[fixind[cond,]]<-1                                                      # puts a catch in the most vulnerable age class
        
        CAA<-array(NA,dim=c(nsim,interval,maxage))                                  # Catch  at age array
        # for(i in 1:nsim)for(j in 1:interval)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CNtemp[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
        for(i in 1:nsim)for(j in 1:interval)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_ESS[i],CN[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
        CAL <- array(NA,dim=c(nsim,interval,nCALbins))                                # the catch at length array
        CNtemp[is.na(CNtemp)]<-0
        for(i in 1:nsim){
          for(j in 1:interval){
            yy<-yind[j]
            tempCN <- ceiling(-0.5+rmultinom(1, size=CAL_ESS[i], prob=CN[i,j,])*CAL_nsamp[i]/CAL_ESS[i])
            #ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5)          # sample expected age
            lens <- unlist(sapply(1:maxage, function (X) rnorm(tempCN[X],  Len_age[i,X,yy+nyears], LatASD[i,X,yy+nyears])))
            lens[lens > (max(Linfarray) + 2 * max(LatASD))|lens>max(CAL_bins)] <- max(Linfarray) + 2 * max(LatASD) # truncate at 2 sd 
            CAL[i,j,] <- hist(lens,CAL_bins,plot=F)$counts                       # assign to bins
          }
        }
        
        I2<-cbind(apply(Biomass,c(1,3),sum),apply(Biomass_P,c(1,3),sum)[,1:(y-1)])*Ierr[,1:(nyears+(y-1))]^betas
        I2[is.na(I2)]<-tiny
        I2<-I2/apply(I2,1,mean)
        
        #Depletion <- apply(Biomass_P[,,y,],1,sum)/apply(Biomass[,,1,],1,sum)
        Depletion <- apply(VBiomass_P[,,y,],1,sum)/apply(VBiomass[,,1,],1,sum)
        Depletion[Depletion < tiny] <- tiny
        A<-apply(VBiomass_P[,,y,],1,sum)
        A[is.na(A)]<-tiny
        OFLreal<-A*FMSY
        
        # assign all the new data
        MSElist[[mm]]@OM$A<-A
        MSElist[[mm]]@Year<-1:(nyears+y-1)
        MSElist[[mm]]@Cat<-cbind(MSElist[[mm]]@Cat,Cobs)
        MSElist[[mm]]@Ind<-I2
        MSElist[[mm]]@Rec<-cbind(MSElist[[mm]]@Rec,Recobs)
        MSElist[[mm]]@t<-rep(nyears+y,nsim)
        MSElist[[mm]]@AvC<-apply(MSElist[[mm]]@Cat,1,mean)
        MSElist[[mm]]@Dt<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
        oldCAA<-MSElist[[mm]]@CAA
        MSElist[[mm]]@CAA<-array(0,dim=c(nsim,nyears+y-1,maxage))
        MSElist[[mm]]@CAA[,1:(nyears+y-interval-1),]<-oldCAA
        MSElist[[mm]]@CAA[,nyears+yind,]<-CAA
        #MSElist[[mm]]@Dep<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
        MSElist[[mm]]@Abun<-A*Abias*rlnorm(nsim,mconv(1,Aerr),sdconv(1,Aerr))
        MSElist[[mm]]@CAL_bins<-CAL_bins
        oldCAL<-MSElist[[mm]]@CAL
        MSElist[[mm]]@CAL<-array(0,dim=c(nsim,nyears+y-1,nCALbins))
        MSElist[[mm]]@CAL[,1:(nyears+y-interval-1),]<-oldCAL
        MSElist[[mm]]@CAL[,nyears+yind,]<-CAL[,1:interval,]
        
        temp<-CAL*rep(MLbin,each=nsim*interval)
        MSElist[[mm]]@ML<-cbind(MSElist[[mm]]@ML,apply(temp,1:2,sum)/apply(CAL,1:2,sum)) 
        MSElist[[mm]]@Lc<-cbind(MSElist[[mm]]@Lc,array(MLbin[apply(CAL,1:2,which.max)],dim=c(nsim,interval)))
        nuCAL<-CAL
        for(i in 1:nsim)for(j in 1:interval)nuCAL[i,j,1:match(max(1,MSElist[[mm]]@Lc[i,j]),MLbin)]<0
        temp<-nuCAL*rep(MLbin,each=nsim*interval)
        MSElist[[mm]]@Lbar<-cbind(MSElist[[mm]]@Lbar,apply(temp,1:2,sum)/apply(nuCAL,1:2,sum))
        
        MSElist[[mm]]@Ref<-OFLreal
        MSElist[[mm]]@Ref_type<-'Simulated OFL'
        MSElist[[mm]]@Misc  <- DLM_data@Misc
        
        # assign("DLM_data",MSElist[[mm]],envir=.GlobalEnv) # for debugging fun
        
        if(class(match.fun(MPs[mm]))=="DLM_output"){
          DLM_data <- Sam(MSElist[[mm]],MPs=MPs[mm],perc=pstar,reps=reps)
          TACused<-apply(DLM_data@TAC,3,quantile,p=pstar,na.rm=TRUE) #
          NAs <- which(is.na(TACused))
          if (length(NAs) >0 ) { # robustifying TAC setting!
            TACused[NAs] <- TACa[NAs,mm,y-1] #
            if (!exists("store")) store <- list()
            store <- append(store, c(MPs[mm], NAs))
          }
          TACa[,mm,y]<-TACused
          MSElist[[mm]]@MPrec<-TACused
          fishdist<-(apply(VBiomass_P[,,y,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass_P[,,y,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          
          CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-V_P[SAYt]*fishdist[SR]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          
          temp<-CB_P[,,y,]/apply(CB_P[,,y,],1,sum)   # how catches are going to be distributed
          CB_P[,,y,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
          
          temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
          temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
          FM_P[SAYR]<--log(1-temp)
          
        }else{
          MSElist[[mm]]@MPeff<-Ei
          runIn <- runInMP(MSElist[[mm]],MPs=MPs[mm], reps=reps) # Apply input control MP
          inc <- runIn[[1]]
          DLM_data <- runIn[[2]]
          Ai<-inc[1,,1]
          Ei<-inc[2,,1]
          Si<-t(inc[3:4,,1])
          newSel<-(inc[5:6,,1])
          
          newUppLim <- inc[7,,1]
          newVmax <- inc[8,,1]
          
          chngSel <- which(colSums(apply(newSel, 2, is.na))==0) # selectivity pattern changed 
          if (length(chngSel) >0) {
            L5[y+nyears,chngSel] <- newSel[1,chngSel]
            LFS[y+nyears,chngSel] <- newSel[2,chngSel]
            if (any(!is.na(inc[7,,1]))) {
              ind <- which(!is.na(inc[7,,1]))
              Vmaxlen[y+nyears,ind] <- inc[7,ind,1]
            }
          }	
          Vi <- t(sapply(1:nsim, SelectFun, L5[y+nyears,], LFS[y+nyears,], 
                         Vmaxlen[y+nyears,], Linfs=Linfarray[,y+nyears], 
                         Lens=Len_age[,,y+nyears]))
          
          # Maximum Size Limit 
          if (!all(is.na(newUppLim))) { # a upper size limit has been set 
            Vi[Len_age[,,y+nyears] >= newUppLim] <- 0   
          }
          # Vuln flag 
          Vchange <- any(!is.na(inc[5:8]))	
          
          if(sum(Si!=1)==0){ # if there is no spatial closure
            if(!Vchange){ # if no vulnerability schedule is specified
              newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
              fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V_P[SAYt]*fishdist[SR]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
              
              QVar <- apply(qvar[,1:y], 1, prod) # qvar from 1:y
              FM_P[SAYR] <- FinF[S1] * qs[S1] * QVar[S1] * (1+qinc[S1]/100)^y  * V_P[SAYt] * Ei[S1] * fishdist[SR]
              
            }else{
              V_P[,,(y+nyears+1):(proyears+nyears)] <- Vi # Update vulnerability schedule for all future years
              newVB<-apply(VBiomass_P[,,y,]*Vi[SA],c(1,3),sum) # vulnerability modified
              fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA]*fishdist[SR]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
              
              QVar <- apply(qvar[,1:y], 1, prod) # qvar from 1:y
              FM_P[SAYR] <- FinF[S1] * qs[S1] * QVar[S1] * (1+qinc[S1]/100)^y * Vi[SA] *Ei[S1] * fishdist[SR]
            }
          }else{  # A spatial closure
            if(!Vchange){ # if no vulnerability schedule is specified
              newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
              fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR]<-FinF[S1]*Ei[S1]*V_P[SAYt]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y 
              
              QVar <- apply(qvar[,1:y], 1, prod) # qvar from 1:y
              FM_P[SAYR] <- FinF[S1] * qs[S1] * QVar[S1] * (1+qinc[S1]/100)^y * V_P[SAYt] *Ei[S1] * fishdist[SR] * Emult[S1]
            }else{
              V_P[,,(y+nyears+1):(proyears+nyears)] <- Vi # Update vulnerability schedule for all future years
              newVB<-apply(VBiomass_P[,,y,]*Vi[SA],c(1,3),sum) # vulnerability modified
              fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y 
              
              QVar <- apply(qvar[,1:y], 1, prod) # qvar from 1:y
              FM_P[SAYR] <- FinF[S1] * qs[S1] * QVar[S1] * (1+qinc[S1]/100)^y * Vi[SA] *Ei[S1] * fishdist[SR] * Emult[S1]
            } #vuln not changed
          }   # spatial closure
          
        }   # input or output control 
        Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt] 
        # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
        CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] *  Biomass_P[SAYR]*(1-exp(-Z_P[SAYR]))
        
        TACused <- apply(CB_P[,,y,], 1, sum) # Set last years TAC to actual catch from last year
        TACa[,mm,y]<-TACused
        MSElist[[mm]]@MPrec<-TACused		
      }else{ # not an update yr 
        vbio <- apply(VBiomass_P[,,y,],c(1,3),sum)
        fishdist <- (vbio^Spat_targ)/apply(vbio^Spat_targ,1,mean) # calculate distribution of effort 	  
        if(class(match.fun(MPs[mm]))=="DLM_output"){
          CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-fishdist[SR]*V_P[SAYt]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          temp<-CB_P[,,y,]/apply(CB_P[,,y,],1,sum)   # how catches are going to be distributed
          CB_P[,,y,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
          temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
          temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
          FM_P[SAYR]<--log(1-temp)
        }else{ #input control
          # FM_P[SAYR] <- FM_P[SAY1R]*qvar[SY] *(1+qinc[S1]/100)^y  # add fishing efficiency changes and variability
          FM_P[SAYR] <- FM_P[SAY1R]*qvar[SY] *(1+qinc[S1]/100)  # add fishing efficiency changes and variability
        }
        Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt]
        # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
        CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] *  Biomass_P[SAYR]*(1-exp(-Z_P[SAYR]))
        
      } # not an update year
      
    } # end of year
    
    
    B_BMSYa[,mm,]<-apply(SSB_P,c(1,3),sum)/SSBMSY # SSB relative to SSBMSY
    
    F_FMSYa[,mm,]<-(-log(1-apply(CB_P,c(1,3),sum)/(apply(CB_P,c(1,3),sum)+apply(VBiomass_P,c(1,3),sum))))/FMSY
    Ba[,mm,]<-apply(Biomass_P,c(1,3),sum)
    FMa[,mm,]<--log(1-apply(CB_P,c(1,3),sum)/(apply(CB_P,c(1,3),sum)+apply(VBiomass_P,c(1,3),sum)))
    Ca[,mm,]<-apply(CB_P,c(1,3),sum)
    cat("\n")
  }    # end of mm methods
  new('MSE',Name=OM@Name,nyears,proyears,nMP,MPs,nsim,OMtable=DLM_data@OM,DLM_data@Obs,B_BMSYa,F_FMSYa,Ba,FMa,Ca,TACa,SSB_hist=SSB,CB_hist=CB,FM_hist=FM)
}
environment(runMSE) <- asNamespace('DLMtool')
sfExport("runMSE")

# DCAC static
DCAC_s<-function(x,DLM_data,reps=100){
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_tot<-sum(DLM_data@Cat[x,1:refyr])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])  
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) 
  depo<-max(0.01,min(0.99,DLM_data@Dep[x,1]))  
  B_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x]))) 
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x])) 
  TACfilter(C_tot/(refyr+((1-B_K)/(BMSY_K*FMSY_M*Mdb))))
}
class(DCAC_s)<-"DLM_output"
environment(DCAC_s)<-asNamespace('DLMtool')
sfExport("DCAC_s")

#DCAC adjusted
DCAC2_s<-function(x,DLM_data,reps=100){
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_tot<-sum(DLM_data@Cat[x,1:refyr])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) 
  depo<-max(0.01,min(0.99,DLM_data@Dep[x,1]))  
  B_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  susC<-C_tot/(refyr+((1-B_K)/(BMSY_K*FMSY_M*Mdb)))
  TAC<-susC*B_K/BMSY_K
  TACfilter(TAC)
}
class(DCAC2_s)<-"DLM_output"
environment(DCAC2_s)<-asNamespace('DLMtool')
sfExport("DCAC2_s")

#DBSRA Static
DBSRA_s<-function(x,DLM_data,reps=100){ 
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_hist<-DLM_data@Cat[x,1:refyr]
  MSYc<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x,1]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x,1])) 
    B_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    B_K<-B_K[B_K>0.00999&B_K<0.99001][1] 
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]   
    if(is.na(Mdb)) Mdb<-0.9 
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1] 
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm=TRUE)
      Max <- max(BMSY_K, na.rm=TRUE)
      if (Max <= 0.05) BMSY_K <- 0.05
      if (Min >= 0.95) BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K)) BMSY_K <- tryBMSY_K
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt2,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
                  FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=B_K,adelay=adelay,tol=0.01)
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc[DBSRAcount]<-Kc*BMSY_K*UMSYc
    DBSRAcount<-DBSRAcount+1
  } 
  TAC<-MSYc*B_K/BMSY_K
  TACfilter(TAC)
}  
class(DBSRA_s)<-"DLM_output"
environment(DBSRA_s)<-asNamespace('DLMtool')
sfExport("DBSRA_s")

# the optimization for B0 given DBSRA assumptions
DBSRAopt2<-function(lnK,C_hist,nys,Mdb,FMSY_M,BMSY_K,Bt_K,adelay){        
  Kc<-exp(lnK)
  n<-getn(BMSY_K)
  g<-gety(n)
  FMSY<-FMSY_M*Mdb
  UMSY<-(FMSY/(FMSY+Mdb))*(1-exp(-(FMSY+Mdb)))
  MSY<-Kc*BMSY_K*UMSY
  # Bjoin rules from Dick & MacCall 2011  ---------
  Bjoin_K<-0.5
  if(BMSY_K<0.3)Bjoin_K<-0.5*BMSY_K
  if(BMSY_K>0.3&BMSY_K<0.5)Bjoin_K<-0.75*BMSY_K-0.075
  Bjoin<-Bjoin_K*Kc
  PBjoin<-prodPTF(Bjoin_K,n,MSY)
  cp<-(1-n)*g*MSY*(Bjoin^(n-2))*Kc^-n
  Bc<-rep(NA,nys)
  Bc[1]<-Kc
  obj<-0
  for(yr in 2:nys){
    yref<-max(1,yr-adelay)
    if(Bc[yref]>Bjoin|BMSY_K>0.5){
      Bc[yr]<-Bc[yr-1]+g*MSY*(Bc[yref]/Kc)-g*MSY*(Bc[yref]/Kc)^n-C_hist[yr-1]
    }else{
      Bc[yr]<-Bc[yr-1]+Bc[yref]*((PBjoin/Bjoin)+cp*(Bc[yref]-Bjoin))-C_hist[yr-1]
    }
    if(Bc[yr]<0)obj<-obj+log(-Bc[yr])
    Bc[yr]<-max(0.000001,Bc[yr])
  }
  obj+((Bc[nys]/Kc)-Bt_K)^2
} 
environment(DBSRAopt2)<-asNamespace('DLMtool')
sfExport("DBSRAopt2")

#DCAC adjusted
DCAC2_adj<-function(x,DLM_data,reps=100){
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_tot<-sum(DLM_data@Cat[x,1:refyr])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) 
  depo<-max(0.01,min(0.99,DLM_data@Dep[x,1]))  
  B_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  susC<-C_tot/(refyr+((1-B_K)/(BMSY_K*FMSY_M*Mdb)))
  depo<-max(0.01,min(0.99,DLM_data@Dt[x]))  
  Bt_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  TAC<-susC*Bt_K/BMSY_K
  TACfilter(TAC)
}
class(DCAC2_adj)<-"DLM_output"
environment(DCAC2_adj)<-asNamespace('DLMtool')
sfExport("DCAC2_adj")


#DBSRA adjusted
DBSRA_adj<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_hist<-DLM_data@Cat[x,1:refyr]
  MSYc<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x,1]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x,1]))  
    B_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    B_K<-B_K[B_K>0.00999&B_K<0.99001][1] 
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]   
    if(is.na(Mdb)) Mdb<-0.9 
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@Dep[x,2],DLM_data@Dep[x,2]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1] 
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm=TRUE)
      Max <- max(BMSY_K, na.rm=TRUE)
      if (Max <= 0.05) BMSY_K <- 0.05
      if (Min >= 0.95) BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K)) BMSY_K <- tryBMSY_K
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt2,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
                  FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=B_K,adelay=adelay,tol=0.01)
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc[DBSRAcount]<-Kc*BMSY_K*UMSYc
    DBSRAcount<-DBSRAcount+1
  } 
  depo<-max(0.01,min(0.99,DLM_data@Dt[x]))  
  Bt_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  TAC<-MSYc*Bt_K/BMSY_K
  TACfilter(TAC)
} 
class(DBSRA_adj)<-"DLM_output"
environment(DBSRA_adj)<-asNamespace('DLMtool')
sfExport("DBSRA_adj")

Fmsy_PI<-function(x,DLM_data,reps=100){
  Z<-DLM_data@OM$FMSY[x] + DLM_data@OM$M[x]
  trlnorm(reps,  (DLM_data@OM$FMSY[x]/Z)*(1 - exp(-Z))*DLM_data@OM$A[x],0.01)
}
class(Fmsy_PI)<-"DLM_output"
environment(Fmsy_PI)<-asNamespace('DLMtool')
sfExport("Fmsy_PI")

SPMSY_s<-function (x, DLM_data, reps = 100) 
{
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_hist<-DLM_data@Cat[x,1:refyr]
  maxage <- DLM_data@MaxAge
  nsamp <- reps*200+200000
  rule <- rep(4, 3)
  if (DLM_data@vbK[x] > 0.3) {
    rule[1] <- 1
  } else if (DLM_data@vbK[x] < 0.3 & DLM_data@vbK[x] > 0.16) {
    rule[1] <- 2
  } else if (DLM_data@vbK[x] < 0.16 & DLM_data@vbK[x] > 0.05) {
    rule[1] <- 3
  }
  AM <-iVB(DLM_data@vbt0[x], DLM_data@vbK[x], DLM_data@vbLinf[x], 
            DLM_data@L50[x])
  if (AM < 1.5) {
    rule[2] <- 1
  } else if (AM < 4.5 & AM > 1.5) {
    rule[2] <- 2
  } else if (AM < 10 & AM > 4.5) {
    rule[2] <- 3
  }
  if (DLM_data@MaxAge < 4) {
    rule[3] <- 1
  } else if (DLM_data@MaxAge < 11 & DLM_data@MaxAge > 3) {
    rule[3] <- 2
  } else if (DLM_data@MaxAge < 31 & DLM_data@MaxAge > 10) {
    rule[3] <- 3
  }
  if (mean(rule) < 1.5) 
    rsamp <- runif(nsamp, 0.6, 1.5)
  if (mean(rule) > 1.5 & mean(rule) < 2.5) 
    rsamp <- runif(nsamp, 0.2, 1)
  if (mean(rule) > 2.5 & mean(rule) < 3.5) 
    rsamp <- runif(nsamp, 0.05, 0.5)
  if (mean(rule) > 3.5) 
    rsamp <- runif(nsamp, 0.015, 0.1)
  Ksamp <- runif(nsamp, mean(C_hist)/rsamp, (10*mean(C_hist))/rsamp)
  nyears <- length(C_hist)
  B <- array(NA, dim = c(nsamp, nyears))
  if (C_hist[1] < (0.5 * max(C_hist))) {
    B[, 1] <- Ksamp * runif(nsamp, 0.5, 0.95)
  } else {
    B[, 1] <- Ksamp * runif(nsamp, 0.3, 0.6)
  }
  if (C_hist[nyears] < (0.5 * max(C_hist))) {
    LB <- 0.01
    UB <- 0.4
  } else {
    LB <- 0.3
    UB <- 0.7
  }
  for (i in 2:nyears) {
    B[, i] <- B[, i - 1] - C_hist[i - 1]
    B[, i] <- B[, i] + rsamp * B[, i] * (1 - B[, i]/Ksamp)
  }
  B <- B/rep(Ksamp, nyears)
  cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  if (sum(cond) < 1) {
    B[B[, nyears] >= UB, nyears] <- UB
    cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  }
  
  if(reps>1){
    dep <- B[cond, nyears][1:reps]
    MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
    Kc <- Ksamp[cond][1:reps]
    rc <- rsamp[cond][1:reps]
  } else {
    dep <- mean(B[cond, nyears])
    MSY <- mean(rsamp[cond] * Ksamp[cond]/4)
    Kc <- mean(Ksamp[cond])
    rc <- mean(rsamp[cond])
  }
  if (sum(!is.na(Kc * dep * rc/2)) < ceiling(reps/10)) {
    cond <- (B[, nyears] >= 0.01) & (B[, nyears] <= 0.7)
    if(reps>1){
      dep <- B[cond, nyears][1:reps]
      MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
      Kc <- Ksamp[cond][1:reps]
      rc <- rsamp[cond][1:reps]
    } else {
      dep <- mean(B[cond, nyears])
      MSY <- mean(rsamp[cond] * Ksamp[cond]/4)
      Kc <- mean(Ksamp[cond])
      rc <- mean(rsamp[cond])
    }
  }
  depo<-max(0.01,min(0.99,DLM_data@Dep[x,1]))  
  B_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  TAC<-MSY*B_K/0.5
  TACfilter(TAC)
}
class(SPMSY_s)<-"DLM_output"
environment(SPMSY_s)<-asNamespace('DLMtool')
sfExport("SPMSY_s")


SPMSY_adj<-function (x, DLM_data, reps = 100) 
{
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_hist<-DLM_data@Cat[x,1:refyr]
  maxage <- DLM_data@MaxAge
  nsamp <- reps*200+200000
  rule <- rep(4, 3)
  if (DLM_data@vbK[x] > 0.3) {
    rule[1] <- 1
  } else if (DLM_data@vbK[x] < 0.3 & DLM_data@vbK[x] > 0.16) {
    rule[1] <- 2
  } else if (DLM_data@vbK[x] < 0.16 & DLM_data@vbK[x] > 0.05) {
    rule[1] <- 3
  }
  AM <-iVB(DLM_data@vbt0[x], DLM_data@vbK[x], DLM_data@vbLinf[x], 
           DLM_data@L50[x])
  if (AM < 1.5) {
    rule[2] <- 1
  } else if (AM < 4.5 & AM > 1.5) {
    rule[2] <- 2
  } else if (AM < 10 & AM > 4.5) {
    rule[2] <- 3
  }
  if (DLM_data@MaxAge < 4) {
    rule[3] <- 1
  } else if (DLM_data@MaxAge < 11 & DLM_data@MaxAge > 3) {
    rule[3] <- 2
  } else if (DLM_data@MaxAge < 31 & DLM_data@MaxAge > 10) {
    rule[3] <- 3
  }
  if (mean(rule) < 1.5) 
    rsamp <- runif(nsamp, 0.6, 1.5)
  if (mean(rule) > 1.5 & mean(rule) < 2.5) 
    rsamp <- runif(nsamp, 0.2, 1)
  if (mean(rule) > 2.5 & mean(rule) < 3.5) 
    rsamp <- runif(nsamp, 0.05, 0.5)
  if (mean(rule) > 3.5) 
    rsamp <- runif(nsamp, 0.015, 0.1)
  Ksamp <- runif(nsamp, mean(C_hist)/rsamp, (10*mean(C_hist))/rsamp)
  nyears <- length(C_hist)
  B <- array(NA, dim = c(nsamp, nyears))
  if (C_hist[1] < (0.5 * max(C_hist))) {
    B[, 1] <- Ksamp * runif(nsamp, 0.5, 0.95)
  } else {
    B[, 1] <- Ksamp * runif(nsamp, 0.3, 0.6)
  }
  if (C_hist[nyears] < (0.5 * max(C_hist))) {
    LB <- 0.01
    UB <- 0.4
  } else {
    LB <- 0.3
    UB <- 0.7
  }
  for (i in 2:nyears) {
    B[, i] <- B[, i - 1] - C_hist[i - 1]
    B[, i] <- B[, i] + rsamp * B[, i] * (1 - B[, i]/Ksamp)
  }
  B <- B/rep(Ksamp, nyears)
  cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  if (sum(cond) < 1) {
    B[B[, nyears] >= UB, nyears] <- UB
    cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  }
  
  if(reps>1){
    dep <- B[cond, nyears][1:reps]
    MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
    Kc <- Ksamp[cond][1:reps]
    rc <- rsamp[cond][1:reps]
  } else {
    dep <- mean(B[cond, nyears])
    MSY <- mean(rsamp[cond] * Ksamp[cond]/4)
    Kc <- mean(Ksamp[cond])
    rc <- mean(rsamp[cond])
  }
  if (sum(!is.na(Kc * dep * rc/2)) < ceiling(reps/10)) {
    cond <- (B[, nyears] >= 0.01) & (B[, nyears] <= 0.7)
    if(reps>1){
      dep <- B[cond, nyears][1:reps]
      MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
      Kc <- Ksamp[cond][1:reps]
      rc <- rsamp[cond][1:reps]
    } else {
      dep <- mean(B[cond, nyears])
      MSY <- mean(rsamp[cond] * Ksamp[cond]/4)
      Kc <- mean(Ksamp[cond])
      rc <- mean(rsamp[cond])
    }
  }
  depo<-max(0.01,min(0.99,DLM_data@Dt[x]))  
  Bt_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  TAC<-MSY*Bt_K/0.5
  TACfilter(TAC)
}
class(SPMSY_adj)<-"DLM_output"
environment(SPMSY_adj)<-asNamespace('DLMtool')
sfExport("SPMSY_adj")

MC_s<-function(x,DLM_data,reps=100){
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  mean(DLM_data@Cat[x,(refyr-19):refyr])
}
class(MC_s)<-"DLM_output"
environment(MC_s)<-asNamespace('DLMtool')
sfExport("MC_s")

MC_adj<-function(x,DLM_data,reps=100){
  refyr<-match(DLM_data@LHYear,DLM_data@Year)
  C_av<-mean(DLM_data@Cat[x,(refyr-19):refyr])
  depo<-max(0.01,min(0.99,DLM_data@Dt[x])) 
  Bt_K<-rbeta(reps,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  TAC<-C_av*2*Bt_K
  TACfilter(TAC)
}
class(MC_adj)<-"DLM_output"
environment(MC_adj)<-asNamespace('DLMtool')
sfExport("MC_adj")


#############################################################
#Run the analysis

sfExportAll()
set.seed(1) 

Herring@Perr<-c(0.7,0.9)
Rockfish@Perr<-c(0.5, 0.7)
Rockfish@h<-c(0.6, 0.95)
Rockfish@M<-c(0.025,0.06)

Stocks<-c("Herring","Bluefin_tuna_WAtl","Rockfish")
Obs<-c("Perfect_Info","Precise_Unbiased","Imprecise_Biased")
Dep<-matrix(c(0.05,0.2,0.2,0.6),nrow=2)

design<-expand.grid(Stocks,Obs,1:2)
nOMs<-nrow(design)

MSEs<-new('list')
OMs<-new('list')
Meths<-c("DCAC_s", "DCAC2_s", "DCAC2_adj", "DBSRA_s", "DBSRA_adj", "SPMSY_s", "SPMSY_adj", "Fmsy_PI", "MC_s", "MC_adj")

#Check convergence sample size
nsim<-500
proyears<-60
interval<-3
reps<-1
i=17
print(paste("=========  ", i," - ", design[i,1],design[i,2],design[i,3], "  ================"))
assign('temp_stock',get(as.character(design[i,1])))
temp_stock@D<-Dep[design[i,3],]
OMs[[i]]<-new('OM',temp_stock,Generic_fleet,get(as.character(design[i,2])))
OMs[[i]]@Name<-paste(as.character(design[i,1]),c("PI","DataRich","DataPoor")[as.integer(design[i,2])],c("5_30","30_70")[design[i,3]],sep="_")
MSEs[[i]]<-runMSE(OMs[[i]],Meths,nsim,proyears,interval,pstar=0.5,maxF=0.8,timelimit=20,reps)
Tplot(MSEs[[i]])

tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Convergence.tiff",width=10,height=8,units="in",res=600)
CheckConverg(MSEs[[i]])
dev.off()


#Run full analysis
nsim<-500
proyears<-60
interval<-3
reps<-1
for(i in 1:nOMs){
  print(paste("=========  ", i," - ", design[i,1],design[i,2],design[i,3], "  ================"))
  assign('temp_stock',get(as.character(design[i,1])))
  temp_stock@D<-Dep[design[i,3],]
  #OMs[[i]]<-new('OM',temp_stock,Generic_fleet,get(as.character(design[i,2])))
  
  OMs[[i]]<-new('OM',temp_stock,Generic_FlatE,get(as.character(design[i,2])))
  
  OMs[[i]]@Name<-paste(as.character(design[i,1]),c("PI","DataRich","DataPoor")[as.integer(design[i,2])],c("5_30","30_70")[design[i,3]],sep="_")
  MSEs[[i]]<-runMSE(OMs[[i]],Meths,nsim,proyears,interval,pstar=0.5,maxF=0.8,timelimit=20,reps)
  Tplot(MSEs[[i]])
}
save(OMs,file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/OMs_revision2_GenericFlatE.R")
save(MSEs,file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/MSEs_revision2_GenericFlatE.R")


#save(OMs,file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/OMs_revision1.R")
#save(MSEs,file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/MSEs_revision1.R")


# Data plots ==============================================================================================
load(file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/OMs_revision1.R")
load(file="C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/simulations/MSEs_revision1.R")

Stocks<-c("Herring","Bluefin tuna","Rockfish")
Obs<-c("Perfect_Info","Precise-unbiased","Imprecise-biased")
design<-expand.grid(Stocks,Obs,c("Depleted","Less depleted"))

library(ggplot2)
library(ggrepel)
library(wordcloud)
makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
coly<-makeTransparent(c("black","black",rep("grey",10)))
VYvY<-function(MSEobj,coly,px=T,py=T){
  VY<-rep(NA,MSEobj@nMPs)
  Y<-rep(NA,MSEobj@nMPs)
  yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  for(mm in 1:MSEobj@nMPs){
    AAVY<-apply(((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)/apply(MSEobj@C[,mm,y2],1,mean)
    VY[mm]<-round(sum(AAVY<0.1,na.rm=T)/MSEobj@nsim,3)
    Y[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)
  }
  return(cbind(VY, Y))
} 
B50vY_S<-function(MSEobj,coly,px=T,py=T){
  Y<-rep(NA,MSEobj@nMPs)
  B50<-rep(NA,MSEobj@nMPs)
  ystart<-1:10
  RefYd<-MSEobj@OM$RefY
  for(mm in 1:MSEobj@nMPs){
    Y[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)
    B50[mm]<-round(sum(MSEobj@B_BMSY[,mm,ystart]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,ystart])),3)
  }
  return(cbind(B50, Y))
} 
B50vY_L<-function(MSEobj,coly,px=T,py=T){
  Y<-rep(NA,MSEobj@nMPs)
  B50<-rep(NA,MSEobj@nMPs)
  yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  for(mm in 1:MSEobj@nMPs){
    Y[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)
    B50[mm]<-round(sum(MSEobj@B_BMSY[,mm,yend]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,yend])),3)
  }
  return(cbind(B50, Y))
} 


#######################
#Fig. Historic patterns plot
#######################
#select design row to plot
#BFT, depleted, imprecise-biased
MSEobj<-MSEs[[8]]

tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Fig 1 historic.tiff",width=3.5,height=8,units="in",res=600)
nf<-layout(matrix(c(1,2,3), 3, 1, byrow=TRUE))

#plot1
y<-apply(MSEobj@SSB_hist, c(1,3), sum)
M<-apply(y/y[,1], 2, median, na.rm=TRUE)
L<-apply(y/y[,1], 2, quantile, probs=0.025, na.rm=TRUE)
U<-apply(y/y[,1], 2, quantile, probs=0.975, na.rm=TRUE)
x<-seq(1,NROW(M),1)
par(mar=c(4,4,3,0))
plot(M, type = 'n', ylim=c(0,1.1), bty="n", las=1, ylab="", xlab="", font=2, font.lab=2)
mtext("SSB / SSB0", side=2, font=2, cex=0.8, line=2.5)
mtext("Historic year", side=1, font=2, cex=0.8, line=2)
polygon(c(rev(x), x), c(rev(U), L), col = 'grey80', border = NA)
lines(x, M, col="black", lwd=2)
mtext("(A)", side=3, adj=1, font=2, cex=0.8)

#plot2
y<-apply(MSEobj@CB_hist, c(1,3), sum)
M<-apply(y/y[,1], 2, median, na.rm=TRUE)
L<-apply(y/y[,1], 2, quantile, probs=0.025, na.rm=TRUE)
U<-apply(y/y[,1], 2, quantile, probs=0.975, na.rm=TRUE)
x<-seq(1,NROW(M),1)
par(mar=c(4,4,3,0))
plot(M, type = 'n', ylim=c(0,75), bty="n", las=1, ylab="", xlab="", font=2, font.lab=2)
mtext("Relative catch", side=2, font=2, cex=0.8, line=2.5)
mtext("Historic year", side=1, font=2, cex=0.8, line=2)
polygon(c(rev(x), x), c(rev(U), L), col = 'grey80', border = NA)
lines(x, M, col="black", lwd=2)
mtext("(B)", side=3, adj=1, font=2, cex=0.8)

#plot3
y<-apply(MSEobj@CB_hist, c(1,3), sum)/apply(MSEobj@SSB_hist, c(1,3), sum)
M<-apply(y/y[,1], 2, median, na.rm=TRUE)
L<-apply(y/y[,1], 2, quantile, probs=0.025, na.rm=TRUE)
U<-apply(y/y[,1], 2, quantile, probs=0.975, na.rm=TRUE)
x<-seq(1,NROW(M),1)
par(mar=c(4,4,3,0))
plot(M, type = 'n', ylim=c(0,160), bty="n", las=1, ylab="", xlab="", font=2, font.lab=2)
mtext("Relative fishing effort", side=2, font=2, cex=0.8, line=2.5)
mtext("Historic year", side=1, font=2, cex=0.8, line=2)
polygon(c(rev(x), x), c(rev(U), L), col = 'grey80', border = NA)
lines(x, M, col="black", lwd=2)
mtext("(C)", side=3, adj=1, font=2, cex=0.8)
dev.off()


#########
#Fig.Short term trade-off plot
########
ys<-c(4,7,13,16)
xs<-c(7:9,16:18)

j<-0
output<-matrix(ncol=6, nrow=0)
for(i in c(4:9,13:18)){
  j<-j+1
  out<-B50vY_S(MSEs[[i]],py=i%in%ys,px=i%in%xs,coly=coly)
  output<-rbind(output, cbind(out[,1], out[,2], design[i,], as.factor(MSEs[[i]]@MPs[1:MSEs[[i]]@nMPs]), deparse.level=0))
}  
colnames(output)<-c("B50", "Y", "Var1", "Var2", "Var3", "Lab")

plot.new()
dep.mod<-c(Depleted="More depleted", "Less depleted"="Less depleted")
ggplot(output, aes(B50, Y)) +
  geom_point(size=1) +
  geom_text_repel(data=output,
                  aes(output$B50, output$Y, label = output$Lab),
                  segment.size = 0.3,
                  size=2, 
                  fontface = 'bold', color="blue",
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0, "lines"),
                  show.legend = FALSE
  ) +
  xlim(0,1.2)+
  ylim(0,1.4) +
  abline(h=0)+
  abline(v=1.25)+
  ylab("Fraction of years yield > 0.5MSY") +
  xlab("Fraction of years biomass > 0.5BMSY") +
  theme_classic() +
  theme(strip.text.x=element_text(colour = "black", size=8, face="bold"), 
        strip.text.y=element_text(colour = "black", size=8, face="bold"),
        strip.background = element_rect(fill ="lightgrey"),
        axis.text=element_text(size=6),
        panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))+
  facet_grid(Var3 + Var2 ~ Var1, labeller=labeller(Var3 = dep.mod)) 
ggsave("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Fig 2 B50 v Y _S.tiff", width=7, height=5, units="in",dpi=300)



##########
#Fig Long term trade-off plot
#########
ys<-c(4,7,13,16)
xs<-c(7:9,16:18)

j<-0
output<-matrix(ncol=6, nrow=0)
for(i in c(4:9,13:18)){
  j<-j+1
  out<-B50vY_L(MSEs[[i]],py=i%in%ys,px=i%in%xs,coly=coly)
  output<-rbind(output, cbind(out[,1], out[,2], design[i,], as.factor(MSEs[[i]]@MPs[1:MSEs[[i]]@nMPs]), deparse.level=0))
}  
colnames(output)<-c("B50", "Y", "Var1", "Var2", "Var3", "Lab")

plot.new()
dep.mod<-c(Depleted="More depleted", "Less depleted"="Less depleted")
ggplot(output, aes(B50, Y)) +
  geom_point(size=1) +
  geom_text_repel(data=output,
                  aes(output$B50, output$Y, label = output$Lab),
                  segment.size = 0.3,
                  size=2, 
                  fontface = 'bold', color="blue",
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0, "lines"),
                  show.legend = FALSE
  ) +
  xlim(0,1.2)+
  ylim(0,1.4) +
  ylab("Fraction of years yield > 0.5MSY") +
  xlab("Fraction of years biomass > 0.5BMSY") +
  theme_classic() +
  theme(strip.text.x=element_text(colour = "black", size=8, face="bold"), 
        strip.text.y=element_text(colour = "black", size=8, face="bold"),
        strip.background = element_rect(fill ="lightgrey"),
        axis.text=element_text(size=6),
        panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))+
  facet_grid(Var3 + Var2 ~ Var1, labeller=labeller(Var3 = dep.mod)) 
ggsave("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Fig 3 B50 v Y _L.tiff", width=7, height=5, units="in",dpi=300)


#############
#Fig. Progression plot biomass
############
#select design row to plot
#BFT, depleted, imprecise-biased
MSEobj<-MSEs[[8]]

name<-c("(I)", "(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(J)", "(G)", "(H)")
tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/New Fig 4 progression.tiff",width=4,height=8,units="in",res=600)
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 5, 4, byrow=TRUE), widths=c(4,1,4,1))
#layout.show(nf)
for (mm in c(2,3,4,5,6,7,9,10,1,8)){ 
  
  M<-apply(MSEobj@B_BMSY[,mm,], 2, median, na.rm=TRUE)
  L<-apply(MSEobj@B_BMSY[,mm,], 2, quantile, probs=0.025, na.rm=TRUE)
  U<-apply(MSEobj@B_BMSY[,mm,], 2, quantile, probs=0.975, na.rm=TRUE)
  x<-seq(1,NROW(M),1)
  par(mar=c(4,4,2,0))
  plot(M, type = 'n', ylim=c(0,5), bty="n", las=1, ylab="", xlab="", font=2, font.lab=2, xlim=c(0,MSEobj@proyears))
  mtext("B relative to BMSY", side=2, font=2, cex=0.8, line=2.5)
  mtext("Projection year", side=1, font=2, cex=0.8, line=2)
  
  polygon(c(rev(x), x), c(rev(U), L), col = 'grey80', border = NA)
  lines(x, M, col="black", lwd=2)
  mtext(MSEobj@MPs[mm], side=3, adj=1, font=2, cex=0.8)
  abline(h=1, lty=3)
  
  par(mar=c(4,0,2,1))
  yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
  d<-density(rowMeans(MSEobj@B_BMSY[,mm,yend]),bw=0.4, na.rm=TRUE)
  plot(d$y, d$x, type="l", bty="n", axes=FALSE, ylim=c(0,5), xlab="", ylab="")
  polygon(d$y, d$x, col="black")
  mtext(name[mm], side=3, adj=1, font=2, cex=0.8)
}
dev.off()

#############
#Fig. Progression plot catch
############
#select design row to plot
#BFT, depleted, imprecise-biased
MSEobj<-MSEs[[8]]

name<-c("(I)", "(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(J)", "(G)", "(H)")
tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/New Fig 4 progression catch.tiff",width=4,height=8,units="in",res=600)
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 5, 4, byrow=TRUE), widths=c(4,1,4,1))
#layout.show(nf)
for (mm in c(2,3,4,5,6,7,9,10,1,8)){ 
  
  M<-apply(MSEobj@C[,mm,]/MSEobj@OM$MSY, 2, median, na.rm=TRUE)
  L<-apply(MSEobj@C[,mm,]/MSEobj@OM$MSY, 2, quantile, probs=0.025, na.rm=TRUE)
  U<-apply(MSEobj@C[,mm,]/MSEobj@OM$MSY, 2, quantile, probs=0.975, na.rm=TRUE)
  x<-seq(1,NROW(M),1)
  par(mar=c(4,4,2,0))
  plot(M, type = 'n', ylim=c(0,2), bty="n", las=1, ylab="", xlab="", font=2, font.lab=2, xlim=c(0,MSEobj@proyears))
  mtext("C relative to MSY", side=2, font=2, cex=0.8, line=2.5)
  mtext("Projection year", side=1, font=2, cex=0.8, line=2)
  
  polygon(c(rev(x), x), c(rev(U), L), col = 'grey80', border = NA)
  lines(x, M, col="black", lwd=2)
  mtext(MSEobj@MPs[mm], side=3, adj=1, font=2, cex=0.8)
  abline(h=1, lty=3)
  
  par(mar=c(4,0,2,1))
  yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
  d<-density(rowMeans(MSEobj@C[,mm,yend]/MSEobj@OM$MSY),bw=0.4, na.rm=TRUE)
  plot(d$y, d$x, type="l", bty="n", axes=FALSE, ylim=c(0,2), xlab="", ylab="")
  polygon(d$y, d$x, col="black")
  mtext(name[mm], side=3, adj=1, font=2, cex=0.8)
}
dev.off()


##############
#Fig. Long term stabilty plot
##############
ys<-c(4,7,13,16)
xs<-c(7:9,16:18)

j<-0
output<-matrix(ncol=6, nrow=0)
for(i in c(4:9,13:18)){
  j<-j+1
  out<-VYvY(MSEs[[i]],py=i%in%ys,px=i%in%xs,coly=coly)
  output<-rbind(output, cbind(out[,1], out[,2], design[i,], as.factor(MSEs[[i]]@MPs[1:MSEs[[i]]@nMPs]), deparse.level=0))
}  
colnames(output)<-c("VY", "Y", "Var1", "Var2", "Var3", "Lab")
plot.new()
dep.mod<-c(Depleted="More depleted", "Less depleted"="Less depleted")
ggplot(output, aes(VY, Y)) +
  geom_point(size=1) +
  geom_text_repel(data=output,
                  aes(output$VY, output$Y, label = output$Lab),
                  segment.size = 0.3,
                  size=2, 
                  fontface = 'bold', color="blue",
                  box.padding = unit(1.0, "lines"),
                  point.padding = unit(0, "lines"),
                  show.legend = FALSE
  ) +
  xlim(0,1.2)+
  ylim(0,1.4) +
  xlab("Stability in yield") +
  ylab("Fraction of years yield > 0.5MSY") +
  theme_classic() +
  theme(strip.text.x=element_text(colour = "black", size=8, face="bold"), 
        strip.text.y=element_text(colour = "black", size=8, face="bold"),
        strip.background = element_rect(fill ="lightgrey"),
        axis.text=element_text(size=6),
        panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))+
  facet_grid(Var3 + Var2 ~ Var1, labeller=labeller(Var3 = dep.mod)) 
ggsave("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Fig 6 VY vs Y.tiff", width=7, height=5, units="in",dpi=300)


##########
#Fig Depletion sensitivity plot
#########
 
tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/Fig 7 Depletion effect.tiff",width=6,height=8,units="in",res=600)

MSEobj<-MSEs[[8]]
Qnam<-c("Observation precision","Observation bias","'True' initial depletion")
qs<-cbind(MSEobj@Obs$Derr,MSEobj@Obs$Dbias,MSEobj@OM$Depletion)
yend<-max(MSEobj@proyears-9,1):MSEobj@proyears

mind<-c(1,2,3,4,5,6,7)
nm<-length(mind)
MPnam<-MSEobj@MPs[mind]
Yds<-apply(MSEobj@C[,mind,yend],1:2,mean)/array(MSEobj@C[,8,yend],c(MSEobj@nsim,nm))

plabs<-c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)","(q)","(r)","(s)", "(t)", "(u)", "(v)", "(w)", "(x)")
plabs<-plabs[c(1,4,7,10,13,16,19,22, 
               2,5,8,11,14,17,20,23,
               3,6,9,12,15,18,21,24)]

par(mfcol=c(nm,3),mai = c(0.15, 0.15, 0.2, 0.2),omi=c(0.2,0.3,0.04,0.04))

nbins=6
OMv<-array(NA,c(nm,3,nbins))
OMp<-apply(qs,2,quantile,p=seq(0,1,length.out=nbins+1))

for(q in 1:3){
  for(j in 1:nbins){
    for(mm in 1:nm){
      cond<-qs[,q]>OMp[j,q]&qs[,q]<OMp[j+1,q]
      OMv[mm,q,j]<-mean(Yds[cond,mm],na.rm=T)
    }
  }
}
OMs<-apply(OMv,1:2,sd,na.rm=T)

ncols<-40
#colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),90)[ncols:1]
minsd<-0
maxsd<-max(OMs,na.rm=T)
coly<-ceiling(OMs/maxsd*ncols)

k<-0
for(q in 1:3){
  for(mm in 1:nm){
   k<-k+1
     x<-qs[,q]
    y<-Yds[,mm]
    plot(x,y,col="white",ylim=c(0,2),axes=F)
    axis(1,pretty(OMp[,q]),pretty(OMp[,q]),cex.axis=0.8,padj=-1.5)
    axis(2,cex.axis=0.8, las=1)
    abline(v=OMp[,q],col="#99999960")
    points(x,y,col=colsse[coly[mm,q]],pch=19,cex=0.8)
    x2<-(OMp[1:nbins,q]+OMp[2:(nbins+1),q])/2
    y2<-OMv[mm,q,]
    lines(x2,y2)
    if(q==1)mtext(MPnam[mm],3,line=0,cex=0.7,adj=-0.1)
    if(mm==nm)mtext(Qnam[q],1,line=1.5,cex=0.75)
    mtext(plabs[k],3,adj=0.98,line=0.02,cex=0.7)
  }
}  
mtext("Yield relative to Fmsy_PI yield",2,outer=T,line=1,cex=0.85)
dev.off()



##############################################################
#Supplementary figures

##########
#Catch sensitivity plot
#########

tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/S1 Catch effect.tiff",width=4,height=8,units="in",res=600)

MSEobj<-MSEs[[8]]
Qnam<-c("Observation precision","Observation bias")
qs<-cbind(MSEobj@Obs$Csd, MSEobj@Obs$Cbias)
yend<-max(MSEobj@proyears-9,1):MSEobj@proyears

mind<-c(1,2,3,4,5,6,7)
nm<-length(mind)
MPnam<-MSEobj@MPs[mind]
Yds<-apply(MSEobj@C[,mind,yend],1:2,mean)/array(MSEobj@C[,8,yend],c(MSEobj@nsim,nm))

plabs<-c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)","(q)","(r)","(s)", "(t)", "(u)", "(v)", "(w)", "(x)")
plabs<-plabs[c(1,3,5,7,9,11,13,15, 
               2,4,6,8,10,12,14,16)]

par(mfcol=c(nm,2),mai = c(0.15, 0.15, 0.2, 0.2),omi=c(0.2,0.3,0.04,0.04))

nbins=6
OMv<-array(NA,c(nm,2,nbins))
OMp<-apply(qs,2,quantile,p=seq(0,1,length.out=nbins+1))

for(q in 1:2){
  for(j in 1:nbins){
    for(mm in 1:nm){
      cond<-qs[,q]>OMp[j,q]&qs[,q]<OMp[j+1,q]
      OMv[mm,q,j]<-mean(Yds[cond,mm],na.rm=T)
    }
  }
}
OMs<-apply(OMv,1:2,sd,na.rm=T)

ncols<-40
#colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),90)[ncols:1]
minsd<-0
maxsd<-max(OMs,na.rm=T)
coly<-ceiling(OMs/maxsd*ncols)

k<-0
for(q in 1:2){
  for(mm in 1:nm){
    k<-k+1
    x<-qs[,q]
    y<-Yds[,mm]
    plot(x,y,col="white",ylim=c(0,2),axes=F)
    axis(1,pretty(OMp[,q]),pretty(OMp[,q]),cex.axis=0.8,padj=-1.5)
    axis(2,cex.axis=0.8, las=1)
    abline(v=OMp[,q],col="#99999960")
    points(x,y,col=colsse[coly[mm,q]],pch=19,cex=0.8)
    x2<-(OMp[1:nbins,q]+OMp[2:(nbins+1),q])/2
    y2<-OMv[mm,q,]
    lines(x2,y2)
    if(q==1)mtext(MPnam[mm],3,line=0,cex=0.7,adj=-0.1)
    if(mm==nm)mtext(Qnam[q],1,line=1.5,cex=0.75)
    mtext(plabs[k],3,adj=0.98,line=0.02,cex=0.7)
  }
}  
mtext("Yield relative to Fmsy_PI yield",2,outer=T,line=1,cex=0.85)
dev.off()


#####################
#Ref point estimation
#####################

tiff("C:/Users/Bill/Documents/Box Sync/Data limited/Carruthers manuscript/Figures/S4 Reference point.tiff",width=8,height=10,units="in",res=600)

dep<-c(1,2,3)
lessDep<-c(10,11,12)
Meths<-c("DCAC_s", "DCAC2_s", "DCAC2_adj", "DBSRA_s", "DBSRA_adj", "SPMSY_s", "SPMSY_adj", "Fmsy_PI", "MC_s", "MC_adj")
mind<-c(1,2,4,6)

par(mfcol=c(4,3))
for (st in 1:3){
  for (nm in mind) {
    MSEobj<-MSEs[[dep[st]]]
    x<-MSEobj@OM$Depletion/MSEobj@OM$BMSY_B0
    y<-MSEobj@TAC[,nm,1]/MSEobj@OM$MSY
    plot(x, y, ylim=c(0,3), xlim=c(0,3), pch=16, col="red", las=1, bty="n", xlab="B / BMSY", ylab="TAC / MSY")
    
    MSEobj<-MSEs[[lessDep[st]]]
    x<-MSEobj@OM$Depletion/MSEobj@OM$BMSY_B0
    y<-MSEobj@TAC[,nm,1]/MSEobj@OM$MSY
    points(x, y, col="blue")
    abline(h=1, lty=3)
    lines(c(1,1), c(0,1), lty=3)
    mtext(Meths[nm], side=3, adj=1)
    if(nm==1) mtext(design$Var1[st], font=2, line=2)
  }
}
dev.off()






