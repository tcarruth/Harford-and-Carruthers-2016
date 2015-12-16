
# ====================== Simulation testing of DDCAC, EDCAC, MCD, DBSRA, AC and perfect information OFL ======================
# 
# Harford and Carruthers 2015. Testing a simple modification to MacCall's DCAC
#
# December 26th 2015
#
# ============================================================================================================================
Drive<-"C"
rm(list=ls())
graphics.off()
library(DLMtool)
library(wordcloud)
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]])
sfInit(parallel=T,cpus=detectCores())


EDCAC<-function (x, DLM_data, reps = 100) 
{
  dependencies = "DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot <- DLM_data@AvC[x] * DLM_data@t[x]
  Mdb <- trlnorm(reps, DLM_data@Mort[x], DLM_data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, DLM_data@FMSY_M[x], DLM_data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, DLM_data@Dt[x], DLM_data@CV_Dt[x])
  BMSY_K <- rbeta(reps, alphaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]), 
                  betaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]))
  dcac<-C_tot/(DLM_data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  TAC<-dcac*Bt_K/BMSY_K
  TACfilter(TAC)
}
class(EDCAC)<-"DLM_output"
environment(EDCAC)<-asNamespace('DLMtool')
sfExport("EDCAC")

HDCAC<-function (x, DLM_data, reps = 100) 
{
  dependencies = "DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot <- DLM_data@AvC[x] * DLM_data@t[x]
  Mdb <- trlnorm(reps, DLM_data@Mort[x], DLM_data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, DLM_data@FMSY_M[x], DLM_data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, DLM_data@Dt[x], DLM_data@CV_Dt[x])
  BMSY_K <- rbeta(reps, alphaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]), 
                  betaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]))
  dcac<-C_tot/(DLM_data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  ddcac<-dcac*Bt_K/BMSY_K
  TAC<-dcac
  TAC[Bt_K<BMSY_K]<-ddcac[Bt_K<BMSY_K]
  TACfilter(TAC)
}
class(HDCAC)<-"DLM_output"
environment(HDCAC)<-asNamespace('DLMtool')
sfExport("HDCAC")

AC<-function(x,DLM_data,reps=100){
  C_tot <- DLM_data@AvC[x] * DLM_data@t[x] 
  trlnorm(reps,C_tot, DLM_data@CV_Cat[x]) 
}
class(AC)<-"DLM_output"
environment(AC)<-asNamespace('DLMtool')
sfExport("AC")

OFL_PI<-function(x,DLM_data,reps=100){
  trlnorm(reps, DLM_data@OM$A[x] * (1 - exp(-DLM_data@OM$FMSY[x])),0.01)
}
class(OFL_PI)<-"DLM_output"
environment(OFL_PI)<-asNamespace('DLMtool')
sfExport("OFL_PI")

sfExportAll()
set.seed(1) 

Stocks<-c("Bluefin_tuna_WAtl","Snapper","Rockfish")
Obs<-c("Perfect_Info","Precise_Unbiased","Imprecise_Biased")
Dep<-matrix(c(0.05,0.3,0.3,0.7),nrow=2)

design<-expand.grid(Stocks,Obs,1:2)
nOMs<-nrow(design)

MSEs<-new('list')
OMs<-new('list')
Meths<-c("OFL_PI","AC","MCD","DCAC","EDCAC","HDCAC","DBSRA")
nsim<-200
proyears<-40
interval<-3
reps<-1


for(i in 1:nOMs){
  print(paste("=========  ", i," - ", design[i,1],design[i,2],design[i,3], "  ================"))
  assign('temp_stock',get(as.character(design[i,1])))
  temp_stock@D<-Dep[design[i,3],]
  OMs[[i]]<-new('OM',temp_stock,Generic_fleet,get(as.character(design[i,2])))
  OMs[[i]]@Name<-paste(as.character(design[i,1]),c("PI","DataRich","DataPoor")[as.integer(design[i,2])],c("5_30","30_70")[design[i,3]],sep="_")
  MSEs[[i]]<-runMSE(OMs[[i]],Meths,nsim,proyears,interval,pstar=0.5,maxF=0.8,timelimit=20,reps)
  Tplot(MSEs[[i]])
}

save(OMs,file="C:/EDCAC/OMs3.R")
save(MSEs,file="C:/EDCAC/MSEs3.R")


# Data Analysis ==============================================================================================
Drive<-"F"
load(file=paste(Drive,":/EDCAC/OMs3.R",sep=""))
load(file=paste(Drive,":/EDCAC/MSEs3.R",sep=""))

Stocks<-c("Bluefin tuna","Snapper","Rockfish")
#Stocks<-c("Herring","Bluefin tuna","Rockfish")
Obs<-c("Perfect_Info","Precise-unbiased","Imprecise-biased")
design<-expand.grid(Stocks,Obs,c("Depleted","Less depleted"))

makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

tradeoffplot2<-function(x,y,xlab,ylab,labs,cex=1,vl,hl,coly,px=T,py=T){
  
  levs<-(0:100)/5
  plot(NA,xlim=c(-0.0,1.),ylim=c(-0.0,1),xlab=xlab,ylab=ylab,axes=F)
  
  if(px)axis(1,levs,levs)
  if(!px)axis(1,levs,rep(NA,length(levs)))
  if(py)axis(2,levs,levs)
  if(!py)axis(2,levs,rep(NA,length(levs)))
  
  abline(v=levs,col="#99999920",lwd=2)
  abline(h=levs,col="#99999920",lwd=2)
  #text(x,y,labs,font=2,col=coly,cex=1)
  textplot(x,y,labs,col=c(rep('black',2),rep('blue',5)),cex=0.85,new=FALSE,show.lines=T,font=2,xlim=c(0,1),ylim=c(0,1))
  
}

tradeoffplot3<-function(x,y,xlab,ylab,labs,cex=1,vl,hl,coly,px=T,py=T){
  
  levs<-(0:100)/10
  plot(x,y,col="white",ylim=c(0.6,1.05),xlim=c(0.5,1.05),xlab=xlab,ylab=ylab,axes=F)
  
  if(px)axis(1,levs,levs)
  if(!px)axis(1,levs,rep(NA,length(levs)))
  if(py)axis(2,levs,levs)
  if(!py)axis(2,levs,rep(NA,length(levs)))
  
  abline(v=levs,col="#99999920",lwd=2)
  abline(h=levs,col="#99999920",lwd=2)
  #text(x,y,labs,font=2,col=coly,cex=1)
  textplot(x,y,labs,col=c('black',rep('blue',5)),cex=0.85,new=FALSE,show.lines=T,font=2,xlim=c(0,1),ylim=c(0,1))
  
}



coly<-makeTransparent(rep(c("black","grey"),10))

LTYvSTY<-function(MSEobj,coly,px=T,py=T){
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
  ystart<-1:10
  RefYd<-MSEobj@OM$RefY
  for(mm in 1:MSEobj@nMPs){
    LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)
    STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)
  }
  tradeoffplot2(STY,LTY,"","",labs<-MSEobj@MPs[1:MSEobj@nMPs],coly=coly,px=px,py=py)
}

B50vVY<-function(MSEobj,coly,px=T,py=T){
  VY<-rep(NA,MSEobj@nMPs)
  B10<-rep(NA,MSEobj@nMPs)
  RefYd<-MSEobj@OM$RefY
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  for(mm in 1:MSEobj@nMPs){
    AAVY<-apply(((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)/apply(MSEobj@C[,mm,y2],1,mean)
    VY[mm]<-round(sum(AAVY<0.1,na.rm=T)/MSEobj@nsim,3)
    B10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,])),3)
  }
  tradeoffplot3(B10,VY,xlab="",ylab="",MSEobj@MPs[1:MSEobj@nMPs],coly=coly,px=px,py=py)
} 




jpeg(paste(Drive,":/EDCAC/Figures/Fig 3 STY vs LTY.jpg",sep=""),width=6,height=8,units="in",res=600)

m<-cbind(matrix(rep(rep(1:3,each=3),3),nrow=9),matrix(rep(rep(4:6,each=3),3),nrow=9))
m<-t(cbind(m,rep(13,9),matrix(rep(rep(7:9,each=3),3),nrow=9),matrix(rep(rep(10:12,each=3),3),nrow=9)))
layout(m)
par(mai = c(0.15, 0.15, 0.02, 0.02),omi=c(0.5,0.5,0.6,0.5))
plabs<-c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)")

ys<-c(4,7,13,16)
xs<-c(7:9,16:18)

j<-0
for(i in c(4:9,13:18)){
  j<-j+1
  MSEs[[i]]@MPs[1]<-"OFL"
  LTYvSTY(MSEs[[i]],py=i%in%ys,px=i%in%xs,coly=coly)
  if(i==4)mtext("- More depleted -",3,adj=-0.1,line=0.5,cex=0.8)
  if(i==13)mtext("- Less depleted -",3,adj=-0.1,line=0.5,cex=0.8)
  if(i<10)mtext(Stocks[i-3],line=2.5,cex=0.8)
  if(design[i,1]=="Rockfish")mtext(design[i,2],4,las=3,line=0.3,cex=0.8)
  mtext(plabs[j],3,adj=0.98,line=-0.1,cex=0.75)
}  
mtext("Short-term yield (STY)",1,line=1.6,outer=T,cex=0.8)
mtext("Long-term yield (LTY)",2,line=1.3,outer=T,cex=0.8)

dev.off()


jpeg(paste(Drive,":/EDCAC/Figures/Fig 4 B50 vs VY.jpg",sep=""),width=6,height=8,units="in",res=600)

m<-cbind(matrix(rep(rep(1:3,each=3),3),nrow=9),matrix(rep(rep(4:6,each=3),3),nrow=9))
m<-t(cbind(m,rep(13,9),matrix(rep(rep(7:9,each=3),3),nrow=9),matrix(rep(rep(10:12,each=3),3),nrow=9)))
layout(m)
par(mai = c(0.15, 0.15, 0.02, 0.02),omi=c(0.5,0.5,0.6,0.5))
plabs<-c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)")

ys<-c(4,7,13,16)
xs<-c(7:9,16:18)
j<-0
for(i in c(4:9,13:18)){
  j<-j+1
  MSEs[[i]]@MPs[1]<-"OFL"
  B50vVY(Sub(MSEs[[i]],MSEs[[i]]@MPs[c(1,3:7)]),py=i%in%ys,px=i%in%xs,coly=coly)
  if(i==4)mtext("- More depleted -",3,adj=-0.1,line=0.5,cex=0.8)
  if(i==13)mtext("- Less depleted -",3,adj=-0.1,line=0.5,cex=0.8)
  if(i<10)mtext(Stocks[i-3],line=2.5,cex=0.8)
  if(design[i,1]=="Rockfish")mtext(design[i,2],4,las=3,line=0.3,cex=0.8)
  mtext(plabs[j],3,adj=0.98,line=-0.1,cex=0.75)
}  
mtext("Biomass levels (B50)",1,line=1.6,outer=T,cex=0.8)
mtext("Stability in yield (SY)",2,line=1.3,outer=T,cex=0.8)

dev.off()

