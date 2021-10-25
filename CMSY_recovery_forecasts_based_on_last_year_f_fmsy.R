

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Estimate recovery proportion over 10 years
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Set working directory

#upload data
#datafile =  "Overall_stocks_outputs.csv"
datafile =  "Out_Ukr_ID.csv"#"Out_March062019_O_Stocks_ID_19_Med_single.csv"
#DATA = read.csv("All_stocks_q.csv") 
DATA = read.csv(datafile) 

#DATA = subset(DATA, Subregion != "Wide ranging")
dir.create("output",showWarnings =FALSE)
dir.create("output/Files",showWarnings =FALSE)
dir.create("output/Figs",showWarnings =FALSE)

# n permutations
mcs = 1000
ci = rep(c(1,rep(0,(round(mcs/50,1)-1))),mcs)

# get Info
Region = DATA$Region
Subregion = DATA$Subregion

longterm = FALSE

# number of projection years 
nyrs =  15
if(longterm==TRUE) nyrs =15 
# start year
styr = 2017
endyr = styr+nyrs
# Period
years = styr:endyr

# Define F as a fraction  of Fmsy for target F: i.e. F = Fx*Fmsy 
Fx = c(0.5,0.6,0.8,0.95)


#########################################################################
# Adjusted Schaefer model for B/Bmsy < 0.5
# after Froese et al. 2016 (CMYS, FaF)
#########################################################################

# test
#B_Bmsy = 0.2; Bmsy = 500; Fmsy= 0.2; Fcur = 0.21; Fx = 0.8
#SPF(B_Bmsy,Bmsy,Fmsy,Fcur,Fx,nyrs)

SPF = function(B_Bmsy,Bmsy,Fmsy,Fcur,Fx,nyrs){ # Ft = Target F
  
  Pmsy = mat.or.vec(length(B_Bmsy),(nyrs+2))
  B = mat.or.vec(length(B_Bmsy),(nyrs+1))
  C =  mat.or.vec(length(B_Bmsy),(nyrs+1))
  Pmsy[,1] = B_Bmsy
  
  FxRed = Fx*Fcur/Fmsy
  for(y in 1:(nyrs+1))
  {
    Fproj = ifelse(FxRed <  rep(0.95,length(Fmsy)) & Pmsy[,y]<0.5,2*Pmsy[,y]*Fmsy*FxRed,Fmsy*FxRed)
    
   # for (ifx in 1:length(FxRed)){
    #  if(FxRed [ifx]== 0.5) Fproj = ifelse(Pmsy[,y][ifx]<0.5,Fmsy*0.01,Fproj)
    #}  
    
    if(y>5){Ft =Fproj } else {Ft=ifelse(Pmsy[,y]<0.2 & Fcur>0.5*Fmsy,2*Pmsy[,y]*Fmsy,Fcur)}
    
    #Ft = ifelse(Ft>Fcur,Fcur,Ft)
    Pmsy[,y] = ifelse(Pmsy[,y]<0.01,0.01,Pmsy[,y])
    Pmsy[,y+1] = ifelse(Pmsy[,y]<0.5,
    Pmsy[,y]+4*Fmsy*Pmsy[,y]^2*(1-0.5*Pmsy[,y])-Ft*Pmsy[,y],                    
    Pmsy[,y]+2*Fmsy*Pmsy[,y]*(1-0.5*Pmsy[,y])-Ft*Pmsy[,y]) 
    B[,y] = Pmsy[,y]*Bmsy 
    C[,y] = Ft*Pmsy[,y]*Bmsy 
  }  
  return(list(BtoBmsy=Pmsy[,1:(nyrs+1)],Catch = C,Bt = B))
}
#########################################################################

# 1) Run By Region
Regs = c(rev(levels(Region)))
nReg = length(Regs)
regs = c("BS")

for(j in 1:nReg){

  
# subset by region
dat=DATA
if(j<3) dat = subset(DATA,Region==Regs[j])  
B = dat$B
Bmsy = dat$Bmsy
lcl.Bmsy =  dat$lcl.11 
sd.Bmsy = (log(Bmsy)-log(lcl.Bmsy))/1.96 #log-space

B_Bmsy = B/Bmsy

Fmsy = dat$F_msy
Fcur = dat$F
lcl.Fmsy= dat$lcl.8 #log-space
sd.Fmsy = (log(Fmsy)-log(lcl.Fmsy))/1.96 #log-space
B_Bmsy = B/Bmsy

MSY = dat$MSY
lcl.MSY = dat$lcl.10

bk = B_Bmsy/2
lcl.bk = dat$lcl.13/2
ucl.bk = dat$ucl.13/2


# recover sd from B_Bmsy...assuming logistic as a function of b/k
lg.bk = log(bk/(1-bk))
lg.bk.lc = log(lcl.bk/(1-lcl.bk))
lg.bk.uc = log(ucl.bk/(1-ucl.bk))
sd.bk = (lg.bk.uc - lg.bk.lc) /(1.96*2)


# Start projecting using Monte-Carlo
cat(paste0("\n","Analysing: ",Regs[j]," Stocks","\n"))
cat(paste0("\n","Projection period: ",nyrs," years until ",years[nyrs+1],"\n"))
cat(paste0("\n","Running ",mcs," Monte-Carlo simulations...","\n"))


# Monte-Carlo Sims over all 3 fishing scenarios
b_bmsy.F = list()
b_bmsy05.F = list()
b_bmsy02.F = list()

C.F = list()
B.F = list()
dB.F = list()

for(f in 1:4){
cat(paste0("\n","Running Scenario F = ",Fx[f],"*Fmsy","\n"))
Fproj = Fx[f]*Fmsy


    
# All combined
Pmsy.MCS = mat.or.vec(mcs,(nyrs+1))
Pmsy05.MCS = mat.or.vec(mcs,(nyrs+1))
Pmsy02.MCS = mat.or.vec(mcs,(nyrs+1))
C.MCS = mat.or.vec(mcs,(nyrs+1))
B.MCS = mat.or.vec(mcs,(nyrs+1))
dB.MCS = mat.or.vec(mcs,(nyrs+1))
#if(j == 3) BBmsy2030.MCS = Catch2030.MCS = NULL

for(i in 1:mcs){
  
if(i==1) cat(paste("\n","|"))
 
if(ci[i]==1) cat("*")
if(i==mcs) cat(paste("|","\n"))
  
rFmsy = rlnorm(length(Fmsy),log(Fmsy),sd.Fmsy) 
rBmsy =rlnorm(length(Bmsy),log(Bmsy),sd.Bmsy) 

rbk = 1/(1+exp(-rnorm(length(sd.bk),lg.bk,sd.bk))) 
rB_Bmsy = rbk*2


get_SPF = SPF(rB_Bmsy,rBmsy,rFmsy,Fcur,Fx[f],nyrs)
rPmsy = get_SPF$BtoBmsy  
rCatch = get_SPF$Catch
rBt = get_SPF$Bt
Pmsy.bin = ifelse(rPmsy<1,0,1)
Pmsy05 = ifelse(rPmsy<0.5,0,1) #ifelse(rPmsy<0.5 & rPmsy >= 0.2,0,1)
Pmsy02 = ifelse(rPmsy<0.2,0,1)

Pmsy.MCS[i,] =  apply(Pmsy.bin,2,sum)/apply(Pmsy.bin,2,length)
Pmsy05.MCS[i,] =  apply(Pmsy05,2,sum)/apply(Pmsy05,2,length)
Pmsy02.MCS[i,] =  apply(Pmsy02,2,sum)/apply(Pmsy02,2,length)
C.MCS[i,] =  apply(rCatch,2,sum,na.rm=T)
B.MCS[i,] =  apply(rBt,2,sum)#/apply(rBt,2,sum)
dB.MCS[i,] = apply(rBt,2,sum)/sum(rBt[,1])

}

mu = apply(Pmsy.MCS,2,mean,na.rm=T)*100
CIs = apply(Pmsy.MCS,2,quantile,c(0.025,0.975),na.rm=T)*100
b_bmsy.F[[f]] = rbind(mu,CIs)


mu = apply(Pmsy05.MCS,2,mean,na.rm=T)*100
CIs = apply(Pmsy05.MCS,2,quantile,c(0.025,0.975),na.rm=T)*100
b_bmsy05.F[[f]] = rbind(mu,CIs)

mu = apply(Pmsy02.MCS,2,mean,na.rm=T)*100
CIs = apply(Pmsy02.MCS,2,quantile,c(0.025,0.975),na.rm=T)*100
b_bmsy02.F[[f]] = rbind(mu,CIs)


mu = apply(C.MCS,2,mean,na.rm=T)
CIs = apply(C.MCS,2,quantile,c(0.025,0.975),na.rm=T)
C.F[[f]] =  rbind(mu,CIs)

mu = apply(B.MCS,2,mean)
CIs = apply(B.MCS,2,quantile,c(0.025,0.975),na.rm=T)
B.F[[f]] =  rbind(mu,CIs)

mu = apply(dB.MCS,2,mean)
CIs = apply(dB.MCS,2,quantile,c(0.025,0.975),na.rm=T)
dB.F[[f]] =  rbind(mu,CIs)

}

ColName = c("Year",paste0(c("mu","lcl","ucl"),rep(Fx,each=3)))

Catches = data.frame(years,t(C.F[[1]]))
for(f in 2:4) Catches=cbind(Catches,t(C.F[[f]]))
colnames(Catches) = ColName

BtoBmsy = data.frame(years,t(b_bmsy.F[[1]]))
for(f in 2:4) BtoBmsy =cbind(BtoBmsy ,t(b_bmsy.F[[f]]))
colnames(BtoBmsy) = ColName

BtoBmsy05 = data.frame(years,t(b_bmsy05.F[[1]]))
for(f in 2:4) BtoBmsy05 =cbind(BtoBmsy05 ,t(b_bmsy05.F[[f]]))
colnames(BtoBmsy05) = ColName

BtoBmsy02 = data.frame(years,t(b_bmsy02.F[[1]]))
for(f in 2:4) BtoBmsy02 =cbind(BtoBmsy02 ,t(b_bmsy02.F[[f]]))
colnames(BtoBmsy02) = ColName

Biomass = data.frame(years,t(B.F[[1]]))
for(f in 2:4) Biomass =cbind(Biomass ,t(B.F[[f]]))
colnames(Biomass) = ColName

Profitability = data.frame(years,t(dB.F[[1]]))
for(f in 2:4) Profitability  =cbind(Profitability  ,t(dB.F[[f]]))
colnames(Profitability) = ColName

#CHECK BIOMASS
# dump results
write.csv(Catches,paste0("output/Files/Catch_",regs[j],endyr,".csv"))
write.csv(BtoBmsy,paste0("output/Files/BtoBmsy_",regs[j],endyr,".csv"))
write.csv(Biomass,paste0("output/Files/Biomass_",regs[j],endyr,".csv"))
write.csv(Profitability,paste0("output/Files/Profitability_",regs[j],endyr,".csv"))
write.csv(BtoBmsy05,paste0("output/Files/BtoBmsy05_",regs[j],endyr,".csv"))
write.csv(BtoBmsy02,paste0("output/Files/BtoBmsy02_",regs[j],endyr,".csv"))
} # End of loop

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Plots by region 
# BtoBmsy

xlab="Projection Years"

Par = list(mfrow=c(1,1),mai=c(0.2,0.15,.2,0),omi = c(0.2,0.25,0,0) + 0.1,mgp=c(3,0.7,0), tck = -0.02,cex=0.8)

png(file = paste0("output/Figs/BtoBmsy",endyr,".png"), width = 12, height = 4.3, 
    res = 500, units = "in")
par(Par)


for(j in 1:length(Regs)){
  dat=DATA
  if(j<3) dat = subset(DATA,Region==Regs[j])  
  
  d = read.csv(paste(paste0("output/Files/BtoBmsy_",regs[j],endyr,".csv")))[,-1]
  years = d[,1]
  sc = seq(2,11,3) # select cols
  
  ylab= "B/Bmsy"
  plot(x = years,y = seq(0,100,by = 100/(length(years)-1)),ylim=c(0,100),xlim=c(years[2],(years[nyrs]+2.3)),xaxs="i",yaxs="i",frame=FALSE,xlab="",ylab="",type="n",lwd=2,main=paste(Regs[j]), yaxt="n")
  
  if(j==1){axis(2,labels = T)} else {axis(2,labels = F)}
  axis(2,seq(-1e6,1e6,2e6))
  
  lines(c(0,years[nyrs]+1),rep(0,2))
  
  # Add colors
  rect(0, -10,years[nyrs]+1,20,col="darkred",border = NA)
  rect(0, 20, years[nyrs]+1, 40,col="red",border = NA)
  rect(0, 40, years[nyrs]+1, 60,col="orange",border = NA)
  rect(0, 60, years[nyrs]+1, 80,col="yellow",border = NA)
  rect(0, 80, years[nyrs]+1, 110,col="green",border = NA)
  legend("topleft",paste("N =",(nrow(dat))),bty="n",cex=1, x.intersp = -0.5)
  lines(c(0,years[nyrs]+1),rep(99.8,2),col="darkgreen",lwd=3)
  
  grey.scale = rev(c(0.6,0.7,0.8,0.9))
  
  for(f in 1:4) polygon(c(years,rev(years)),c(d[,sc[f]+1],rev(d[,sc[f]+2])),col=adjustcolor(grey(grey.scale[f] ,1),alpha=0.3),border = NA)
  
  for(f in 1:4){
    lines(years,d[,sc[f]],lwd=2,lty=f)
    text(years[nyrs]+1.2, d[nyrs*1.0,sc[f]],paste(Fx[f]),cex=1)  
  }
}
mtext(paste(xlab), side=1, outer=T, at=0.5,line=1,cex=0.9)
mtext(paste(ylab), side=2, outer=T, at=0.5,line=1,cex=0.9)

dev.off()

# CATCHES
ylab= "Catch ('000 t)"
Par = list(mfrow=c(1,1),mai=c(0.2,0.28,.2,0),omi = c(0.2,0.25,0,0) + 0.1,mgp=c(3,0.7,0), tck = -0.02,cex=0.8)
png(file = paste0("output/Figs/Catch",endyr,".png"), width = 12, height = 4.3, 
    res = 500, units = "in")
par(Par)

for(j in 1:length(Regs)){
  dat=DATA
  if(j<3) dat = subset(DATA,Region==Regs[j])  
  MSY = dat$MSY
  lcl.MSY = dat$lcl.10
  
  d = read.csv(paste(paste0("output/Files/Catch_",regs[j],endyr,".csv")))[,-1]
  years = d[,1]
  sc = seq(2,11,3) # select cols

  # CATCH
plot(years,years,ylim=c(0,sum(MSY)*1.2),xlim=c(years[2],years[nyrs]+1),xaxs="i",yaxs="i",frame=FALSE,xlab="Projection Years",ylab="",type="n",lwd=2,main=paste(Regs[j]))
axis(1,seq(0,1e6,1e6))
axis(2,seq(-1e6,1e6,2e6))

lines(c(0,years[nyrs]+1),rep(sum(lcl.MSY),2),lty=2)
text(years[1]+2.3,sum(lcl.MSY)*1.03,"lcl MSY")
lines(c(0,years[nyrs]+1),rep(sum(MSY),2),lty=2)
text(years[1]+2.1,sum(MSY)*1.03,"MSY")
lines(c(0,years[nyrs]+1),rep(0,2),lwd=1.2)

legend("topleft",paste("N =",(nrow(dat))),bty="n",cex=1, x.intersp = -0.5)
for(f in 4:1) polygon(c(years,rev(years)),c(d[,sc[f]+1],rev(d[,sc[f]+2])),col=grey(rev(grey.scale)[f] ,0.8),border = NA)

for(f in 1:4){
  lines(years,d[,sc[f]],lwd=2,lty=f)
  text(years[6]+0.3, d[6,sc[f]]*ifelse(f<3,0.95,1.05),paste(Fx[f]),cex=1.)  
  
}
} 
mtext(paste(xlab), side=1, outer=T, at=0.5,line=1,cex=0.9)
mtext(paste(ylab), side=2, outer=T, at=0.5,line=1,cex=0.9)

dev.off()

# Profitability #B/Bstartyear

Par = list(mfrow=c(1,1),mai=c(0.2,0.15,.2,0),omi = c(0.2,0.25,0,0) + 0.1,mgp=c(3,0.7,0), tck = -0.02,cex=0.8)
png(file = paste0("output/Figs/Profitability",endyr,".png"), width = 12, height = 4.3, 
    res = 500, units = "in")
par(Par)

for(j in 1:length(Regs)){
  dat=DATA
  if(j<3) dat = subset(DATA,Region==Regs[j])  
  d = read.csv(paste(paste0("output/Files/Profitability_",regs[j],endyr,".csv")))[,-1]
  years = d[,1]
  sc = seq(2,11,3) # select cols
  
  ylab=expression('B /  B'[forecast_initial_year]*' (%)')
# B/Bstartyear
plot(years,years,ylim=range(80,300),xlim=c(years[2],years[nyrs]+2.3),xaxs="i",yaxs="i",frame=FALSE,axes=F,ylab="",type="n",lwd=2,main=paste(Regs[j]))
if(j==1){axis(2,labels = T)} else {axis(2,labels = F)}
axis(1,labels = T)
axis(2,seq(-1e6,1e6,2e6))
lines(c(0,years[nyrs]+1),rep(80,2))


legend("topleft",paste("N =",(nrow(dat))),bty="n",cex=1, x.intersp = -0.5)
for(f in 1:4) polygon(c(years,rev(years)),c(d[,sc[f]+1]*100,rev(d[,sc[f]+2]*100)),col=grey(grey.scale[f] ,1),border = NA)

for(f in 1:4){
  lines(years,d[,sc[f]]*100,lwd=2,lty=f)
  #text(years[nyrs]+0., d[nyrs,sc[f]]*100*ifelse(f>2,0.95,1.03),paste(Fx[f],"Fmsy"))  
  text(years[nyrs]+1.7, d[nyrs,sc[f]]*100,paste(Fx[f]),cex=1)
}
lines(c(0,years[nyrs+1]),rep(100,2),lty=2)
}

mtext(paste(xlab), side=1, outer=T, at=0.5,line=1,cex=0.9)
mtext(ylab, side=2, outer=T, at=0.5,line=1,cex=0.9)

dev.off()


# Deterministic projections

#upload data
dat = read.csv(datafile) 

B = dat$B
Bmsy = dat$Bmsy
B_Bmsy = B/Bmsy
Fmsy = dat$F_msy
Fcur = dat$F
MSY = dat$MSY


Par = list(mfrow=c(1,1),mar = c(5, 5, 2, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.5, 0.1),mex=0.8, tck = -0.02,cex=0.6)
png(file = paste0("output/Files/checkTop20.png"), width = 4, height = 4, 
    res = 500, units = "in")
par(Par)
# run all three scenarios
for(f in 4:4){
 
 
get_SPF = SPF(B_Bmsy,Bmsy,Fmsy,Fcur,Fx[f],nyrs)
Catches = get_SPF$Catch 
dC = Catches[,2:ncol(Catches)]-Catches[,1:(ncol(Catches)-1)]
plot(years,years,type="n",ylim=range(dC),xlim=c(2017,2022),ylab=expression(Delta*"Catch (mt)"))
for(i in 1:nrow(Catches)){ lines(years[-1],dC[i,])}

# Find top ten
max.dC = dat[order(-dC[,5]),][1:15,]
legend("topright",c("Top 20 Stocks",paste(1:20,max.dC$Stock,", B/Bmsy =",round(max.dC$B_Bmsy,2),", F/Fmsy =",round(max.dC$F_msy,2))),cex=0.8,bty="n")

}

dev.off()



### DEPLETION
        
      
Par = list(mfrow=c(1,1),mai=c(0.2,0.15,.2,0),omi = c(0.2,0.25,0,0) + 0.1,mgp=c(3,0.7,0), tck = -0.02,cex=0.8)
png(file = paste0("output/Figs/Depleted_",endyr,".png"), width = 12, height = 4.3, 
    res = 500, units = "in")
par(Par)

# Depleted
for(j in 1:length(Regs)){
  dat=DATA
  if(j<3) dat = subset(DATA,Region==Regs[j])  
  
  d = read.csv(paste(paste0("output/Files/BtoBmsy05_",regs[j],endyr,".csv")))[,-1]
  years = d[,1]
  ylab="Depleted B/Bmsy"
  d[,-1] =  100-d[,-1] 
  sc = seq(2,11,3) # select cols
  plot(years,years,ylim=c(0,65),xlim=c(years[2],years[nyrs]+2.3),xaxs="i",yaxs="i",frame=FALSE,axes=F,ylab="",type="n",lwd=2,main=paste(Regs[j]))
  
  if(j==1){axis(2,labels = T)} else {axis(2,labels = F)}
  #axis(2,labels = T)
  axis(1,labels = T)
  axis(2,seq(-1e6,1e6,2e6))
  lines(c(0,years[nyrs]+1),rep(0,2))
  
  # Add colors
  
  #legend("topleft",paste("N =",(nrow(dat))),bty="n",cex=1, x.intersp = -0.5)
  
  grey.scale = (c(0.6,0.7,0.8,0.9))
  
  for(f in 1:4) polygon(c(years,rev(years)),c(d[,sc[f]+1],rev(d[,sc[f]+2])),col=grey(grey.scale[f] ,1),border = NA)
  
  for(f in 1:4){
    lines(years,d[,sc[f]],lwd=2,lty=f)
    #text(years[nyrs]+0.,d[nyrs*1.0,sc[f]]+1,paste(Fx[f],"Fmsy"))  
    text(years[nyrs]+1.7, d[nyrs,sc[f]],paste(Fx[f]),cex=1)
  }
}

mtext(paste(xlab), side=1, outer=T, at=0.5,line=1,cex=0.9)
mtext(ylab, side=2, outer=T, at=0.5,line=1,cex=0.9)

dev.off()


# Collapsed
Par = list(mfrow=c(1,1),mai=c(0.2,0.15,.2,0),omi = c(0.2,0.25,0,0) + 0.1,mgp=c(3,0.7,0), tck = -0.02,cex=0.8)
png(file = paste0("output/Figs/Collapsed_",endyr,".png"), width = 12, height = 4.3, 
    res = 500, units = "in")
par(Par)



for(j in 1:length(Regs)){
  dat=DATA
  if(j<3) dat = subset(DATA,Region==Regs[j])  
  
  d = read.csv(paste(paste0("output/Files/BtoBmsy02_",regs[j],endyr,".csv")))[,-1]
  years = d[,1]
  ylab="Collapsed B/Bmsy"
  d[,-1] =  100-d[,-1] 
  sc = seq(2,11,3) # select cols
  plot(years,years,ylim=c(0,50),xlim=c(years[2],years[nyrs]+2.3),xaxs="i",yaxs="i",frame=FALSE,axes=F,ylab="",type="n",lwd=2,main=paste(Regs[j]))
  
  if(j==1){axis(2,labels = T)} else {axis(2,labels = F)}
  #axis(2,labels = T)
  axis(1,labels = T)
  axis(2,seq(-1e6,1e6,2e6))
  lines(c(0,years[nyrs]+1),rep(0,2))
  
  # Add colors
  
  #legend("topleft",paste("N =",(nrow(dat))),bty="n",cex=1, x.intersp = -0.5)
  
  grey.scale = (c(0.6,0.7,0.8,0.9))
  
  for(f in 1:4) polygon(c(years,rev(years)),c(d[,sc[f]+1],rev(d[,sc[f]+2])),col=grey(grey.scale[f] ,1),border = NA)
  
  for(f in 1:4){
    lines(years,d[,sc[f]],lwd=2,lty=f)
    #text(years[nyrs]+0.,d[nyrs*1.0,sc[f]]+1,paste(Fx[f],"Fmsy"))  
    text(years[nyrs]+1.7, d[nyrs,sc[f]],paste(Fx[f]),cex=1)
  }
}

mtext(paste(xlab), side=1, outer=T, at=0.5,line=1,cex=0.9)
mtext(ylab, side=2, outer=T, at=0.5,line=1,cex=0.9)

dev.off()

