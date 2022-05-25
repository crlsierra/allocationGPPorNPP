library(SoilR)
library(expm)
library(RColorBrewer)
library(ggplot2)

setwd("~/Publications/allocationGPPorNPP/")

pal=palette.colors(4, palette="Classic Tableau",alpha=0.3)

#ages=read.csv("CO2_age.csv")
ages=read.csv2("CO2_Age_healthy_trees.csv")

Stem_chamber=ages[ages$organ=="stem",]
Stem_core=ages[ages$organ=="Stem-core",]
Coarse_roots=ages[ages$organ=="Coarse roots",]
Fine_roots=ages[ages$organ=="fine roots",]
Roots=ages[ages$organ=="Coarse roots" | ages$organ=="fine roots",]

pdf("ageRespiredC.pdf")
par(mar=c(4,4,1,1))
boxplot(list(Roots$CO2_age, Stem_chamber$CO2_age, Stem_core$CO2_age), col=pal, frame.plot=FALSE,xlab="",ylim=c(0,20),
        ylab=expression(paste("Age of respired ", CO[2], " (yr)")),pch=20,
        names=c("Roots", "Stem-chambers", "Stem-cores"))
text(x=c(1,2,3,4), y=20, labels=c(nrow(Coarse_roots)+nrow(Fine_roots),nrow(Stem_chamber),nrow(Stem_core)))
abline(h=1,lty=2)
dev.off()


#### Emanuel model
tm=seq(0,600,by=0.1)
pools_emanuel=c("Non-woody tree parts", "Woody tree parts", "Ground vegetation", "Detritus/decomposers", "Active soil carbon")

xssG=matrix(c(37, 452, 69, 81, 1121),5,5, byrow=TRUE)
FlG=matrix(c(-(21+31+25), 0, 0, 0, 0,
            31, -(14+15+2), 0, 0, 0,
            0, 0, -(18+12+6),0, 0,
            21, 15, 12, -(45+3), 0,
            0, 2, 6, 3, -(11)), 5,5, byrow=TRUE)
B_G=FlG/xssG
u_G=matrix(c(77, 0, 36, 0, 0), ncol=1)

x_G=-1*solve(B_G)%*%u_G
Rt=function(t, B, u){diag(colSums(-B))%*%expm(t*B)%*%u}

z=colSums(B_G)
R=-z*x_G
TT_G=transitTime(A=B_G, u=u_G, a=tm)
SA_G=systemAge(A=B_G, u=u_G, a=tm)
# M_emanuel=t(sapply(tm,FUN=M, B=B_emanuel, u=u_emanuel))
# M1_emanuel=t(sapply(tm,FUN=M1, B=B_emanuel, u=u_emanuel))
Rt_G=t(sapply(tm, Rt, B=B_G, u=u_G))

xssN=matrix(c(18, 248, 35, 81, 1121),5,5, byrow=TRUE)
FlN=matrix(c(-(21+17), 0, 0, 0, 0,
             17, -(15+2), 0, 0, 0,
             0, 0, -(12+6),0, 0,
             21, 15, 12, -(45+3), 0,
             0, 2, 6, 3, -(11)), 5,5, byrow=TRUE)
B_N=FlN/xssN
u_N=matrix(c(38, 0, 18, 0, 0), ncol=1)

x_N=-1*solve(B_N)%*%u_N
TT_N=transitTime(A=B_N, u=u_N, a=tm)
SA_N=systemAge(A=B_N, u=u_N, a=tm)

plot(tm,TT_G$transitTimeDensity,type="l",xlim=c(0,50),xlab="Transit time (year)",ylab="Density", bty="n")
abline(v=TT_G$quantiles[2],lty=2)
lines(tm,TT_N$transitTimeDensity,col=2)
abline(v=TT_N$quantiles[2],lty=2,col=2)

pdf("logTT.pdf")
par(mar=c(4,4,0,0))
plot(tm,TT_G$transitTimeDensity,type="l", log="y",xlab="Transit time (year)",
     ylab="Proportion of respired C", yaxt="n", bty="n")
aty <- axTicks(2)
labels <- sapply(5:1,function(i)
  as.expression(bquote(10^ .(-i)))
)
axis(2,at=aty,labels=labels)
dev.off()

GPP=sum(u_G)
NPP=sum(u_N)

Ra=GPP-NPP

f_TN=NPP*TT_N$transitTimeDensity
f_TN[1]=Ra

f_TG=GPP*TT_G$transitTimeDensity

palB=brewer.pal(5,"Dark2")

pdf('GPPversusNPP_TT.pdf')
plot(tm,f_TG,type="l",col=4,xlim=c(0,20), xlab="Age of respired carbon (transit time) in yr", 
     ylab=expression(paste("Mass density (PgC y",r^-1,")")), bty="n")
lines(tm,f_TN,col=2)
abline(v=TT_G$quantiles[2],lty=2,col=4)
legend("topright", c("GPP-based model", "NPP-based model"),lty=1,col=c(4,2),bty="n")
dev.off()

pdf('TT_contributions.pdf')
plot(tm,f_TG,type="l",xlim=c(0,10), xlab="Age of respired carbon (transit time), in yr",lty=3, 
     ylab=expression(paste("Mass density (PgC y",r^-1,")")), bty="n")
lines(tm, rowSums(Rt_G[,1:3]), col=palB[1])
lines(tm, rowSums(Rt_G[,4:5]), col=palB[2])
legend("topright",c("Ecosystem respiration", "Autotrophic respiration", "Heterotrophic respiration"),lty=c(3,1,1),col=c(1,palB[1:2]),bty="n")
#matlines(tm,Rt_G,col=palB,lty=c(1,1,1,2,2))
#legend("topright",c("Ecosystem respiration", pools_emanuel),lty=c(3,1,1,1,2,2),col=c(1,palB),bty="n")
dev.off()

# Radiocarbon simulations
yr=seq(1851,2010,by=1/12)
AtmF14=BoundFc(map = Graven2017[,1:2], lag = 0, format = "Delta14C")
modG=GeneralModel_14(t=yr,A=B_G,ivList = as.numeric(x_G), initialValF = ConstFc(rep(0,5),"Delta14C"), 
                     inputFluxes = as.numeric(u_G), inputFc = AtmF14)
Rt14G=getF14R(modG)

modN=GeneralModel_14(t=yr,A=B_N,ivList = as.numeric(x_N), initialValF = ConstFc(rep(0,5),"Delta14C"), 
                     inputFluxes = as.numeric(u_N), inputFc = AtmF14, pass=TRUE)
Rh14N=getF14R(modN) # This only includes respiration from heterotrophic pools
# The radiocarbon signal of the autotrophic pools is equal to the atmosphere
C14spline=splinefun(x=Graven2017[,1], y=Graven2017[,2])
Ra14N=C14spline(yr)
Rh=rowSums(getReleaseFlux(modN))

Re14N=((Rh14N*Rh)+(Ra14N*Ra))/(Rh+Ra)

pdf('radiocarbon.pdf', encoding = 'WinAnsi.enc')
par(mar=c(4,4.5,1,1))
plot(yr,Rt14G, col=4, type="l", xlim=c(1950,2010), ylim=c(0,600), xlab="Calendar year", 
     ylab=expression(paste(Delta^14,"C (\u2030)")), bty="n")
#lines(yr,Rt14G, col=4)
lines(yr,Re14N,col=2)
legend("topright", c("GPP-based model", "NPP-based model"),lty=1,col=c(4,2),bty="n")
dev.off()


