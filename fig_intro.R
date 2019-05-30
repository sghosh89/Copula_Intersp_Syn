library(VineCopula)

pdf("./Results/fig_intro_A.pdf",height=8,width=8)
op<-par(mar=c(5,5,1,2),mgp=c(3,0.5,0))
cc<-BiCopSim(N=50,family=3,par=2)
plot(cc[,1],cc[,2],xlim=c(0,1),ylim=c(0,1),xlab="u",ylab="v",
     asp=1,col=rgb(0,0,0,0.3),pch=16,cex=1.5,cex.lab=3,cex.axis=2)
rect(0,0,1,1)

abline(a=0,b=1,lty="dashed")
abline(a=0,b=-1,col="darkgrey")
abline(a=2,b=-1,col="darkgrey")
#segments(0,0.5,0.5,0)
segments(0,1,1,0)
segments(0,2/3,2/3,0,col="red")
#segments(1,0.5,0.5,1)
segments(1,1/3,1/3,1,col="blue")
legend("topleft","A",bty="n",cex=3,x.intersp=-0.1,y.intersp = -0.4)
par(op)
dev.off()


pdf("./Results/fig_intro_BCDE.pdf",height=8,width=8)
op<-par(mfrow=c(2,2),mar=c(6,8,1,2),mgp=c(3,0.5,0))
cc1<-BiCopSim(N=50,family=3,par=3)
cc2<-BiCopSim(N=50,family=3,par=2.5)
scc1<-BiCopSim(N=50,family=13,par=3)
scc2<-BiCopSim(N=50,family=13,par=2.5)

plot(cc1[,1],cc1[,2],xlim=c(0,1),ylim=c(0,1),xlab="-Temp.",ylab=expression(Sp[i]),
     col=rgb(0,0,0,0.3),pch=16,cex=1.5,cex.lab=3,cex.axis=2)
segments(0,1,1,0)
legend("topleft","B",bty="n",cex=3,x.intersp=0.1,y.intersp = -0.3)

plot(cc2[,1],cc2[,2],xlim=c(0,1),ylim=c(0,1),xlab=expression(Sp[i]),ylab=expression(Sp[j!=i]),
     col=rgb(0,0,0,0.3),pch=16,cex=1.5,cex.lab=3,cex.axis=2)
segments(0,1,1,0)
legend("topleft","C",bty="n",cex=3,x.intersp=0.1,y.intersp = -0.3)

plot(scc1[,1],scc1[,2],xlim=c(0,1),ylim=c(0,1),xlab="-Temp.",ylab=expression(Sp[i]),
     col=rgb(0,0,0,0.3),pch=16,cex=1.5,cex.lab=3,cex.axis=2)
segments(0,1,1,0)
legend("topleft","D",bty="n",cex=3,x.intersp=0.1,y.intersp = -0.3)

plot(scc2[,1],scc2[,2],xlim=c(0,1),ylim=c(0,1),xlab=expression(Sp[i]),ylab=expression(Sp[j!=i]),
     col=rgb(0,0,0,0.3),pch=16,cex=1.5,cex.lab=3,cex.axis=2)
segments(0,1,1,0)
legend("topleft","E",bty="n",cex=3,x.intersp=0.1,y.intersp = -0.3)

par(op)
dev.off()




