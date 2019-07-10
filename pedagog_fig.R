library(copula)

#clayton data
cc5<-claytonCopula(5)
dcc5<-rCopula(250,cc5)

#survival clayton data
sc5<-rotCopula(cc5)
dsc5<-(-dcc5+1)

#make marginals normal
dcc5_nm<-qnorm(dcc5)
dsc5_nm<-qnorm(dsc5)

cor(dcc5[,1],dcc5[,2],method="pearson")
cor(dsc5[,1],dsc5[,2],method="pearson")

cor(dcc5[,1],dcc5[,2],method="spearman")
cor(dsc5[,1],dsc5[,2],method="spearman")

cor(dcc5_nm[,1],dcc5_nm[,2],method="spearman")
cor(dsc5_nm[,1],dsc5_nm[,2],method="spearman")

pdf("./Results/pedagog_fig.pdf",height=9,width=6)
op<-par(mfrow=c(3,2),mar=c(5,5,2,1),pty="s")

rg<-max(abs(dcc5_nm))
plot(dcc5_nm[,1],dcc5_nm[,2],type="p",col=rgb(0,0,0,0.3),pch=16,ylim=c(-rg,rg),xlim=c(-rg,rg),
     xlab=expression(x[t]),ylab=expression(y[t]),cex.lab=2,cex.axis=2)
legend("topleft","A",bty="n",cex=2)

rg<-max(abs(dsc5_nm))
plot(dsc5_nm[,1],dsc5_nm[,2],type="p",col=rgb(0,0,0,0.3),pch=16,ylim=c(-rg,rg),xlim=c(-rg,rg),
     xlab=expression(x[t]),ylab=expression(y[t]),cex.lab=2,cex.axis=2)
legend("topleft","B",bty="n",cex=2)

plot(dcc5[,1],dcc5[,2],type="p",col=rgb(0,0,0,0.3),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend("topleft","C",bty="n",cex=2)

plot(dsc5[,1],dsc5[,2],type="p",col=rgb(0,0,0,0.3),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend("topleft","D",bty="n",cex=2)

plot(NA,NA,type="p",col=rgb(0,0,0,0.3),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab="?",ylab="?",cex.lab=2,cex.axis=2)
legend("topleft","E",bty="n",cex=2)

par(op)
dev.off()
