library(copula)
#set.seed(101)
#clayton data
cc5<-claytonCopula(5)
dcc5<-rCopula(250,cc5)

#survival clayton data
sc5<-rotCopula(cc5)
dsc5<-(-dcc5+1)

#make marginals normal
dcc5_nm<-qnorm(dcc5)
dsc5_nm<-qnorm(dsc5)


p1<-cor(dcc5_nm[,1],dcc5_nm[,2],method="pearson")
p2<-cor(dsc5_nm[,1],dsc5_nm[,2],method="pearson")

s1<-cor(dcc5_nm[,1],dcc5_nm[,2],method="spearman")
s2<-cor(dsc5_nm[,1],dsc5_nm[,2],method="spearman")

pdf("./Results/pedagog_fig.pdf",height=9,width=6)
op<-par(mfrow=c(3,2),mar=c(5.2,5,2.5,1),pty="s")

rg<-max(abs(dcc5_nm))
plot(dcc5_nm[,1],dcc5_nm[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(-rg,rg),xlim=c(-rg,rg),
     xlab=expression(x[t]),ylab=expression(y[t]),cex.lab=2,cex.axis=2)
legend(x=-5,y=4,"A",bty="n",cex=3)
text(3,-2,paste0("P=",round(p1,2)),adj=c(1,0),cex=1.5)
text(3,-2.5,paste0("S=",round(s1,2)),adj=c(1,0),cex=1.5)

rg<-max(abs(dsc5_nm))
plot(dsc5_nm[,1],dsc5_nm[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(-rg,rg),xlim=c(-rg,rg),
     xlab=expression(x[t]),ylab=expression(y[t]),cex.lab=2,cex.axis=2)
legend(x=-5,y=4,"B",bty="n",cex=3)
text(3,-2,paste0("P=",round(p2,2)),adj=c(1,0),cex=1.5)
text(3,-2.5,paste0("S=",round(s2,2)),adj=c(1,0),cex=1.5)

plot(dcc5[,1],dcc5[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend(x=-0.3,y=1.1,"C",bty="n",cex=3)

plot(dsc5[,1],dsc5[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend(x=-0.3,y=1.1,"D",bty="n",cex=3)

plot(dcc5[,1],dcc5[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend(x=-0.3,y=1.1,"E",bty="n",cex=3)
abline(a=0.4,b=-1)
abline(a=0.8,b=-1)
arrows(x0=0,y0=0,x1=0.8,y1=0,code=3,lwd=1.5,length=0.2)
arrows(x0=0,y0=0,x1=0,y1=0.4,code=3,lwd=1.5,length=0.2)
text(0.6,0.01,expression(2*u[b]),adj=c(1,0),cex=1.5)
text(0.12,0.16,expression(2*l[b]),adj=c(1,0),cex=1.5)
text(0.12,0.4,expression(u[t]+v[t]==2*l[b]),cex=1.5,srt=-45)
text(0.65,0.27,expression(u[t]+v[t]==2*u[b]),cex=1.5,srt=-45)

plot(dcc5[,1],dcc5[,2],type="p",col=rgb(0,0,0,0.2),pch=16,ylim=c(0,1),xlim=c(0,1),
     xlab=expression(u[t]),ylab=expression(v[t]),cex.lab=2,cex.axis=2)
legend(x=-0.3,y=1.1,"F",bty="n",cex=3)
abline(a=0.4,b=-1)
abline(a=1.6,b=-1)
abline(a=0,b=-1)
abline(a=2,b=-1)
arrows(x0=0,y0=0,x1=0.4,y1=0,code=3,lwd=1.5,length=0.2)
text(0.22,0.01,expression(2*b),adj=c(1,0),cex=1.5)
arrows(x0=1,y0=1,x1=1,y1=0.6,code=3,lwd=1.5,length=0.2)
text(0.98,0.8,expression(2*b),adj=c(1,0),cex=1.5)

text(0.12,0.4,expression(u[t]+v[t]==2*b),cex=1.5,srt=-45)
text(-0.2,0.18,expression(u[t]+v[t]==0),cex=1.5,srt=-40,xpd=NA)
text(0.64,0.83,expression(u[t]+v[t]==2-2*b),cex=1.5,srt=-45)
text(0.82,1.18,expression(u[t]+v[t]==2),cex=1.5,srt=-45,xpd=NA)

par(op)
dev.off()
