source("tailsignif.R")

ub<-.5
numpts<-30
spcors<-seq(from=-.9,to=.9,by=.1)
numsims<-500
ploton<-TRUE
res<-tailsignif(ub=ub,numpts=numpts,spcors=spcors,numsims=numsims,ploton=TRUE)

dim(res)

