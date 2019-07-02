library(copula)
library(mvtnorm)
source("CopulaFunctions_flexible.R")

#' A function for looking at statistical significance of a measure of 
#' asymmetry of tail dependence, against a normal-copula null
#' hypothesis. The measure is corl-coru for positively associated 
#' variables and the absolute value of this quantity, computed after
#' flipping one of the variables, for negatively associated variables.
#' 
#' @param ub The measure of asymmetry of tail dependence is cor_{0,ub}-
#' cor_{1-ub,1} for positively associated variables, and the absolute 
#' value of this quantity, computed after flipping one of the variables, 
#' for negatively associated variables
#' @param numpts This many points from the null generated for each run
#' @param spcors The distribution of the statistic for these values of
#' the Spearman correlation is computed
#' @param numsims For each value of spcors, this many runs done
#' @param xy An n by 2 matrix of values to be added as points to the 
#' plot. Ignored if ploton is FALSE. These are the values of the 
#' statistics for real data, for comparison. The first column (x)
#' is actual Spearman correlations, and the second is the asymmetry
#' statistic (corl-coru for positive values of x and the absolute value
#' of that for negative values of x).
#' @param ploton TRUE if you want a plot on the default plotting device,
#' a character string if you want to save a pdf with that name, FALSE
#' for no plot
#' 

tailsignif<-function(ub,numpts,spcors,numsims,ploton=TRUE)
{
  #a small amount of error checking
  if (!is.numeric(ub) || length(ub)!=1 || ub>.5 || ub<=0)
  {
    stop("Error in tailsignif: bad value for ub")  
  }
  if (any(spcors<=-1) || any(spcors>=1))
  {
    stop("Error in tailsignif: bad value for spcors")
  }
  
  #get the numeric results
  allres<-matrix(NA,numsims,length(spcors))
  dumnorm<-normalCopula(.1,2)
  for (ct_spcors in 1:length(spcors))
  {
    #find the parameter for which the normal copula has this Spearman
    thisparm<-iRho(dumnorm,spcors[ct_spcors])
    thissig<-matrix(c(1,thisparm,thisparm,1),2,2)
    
    for (ct_sims in 1:numsims)
    {
      #generate numpts data from the normal copula
      d<-rmvnorm(numpts,mean=c(0,0),sig=thissig)
      d<-pobs(d) #note are doing it this way instead of using normalCopula and rCopula
      #because real data will also have had pobs applied to it
      
      #compute the statistic and store
      if (spcors[ct_spcors]>=0)
      {
        allres[ct_sims,ct_spcors]<-Corbds(d[,1],d[,2],0,ub)-Corbds(d[,1],d[,2],1-ub,1)
      }
      if (spcors[ct_spcors]<0)
      {
        allres[ct_sims,ct_spcors]<-abs(Corbds(d[,1],-d[,2]+1,0,ub)-Corbds(d[,1],-d[,2]+1,1-ub,1))
      }
    }
  }
  
  #do the plot if desired
  if (is.character(ploton))
  {
    pdf(file=paste0(ploton,".pdf"))
  }
  qtl<-apply(FUN=quantile,X=allres,MARGIN=2,prob=c(.025,.25,.5,.75,.975))
  plot(c(-1,1),c(0,0),ylim=c(-1,1),xlab="Spearman",ylab="Statistic",type='l',col='red')
  lines(c(0,0),c(-1,1),type="l",col="red")
  lines(spcors[spcors>=0],qtl[3,spcors>=0],type='l')
  lines(spcors[spcors>=0],qtl[2,spcors>=0],type='l',lty="dashed")
  lines(spcors[spcors>=0],qtl[4,spcors>=0],type='l',lty="dashed")
  lines(spcors[spcors>=0],qtl[1,spcors>=0],type='l',lty="dotted")
  lines(spcors[spcors>=0],qtl[5,spcors>=0],type='l',lty="dotted")
  lines(spcors[spcors<0],qtl[3,spcors<0],type='l')
  lines(spcors[spcors<0],qtl[2,spcors<0],type='l',lty="dashed")
  lines(spcors[spcors<0],qtl[4,spcors<0],type='l',lty="dashed")
  lines(spcors[spcors<0],qtl[1,spcors<0],type='l',lty="dotted")
  lines(spcors[spcors<0],qtl[5,spcors<0],type='l',lty="dotted")
  if (is.character(ploton))
  {
    dev.off()
  }
  
  return(allres)
}