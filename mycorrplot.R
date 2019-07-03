# This plotter function is to visualize a matrix
# Input:
#     z : a matrix
#     posnI_ind : a matrix containing the row and col indices of z for which z indicates indep values
#     posnN_ind : a matrix containing the row and col indices of z for which z indicates sig. neg cor values
#     colrange : a vector containing min and max value of the color range
#     nvar : an integer showing number of external biotic/abiotic variables affecting cross-sp. synchrony
# Output :
#     A matrix plot comes with green spot for the cell with -vely correlated pair
library(corrplot)

mycorrplot<-function(z,posnI_ind,posnN_ind,colrange,nvar){
  
  col1 <- colorRampPalette(c("blue","white","red")) 
  
  z[is.na(z)]<-mean(colrange)
  diag(z)[1]<-colrange[1] # just to ensure that plot always have specific colorbar range even 
  diag(z)[2]<-colrange[2]     # though all entries are either +ve or -ve
  
  dimrs<-nrow(z)-nvar
  zrs<- z[1:dimrs,]# square matrix reshaped into rectangular display
  
  corrplot(zrs,is.corr = F,col=col1(100),method="color",addgrid.col = "black",
           diag=F,bg = "white",tl.cex=4.5,tl.col = "black",
           cl.cex = 4,cl.lim = colrange,
           cl.align.text = "l",cl.ratio = 0.2)
  
  
  
  # colorize as black for diagonal indices
  Dg <- matrix(NA,nrow(zrs),nrow(zrs))
  diag(Dg)<- 1 
  
  corrplot(Dg, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",
           bg = "transparent", tl.col = "transparent",col="black",method="color")
  
  #colorize as yellow for indep posn indices
  if(dim(posnI_ind)[1]!=0){
    
    I <- matrix(NA,nrow(z),ncol(z))
    I[posnI_ind]<- 1 
    I<- I[1:dimrs,]# square matrix reshaped into rectangular display
    
    allNA<-is.na(I)
    if(all(allNA)==F){
      corrplot(I, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black", 
               bg = "transparent", tl.col = "transparent",col="yellow",method="color")
    }
    
  } 
  
  #colorize as green for -ve correlated (siginificantly) posn indices
  if(dim(posnN_ind)[1]!=0){
    
    N <- matrix(NA,nrow(z),ncol(z))
    N[posnN_ind]<- -1 
    N<- N[1:dimrs,]# square matrix reshaped into rectangular display
    
    allNA<-is.na(N)
    
    if(all(allNA)==F){
      corrplot(N, cl.pos = "n", na.label = " ", add = T,addgrid.col = "transparent",  
               bg = "transparent", tl.col = "transparent",p.mat = N,sig.level = -2,col="transparent",
               pch=20,pch.col="green",pch.cex = 5,number.cex = 2)
      
      
      # If two sp. are negatively correlated then color cells as full green 
      sp_negcor<-which(N[,1:dimrs]==-1,arr.ind=T)
      
      if(nrow(sp_negcor)!=0){
        Nsp<-matrix(NA,nrow(z),ncol(z))
        Nsp[sp_negcor]<- -1 
        Nsp<- Nsp[1:dimrs,]# square matrix reshaped into rectangular display
        
        corrplot(Nsp, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black", 
                 bg = "transparent", tl.col = "transparent",col="green",method="color")
      }
    }
    
  }
  
}

