# This function is written to generate matrix plot for non-parametric stat results for all locations
# Input : 
#   data_ln_all :        a list of length of all locations (output from NonParamStat.R for all locations)
#   resloc : location to save the results
#   nvar : an integer : number of factors which influence css between considered species
#   nvar_names : a vector of characters containing names of n-variables
#   tagon : logical (argument for vivj_matrix fn call)
#   npa_stats : a character tag to choose any one of 3 npa stats
#---------------------------
source("mycorrplot.R")
#---------------------------

NonParamStat_matrixplot<-function(data_ln_all,resloc,nvar,nvar_names,tagon,npa_stats){
  
  #------------------------------------
  # initialize
  CorlmCoru<-NA
  PlmPu<-NA
  D2umD2l<-NA
  summary_nLU_CorlmCoru<-NA
  summary_nLU_PlmPu<-NA
  summary_nLU_D2umD2l<-NA
  summary_LU_CorlmCoru<-NA
  summary_LU_PlmPu<-NA
  summary_LU_D2umD2l<-NA
  #-------------------------------------
  
  
  numloc<-length(data_ln_all)
  selected_loc<-names(data_ln_all)
  
  #--------------------------Spearman plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$spear
    indI<-data_ln_all[[i]]$posnI
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Spearman.pdf",sep=''),width=24, height=24)
    z<-tempo[[loc]]
    
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=data_ln_all[[loc]]$posnI,
               posnN_ind=data_ln_all[[loc]]$posnN,
               colrange=c(-cr,cr))
    segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
    segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
    dev.off()
  }
  
  #--------------------------Kendall plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$kend
    indI<-data_ln_all[[i]]$posnI
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Kendall.pdf",sep=''),width=24, height=24)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=data_ln_all[[loc]]$posnI,
               posnN_ind=data_ln_all[[loc]]$posnN,
               colrange=c(-cr,cr))
    segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
    segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
    dev.off()
  }
  
  #========================================= For cor npa stats ===============================================
  
  if(npa_stats=="cor"){
    #--------------------------Corl plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Corl
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Corl.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    #--------------------------Coru plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Coru
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Coru.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    #--------------------------Corl-Coru plot---------------------------
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Corl-data_ln_all[[i]]$Coru
      indI<-data_ln_all[[i]]$posnI
      indN<-data_ln_all[[i]]$posnN
      tempo[[i]][indI]<-NA
      diag(tempo[[i]])<-NA
    }
    
    CorlmCoru<-tempo
    names(CorlmCoru)<-selected_loc
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    summary_nLU_CorlmCoru<-matrix(NA,2,numloc) # to keep total count on CorlmCoru for +ve or -ve numbers
    colnames(summary_nLU_CorlmCoru)<-selected_loc
    rownames(summary_nLU_CorlmCoru)<-c("nL","nU")
    
    summary_LU_CorlmCoru<- summary_nLU_CorlmCoru  # to keep sum on CorlmCoru values only for +ve or -ve numbers
    rownames(summary_LU_CorlmCoru)<-c("L","U")
    
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Corl-Coru.pdf",sep=''),width=24, height=24)
      #op<-par(mar=c(5.1, 2, 4.1, 8), mgp=c(3, 1, 0), las=0)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      z[data_ln_all[[loc]]$posnN]<-NA
      dl2<-nrow(z)-nvar
      z1<-z[1:dl2,1:dl2] # only select sp-sp interaction matrix to calculate nL,nU,L,U
      nL<-sum(z1>0,na.rm = T)
      nU<-sum(z1<0,na.rm = T)
      L<-sum(z1[which(z1>0,arr.ind=T)])
      U<-sum(z1[which(z1<0,arr.ind=T)])
      summary_nLU_CorlmCoru[1,loc]<-nL
      summary_nLU_CorlmCoru[2,loc]<-nU
      summary_LU_CorlmCoru[1,loc]<-L
      summary_LU_CorlmCoru[2,loc]<-U
      if(tagon == T){
        mtext(paste0(selected_loc[loc]," : "),cex=6,side=1,adj=0.3,line=-1)
      }
      #mtext(paste0("nL =",nL,", nU =",nU),cex=5,side=1,adj=0.7)
      mtext((as.expression(bquote('N'['L']*' = '*.(nL)))),cex=6,side=1,adj=0.5,col="red")
      mtext((as.expression(bquote('N'['U']*' = '*.(nU)))),cex=6,side=1,adj=0.8,col="blue")
      #par(op)
      dev.off()
    }
    
    res<-list(CorlmCoru_all_ln_list=CorlmCoru,
         summary_nLU_CorlmCoru=summary_nLU_CorlmCoru,
         summary_LU_CorlmCoru=summary_LU_CorlmCoru)
  }
 
  #=========================================== For P npa stats ==========================================
  
  if(npa_stats=="P"){
    #--------------------------Pl plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Pl
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Pl.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    
    #--------------------------Pu plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Pu
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Pu.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    
    #--------------------------Pl-Pu plot---------------------------
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$Pl-data_ln_all[[i]]$Pu
      indI<-data_ln_all[[i]]$posnI
      indN<-data_ln_all[[i]]$posnN
      tempo[[i]][indI]<-NA
      diag(tempo[[i]])<-NA
    }
    
    PlmPu<-tempo
    names(PlmPu)<-selected_loc
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    summary_nLU_PlmPu<-matrix(NA,2,numloc)
    colnames(summary_nLU_PlmPu)<-selected_loc
    rownames(summary_nLU_PlmPu)<-c("nL","nU")
    
    summary_LU_PlmPu<- summary_nLU_PlmPu  # to keep sum on PlmPu values only for +ve or -ve numbers
    rownames(summary_LU_PlmPu)<-c("L","U")
    
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_Pl-Pu.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      z[data_ln_all[[loc]]$posnN]<-NA
      dl2<-nrow(z)-nvar
      z1<-z[1:dl2,1:dl2] # only select sp-sp interaction matrix to calculate nL,nU,L,U
      nL<-sum(z1>0,na.rm = T)
      nU<-sum(z1<0,na.rm = T)
      L<-sum(z1[which(z1>0,arr.ind=T)])
      U<-sum(z1[which(z1<0,arr.ind=T)])
      summary_nLU_PlmPu[1,loc]<-nL
      summary_nLU_PlmPu[2,loc]<-nU
      summary_LU_PlmPu[1,loc]<-L
      summary_LU_PlmPu[2,loc]<-U
      if(tagon == T){
        mtext(paste0(selected_loc[loc]," : "),cex=6,side=1,adj=0.3,line=-1)
      }
      mtext((as.expression(bquote('N'['L']*' = '*.(nL)))),cex=6,side=1,adj=0.5,col="red")
      mtext((as.expression(bquote('N'['U']*' = '*.(nU)))),cex=6,side=1,adj=0.8,col="blue")
      dev.off()
    }
    
    res<-list(PlmPu_all_ln_list=PlmPu,
         summary_nLU_PlmPu=summary_nLU_PlmPu,
         summary_LU_PlmPu=summary_LU_PlmPu)
  }
  
  #=========================================== For D2 npa stats ==========================================
  
  if(npa_stats=="D2"){
    
    #--------------------------D2l plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$D2l
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_D2l.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    #--------------------------D2u plot---------------------------
    
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$D2u
      indI<-data_ln_all[[i]]$posnI
      tempo[[i]][indI]<-NA
    }
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_D2u.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      dev.off()
    }
    
    #--------------------------D2u-D2l plot---------------------------
    tempo<-vector("list",numloc)
    for(i in 1:numloc){
      tempo[[i]]<-data_ln_all[[i]]$D2u-data_ln_all[[i]]$D2l
      indI<-data_ln_all[[i]]$posnI
      indN<-data_ln_all[[i]]$posnN
      tempo[[i]][indI]<-NA
      diag(tempo[[i]])<-NA
    }
    
    D2umD2l<-tempo
    names(D2umD2l)<-selected_loc
    
    minval<-min(unlist(lapply(tempo,min,na.rm=T)))
    maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
    
    cr<-max(abs(minval),abs(maxval))
    
    summary_nLU_D2umD2l<-matrix(NA,2,numloc)
    colnames(summary_nLU_D2umD2l)<-selected_loc
    rownames(summary_nLU_D2umD2l)<-c("nL","nU")
    
    summary_LU_D2umD2l<- summary_nLU_D2umD2l  # to keep sum on D2umD2l values only for +ve or -ve numbers
    rownames(summary_LU_D2umD2l)<-c("L","U")
    
    for(loc in 1:numloc){
      resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
      pdf(paste(resloc2,selected_loc[loc],file="_D2u-D2l.pdf",sep=''),width=24, height=24)
      z<-tempo[[loc]]
      if(nvar!=0){
        dl<-nrow(z)-nvar+1
        rownames(z)[c(dl:nrow(z))]<-nvar_names
      }
      colnames(z)<-rownames(z)
      mycorrplot(z=z,
                 posnI_ind=data_ln_all[[loc]]$posnI,
                 posnN_ind=data_ln_all[[loc]]$posnN,
                 colrange=c(-cr,cr))
      segments(x0=c(0.5, 0.5), y0=c(0.5+nvar, 0.5+nrow(z)), x1=c(0.5+nrow(z), 0.5+nrow(z)), y1=c(0.5+nvar,0.5+nrow(z)), lwd=6,col="green")
      segments(c(0.5,0.5+nrow(z)), rep(0.5+nvar,1), c(0.5,0.5+nrow(z)), rep(0.5+nrow(z),1), lwd=6, col="green")
      z[data_ln_all[[loc]]$posnN]<-NA
      dl2<-nrow(z)-nvar
      z1<-z[1:dl2,1:dl2] # only select sp-sp interaction matrix to calculate nL,nU,L,U
      nL<-sum(z1>0,na.rm = T)
      nU<-sum(z1<0,na.rm = T)
      L<-sum(z1[which(z1>0,arr.ind=T)])
      U<-sum(z1[which(z1<0,arr.ind=T)])
      summary_nLU_D2umD2l[1,loc]<-nL
      summary_nLU_D2umD2l[2,loc]<-nU
      summary_LU_D2umD2l[1,loc]<-L
      summary_LU_D2umD2l[2,loc]<-U
      if(tagon == T){
        mtext(paste0(selected_loc[loc]," : "),cex=6,side=1,adj=0.3,line=-1)
      }
      mtext((as.expression(bquote('N'['L']*' = '*.(nL)))),cex=6,side=1,adj=0.5,col="red")
      mtext((as.expression(bquote('N'['U']*' = '*.(nU)))),cex=6,side=1,adj=0.8,col="blue")
      dev.off()
    }
    
    res<-list(D2umD2l_all_ln_list=D2umD2l,
              summary_nLU_D2umD2l=summary_nLU_D2umD2l,
              summary_LU_D2umD2l=summary_LU_D2umD2l)
  }

  return(res)
}

