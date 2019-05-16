#----------------------------------------------------------------------------------------------------------------
# This function makes data shaped as data[[loc]][[sp]] for cross species synchrony calculation data for aphid
# Input : 
#       d0 : a raw data matrix given for count, firstflight or flightduration data set

aphid_css_datacall<-function(d0){
  
  d0[d0==-999]<-NA
  
  name.aphid<-c("Apple-grass aphid","Bird cherry-oat aphid","Black bean aphid", "Blackberry-cereal aphid","Blackcurrant-sowthistle aphid",
                "Corn leaf aphid","Currant-lettuce aphid","Damson-hop aphid","Grain aphid","Green spruce aphid",
                "Leaf-curling plum aphid","Mealy cabbage aphid","Mealy plum aphid","Pea aphid","Peach-potato aphid",
                "Potato aphid","Rose-grain aphid","Shallot aphid","Sycamore aphid","Willow-carrot aphid")
  
  ds_c<-vector("list",dim(d0)[1])
  names(ds_c)<-paste("loc",c(1:dim(d0)[1]),sep="")
  
  for(il in 1:length(ds_c)){
    ds1<-vector("list", 20)
    names(ds1)<-name.aphid
    
    Years<-1976:2010
    
    for(is in 1:length(ds1)){
      temp<-d0[il,((is-1)*35+1):(is*35)]
      ds1[[is]]<-data.frame(Year=Years,Dat=temp)
    }
    
    ds_c[[il]]<-ds1
  }
  
  return(ds_c)
}
#----------------------------------------------------------



