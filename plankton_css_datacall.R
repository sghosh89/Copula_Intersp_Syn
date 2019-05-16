#----------------------------------------------------------------------------------------------------------------
# This function makes data shaped as data[[loc]][[sp]] for cross species synchrony calculation data for plankton
# Input : 
#       d0 : a raw data matrix given for plankton

plankton_css_datacall<-function(d0){
  
  d0[d0==-999]<-NA
  
  name.plankton<-c('Calanus_I-IV','Para-Pseudocalanus_spp','Temora_longicornis','Acartia_spp_(unidentified)','Centropages_typicus',
                   'Oithona_spp','Echinoderm_larvae','Calanus_finmarchicus','Calanus_helgolandicus','Metridia_lucens','Decapoda_larvae_(Total)',
                   'Euphausiacea_Total','Thalassiosira_spp','Rhizosolenia_styliformis','Ceratium_fusus','Ceratium_furca','Ceratium_tripos',
                   'Ceratium_macroceros','Rhizosolenia_alata_alata','Pseudocalanus_elongatus_Adult','Nitzschia_delicatissima','Nitzschia_seriata')
  
  name.plankton.selected<-c('Ceratium_fusus','Ceratium_furca','Ceratium_tripos','Ceratium_macroceros',
                            'Calanus_finmarchicus') # last one is zooplankton
  
  ds_c<-vector("list",26)
  names(ds_c)<-paste("loc",c(1:26),sep="")
  
  for(il in 1:length(ds_c)){
    ds1<-vector("list", 22)
    names(ds1)<-name.plankton
    
    Years<-1958:2013
    
    for(is in 1:length(ds1)){
      ind<-il+((is-1)*26)
      temp<-d0[ind,]
      ds1[[is]]<-data.frame(Year=Years,Dat=temp)
    }
    
    ds1<-ds1[name.plankton.selected]
    
    ds_c[[il]]<-ds1
    
  }
  
  return(ds_c)
}