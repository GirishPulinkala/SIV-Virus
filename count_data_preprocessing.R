


count_data_pro<-function(count_data,timepoint){
  
  if(timepoint=='D0')(count_data<-subset(count_data,select=c(V1,grep("W0_D0|_[0-9]+_D0",colnames(count_data)))))
  if(timepoint=='D3')(count_data<-subset(count_data,select=c(V1,grep("W0_D3|_[0-9]+_D3",colnames(count_data)))))
  if(timepoint=='D7')(count_data<-subset(count_data,select=c(V1,grep("W0_D7|_[0-9]+_D7",colnames(count_data)))))
  if(timepoint=='D126')(count_data<-subset(count_data,select=c(V1,grep("W18_D0|_[0-9]+_D98|_[0-9]+_D126",colnames(count_data)))))
  if(timepoint=='D129')(count_data<-subset(count_data,select=c(V1,grep("W18_D3|_[0-9]+_D101|_[0-9]+_D129",colnames(count_data)))))
  if(timepoint=='D133')(count_data<-subset(count_data,select=c(V1,grep("W18_D7|_[0-9]+_D105|_[0-9]+_D133",colnames(count_data)))))
  if(timepoint=='PreChal')(count_data<-subset(count_data,select=c(V1,grep("W88_D0|_[0-9]+_D462|_[0-9]+_D455|_[0-9]+_D538|_[0-9]+_D448|_[0-9]+_D540",colnames(count_data)))))
  

  count_data<-count_data %>% column_to_rownames(.,var='V1')
  
  #logarithmic conversion to base 10
  count_data<-log(count_data,10)
  
  #convert infinitely negative values to zero
  invisible(lapply(names(count_data),function(.name) set(count_data, which(is.infinite(count_data[[.name]])), j = .name,value =0)))
  
  
  count_data<-as.data.frame(t(count_data))
  
  return(count_data)
  
  
}


