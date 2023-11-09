
data_pre_td<-function(D,data){
  
  D <- D %>% column_to_rownames(., var = 'V1')
  
  #logarithmic conversion to base 10
  D<-log(D,10)
  
  #convert infinitely negative values to zero
  invisible(lapply(names(D),function(.name) set(D, which(is.infinite(D[[.name]])), j = .name,value =0)))
  
  
  #Transpose
  df<-t(D)
  
  #replace NAN values with zero
  df[is.na(df)]<-0
  
  df<-as.data.frame(df)
  #merge data
  df <- df %>% rownames_to_column(., var = 'V1')
  
  #join Animal ID from the metadata
  df<- df %>% left_join(select(data,V1,AnimalID),by='V1')
  df<- df %>% column_to_rownames(.,'AnimalID')
  
  #remove column with samples named w.r.t week  
  df<-subset(df,select=-(V1))
  
  
  df<- df%>% rownames_to_column(.,'AnimalID')
  
  return(df)
  
}


