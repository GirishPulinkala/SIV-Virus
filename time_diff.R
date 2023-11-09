library(import)
setwd('C://Users//giris//Downloads//SIV-Virus//')
import::from(data_pre_td.R,data_pre_td)
library(data.table)
library(dplyr)
library(tidyverse)
#read metadata for the samples
data<-fread('C://Users//giris//Downloads//target_renamed.csv')



#read Count data for gene expressions
count_data<-fread('C://Users//giris//Downloads//count_matrix_renamed.csv')



D0<-subset(count_data,select=c(V1,grep("W0_D0|_[0-9]+_D0",colnames(count_data))))
D3<-subset(count_data,select=c(V1,grep("W0_D3|_[0-9]+_D3",colnames(count_data))))
D7<-subset(count_data,select=c(V1,grep("W0_D7|_[0-9]+_D7",colnames(count_data))))
D126<-subset(count_data,select=c(V1,grep("W18_D0|_[0-9]+_D98|_[0-9]+_D126",colnames(count_data))))
D129<-subset(count_data,select=c(V1,grep("W18_D3|_[0-9]+_D101|_[0-9]+_D129",colnames(count_data))))
D133<-subset(count_data,select=c(V1,grep("W18_D7|_[0-9]+_D105|_[0-9]+_D133",colnames(count_data))))
prechal<-subset(count_data,select=c(V1,grep("W88_D0|_[0-9]+_D462|_[0-9]+_D455|_[0-9]+_D538|_[0-9]+_D448|_[0-9]+_D540",colnames(count_data))))





time_points<-list('D0'=D0,'D3'=D3,'D7'=D7,'D126'=D126,'D129'=D129,'D133'=D133,'PreChal'=prechal)


for(i in 1:length(time_points)){
  
  for(j in (i+1):length(time_points)){
    s1<-names(time_points[i])
    s2<-names(time_points[j])
    print(s1)
    print(s2)
    
    df1<-data_pre_td(as.data.frame(time_points[[i]]),data)
    df2<-data_pre_td(as.data.frame(time_points[[j]]),data)
    
    
    
    df_diff <- merge(df1,df2, by = "AnimalID", suffix = c(paste0('_',s1), paste0('_',s2)))
    
    res <- df_diff[, grepl(s1, colnames(df_diff))] - df_diff[, grepl(s2, colnames(df_diff))]
    
    colnames(res) <- paste(colnames(df_diff[, grepl(s1, colnames(df_diff))]),s2, sep = "-")
    
    res<-cbind('AnimalID'=df_diff[,1],res)
    
    # to make consolidated matrix of all timepoint differentials
    #if(j==2)(fin_res<-res)
    #else(fin_res<-inner_join(fin_res,res,by='AnimalID'))
    
    
    #to write time differentials for al combinations
    write.csv(res,paste('C://Users//giris//Downloads//SIV-Virus//time_diff//','count_data',s1,s2,'.csv',sep='_'))
    
  }  
  
  if(i==length(time_points)-1)(break)
  
}





write.csv(fin_res,'/proj/snic2022-22-995/SIV_immune_response/girish_codes/time_diff_counts.csv')








