####ACTUAL rosetta ON SUBSTE GENES
library(purrr)
library(tidyverse)
library(dplyr)
library(data.table)
library(infer)
library(R.ROSETTA)
library(VisuNet)
library(import)

setwd('C://Users//giris//Downloads//SIV-Virus//')
import::from(count_data_preprocessing.R,count_data_pro)
import::from(after_cv.R,after_cv)
import::from(discretize_data.R,discretize_quartile)
import::from(annotation_rmcfs_results.R,annotation_man)
import::from(cluster_rules.R,heatmap_rules)

#read metadata for the samples
data<-fread('C://Users//giris//Downloads//target_renamed.csv')

#read Count data for gene expressions
count_data<-fread('C://Users//giris//Downloads//count_matrix_renamed.csv')

timepoint<-'D3'

file_path<-'C://Users//giris//Downloads//test//rmcfs_D3.zip'

count_data<-count_data_pro(count_data,timepoint)

after_cv<-after_cv(timepoint)

#metadata for quartile range
quartiles_df <-after_cv

# Column names for columns to be discretized
cut_columns <- rownames(after_cv)


quality<-vector()


annot<-annotation_man(file_path)


for(i in 1:(ceiling((nrow(annot)/10)))){
  
 if(i==ceiling((nrow(annot)/10)))( sig_genes<-annot$ENSEMBL[1:nrow(annot)])  
 else(sig_genes<-annot$ENSEMBL[1:(i*10)])
 
  #sig_genes<-annot$ENSEMBL[1:n]
  
  #data to be discretized
  data_df <- subset(count_data,select=sig_genes)

  # Discretize the data based on quartiles
  
  discretized_data_df <- discretize_quartile(data_df, cut_columns, quartiles_df)
  
  
  
  discretized_data_df <- discretized_data_df %>% rownames_to_column(., var = 'V1')
  discretized_data_df <- discretized_data_df %>% left_join(dplyr::select(data,V1,protectionStatus), by="V1")
  discretized_data_df<- discretized_data_df %>% column_to_rownames(.,var = 'V1')
  discretized_data_df<-discretized_data_df[discretized_data_df$protectionStatus=='Prot'| discretized_data_df$protectionStatus=='NonProt', ]
  
  
  ind<-match(colnames(discretized_data_df),annot$ENSEMBL)
  new_gname<-annot$SYMBOL[ind]
  print(dim(discretized_data_df))
  print(new_gname)
  colnames(discretized_data_df)<-new_gname
  ros<-rosetta(discretized_data_df,reducer = 'Johnson',clroc='Prot',discrete = TRUE,roc=TRUE,JohnsonParam = list(Modulo=TRUE, BRT=TRUE, BRTprec=0.1, Precompute=FALSE, Approximate=TRUE, Fraction=0.99))
  
  if(ros$quality$ROC.AUC.MEAN > max(quality)){
    ros_max<-ros;
    discretized_data_df_max<-discretized_data_df;
    recal_max<-recalculateRules(discretized_data_df_max,ros_max$main,discrete = TRUE)
  }
  quality<-append(quality,ros$quality$ROC.AUC.MEAN)
  print(quality)
 
  
  if(i==(ceiling((nrow(annot)/10))))(plot(quality,type='o',xlab='NumberofFeatures/10',ylab='Accuracy',main='Boosting Performance'))
  
  
  
}

recal<-recalculateRules(discretized_data_df,ros$main,discrete = TRUE)


#heatmap
heatmap_rules(recal_max,discretized_data_df_max)
heatmap_rules(recal,discretized_data_df)

#model accuracy
ros_max$quality
ros$quality

#rules
View(viewRules(ros_max$main))
View(viewRules(ros$main))

#visunet
visunet(ros_max$main)
visunet(ros$main)

#save
#saveRDS(ros,'proj/snic2022-22-995/SIV_immune_response/ros.RData')

