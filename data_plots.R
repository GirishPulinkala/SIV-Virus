rm=(list=ls()) #CLEAR THE ENVIRONMENT

library(gapminder)
library(tidyverse)
library(infer)
library(data.table)
library(treemap)
library(data.tree)
library(networkD3)
library(MVN)
library(dplyr)
library(org.Mmu.eg.db)
#read data target_renamed.csv to create metadata
data<-fread('C://Users//giris//Downloads//target_renamed.csv')

########Understanding Data#############

#create root node and its subnode
data$pathString <- paste("PS",
                            data$protectionStatus,
                            data$Group,
                            data$Sex,
                            data$AnimalID,
                            
                            sep = "/")

#set node
population <- as.Node(data)
#print(population,'population$Prot$T')
#plot(population)
#plot radial network
useRtreeList <- ToListExplicit(population, unname = TRUE)
radialNetwork( useRtreeList)
#plot simple network
acmeNetwork <- ToDataFrameNetwork(population, "name")
simpleNetwork(acmeNetwork[-3], fontSize = 12)

#barplot for number of indivduals in each timepont w.r.t protection status
tbl <- with(data, table(Time_Point,protectionStatus))
pal <- colorRampPalette(colors = c('lightblue',"lightgreen", "lightred"))(7)
b<-barplot(tbl, beside = TRUE, legend = TRUE, ylim = c(0,70),ylab='Number of Animals',col = pal)
text(b, tbl + 5, tbl, font=2)


#########################################
#read count matrix
count_data<-fread('C://Users//giris//Downloads//count_matrix_renamed.csv')
View(count_data)
k<-1

#Subsetting data based on timepoints for normalization and further steps
D0<-subset(count_data,select=c(V1,grep("W0_D0|_[0-9]+_D0",colnames(count_data))))
D3<-subset(count_data,select=c(V1,grep("W0_D3|_[0-9]+_D3",colnames(count_data))))
D7<-subset(count_data,select=c(V1,grep("W0_D7|_[0-9]+_D7",colnames(count_data))))
D126<-subset(count_data,select=c(V1,grep("W18_D0|_[0-9]+_D98|_[0-9]+_D126",colnames(count_data))))
D129<-subset(count_data,select=c(V1,grep("W18_D3|_[0-9]+_D101|_[0-9]+_D129",colnames(count_data))))
D133<-subset(count_data,select=c(V1,grep("W18_D7|_[0-9]+_D105|_[0-9]+_D133",colnames(count_data))))
prechal<-subset(count_data,select=c(V1,grep("W88_D0|_[0-9]+_D462|_[0-9]+_D455|_[0-9]+_D538|_[0-9]+_D448|_[0-9]+_D540",colnames(count_data))))

#D<-D129

D<-count_data

time_points<-list('D0'=D0,'D3'=D3,'D7'=D7,'D126'=D126,'D129'=D129,'D133'=D133,'PreChal'=prechal)

for(i in 1:length(time_points)){
  D<-as.data.frame(time_points[[i]])
  #print(D[1:2,1:2])
  t<-names(time_points[i])
  c_v<-vector()
  #print(t)  
  
  #CONVERT gene names as index or rownames
  D <- D %>% column_to_rownames(., var = 'V1')
  
  #logarithmic conversion to base 10
  D<-log(D,10)
  
  #convert infinitely negative values to zero
  invisible(lapply(names(D),function(.name) set(D, which(is.infinite(D[[.name]])), j = .name,value =0)))
  
  #remove rows with zero as read count for all animals
  D<-D[apply(D[,-1], 1, function(x) !all(x==0)),]
  
  #Transpose
  df<-t(D)
  
  #replace NAN values with zero
  df[is.na(df)]<-0
  
  df<-df[apply(df[,-1], 1, function(x) !all(x==0)),]
  
  df<-as.data.frame(df)
  #merge data
  df <- df %>% rownames_to_column(., var = 'V1')
  data_D<-data[grep(t, data$Time_Point), ]
  
  df_merge <- df %>% left_join(select(data_D,V1,protectionStatus), by="V1")
  df_merge<- df_merge %>% column_to_rownames(.,var = 'V1')
  #table(df_merge$protectionStatus)

  
  df <- df %>% column_to_rownames(., var = 'V1')
  
  
  
  #BOTH Classes Calculate CV and subset genes
  results<-mvn(df[,1:20], univariateTest ="SW")
  results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
  results$Descriptives$CVround<-round(results$Descriptives$CV)
  after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
  print(i)
  df_merge<-fix.data(df_merge)
  mcfs_data<-subset(df_merge,select=c(rownames(after_cv),'protectionStatus'))
  mcfs_result <- mcfs(protectionStatus~., mcfs_data, featureFreq = 100, cutoffPermutations = 2, threadsNumber = 8)
  print(i)
  export.result(mcfs_result, path = "C://Users//giris//Downloads//", label = paste0("rmcfs", names(time_points[i])), zip = TRUE)
  

  if(i==1) rm(D0,results)
  if(i==2) rm(D3,results)
  if(i==3) rm(D7,results)
  if(i==4) rm(D126,results)
  if(i==5) rm(D129,results)
  if(i==6) rm(D133,results)
  if(i==7) rm(prechal,results)
  
  
  }

rm(list=ls()) #CLEAR THE ENVIRONMENT












#STANDARD DEVIATION

#CONVERT gene names as index or rownames
D <- D %>% column_to_rownames(., var = 'V1')
D0 <- D0 %>% column_to_rownames(., var = 'V1')
D3 <- D3 %>% column_to_rownames(., var = 'V1')

#animal<-animal %>% column_to_rownames(.,var='V1')

#logarithmic conversion to base 10
D<-log(D,10)
D0<-log(D0,10)
D3<-log(D3,10)

#convert infinitely negative values to zero
invisible(lapply(names(D),function(.name) set(D, which(is.infinite(D[[.name]])), j = .name,value =0)))
invisible(lapply(names(D0),function(.name) set(D0, which(is.infinite(D0[[.name]])), j = .name,value =0)))
invisible(lapply(names(D3),function(.name) set(D3, which(is.infinite(D3[[.name]])), j = .name,value =0)))

#remove rows with zero as read count for all animals
D<-D[apply(D[,-1], 1, function(x) !all(x==0)),]
D0<-D0[apply(D0[,-1], 1, function(x) !all(x==0)),]
D3<-D3[apply(D3[,-1], 1, function(x) !all(x==0)),]

df<-t(D)
df<-t(D0)
df<-t(D3)

#replace NAN values with zero
df[is.na(df)]<-0


df<-df[apply(df[,-1], 1, function(x) !all(x==0)),]

df<-as.data.frame(df)



#merge data
df <- df %>% rownames_to_column(., var = 'V1')
data_D<-data[grep('D129', data$Time_Point), ]

data_D0<-data[grep("D0", data$Time_Point), ]
data_D3<-data[grep("D3", data$Time_Point), ]


df_merge <- df[,1:200] %>% left_join(select(data_D0,V1,protectionStatus), by="V1")
df_merge<- df_merge %>% column_to_rownames(.,var = 'V1')
table(df_merge$protectionStatus)

df_merge_prot<-df_merge[grep('^Prot',df_merge$protectionStatus), ]
df_merge_nprot<-df_merge[grep('NonProt',df_merge$protectionStatus), ]

df <- df %>% column_to_rownames(., var = 'V1')


library(MVN)
#BOTH
mvn(df[,c('ENSMMUG00000007658','ENSMMUG00000000021')], mvnTest = "energy", univariatePlot = "histogram")
results<-mvn(df[,,c('ENSMMUG00000007658','ENSMMUG00000000015')], mvnTest = "energy")
results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
results$Descriptives$CVround<-round(results$Descriptives$CV)
after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
after_cv

#NONPROTECTED
mvn(df_merge_nprot[,annot$ENSEMBL], mvnTest = "energy", univariatePlot = "histogram")
results<-mvn(df_merge_nprot[,1:20], mvnTest = "energy")
results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
results$Descriptives$CVround<-round(results$Descriptives$CV)
after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
after_cv

#PROTECTED
mvn(df_merge_prot[,annot$ENSEMBL], mvnTest = "energy", univariatePlot = "histogram") 
results<-mvn(df_merge_prot[,1:20], mvnTest = "energy")
results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
results$Descriptives$CVround<-round(results$Descriptives$CV)
after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
after_cv


#ANNOTATION
genes <- c('CD276','GDF7','CCL22','IL17F','TMEM132B','CA2','KCNA2','CD300LG','TXNRD3','COL21A1')
ourCols <- c("SYMBOL",'ENSEMBL')
annot <- AnnotationDbi::select(org.Mmu.eg.db, 
                               keys=genes, 
                               columns=ourCols, 
                               keytype="SYMBOL")

###ANIMAL ANALYSIS
animal<-subset(count_data[,],select =c(V1,grep(data$AnimalID[8],colnames(count_data))))
animal<- animal %>% column_to_rownames(.,var='V1')
animal<-log(animal,10)
invisible(lapply(names(animal),function(.name) set(animal, which(is.infinite(animal[[.name]])), j = .name,value =0)))
animal<-animal[apply(animal[,-1], 1, function(x) !all(x==0)),]
animal<-t(animal)
mvn(animal[,1:20], univariateTest = 'SW', univariatePlot = "histogram")
results<-mvn(animal[,1:20], univariateTest = 'SW')
results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
results$Descriptives$CVround<-round(results$Descriptives$CV)
after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
after_cv




#######DESEQ2
library(DESeq2)
metadata<- data %>% column_to_rownames(.,var='V1')
metadata$new<-paste(metadata$protectionStatus,metadata$Time_Point,sep='_')
count_data<- count_data %>% column_to_rownames(., var='V1')
count_data[is.na(count_data)]<-0
dds <- DESeqDataSetFromMatrix(countData=count_data[1:100,], colData=metadata, design= ~ new  )
dds<-DESeq(dds)
resultsNames(dds)
results(dds)
contrast_oe <- c("new",'Prot_PreChal',"NonProt_PreChal")
result_d<-results(dds, contrast=contrast_oe, alpha = 0.05)

otop2Counts <- plotCounts(dds, gene="ENSMMUG00000000015", intgroup=c("Time_Point", "new"), returnData=TRUE)
otop2Counts$z_discrete <- cut(otop2Counts$count, breaks = 2, labels = c("Low", "High"))
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48")
ggplot(otop2Counts, aes(x=Time_Point, y=count, colour=new, group=new)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_gradient(low = 'blue', high = 'red') + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
color<-c('blue','red','yellow')

ggplot(otop2Counts, aes(x=Time_Point, y=count, colour=new, group=new)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")

dds_s<-as.data.frame(result_d)

sigGenes <- rownames(dds_s[dds_s$padj <= 1 & abs(dds_s$log2FoldChange) < 0.08,])
library(ggplot2)
heatmap <- ggplot(dds_s, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


metadata <- data.frame(row.names=colnames(animal[1:2]), )

d <- density(results$Descriptives$CVround) 
d$y <- d$y * length(results$Descriptives$CVround)  ; plot(d, ylim=c(0,4) )

hist(results$Descriptives$CVround, breaks=trunc(min(results$Descriptives$CVround)):(1+trunc(max(results$Descriptives$CVround))), add=TRUE)








mvn(df[,1:10],multivariatePlot = "qq")
mvn(df[,c('ENSMMUG00000002802','ENSMMUG00000013263')], mvnTest = "hz", multivariatePlot = "persp")
mvn(df[,16:30], mvnTest = "royston", multivariatePlot = "qqplot", multivariateOutlierMethod="quan")
mvn(df[,16:40], mvnTest = "energy", univariatePlot = "scatter")
an<-t(animal)
mvn(an[,124:125], mvnTest = "royston", univariatePlot = "histogram")

results
#plot(x, y, type = "l", lwd = 2, col = "blue", ylab = "probability", xlab = "x")

#sd_m<-t(sapply(split.default(D0, rownames(D0)), function(x)  {
# x1 <- unlist(x)
#data.frame(mean = mean(x1, na.rm = TRUE), sd = sd(x1, na.rm = TRUE))}))
#tail(sd_m)


View(ros$main)

try<-ros$main[grep('ACTN3|LIN52|GPR65|RALGDS|LY75|RHOQ|GSAP|HIGD1A|CA11|LCK|PMAIP1|TBK1|CLPX|SHOC2|AGTPBP1|UGCG|BACH1|CLPX|ZNF703|CREBRF|CTSF|SCYL2|LOC704637|GHDC|ALG1|PCGF5|ELL2|PIK3AP1|MYO1E',ros$main$features),]

visunet(try)
vis_out<-visunet(ros$main)

#create a new variable that contains node information for the "all" decision
nodes_RNO <- vis_out$all$nodes

#create a new vector of variables: shape. "dot" is the default shape of nodes
nodes_RNO$shape <- rep("dot", length(nodes_RNO$label))

#mark selected genes as stars using the label attribute 
nodes_RNO$shape[which(as.character(nodes_RNO$label) %in% l)] <- "star"

#create the node object list
nodesL <- list(nodes = nodes_RNO,CustCol =  c("shape"))

#rerun VisuNet with the new shape for nodes
vis_out2 <- visunet(ros$main, CustObjectNodes = nodesL)
 





discretize_fun<- function(x,val){
  

  first_qr<- after_cv[val,]$`25th`
  third_qr<- after_cv[val,]$`75th`
  if(x<first_qr)(x<-'L')
  if (x>third_qr)(x<-'H')
  else(x<-'M')
  
}



for(val in row.names(after_cv)){
  
  print(val)
  first_qr<- after_cv[val,]$`25th`
  third_qr<- after_cv[val,]$`75th`
  new_df<-discretize_fun(df[,val],first_qr,third_qr)

  #cbind(new_df,new_df)
}
colnames(df)[,i]






new_df<-df

discretize_fun<-function(x,row_val){
  first_qr<- after_cv[val,]$`25th`
  third_qr<- after_cv[val,]$`75th`
  if(x<first_qr)(x<-'L')
  if (x>third_qr)(x<-'H')
  else(x<-'M')
  
}
val<-rownames(after_cv)
new_df<-subset(df,select=val)

lapply(new_df, function(x){
  first_qr<- after_cv[val,]$`25th`
  third_qr<- after_cv[val,]$`75th`
  if(x<first_qr)(x<-'L')
  if (x>third_qr)(x<-'H')
  else(x<-'M')
  
})













discretize_quartile <- function(data_df, cuts_df, quartiles_df) {
  
  discretized_data <- data.frame(rownames(data_df))
  
  for (col_name in colnames(data_df)) {
    if (col_name %in% rownames(quartiles_df)) {
      q25 <- as.numeric(quartiles_df[col_name, "25th"])
      q75 <- as.numeric(quartiles_df[col_name, "75th"])
      
      # Discretize the column based on quartiles
      discretized_col <- cut(data_df[[col_name]], breaks = c(-Inf, q25, q75, Inf))
      
      # Rename levels to labels
      levels(discretized_col) <- c("Low", "Medium", "High")
      
      # Add the discretized column to the resulting dataframe
      discretized_data[,col_name]<-discretized_col
    } 
  }
  discretized_data<- discretized_data %>% column_to_rownames(., var='rownames.data_df.')
  return(discretized_data)
}


#metadat for quartile range
quartiles_df <-after_cv

# Column names for columns to be discretized
cut_columns <- rownames(after_cv)

#data to be discretized
data_df <- subset(df,select=cut_columns)

# Discretize the data based on quartiles

discretized_data_df <- discretize_quartile(data_df, cut_columns, quartiles_df)

# Display the resulting dataframe
print(discretized_data_df)













# First let's create a mapped data frame we can join to the gene expression values
_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Mmu.eg.db,
    keys = rmcfs$attribute,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("SYMBOL") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(rmcfs, by = c("SYMBOL" = "attribute"))
 



library(rmcfs)
library(org.Mmu.eg.db)
result_mcfs<-import.result(path='C://Users//giris//Downloads//test//rmcfs_D3.zip')

#genes<-'TBK1'
genes<-result_mcfs$RI[1:result_mcfs$cutoff_value,]

genes<-c('ENSMMUG00000017584',' ENSMMUG00000015702','ENSMMUG00000007419','ENSMMUG00000007457','ENSMMUG00000013416','ENSMMUG00000018394','ENSMMUG00000019862','ENSMMUG00000054945','ENSMMUG00000006692','ENSMMUG00000011577','ENSMMUG00000019851','ENSMMUG00000043581','ENSMMUG00000056626',' ENSMMUG00000008574','ENSMMUG00000009513','ENSMMUG00000023767','ENSMMUG00000014800')
genes<-genes$attribute
ourCols <- c("SYMBOL",'ENSEMBL','GENETYPE')
annot <- AnnotationDbi::select(org.Mmu.eg.db, 
                               keys=genes, 
                               columns=ourCols, 
                               keytype="ENSEMBL")

annot<-na.omit(annot[annot$GENETYPE=='protein-coding',])

list_genes<-annot[1:50,1]
new_names<-annot[1:50,2]
saveRDS(annot,'C://Users//giris//Downloads//test//annot_D3.rds')

saveRDS(list_genes,'C://Users//giris//Downloads//test//list_genes_24.rds')

saveRDS(new_names,'C://Users//giris//Downloads//test//new_genes_24.rds')

boots<-sample(annot$ENSEMBL,size=100,replace=TRUE)
ind<-match(boots,annot$ENSEMBL)
annot$SYMBOL[ind]





###time differentials


  
td<-fread('C://Users//giris//Downloads//SIV-Virus//time_diff//_count_data_D0_D3_.csv')

td<-subset(td,select=-(V1))
td<- td %>% column_to_rownames(.,var='AnimalID')

td<-td[, colSums(td != 0) > 0]

library(MVN)

mvn(td[,1:20], mvnTest = "energy", univariatePlot = "histogram")
results<-mvn(td[,1:50], mvnTest = "energy")
results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
results$Descriptives$CVround<-round(results$Descriptives$CV)
after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
after_cv






stor_file='C://Users//giris//Downloads//SIV-Virus//'

file_path='C://Users//giris//Downloads//SIV-Virus//time_diff//_count_data_D0_D3_.csv'

startTime <- Sys.time() 
for(i in seq(1,ncol(td),50)){
  
  # MVN results 
  results<-mvn(td[,i:(i+49)], mvnTest = "energy")
  
  #write mvn results
  results_name <- paste0(unlist(strsplit(basename(file_path),'[.]'))[1],'result',(i+49),'.csv')
  write.csv(results,paste(stor_file,results_name))
  
  
  #calculate after_cv
  results$Descriptives$CV<-(results$Descriptives$Std.Dev/results$Descriptives$Mean)*100
  results$Descriptives$CVround<-round(results$Descriptives$CV)
  after_cv<-subset(results$Descriptives,results$Descriptives$CV>50)
 
  
  #write after_cv
  aftercv_name <- paste0(unlist(strsplit(basename(file_path),'[.]'))[1],'after_cv',(i+49),'.csv')
  write.csv(after_cv,paste(stor_file,aftercv_name))
  
  if(i==101)(break)
  
  
  
}
endTime <- Sys.time() 



for(i in seq(1,ncol(td),1000)){
  print(i)
  print(i+999)
  
}


paste0(unlist(strsplit(basename('C://Users//giris//Downloads//SIV-Virus//time_diff//_count_data_D0_D3_.csv'),'[.]'))[1],'result',i)


