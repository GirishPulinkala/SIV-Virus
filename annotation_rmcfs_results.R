
annotation_man<-function(file_path){
  library(rmcfs)
  library(org.Mmu.eg.db)
  result_mcfs<-import.result(path=file_path)
  
  #genes<-'TBK1'
  genes<-result_mcfs$RI[1:result_mcfs$cutoff_value,]
  
  #genes<-c('ENSMMUG00000017584',' ENSMMUG00000015702','ENSMMUG00000007419','ENSMMUG00000007457','ENSMMUG00000013416','ENSMMUG00000018394','ENSMMUG00000019862','ENSMMUG00000054945','ENSMMUG00000006692','ENSMMUG00000011577','ENSMMUG00000019851','ENSMMUG00000043581','ENSMMUG00000056626',' ENSMMUG00000008574','ENSMMUG00000009513','ENSMMUG00000023767','ENSMMUG00000014800')
  genes<-genes$attribute
  ourCols <- c("SYMBOL",'ENSEMBL','GENETYPE')
  annot <- AnnotationDbi::select(org.Mmu.eg.db, 
                                 keys=genes, 
                                 columns=ourCols, 
                                 keytype="ENSEMBL")
  
  annot<-na.omit(annot[annot$GENETYPE=='protein-coding',])

  return(annot)

}




