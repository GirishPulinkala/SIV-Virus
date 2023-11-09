heatmap_rules<-function(recal,discretized_data_df){
  library(ComplexHeatmap)
  dataM <- data.frame(matrix(ncol = length(recal$features[1:15]), nrow = length(row.names(discretized_data_df))))
  
  
  discretized_data_df<-discretized_data_df[,1:(ncol(discretized_data_df)-1)]
  discretized_data_df<- discretized_data_df %>% rownames_to_column(.,var = 'V1')
  
  
  discretized_data_df <- discretized_data_df %>% left_join(dplyr::select(data,V1,Group,protectionStatus), by="V1")
  discretized_data_df<- discretized_data_df %>% column_to_rownames(.,var = 'V1')
  discretized_data_df<-discretized_data_df[discretized_data_df$protectionStatus=='Prot'| discretized_data_df$protectionStatus=='NonProt', ]
  
  
  
  rownames(dataM)<-rownames(discretized_data_df)
  
  
  for(i in 1:length(recal$features[1:15])){
    
    list_of_features<-unlist(strsplit(recal$supportSetLHS[i],split=','))
    
    for(j in 1:length(list_of_features)){
      if(list_of_features[j] %in% rownames(dataM))(dataM[list_of_features[j],i]<-1)
      
    }
    
  }
  
  
  dataM[is.na(dataM)] <- 0
  
  ann <- data.frame(discretized_data_df[,(ncol(discretized_data_df)-1)], (discretized_data_df[,(ncol(discretized_data_df))]))
  
  colnames(ann) <- c('Group', 'PS')
  
  colAnn <- HeatmapAnnotation(df = ann,
                              which = 'row',
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))
  
  
  boxplotCol <- HeatmapAnnotation(
    boxplot=anno_boxplot(data.matrix(data),
                         border=TRUE,
                         gp=gpar(fill="#CCCCCC"),
                         pch=".",
                         size=unit(2, "mm"),
                         axis=TRUE,
                         
                         annotation_width=unit(c(1, 5.0), "cm"),
                         which="row")    )                        
  
  
  
  hmap <- Heatmap(
    dataM,
    name = "VALUE",
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    width = unit(350, "mm"),
    right_annotation=colAnn)
  
  
  draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
  
  
  
}


