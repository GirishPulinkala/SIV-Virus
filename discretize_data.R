discretize_quartile <- function(data_df, cuts_df, quartiles_df) {
  
  discretized_data <- data.frame(rownames(data_df))
  
  for (col_name in colnames(data_df)) {
    if (col_name %in% rownames(quartiles_df)) {
      q25 <- as.numeric(quartiles_df[col_name, "X25th"])
      q75 <- as.numeric(quartiles_df[col_name, "X75th"])
      
      if(q25==q75)(q75<-q75+0.05)
      # Discretize the column based on quartiles
      discretized_col <- cut(data_df[[col_name]], breaks = c(-Inf, q25, q75, Inf))
      
      # Rename levels to labels
      levels(discretized_col) <- c(1, 2, 3)
      
      # Add the discretized column to the resulting dataframe
      discretized_data[,col_name]<-discretized_col
    } 
  }
  discretized_data<- discretized_data %>% column_to_rownames(., var='rownames.data_df.')
  return(discretized_data)
}
