



after_cv<-function(timepoint){

library(tidyverse)


#READ DATA about quartile after calculating coefficient of variation
setwd("C://Users//giris//Downloads//subset_run")
filelist <- list.files(recursive=TRUE)



list_of_all_sites <- paste0('after_cv',timepoint)

list_of_all_dataframes <- map(
  list_of_all_sites,
  ~ list.files(pattern = .x, recursive = TRUE) %>%
    map_df(read.csv)
)

after_cv<- as.data.frame(list_of_all_dataframes)

after_cv<- after_cv %>% column_to_rownames(.,var = 'X')

return(after_cv)

}
