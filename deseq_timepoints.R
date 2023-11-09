data<-fread('C://Users//giris//Downloads//target_renamed.csv')
count_data<-fread('C://Users//giris//Downloads//count_matrix_renamed.csv')
count_data<- count_data %>% column_to_rownames(.,var='V1')


library(DESeq2)


deseq2Data <- DESeqDataSetFromMatrix(countData=count_data, colData=data, design= ~ Time_Point+protectionStatus+ Time_Point:protectionStatus)
dim(deseq2Data)
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
deseq2Data <- DESeq(deseq2Data, test = 'LRT',reduced = ~Time_Point+protectionStatus)

deseq2Results <- results(deseq2Data, contrast=c('protectionStatus','Prot','NonProt'))
deseq2ResDF <- as.data.frame(deseq2Results)

deseq2VST <- vst(deseq2Data)
deseq2VST <- assay(deseq2VST)

deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)


sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .1 & abs(deseq2ResDF$log2FoldChange) > 1,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]








otop2Counts <- plotCounts(deseq2Data, gene='ENSMMUG00000010291', intgroup=c("Time_Point", "protectionStatus"), returnData=TRUE)
ggplot(otop2Counts, aes(x=Time_Point, y=count, colour=protectionStatus, group=protectionStatus)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")



