library(gapminder)
library(tidyverse)
library(infer)
library(data.table)
library(treemap)
library(data.tree)
library(networkD3)

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

#Subsetting data based on timepoints for normalization and further steps
D0<-subset(count_data,select=c(V1,grep("W0_D0|_[0-9]+_D0",colnames(count_data))))
D3<-subset(count_data,select=c(V1,grep("W0_D3|_[0-9]+_D3",colnames(count_data))))
D7<-subset(count_data,select=c(V1,grep("W0_D7|_[0-9]+_D7",colnames(count_data))))
D126<-subset(count_data,select=c(V1,grep("W18_D0|_[0-9]+_D98|_[0-9]+_D126",colnames(count_data))))
D129<-subset(count_data,select=c(V1,grep("W18_D3|_[0-9]+_D98|_[0-9]+_D129",colnames(count_data))))
D133<-subset(count_data,select=c(V1,grep("W18_D7|_[0-9]+_D98|_[0-9]+_D133",colnames(count_data))))
prechal<-subset(count_data,select=c(V1,grep("W88_D0|_[0-9]+_D462|_[0-9]+_D455|_[0-9]+_D538|_[0-9]+_D448|_[0-9]+_D540",colnames(count_data))))



#STANDARD DEVIATION

animal<-subset(count_data,select=(grep('29523',colnames(count_data))))
animal<-subset(count_data,select=-c(v1))
animal

x<-as.numeric(animal[36,])
sd(x)
mean(x)
y<-dnorm(x, mean(x), sd(x))

D0<<-D0[,1] 
df<-t(D0)
View(df)

plot(x, y, type = "l", lwd = 2, col = "blue", ylab = "probability", xlab = "x")

t(sapply(split.default(df, names(df)), function(x)  {
  x1 <- unlist(x)
  data.frame(mean = mean(x1, na.rm = TRUE), sd = sd(x1, na.rm = TRUE))}))





