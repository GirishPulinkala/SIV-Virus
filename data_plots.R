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
