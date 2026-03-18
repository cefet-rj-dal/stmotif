source("https://raw.githubusercontent.com/eogasawara/CSA/master/code/stmotifv2.R")

library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

setwd("/Users/heraldoborges/Dropbox/cefet/LIRMM/R-Journal")

load("toydataset.RData")


#input toydataset
w  <- 4 #tamanho da palavra
a  <- 5 #alfabeto
si <- 3 #Minimum number of occurrences inside each block
ka <- 2 #Minimum number of spatial-time series with occurrences inside each block
tb <- 10 #temporal blocks
sb <- 9 # spatial blocks


#CSA Mining Process 
D <- STSADatasetAdjust(vectorMatrix, tb, sb)
#normalization and SAX indexing 
DS <- NormSAX(D,a)
#Check and filter the stmotifs, grouping the motifs from the neighboring block.
stmotifs <- SearchSTMotifs(D,DS,w,a,sb,tb,si,ka)
#Rank the stmotifs list
rstmotifs <- RankSTMotifs(stmotifs)

#CSA Mining Process - All this process can be summarized in this function 
rstmotifs <- CSAMiningProcess(D,DS,w,a,sb,tb,si,ka)

#### Visualization
#grafico exportado em 5" x 4.25"
display_motifsDataset(dataset = D, rstmotifs,  5)

#grafico exportado em 5" x 4.25"
display_motifsSTSeries(dataset = D,rstmotifs[c(1:4)],space = c(1:4,10:12))


## Plots  Paper 

######## 1.SAX 
dataSTS <- as.data.frame(D)
vector <- as.matrix(D)
vector <- as.vector(vector)
vectorNorm <- (vector-mean(vector, na.rm = T))/stats::sd(vector, na.rm = T)
mybin <- binning(vector, a)
yLabs <- c("<a", "b", "c", "d", ">e")
ggplot(data = dataSTS[1], aes(x = 1:40, y =dataSTS[1]$V1 ), fill=supp)+
  geom_line(color = "#E7B800", size = 0.5) + 
  scale_y_continuous(breaks = c(82,324,598,832,1000), labels=yLabs) + 
  theme_minimal()+
  theme(panel.grid.major.y =  element_line(colour = "#808080", size=0.2))+
  labs( x = "", y = "SAX Values")+ 
  geom_point(size=2, shape=21,color="#E7B800", fill="#E7B800") +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold",  size=14),
        axis.line = element_line( size = 0.5, linetype = "solid")) 



######## 2.display_motifsSTSeries
display_motifsSTSeries(dataset = D, rstmotifs = NULL)
###### 2.1display_motifsSTSeriesWith3BestRankedMotifs 
display_motifsSTSeries(dataset = D, rstmotifs[c(1:3)],space = c(1:36))
###### 2.2display_motifsSTSeriesWithFirstBestRanked
display_motifsSTSeries(dataset = D, rstmotifs[1],space = c(1:36),"#FC8D59")
###### 2.3display_motifsSTSeriesWithSecondBestRanked
display_motifsSTSeries(dataset = D, rstmotifs[2],space = c(1:15),"#ffff00")
###### 2.4display_motifsSTSeriesWithThirdBestRanked
display_motifsSTSeries(dataset = D, rstmotifs[3],space = c(1:18),"#99D594")


####### 3.display_motifsDatasetWith3BestRankedMotifs 
display_motifsDataset(dataset = D,rstmotifs[1:3],  5)
##### 3.1display_motifsDatasetWithFirstBestRanked
display_motifsDataset(dataset = D,rstmotifs[1],  5, "#FC8D59")
##### 3.2display_motifsDatasetWithSecondBestRanked
display_motifsDataset(dataset = D,rstmotifs[2],  5, "#FFFFBF")
##### 3.3display_motifsDatasetWithThirdBestRanked
display_motifsDataset(dataset = D,rstmotifs[3],  5, color ="#99D594")


