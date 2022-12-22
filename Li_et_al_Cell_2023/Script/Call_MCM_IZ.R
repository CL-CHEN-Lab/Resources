library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(scales)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

MCM_R1 = read.table("/Volumes/WWT/MCM_project/BedGraph/rep1.bedgraph")
MCM_R2 = read.table("/Volumes/WWT/MCM_project/BedGraph/rep2.bedgraph")

MCM_R1_Filter = MCM_R1[which(MCM_R1$V4 > 0.6),]
MCM_R2_Filter = MCM_R2[which(MCM_R2$V4 > 0.8),]
colnames(MCM_R1_Filter) <- c("chr","start","end","value")
colnames(MCM_R2_Filter) <- c("chr","start","end","value")


IZ_R1 <- MCM_R1_Filter%>% makeGRangesFromDataFrame() %>% GenomicRanges::reduce(min.gapwidth=6000)%>%as.data.frame()
Remove <- which((IZ_R1$end-IZ_R1$start) <= 10499)
IZ <- IZ_R1[-Remove,]
write.table(IZ,"/Volumes/WWT/MCM_project/IZ_R1.bed",col.names = F,row.names = F,quote = F, sep = "\t")



IZ_R2 <- MCM_R2_Filter%>% makeGRangesFromDataFrame() %>% GenomicRanges::reduce(min.gapwidth=8000)%>%as.data.frame()
Remove <- which((IZ_R2$end-IZ_R2$start) <= 13000)
IZ <- IZ_R2[-Remove,]
write.table(IZ,"/Volumes/WWT/MCM_project/IZ_R2.bed",col.names = F,row.names = F,quote = F, sep = "\t")

