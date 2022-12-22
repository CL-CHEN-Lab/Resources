################### Change R1 or R2 in line 28(Coverage Input), line 92 cutoff value and line 140 (plot tittle) ###################

Path = '/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/'
Intergenic.Ensemble.ActiveGene.Hela.RPKM_1 = read.table(paste0(Path,'Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.tab'),header =T)
Intergenic.Ensemble.ActiveGene.Hela.RPKM_1 = Intergenic.Ensemble.ActiveGene.Hela.RPKM_1[-which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Start > Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$End),]
Intergenic = data.frame(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1[,1:3])
Intergenic$Chr = paste0('chr',Intergenic$Chr)
write.table(Intergenic,paste0(Path,"Intergenic.bed"),quote = F,col.names = F,row.names = F,sep="\t")
# run command below
# cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot
# bedtools coverage -counts -a Intergenic.bed -b MCM_R1_IZ.bed > MCM_R1_Coverage.bedgraph
# bedtools coverage -counts -a Intergenic.bed -b MCM_R2_IZ.bed > MCM_R2_Coverage.bedgraph
# bedtools coverage -counts -a Intergenic.bed -b MCM_IZ_R1_early.bed > MCM_R1_Coverage_Early.bedgraph
# bedtools coverage -counts -a Intergenic.bed -b MCM_IZ_R2_early.bed > MCM_R2_Coverage_Early.bedgraph

Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Chr = paste0("chr",Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Chr)

Intergenic.Ensemble.ActiveGene.Hela.RPKM_1 = data.frame(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1,ID = paste(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Chr,Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Start,sep="_"))

######  The Intergenic_Early.bed here means S50 < 0.5 ######
Intergenic_Early = read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_Early.bed")
Intergenic_Early = data.frame(Intergenic_Early , ID = paste(Intergenic_Early$V1,Intergenic_Early$V2,sep="_"))
######  Select the Intergenic early regions ######
Index = which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$ID%in%Intergenic_Early$ID)

Intergenic.Ensemble.ActiveGene.Hela.RPKM_1 = Intergenic.Ensemble.ActiveGene.Hela.RPKM_1[Index,]

Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage = read.table(paste0(Path,'MCM_R2_Coverage.bedgraph'))
Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage=Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[,-4]


##################  Intergenic Region Length distribution  ################
Coverage = Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage
Index_Coverage = which(Coverage[,4]>0)
Region = Intergenic.Ensemble.ActiveGene.Hela.RPKM_1
Region_WithCoverage = Intergenic.Ensemble.ActiveGene.Hela.RPKM_1[Index_Coverage,]

round(quantile(Region$Length/1000),3)
#0%       25%       50%       75%      100% 
#0.001     8.290    37.775   130.654 30583.880 

round(quantile(Region_WithCoverage$Length/1000),3)
#0%       25%       50%       75%      100% 
#0.016    46.916   100.709   241.097 30583.880 

round(mean(Region$Length/1000),3)
round(sd(Region$Length/1000),3)
#142kb +/- 590kb

round(mean(Region_WithCoverage$Length/1000),3)
round(sd(Region_WithCoverage$Length/1000),3)
#233kb +/- 758kb

##################  MCM Binding Region Length distribution  ################
MCM_R1 <- read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData/MCM/MCM_R1_IZ.bed")
MCM_R2 <- read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData/MCM/MCM_R2_IZ.bed")

Length_R1 = MCM_R1$V3-MCM_R1$V2
Length_R2 = MCM_R2$V3-MCM_R2$V2
round(quantile(Length_R1/1000),1)
round(quantile(Length_R2/1000),1)
round(mean(Length_R1/1000),1)
round(mean(Length_R2/1000),1)
round(sd(Length_R1/1000),1)
round(sd(Length_R2/1000),1)

##################  MCM Early Binding Region Length distribution  ################
MCM_R1 <- read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/2D/MCM_IZ_R1_early0-0.25.bed")
MCM_R2 <- read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/2D/MCM_IZ_R2_early0-0.25.bed")

Length_R1 = MCM_R1$V3-MCM_R1$V2
Length_R2 = MCM_R2$V3-MCM_R2$V2
round(quantile(Length_R1/1000),1)
round(quantile(Length_R2/1000),1)
round(mean(Length_R1/1000),1)
round(mean(Length_R2/1000),1)
round(sd(Length_R1/1000),1)
round(sd(Length_R2/1000),1)


##################  OKseq Length distribution  ################
OKseq <- read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData/OK-seq/OKseq.bed")

Length_OKseq = OKseq$V4
round(quantile(Length_OKseq/1000),1)
round(mean(Length_OKseq/1000),1)
round(sd(Length_OKseq/1000),1)

##################  Draw bar plot take 25%~75% early MCM segment length as cutoff (20~55kb)  ################
#Cutoff <- seq(from=0,to=30,by=5) 
Cutoff <- c(20,55) 
#Cutoff <- c(15,40) 
#Cutoff <- c(0,10,15,50,Inf) 
#Cutoff <- c(0,15,20,55,Inf) 
Counts <- matrix(NA,nrow=3,ncol=length(Cutoff)-1)


for (i in 2:length(Cutoff)) 
{
  ##Divergent
  print("Divergent")
  Chosen1 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Up==c(-1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Down==c(1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length<=Cutoff[i]*1000 & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length>Cutoff[i-1]*1000 & !is.na(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[,4]))
  Chosen2 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[Chosen1,4]>0)
  Counts[1,i-1] <- 100*length(Chosen2)/length(Chosen1)
  print(length(Chosen2))
  print(length(Chosen1))
  
  ##Tandem, 5End
  print("Tandem")
  Chosen1.1 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Up==c(1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Down==c(1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length<=Cutoff[i]*1000 & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length>Cutoff[i-1]*1000 & !is.na(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[,4]))
  Chosen2.1 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[Chosen1.1,4]>0)
  
  Chosen1.2 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Up==c(-1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Down==c(-1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length<=Cutoff[i]*1000 & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length>Cutoff[i-1]*1000 & !is.na(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[,4]))
  Chosen2.2 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[Chosen1.2,4]>0)
  
  Chosen1 <- c(Chosen1.1,Chosen1.2)
  Chosen2 <- c(Chosen2.1,Chosen2.2)
  Counts[2,i-1] <- 100*length(Chosen2)/length(Chosen1)
  print(length(Chosen2))
  print(length(Chosen1))
  
  
  ##Convergent
  print("Convergent")
  Chosen1 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Up==c(1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Strand_Down==c(-1) & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length<=Cutoff[i]*1000 & Intergenic.Ensemble.ActiveGene.Hela.RPKM_1$Length>Cutoff[i-1]*1000 & !is.na(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[,4]))
  Chosen2 <- which(Intergenic.Ensemble.ActiveGene.Hela.RPKM_1.Segment.HMM.Positive.15kb.Hela.Coverage[Chosen1,4]>0)
  Counts[3,i-1] <- 100*length(Chosen2)/length(Chosen1)
  print(length(Chosen2))
  print(length(Chosen1))
}


#Draw <- Counts
#colnames(Draw) <- c("<5kb","5-20kb","20-50kb",">50kb")
#colnames(Draw) <- c("<10kb","10-15kb","15-50kb",">50kb")
#colnames(Draw) <- c("<15kb","15~20kb","20~55kb",">55kb")
#colnames(Draw) <- c("<10kb","10-20kb","20-30kb","30-40kb","40-50kb","50-60kb","60-70kb","70-80kb","80-90kb",">90kb")

Draw <- Counts[,1]
# The title is too long to show in picture, add it in ppt slide "Proportions of all MCM segments (R1) within early (average S50 < 0.5) intergenic regions 
#        between 20~55kb (intergenic regions within 2nd~3rd quantile of early MCM seglength range)"
pdf(file = "/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Final_Result/R2.pdf",
    width = 6, # The width of the plot in inches
    height = 8)
#barplot(Draw,col=c("red","purple","blue"),main="Proportions of all MCM segments (R1) within early (average S50 < 0.5) intergenic regions 
#        between 20~55kb (intergenic regions within 2nd~3rd quantile of early MCM seglength range)",beside=T,ylim=c(0,100),axes=F)

barplot(Draw,col=c("red","purple","blue"),beside=T,ylim=c(0,100),axes=F)
axis(2,las=1,lwd=2)
dev.off()