library(dplyr)

#MCM_R1 = read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/MCM_R1.bed")
MCM_R2 = read.table("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/MCM_R2.bed")

R2_Cutoff = seq(from=35,to=80,by=5)

All_fragment = list()
colnames(MCM_R2)<-c("chr","start","end","Length")
Tmp = MCM_R2

for(i in 1:9)
{
  Index = between(Tmp$Length,R2_Cutoff[i],R2_Cutoff[i+1])
  All_fragment[[i]] <- Tmp[Index,]
  print(nrow(All_fragment[[i]]))
  write.table(All_fragment[[i]],paste0("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/Classify/",R2_Cutoff[i],"~",R2_Cutoff[i+1],"bp.bed"),quote = F,col.names = F,row.names = F,sep="\t")
}