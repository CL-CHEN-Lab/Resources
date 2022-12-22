library(ggplot2)

Path = "/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Box_plot/"

Folder = c("R1","R2")
Type = c("early0-0.25","earlyMid0.25-0.5","midLate0.5-0.75","Late0.75-1")
Name = c("Early","MidEarly","MidLate","Late")
R1 = data.frame()
R2 = data.frame()

All_Data = list(R1,R2)
for(i in 1:2)
{
  for (j in 1:4) 
  {
    Input = paste0(Path,Folder[i],"/MCM_IZ_",Folder[i],"_",Type[j],".bed")
    read = read.table(Input)
    tmp = data.frame(read,log2(read[,3]-read[,2]),Name[j])
    colnames(tmp) = c("Chr","Start","End","S50","Length","Group")
    All_Data[[i]] = rbind(All_Data[[i]],tmp)
  }
}

library(gridExtra)
library(ggplot2)
library(grid)
library(forcats)
library(dplyr)

Plot= list()
Path = "/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Box_plot/Boxplot"
for(i in 1:2)
{
  Plot[[i]] = All_Data[[i]] %>% mutate(Group = fct_relevel(Group,"Early","MidEarly","MidLate","Late"))%>% ggplot(aes(x=Group, y=Length,fill = Group)) + 
    geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=T,notchwidth=0.8,fill=c("red","cyan","blue","orange")) + 
    ggtitle(paste("MCM binding region length in ", Folder[i]))+ 
    #scale_y_continuous(limits = c(0, 3000))+ theme_bw()+
    theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ylab("Log2 Length")
    #ylim(0,3000)
  ggsave(paste0(Path,"_",Folder[i],".pdf"),Plot[[i]],width = 8,height = 8,)
}


#Final_Plot = grid.arrange(Plot[[1]],Plot[[2]],nrow=1,top=textGrob(expression(bold("MCM binding region length in different timing group")), gp=gpar(fontsize=25, face= 'bold')))
#ggsave("/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Box_plot/Boxplot.pdf",Final_Plot,width = 16,height = 8,)