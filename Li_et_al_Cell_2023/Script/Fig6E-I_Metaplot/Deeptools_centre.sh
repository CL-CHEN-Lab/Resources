#!/bin/bash


cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/R1

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed ./MCM_IZ_R1_earlyMid0.25-0.5.bed ./MCM_IZ_R1_midLate0.5-0.75.bed ./MCM_IZ_R1_Late0.75-1.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep1.bw --outFileName ./R1.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros 


   plotProfile -m ./R1.tab.gz -out ./R1.pdf --refPointLabel MCM_Binding_Zone --regionsLabel E ME ML L --colors red cyan blue orange 



   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep1.bw --outFileName ./R1_early.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros 

   plotHeatmap -m ./R1_early.tab.gz -out ./R1_Early.pdf --heatmapHeight 10 --refPointLabel MCM_R1.center --colorMap RdBu --colorList mediumblue,white,firebrick --regionsLabel "" --sortUsing region_length --sortRegions ascend




   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed ./MCM_IZ_R1_earlyMid0.25-0.5.bed ./MCM_IZ_R1_midLate0.5-0.75.bed ./MCM_IZ_R1_Late0.75-1.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep1.bw --outFileName ./R1_region.tab.gz --beforeRegionStartLength 60000 --afterRegionStartLength 60000 --binSize 1000 --missingDataAsZero --skipZeros --regionBodyLength 60000


   plotProfile -m ./R1_region.tab.gz -out ./R1_region.pdf --refPointLabel MCM_Binding_Zone --regionsLabel E ME ML L --colors red cyan blue orange --startLabel "MCM_S" --endLabel "MCM_E" 






  
   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/ORM/ORM_FireEfficiency.bw --outFileName ../ORM/R1_ORM.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../ORM/R1_ORM.tab.gz -out ../ORM/R1_ORM_ByLength.pdf --heatmapHeight 8 --refPointLabel MCM_R1.center --colorMap RdBu --colorList mediumblue,white,firebrick,maroon --regionsLabel "" --sortUsing region_length --sortRegions ascend --linesAtTickMarks

   #plotHeatmap -m ../ORM/R1_ORM.tab.gz -out ../ORM/R1_ORM_Ascend_BySignal.pdf --heatmapHeight 8 --refPointLabel MCM_R1.center --colorMap RdBu --colorList mediumblue,white,firebrick,maroon --regionsLabel "" --sortRegions ascend




   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/RFD/OKseq_RFD.bw --outFileName ../RFD/R1_RFD.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../RFD/R1_RFD.tab.gz -out ../RFD/R1_RFD.pdf --heatmapHeight 8 --refPointLabel MCM_R1.center --regionsLabel "" --colorMap RdBu --sortRegions ascend --sortUsing region_length --colorList mediumblue,white,firebrick



#--colorList lavender,royalblue,gold,darkorange,firebrick
   


########  Chromatin markers  ########  

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData


########  GRO-seq, H2AZ, Pol2-pho  ########  

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed --scoreFileName H2AZ.bw GRO-seq.bw POLR2AphosphoS2.bw --outFileName ../Deeptool/Epigenetic/R1_Epigenetic_1.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/R1_Epigenetic_1.tab.gz -out ../Deeptool/Epigenetic/R1_Epigenetic_1.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel MCM_R1_Center --regionsLabel "" --colorMap RdBu


########  ORC number ########  

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R1_early0-0.25.bed --scoreFileName ORC_SegNumber.bw --outFileName ../Deeptool/Epigenetic/R1_Epigenetic_1_ORC.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/R1_Epigenetic_1_ORC.tab.gz -out ../Deeptool/Epigenetic/R1_Epigenetic_1_ORC_remove2.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel MCM_R1_Center --regionsLabel "" --colorMap RdBu --colorList white,firebrick,darkorange



#######.  Test as postive control  ###########
########  ORC regions MCM signal, ORM signal, RFD signal enrichment ########  

   #computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./ORC-chip.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep1.bw --outFileName ../Deeptool/Epigenetic/ORC_MCM.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./ORC-chip.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep1.bw --outFileName ../Deeptool/Epigenetic/ORC_MCM.tab.gz --beforeRegionStartLength 600 --afterRegionStartLength 600 --binSize 100 --regionBodyLength 600 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/ORC_MCM.tab.gz -out ../Deeptool/Epigenetic/ORC_MCM.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel ORC_Center --regionsLabel MCM_Early_R1 --colorMap RdBu --colorList lavender,royalblue,gold,darkorange,firebrick




   #computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./ORC-chip.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/ORM/ORM_FireEfficiency.bw --outFileName ../Deeptool/Epigenetic/ORC_ORM.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./ORC-chip.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/ORM/ORM_FireEfficiency.bw --outFileName ../Deeptool/Epigenetic/ORC_ORM.tab.gz --beforeRegionStartLength 600 --afterRegionStartLength 600 --binSize 100 --regionBodyLength 600 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/ORC_ORM.tab.gz -out ../Deeptool/Epigenetic/ORC_ORM.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel ORC_Center --regionsLabel ORM_fire_efficiency --colorMap RdBu --colorList lavender,royalblue,gold,darkorange,firebrick



   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./ORC-chip.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/RFD/OKseq_RFD.bw --outFileName ../Deeptool/Epigenetic/ORC_RFD.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/ORC_RFD.tab.gz -out ../Deeptool/Epigenetic/ORC_RFD.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel ORC_Center --regionsLabel RFD --colorMap RdBu --colorList lavender,royalblue,gold,darkorange,firebrick






############    R2    #############

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/R2

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed ./MCM_IZ_R2_earlyMid0.25-0.5.bed ./MCM_IZ_R2_midLate0.5-0.75.bed ./MCM_IZ_R2_Late0.75-1.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep2.bw --outFileName ./R2.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros 

   plotProfile -m ./R2.tab.gz -out ./R2.pdf --refPointLabel MCM_Binding_Zone --regionsLabel E ME ML L --colors red cyan blue orange




   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed ./MCM_IZ_R2_earlyMid0.25-0.5.bed ./MCM_IZ_R2_midLate0.5-0.75.bed ./MCM_IZ_R2_Late0.75-1.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bw/ChIP_MCM_rep2.bw --outFileName ./R2_region.tab.gz --beforeRegionStartLength 60000 --afterRegionStartLength 60000 --binSize 1000 --missingDataAsZero --skipZeros --regionBodyLength 60000

   plotProfile -m ./R2_region.tab.gz -out ./R2_region.pdf --refPointLabel MCM_Binding_Zone --regionsLabel E ME ML L --colors red cyan blue orange --startLabel "MCM_S" --endLabel "MCM_E"






   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/ORM/ORM_FireEfficiency.bw --outFileName ../ORM/R2_ORM.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortRegions keep

   plotHeatmap -m ../ORM/R2_ORM.tab.gz -out ../ORM/R2_ORM_ByLength.pdf --heatmapHeight 8 --refPointLabel MCM_R2.center --colorMap RdBu --colorList mediumblue,white,firebrick,maroon --regionsLabel "" --sortUsing region_length --sortRegions ascend --linesAtTickMarks

   #plotHeatmap -m ../ORM/R2_ORM.tab.gz -out ../ORM/R2_ORM_Ascend_BySignal.pdf --heatmapHeight 8 --refPointLabel MCM_R2.center --colorMap RdBu --colorList mediumblue,white,firebrick,maroon --regionsLabel "" --sortRegions ascend





   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed --scoreFileName /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/Deeptool/RFD/OKseq_RFD.bw --outFileName ../RFD/R2_RFD.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortRegions keep 

   plotHeatmap -m ../RFD/R2_RFD.tab.gz -out ../RFD/R2_RFD.pdf --heatmapHeight 8 --refPointLabel MCM_R2.center --regionsLabel "" --colorMap RdBu --sortRegions ascend --sortUsing region_length --colorList mediumblue,white,firebrick



cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData

computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed --scoreFileName DNAseI.bw GRO-seq_Minus.bw GRO-seq_Plus.bw H2AZ.bw POLR2AphosphoS2.bw ORC.bw --outFileName ../Deeptool/Epigenetic/R2_Epigenetic.tab.gz --referencePoint center -a 60000 -b 60000 --binSize 1000 --missingDataAsZero --skipZeros 

   plotHeatmap -m ../Deeptool/Epigenetic/R2_Epigenetic.tab.gz -out ../Deeptool/Epigenetic/R2_Epigenetic.pdf --heatmapHeight 16 --refPointLabel R2_Binding_Regions_Center --regionsLabel MCM_Early_Binding_Regions --colorMap RdBu --colorList lavender,royalblue,gold,darkorange,firebrick




########  Chromatin markers  ########  

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/RawData

########  GRO-seq, H2AZ, Pol2-pho  ########  

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed --scoreFileName H2AZ.bw GRO-seq.bw POLR2AphosphoS2.bw --outFileName ../Deeptool/Epigenetic/R2_Epigenetic_2.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/R2_Epigenetic_2.tab.gz -out ../Deeptool/Epigenetic/R2_Epigenetic_2.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel MCM_R2_Center --regionsLabel "" --colorMap RdBu


########  ORC number ########  

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./MCM_IZ_R2_early0-0.25.bed --scoreFileName ORC_SegNumber.bw --outFileName ../Deeptool/Epigenetic/R2_Epigenetic_2_ORC.tab.gz --referencePoint center -a 100000 -b 100000 --binSize 1000 --missingDataAsZero --skipZeros --sortUsing region_length --sortRegions ascend

   plotHeatmap -m ../Deeptool/Epigenetic/R2_Epigenetic_2_ORC.tab.gz -out ../Deeptool/Epigenetic/R2_Epigenetic_2_ORC.pdf --heatmapWidth 6 --heatmapHeight 10 --refPointLabel MCM_R2_Center --regionsLabel "" --colorMap RdBu --colorList white,firebrick,darkorange



#--colorList lavender,royalblue,gold,darkorange,firebrick




