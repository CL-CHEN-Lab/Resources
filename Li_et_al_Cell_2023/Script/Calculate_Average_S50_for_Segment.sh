#!/bin/bash

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/VennPlot/S50_1kbBin

#bedtools makewindows -g hg19.txt -w 1000 > hg19_window.bed

#Run Java script GetS50Timing to get the 1kb bin S50 values (hg19_window_S50.bed)


#bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/MCM_R2_IZ.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.25 && $4 != 0)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/MCM_IZ_R2_early.bed

#bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/MCM_R1_IZ.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.25 && $4 != 0)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/MCM_IZ_R1_early.bed



bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$7;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.5 && $4 != 0)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_Early.bed


bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$7;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.25 && $4 != 0)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_VeryEarly.bed


bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$7;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.5 && $4 > 0.25)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_MidEarly.bed





cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot

bedtools coverage -counts -a Intergenic_Early.bed -b MCM_R1_IZ.bed > MCM_R1_Coverage.bedgraph
bedtools coverage -counts -a Intergenic_Early.bed -b MCM_R2_IZ.bed > MCM_R2_Coverage.bedgraph
bedtools coverage -counts -a Intergenic_Early.bed -b OKseq.bed > OKseq_Coverage.bedgraph


bedtools coverage -counts -a Intergenic_VeryEarly.bed -b MCM_R1_IZ.bed > MCM_R1_Coverage_VeryEarly.bedgraph
bedtools coverage -counts -a Intergenic_VeryEarly.bed -b MCM_R2_IZ.bed > MCM_R2_Coverage_VeryEarly.bedgraph
bedtools coverage -counts -a Intergenic_VeryEarly.bed -b OKseq.bed > OKseq_Coverage_VeryEarly.bedgraph


bedtools coverage -counts -a Intergenic_MidEarly.bed -b MCM_R1_IZ.bed > MCM_R1_Coverage_MidEarly.bedgraph
bedtools coverage -counts -a Intergenic_MidEarly.bed -b MCM_R2_IZ.bed > MCM_R2_Coverage_MidEarly.bedgraph
bedtools coverage -counts -a Intergenic_MidEarly.bed -b OKseq.bed > OKseq_Coverage_MidEarly.bedgraph























bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$7;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.25 && $4 != 0)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_Smaller0.25.bed


bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$7;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '($4 <= 0.5 && $4 > 0.25)' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot/Intergenic_0.25_0.5.bed






cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/CO_AO_BarPlot

#bedtools coverage -counts -a Intergenic_early.bed -b MCM_R1_IZ.bed > MCM_R1_Coverage.bedgraph
#bedtools coverage -counts -a Intergenic_early.bed -b MCM_R2_IZ.bed > MCM_R2_Coverage.bedgraph

#bedtools coverage -counts -a Intergenic_early.bed -b MCM_IZ_R1_early.bed > MCM_R1_Coverage_Early.bedgraph
#bedtools coverage -counts -a Intergenic_early.bed -b MCM_IZ_R2_early.bed > MCM_R2_Coverage_Early.bedgraph




bedtools coverage -counts -a Intergenic_Smaller0.25.bed -b OKseq.bed > OK_Coverage_Early.bedgraph
bedtools coverage -counts -a Intergenic_0.25_0.5.bed -b OKseq.bed > OK_Coverage_MidEarly.bedgraph

bedtools coverage -counts -a Intergenic_Smaller0.25.bed -b MCM_R2_IZ.bed > R2_Coverage_Early.bedgraph
bedtools coverage -counts -a Intergenic_0.25_0.5.bed -b MCM_R2_IZ.bed > R2_Coverage_MidEarly.bedgraph








