#!/bin/bash

#Input is copied from /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810/Input_Control.bed

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd/Coverage

bedtools makewindows -g ChrLength.txt -w 1000 | awk '{print $0,NR}' | sed -e 's/ /\t/g' > hg19_1kb_Bin.bed 

bedtools makewindows -g ChrLength.txt -w 5000 | awk '{print $0,NR}' | sed -e 's/ /\t/g' > hg19_5kb_Bin.bed 

bedtools makewindows -g ChrLength.txt -w 10000 | awk '{print $0,NR}' | sed -e 's/ /\t/g' > hg19_10kb_Bin.bed 


grep -v "_" ../short/R1.bed | grep -v "chrM" | sort -k1,1 -k2,2n > ./R1_filter.bed

grep -v "_" ../short/R2.bed | grep -v "chrM" | sort -k1,1 -k2,2n > ./R2_filter.bed

grep -v "-" Input_Control.bed | grep -v "chrM" | sort -k1,1 -k2,2n > ./Input_Control_filter.bed


sort -k 1.4 /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd/Coverage/ChrLength.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd/Coverage/ChrLength_sort.txt 




sort -k1,1 -k2,2n hg19_1kb_Bin.bed > hg19_1kb_Bin_sort.bed
sort -k1,1 -k2,2n hg19_5kb_Bin.bed > hg19_5kb_Bin_sort.bed
sort -k1,1 -k2,2n hg19_10kb_Bin.bed > hg19_10kb_Bin_sort.bed




bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_10kb_Bin_sort.bed \
                     -b ./Input_Control_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./Input_10kb_Coverage.bedgraph

bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_1kb_Bin_sort.bed \
                     -b ./Input_Control_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./Input_1kb_Coverage.bedgraph



bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_10kb_Bin_sort.bed \
                     -b ./R1_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R1_10kb_Coverage.bedgraph

bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_10kb_Bin_sort.bed \
                     -b ./R2_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R2_10kb_Coverage.bedgraph


bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_5kb_Bin_sort.bed \
                     -b ./R1_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R1_5kb_Coverage.bedgraph

bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_5kb_Bin_sort.bed \
                     -b ./R2_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R2_5kb_Coverage.bedgraph


bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_1kb_Bin_sort.bed \
                     -b ./R1_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R1_1kb_Coverage.bedgraph

bedtools coverage -counts \
                     -sorted \
                     -a ./hg19_1kb_Bin_sort.bed \
                     -b ./R2_filter.bed \
                     -g ./ChrLength_sort.txt | cut -f 1-3,5 > ./R2_1kb_Coverage.bedgraph


BedGraphToBigWig R1_10kb_Coverage.bedgraph ChrLength_sort.txt R1_10kb_Coverage.bw
BedGraphToBigWig R2_10kb_Coverage.bedgraph ChrLength_sort.txt R2_10kb_Coverage.bw
BedGraphToBigWig R1_5kb_Coverage.bedgraph ChrLength_sort.txt R1_5kb_Coverage.bw
BedGraphToBigWig R2_5kb_Coverage.bedgraph ChrLength_sort.txt R2_5kb_Coverage.bw
BedGraphToBigWig R1_1kb_Coverage.bedgraph ChrLength_sort.txt R1_1kb_Coverage.bw
BedGraphToBigWig R2_1kb_Coverage.bedgraph ChrLength_sort.txt R2_1kb_Coverage.bw


BedGraphToBigWig Input_10kb_Coverage.bedgraph ChrLength_sort.txt Input_10kb_Coverage.bw






multiBigwigSummary bins -b R1_10kb_Coverage.bw R2_10kb_Coverage.bw Input_10kb_Coverage.bw -o ./results_add_MCM.npz

plotCorrelation -in ./results_add_MCM.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of MCM density of R1, R2 and Input_Control" --whatToPlot scatterplot -o ./Pearson_10kb_MCM.pdf -c pearson --labels R1_10kb R2_10kb Input_Control --plotHeight 30 --plotWidth 30 --removeOutliers

plotCorrelation -in ./results_add_MCM.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of MCM density of R1, R2 and Input_Control" --whatToPlot scatterplot -o ./Spearman_10kb_MCM.pdf -c pearson --labels R1_10kb R2_10kb Input_Control --plotHeight 30 --plotWidth 30 --removeOutliers














multiBigwigSummary bins -b R1_10kb_Coverage.bw R2_10kb_Coverage.bw -o ./results_add_MCM.npz

 plotCorrelation -in ./results_add_MCM.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of MCM density of R1 and R2" --whatToPlot scatterplot -o ./Pearson_${array[$j-1]}_MCM.pdf -c pearson --labels R1_${array[$j-1]} R2_${array[$j-1]} --plotHeight 30 --plotWidth 30 --removeOutliers



array=(
10kb
5kb
1kb
)

for j in $(seq 1 3)

do

   echo ${array[$j-1]}

   #multiBigwigSummary bins -b R1_${array[$j-1]}_Coverage.bw R2_${array[$j-1]}_Coverage.bw -o ./results_${array[$j-1]}_MCM.npz

   plotCorrelation -in ./results_${array[$j-1]}_MCM.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of MCM density of R1 and R2" --whatToPlot scatterplot -o ./Pearson_${array[$j-1]}_MCM.pdf -c pearson --labels R1_${array[$j-1]} R2_${array[$j-1]} --plotHeight 30 --plotWidth 30 --removeOutliers --log1p

   plotCorrelation -in ./results_${array[$j-1]}_MCM.npz --corMethod pearson --skipZeros --plotTitle "Spearman Correlation of MCM density of R1 and R2" --whatToPlot scatterplot -o ./Spearman_${array[$j-1]}_MCM.pdf -c spearman --labels R1_${array[$j-1]} R2_${array[$j-1]} --plotHeight 30 --plotWidth 30 --removeOutliers --log1p

done
