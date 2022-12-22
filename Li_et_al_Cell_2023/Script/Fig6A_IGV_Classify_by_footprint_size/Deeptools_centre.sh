#!/bin/bash

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd 

bamToBed -bedpe -i HChIP_hMCM_rep1.PEASXSuniq.nodup.bam -bedpe -ed | cut -f 1-3 | awk '{ $4 = $3 - $2 } 1' | sed -e 's/ /\t/g' > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/MCM_R1.bed

bamToBed -bedpe -i HChIP_hMCM_rep2.PEASXSuniq.nodup.bam -bedpe -ed | cut -f 1-3 | awk '{ $4 = $3 - $2 } 1' | sed -e 's/ /\t/g' > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/MCM_R2.bed



cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/IGV/RelatedData/MCM_Fragment/


bedtools makewindows -g ChrLength.txt -w 1000 | awk '{print $0,NR}' | sed -e 's/ /\t/g' > hg19_Bin.bed 

cd Classify

for i in $(ls | grep ".bed")
do

  bedtools coverage -counts -a ../hg19_Bin.bed -b $i -F 0.5 | cut -f 1-3,5 > ../Classify_Count/${i/%.bed/.bedgraph}

done