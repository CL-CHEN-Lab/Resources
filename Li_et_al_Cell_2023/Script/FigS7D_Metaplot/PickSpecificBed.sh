#! /bin/bash

cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810

bedtools bamtobed -i input_hMCM.fastq.nodup.bam > Input_Control.bed

cut -f 1-3 Input_Control.bed | awk '{ $4 = $3 - $2 } 1'|awk '($4 == 55)' | sed -e 's/ /\t/g' | sort -u > Input_Control_55bp.bed

bedtools getfasta -fi /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd/Hg19/Hg19.fa -bed Input_Control_55bp.bed | grep -v "chr" > Input_Control_55bp.fa



cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/NewTasks/All_plot/VennPlot/S50_1kbBin


bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810/Input_Control_55bp.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '( $4 <=0.25  && $4 > 0 )' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Early.bed

bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810/Input_Control_55bp.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '( $4 <=0.5  && $4 > 0.25 )' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidEarly.bed

bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810/Input_Control_55bp.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '( $4 <=0.75  && $4 > 0.5 )' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidLate.bed

bedtools intersect -a /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/SingleEnd/hMCM_input_data_20210810/Input_Control_55bp.bed -b ./hg19_window_S50.bed -wao -F 0.5 -f 0.5 -e | awk '{Sum_S50[$1"\t"$2"\t"$3]+=$8;Amount_windows[$1"\t"$2"\t"$3]++} END { for (i in Sum_S50) { print i"\t"(Sum_S50[i]/Amount_windows[i])}}' | awk '( $4 <=1  && $4 > 0.75 )' | sort -u | sortBed -faidx ./Hg19.txt > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Late.bed




cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd 

#### For gglogo plot ####
bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Early.bed | grep -v ">chr"| sort -u > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Early.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidEarly.bed | grep -v ">chr"| sort -u > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidEarly.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidLate.bed | grep -v ">chr"| sort -u > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidLate.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Late.bed | grep -v ">chr"| sort -u > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Late.fa 




#### For oddlogo plot ####
bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Early.bed > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/WithChr/55bp_InputControl_Early.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidEarly.bed > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/WithChr/55bp_InputControl_MidEarly.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_MidLate.bed > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/WithChr/55bp_InputControl_MidLate.fa 

bedtools getfasta -fi ./Hg19/Hg19.fa -bed /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/55bp_InputControl_Late.bed > /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control/WithChr/55bp_InputControl_Late.fa 



logoddslogo=/Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/logoddslogo-1.1.1/logoddslogo

AT=(27 31 29 30.5 27 31 29 30.5)
GC=(23 19 21 19.5 23 19 21 19.5)


cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/LogoWithGC/RawData/Input_Control
j=0
for i in $(ls | grep ".fa");
do
      Echo $i
      Echo $j
   
      $logoddslogo --format PDF -f ./WithChr/$i -D fasta -U 'bits' --composition "{'A':${AT[$j]}, 'C':${GC[$j]}, 'G':${GC[$j]}, 'T':${AT[$j]}}" -n 55 --aspect-ratio 100 -S 0.3 -W 80 -o ../../Plot/Input_Control/${i/%.fa/_Bits.pdf}; 

      #$logoddslogo --format PDF -f ./WithChr/$i -D fasta -U 'probability' --composition "{'A':${AT[$j]}, 'C':${GC[$j]}, 'G':${GC[$j]}, 'T':${AT[$j]}}" -n 55 --aspect-ratio 56 -S 0.3 -W 30 -o ../../Plot/${All_SeqLength[$k]}/${i/%.fa/_Probability.pdf}; 

      j=$(($j+1))

done









