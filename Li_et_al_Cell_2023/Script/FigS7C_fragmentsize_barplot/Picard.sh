#! /bin/bash

cd /Volumes/Seagate\ Expansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd

java -jar picard.jar SortSam I=R1_Bigger_75bp.bam O=R1_Bigger_75bp_sorted.bam SORT_ORDER=coordinate

java -jar picard.jar SortSam I=R2_Bigger_75bp.bam O=R2_Bigger_75bp_sorted.bam SORT_ORDER=coordinate use_jdk_deflater=true use_jdk_inflater=true




java -jar picard.jar BuildBamIndex I=R1_Bigger_75bp_sorted.bam

java -jar picard.jar BuildBamIndex I=R2_Bigger_75bp_sorted.bam