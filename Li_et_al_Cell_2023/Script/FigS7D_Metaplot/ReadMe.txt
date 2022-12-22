Here the raw data are 51~59.bed merged from R1 and R2 files at /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/Rawdata/bam/PairEnd/All_Bed/Odd_51-59bp/Seperate 


cd /Volumes/SeagateExpansi/Chunlong_lab/Ex-Analysis/Old_Project/MCM_project/LastTask/Logo_Sequence/ClassifyByS50/AT_GC_Percentage/RawData

SeqLeg=(53bp 55bp 57bp 59bp)

for i in "${!SeqLeg[@]}";
do

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./Merge_In_Oddbp/${SeqLeg[$i]}.bed ./Merge_In_Oddbp/Input_Control.bed --scoreFileName ./GC_Ratio.bw ./AT_Ratio.bw --outFileName ../Seperate/With_Duplicate/AT_GC_Ratio_${SeqLeg[$i]}.tab.gz --referencePoint center -a 30 -b 30 --binSize 5 --missingDataAsZero --skipZeros 

   plotHeatmap -m ../Seperate/With_Duplicate/AT_GC_Ratio_${SeqLeg[$i]}.tab.gz -out ../Seperate/With_Duplicate/AT_GC_Ratio_${SeqLeg[$i]}.pdf --heatmapHeight 15 --colorMap RdBu --colorList white,lavender,firebrick -x="MCM binding sites +/-30bp" -z "MCM_Binding_Sites" "Input_Control" --yMin 30 --yMax 65
   
done







#############  Then remove duplicate  ##########


for i in $(ls | grep ".bed")
do
     echo $i;
     cut -f 1-3 $i | sort -u > ${i/%.bed/_NoDuplicate.bed}
done


for i in "${!SeqLeg[@]}";
do

   computeMatrix reference-point --numberOfProcessors 8 --regionsFileName ./Merge_In_Oddbp/${SeqLeg[$i]}_NoDuplicate.bed --scoreFileName ./GC_Ratio.bw ./AT_Ratio.bw --outFileName ../Seperate/AT_GC_Ratio_${SeqLeg[$i]}.tab.gz --referencePoint center -a 30 -b 30 --binSize 5 --missingDataAsZero --skipZeros 

   plotHeatmap -m ../Seperate/AT_GC_Ratio_${SeqLeg[$i]}.tab.gz -out ../Seperate/AT_GC_Ratio_${SeqLeg[$i]}.pdf --heatmapHeight 15 --colorMap RdBu --colorList white,lavender,firebrick 
   

done

