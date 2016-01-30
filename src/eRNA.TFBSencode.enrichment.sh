# script to calculate the occurance of each feature (and its counterpart) of featuresA in and beyound intervalA
# Usage: 
# bash $pipeline_path/src/eRNA.TFBSencode.enrichment.sh

#
#       |  A  |  nA  | 
# ------|-----|------|
#   B   |  AB |  nAB | B 
# ------|-----|------|
#   nB  | AnB | nAnB | nB
# ------|-----|------|
#       |  A  |  nA  | 

#intervalB=$1  # a bed file: chr, start, end
featuresA=$1  # a bed file with name column: chr, start, end, name
output=$2

featuresA=~/eRNAseq/externalData/TFBS/factorbookMotifPos.v3.mergebyTF.bed
output=~/eRNAseq/HCILB_SNDA/eRNA.factorbookMotifPos.txt

featuresA=~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12
output=~/eRNAseq/HCILB_SNDA/eRNA.wgEncodeRegTfbsClusteredV3.txt

> $output;
#echo -e "ID\tAB\tAnB\tnAB\tnAnB" 

for i in `cut -f4 $featuresA | sort -u`
do
	nA=`mktemp`
	awk -vNAME=$i '$4==NAME' $featuresA | intersectBed -a $featuresA -b - -v | sortBed | mergeBed > $nA;
	#for j in eRNA eRNA-classI promoter exon random
	for j in  eRNA-classI
	do
	  [ "$j" == "eRNA" ] && intervalB=~/eRNAseq/HCILB_SNDA/eRNA.bed;
	  [ "$j" == "eRNA-classI" ] && intervalB=~/eRNAseq/HCILB_SNDA/eRNA.classI.bed;
	  [ "$j" == "promoter" ] && intervalB=$GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed;
	  [ "$j" == "exon" ] && intervalB=$GENOME/Annotation/Genes/mRNA.innner.exon.bed;
	  [ "$j" == "random" ] && intervalB=~/eRNAseq/HCILB_SNDA/eRNA.random.bed;
	  AB=`awk -vNAME=$i '$4==NAME' $featuresA | intersectBed -a - -b $intervalB -u | wc -l`
	  AnB=`awk -vNAME=$i '$4==NAME' $featuresA | intersectBed -a - -b $intervalB -v | wc -l`
	  nAB=`intersectBed -a $nA -b $intervalB -u | wc -l`
	  nAnB=`intersectBed -a $nA -b $intervalB -v | wc -l`
	  echo -e "$i\t$AB\t$AnB\t$nAB\t$nAnB\t$j";
	  echo -e "$i\t$AB\t$AnB\t$nAB\t$nAnB\t$j" >> $output;	
	done
done

echo "## Fisher test and make plot"
### ##################
Rscript $pipeline_path/src/eRNA.TFBSencode.enrichment.R $output
