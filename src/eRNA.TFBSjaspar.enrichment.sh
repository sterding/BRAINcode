#!/usr/bin/env bash
# script to calculate the occurance of each feature (and its counterpart) of featuresA in and beyound intervalA
# Usage: bash $pipeline_path/src/eRNA.TFBSjaspar.enrichment.sh eRNA.JASPARmotifscan.txt
# Update: re-scan the genome with fimo and all non-redudant vertebrate JASPAR motif. See ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/README.txt (3/18/2018) 
#
#       |  A  |  nA  | 
# ------|-----|------|
#   B   |  AB |  nAB | B 
# ------|-----|------|
#   nB  | AnB | nAnB | nB
# ------|-----|------|
#       |  A  |  nA  | 

output=$1 # eRNA.JASPARmotifscan.txt

> $output;

while read JASPARID TF;
do
	nA=`mktemp`

	A=~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/$JASPARID.meme_in_genome.fa.bed
	#cat `ls  ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/MA*.bed | grep -v $JASPARID` > $nA;
	## cat MA*.bed | LC_ALL=C sort --parallel=24 --buffer-size=5G -k1,1 -k2,2n > ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.in_genome.fa.bed
	sortBed -i $A | intersectBed -a ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.in_genome.fa.bed -b - -v -sorted > $nA;
	
	for j in  eRNA-private eRNA eRNA-classI;
	do
	  [ "$j" == "eRNA" ] && intervalB=eRNA.bed;
	  [ "$j" == "eRNA-classI" ] && intervalB=eRNA.classI.bed;
	  [ "$j" == "eRNA-private" ] && intervalB=eRNA.private.major.bed;
	  [ "$j" == "promoter" ] && intervalB=$GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed;
	  AB=`intersectBed -a $A -b $intervalB -u | wc -l`
	  AnB=`intersectBed -a $A -b $intervalB -v | wc -l`
	  nAB=`intersectBed -a $nA -b $intervalB -u -sorted  | wc -l`
	  nAnB=`intersectBed -a $nA -b $intervalB -v -sorted  | wc -l`
	  echo -e "$TF\t$AB\t$AnB\t$nAB\t$nAnB\t$j";
	  echo -e "$TF\t$AB\t$AnB\t$nAB\t$nAnB\t$j" >> $output;	
	done
done < ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/JASPARmotifsID.TFsymbol.txt

#echo "## Fisher test and make plot"
### ##################
#Rscript $pipeline_path/src/eRNA.TFBSencode.enrichment.R $output
