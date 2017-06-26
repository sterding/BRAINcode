# script to calculate the occurance of each feature (and its counterpart) of featuresA in and beyound intervalA
# Usage: 
# bash $pipeline_path/src/eRNA.TFBSjaspar.enrichment.sh ~/eRNAseq/externalData/TFBS/JASPARmotifscan.hg19.bed eRNA.JASPARmotifscan.txt

#
#       |  A  |  nA  | 
# ------|-----|------|
#   B   |  AB |  nAB | B 
# ------|-----|------|
#   nB  | AnB | nAnB | nB
# ------|-----|------|
#       |  A  |  nA  | 

featuresA=$1 # ~/eRNAseq/externalData/TFBS/JASPARmotifscan.hg19.bed # N=136,626,803
output=$2 # eRNA.JASPARmotifscan.txt

> $output;
#echo -e "ID\tAB\tAnB\tnAB\tnAnB" 

while read JASPARID tf;
do
	nA=`mktemp`
	# change TF to upper case, as all $TF_in_hg19.bed was named in UPPER CASE
	TF=${tf^^} # http://stackoverflow.com/a/2265268
	
	#A=~/neurogen/TF_scan/all_TF_output_test/$TF/${TF}_in_hg19.bed
	#cat `ls  ~/neurogen/TF_scan/all_TF_output_test/*/*_in_hg19.bed | grep -v "/"$TF"/"` > $nA;
	
	A=~/neurogen/TF_scan/all_TF_output_test_by_JASPARID/${JASPARID}_in_hg19.bed
	cat `ls  ~/neurogen/TF_scan/all_TF_output_test_by_JASPARID/MA*.bed | grep -v "/"$JASPARID"/"` > $nA;
	
	#for j in eRNA eRNA-classI promoter exon random
	for j in  eRNA-classI
	do
	  [ "$j" == "eRNA" ] && intervalB=eRNA.bed;
	  [ "$j" == "eRNA-classI" ] && intervalB=eRNA.classI.bed;
	  [ "$j" == "promoter" ] && intervalB=$GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed;
	  [ "$j" == "exon" ] && intervalB=$GENOME/Annotation/Genes/mRNA.innner.exon.bed;
	  [ "$j" == "random" ] && intervalB=eRNA.random.bed;
	  AB=`intersectBed -a $A -b $intervalB -u | wc -l`
	  AnB=`intersectBed -a $A -b $intervalB -v | wc -l`
	  nAB=`intersectBed -a $nA -b $intervalB -u | wc -l`
	  nAnB=`intersectBed -a $nA -b $intervalB -v | wc -l`
	  echo -e "$TF\t$AB\t$AnB\t$nAB\t$nAnB\t$j";
	  echo -e "$TF\t$AB\t$AnB\t$nAB\t$nAnB\t$j" >> $output;	
	done
done < ~/neurogen/TF_scan/motif_contingency_table/JASPARmotifsbyID.txt

#echo "## Fisher test and make plot"
### ##################
#Rscript $pipeline_path/src/eRNA.TFBSencode.enrichment.R $output
