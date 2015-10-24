## script to define HTNE with the following definition
# Usage: $pipeline_path/src/eRNA.define.sh HCILB_SNDA
# for i in HCILB_SNDA HC_nonNeuron HC_PY; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.define.sh $i; done
# Author: Xianjun Dong
# Version: 2.0
# Date: Oct 23, 2015
# preqisition: 
# 1. run bash /data/neurogen/pipeline/RNAseq/src/get.background.sh to generate the background regions (../blacklist.bed and ../toExclude.bed)

################################################
# HTNE definition:
# 1) density higher than the basal level,  
# 2) summit >0.05 RPM, --> p<0.05 comparing to the transcriptional noise
# 3) located in non-generic regions (e.g. 500bp away from any annotated exons),
# 4) at least 100bp in length,
# 5) don't contain any splicing sites (>10 reads in at least 5 samples) from Tophat
# 6) consistenly strong (q-value<0.05 in at least 3 samples) when comparing with random non-functional background
################################################

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

cd ~/projects/PD/results/eRNA

SAMPLE_GROUP=$1 # SAMPLE_GROUP="HC_nonNeuron" #HC_Neuron
inputBG=~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph

mkdir $SAMPLE_GROUP
cd $SAMPLE_GROUP

echo "# step0: measure transcriptional noise in background genomic regions"
# =================

# make merged signal if not existed
[ -e $inputBG ] || _combine_bigwig.sh $SAMPLE_GROUP

# RNAseq signal distribution in the background region
bedtools random -seed 3 -g $GENOME/Annotation/Genes/ChromInfo.txt -l 1 -n 1000000 | sortBed | intersectBed -a - -b ../blacklist.bed -v -sorted | intersectBed -a $inputBG -b - -sorted -u | cut -f4 > transcriptional.noise.rpm.txt
#intersectBed -a $inputBG -b ../blacklist.bed -sorted -v | awk '{OFS="\t"; print $3-$2, $4}' | shuf -n 1000000 > transcriptional.noise.rpm.txt # TOO SLOW!
Rscript /data/neurogen/pipeline/RNAseq/src/_fit.Tx.noise.R
#for HCILB_SNDA: Dsig: 10**-1.105 == 0.079
#for HC_nonNeuron: 0.124
Dsig=`tail -n1 transcriptional.noise.rpm.pvalues.txt`  

echo "# step1: any regions with summit RPM > peakLevel and border > baseLevel"
# =================
basalLevel=`tail -n1 $inputBG | cut -f2 -d'=' | cut -f1`
awk -vmin=$basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' $inputBG | mergeBed -c 5 -o max > eRNA.tmp1

echo "# step2: summit RPM >=Dsig (density with p<0.05)"
# =================
awk -vD=$Dsig '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' eRNA.tmp1 | mergeBed -d 100 -c 5 -o max > eRNA.tmp2

echo "# step3: located in non-generic regions (e.g. 500bp away from any annotated exons)"
# =================
intersectBed -a eRNA.tmp2 -b ../toExclude.bed -v > eRNA.tmp3

echo "# step4: length > 100nt"
# =================
awk '{OFS="\t"; if(($3-$2)>100) print $1,$2,$3,$1"_"$2"_"$3}' eRNA.tmp3 > eRNA.tmp4

echo "# step5: don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly)"
# =================
# cd ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/
# cat cufflinks-cuffmerge/merged.bed trinity-cuffmerge/all_strand_spliced.chr.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' | bed12ToBed6 | awk '{OFS="\t"; print $1, $2-10,$2+10; print $1,$3-10,$3+10;}' | sortBed | uniq > trinitycufflinks.merged.splicingsites.flanking20nt.bed

# more than 10 splicing reads in at least 5 samples
# for i in  ~/neurogen/rnaseq_PD/run_output/*/junctions.bed; do awk '{OFS="\t"; if($5>10) { split($11,a,","); split($12,b,","); print $1,$2+a[1]-10,$2+a[1]+10; print $1,$2+b[2]-10,$2+b[2]+10}}' $i | sortBed | uniq; done | sort | uniq -c | awk '{OFS="\t"; if($1>5) print $2,$3,$4}' > ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed
 
intersectBed -a eRNA.tmp4 -b ~/neurogen/rnaseq_PD/results/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed -v > eRNA.tmp5

echo "# step6: calculate the significance of eRNA"
# =================

#1: create 100,000 random regions (400bp each) as background and calculate their signals
i=~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$SAMPLE_GROUP.bw
bedtools shuffle -excl ../toExclude.bed -noOverlapping -i eRNA.tmp5 -g $GENOME/Annotation/Genes/ChromInfo.txt | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bigWigAverageOverBed $i stdin stdout | cut -f1,5 > $i.rdbg
#bigWigAverageOverBed $i eRNA.tmp5 stdout | cut -f1,5 | sort -k1,1 > $i.eRNA.meanRPM
while read line
do
    echo " - sample:", $line;
    i=~/neurogen/rnaseq_PD/run_output/$line/uniq/accepted_hits.normalized2.bw
    bedtools shuffle -excl ../toExclude.bed -noOverlapping -i eRNA.tmp5 -g $GENOME/Annotation/Genes/ChromInfo.txt | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bigWigAverageOverBed $i stdin stdout | cut -f1,5 > $i.rdbg
    #bigWigAverageOverBed $i eRNA.tmp5 stdout | cut -f1,5 | sort -k1,1 > $i.eRNA.meanRPM
done < ~/neurogen/rnaseq_PD/results/merged/samplelist.$SAMPLE_GROUP

### 2: distribution of random background, in order to define the cutoff with p=0.0001 significance
Rscript /data/neurogen/pipeline/RNAseq/src/_HTNE.consistency.R $SAMPLE_GROUP

# consistence definition: q<0.05 in at least 3 samples (NEED TO CONFIRM WITH JOHN)
awk '{OFS="\t"; split($1,a,"_"); S=0; if($1~/^chr/) {for(i=2;i<=NF;i++) if($i<=0.05) S++; if(S>=3) print a[1],a[2],a[3],$1}}' eRNA.$SAMPLE_GROUP.qvalue.xls > eRNA.bed