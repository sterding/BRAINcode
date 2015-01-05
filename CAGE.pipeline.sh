####################################
# Pipeline for CAGE data Analysis
# Authors: Xianjun DOng
# Email: xdong@rics.bwh.harvard.edu
# Date: 11/21/2014
# Version: 0.0
####################################
#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: $HOME/neurogen/pipeline/RNAseq/CAGE.pipeline.sh /data/neurogen/CAGE_PDBrainMap/rawfiles"
  exit
fi


########################
## 0. setting 
########################
pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

# project folders
input_dir=$1  # input_dir=/data/neurogen/CAGE_PDBrainMap/rawfiles

# create the subfolders (e.g. processed, for_display, results)
processed=$input_dir/../processed
[ -d $processed ] || mkdir $processed

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results 
[ -d $result_dir ] || mkdir $result_dir

########################
# 1. convert bam to bigwig
########################

cd $processed
for i in *bam; do
    bsub -J _bam2bigwig -oo _bam2bw.$i.log -eo _bam2bw.$i.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] -u $EMAIL -N _CAGE_bam2bigwig.sh $i;
done

# continue until the jobs are done

########################
# 2. computing trimmed mean of bedGraph
########################

unionBedGraphs -i `ls *plus.normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); sum=0; for(i=a+1;i<=(c-a);i++) sum+=j[i];return sum/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > trimmedmean.plus.normalized.bg
bedGraphToBigWig trimmedmean.plus.normalized.bg $GENOME/Annotation/Genes/ChromInfo.txt trimmedmean.plus.normalized.bw

unionBedGraphs -i `ls *minus.normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); sum=0; for(i=a+1;i<=(c-a);i++) sum+=j[i];return sum/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > trimmedmean.minus.normalized.bg
bedGraphToBigWig trimmedmean.minus.normalized.bg $GENOME/Annotation/Genes/ChromInfo.txt trimmedmean.minus.normalized.bw




