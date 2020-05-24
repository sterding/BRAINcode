#!/bin/bash
###########################################
# script to combine bigwig, using trimmed mean of unionBedGraphs for samples from the same group
# Usage: $0 HC_SNDA
# bsub -J HC_SN -q normal -n 4 -M 6000 _combine_bigwig.sh HC_SN
# for i in HCILB_SNDA HC_PY HC_nonNeuron HC_Neuron HC_MCPY HC_TCPY HC_SNDA ILB_SNDA HC_PBMC HC_FB PD_SNDA HC_SNDAstranded; do echo $i; bsub -J $i -q normal -n 4 -M 6000 _combine_bigwig.sh $i; done
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
#!/bin/bash

module load ucsc/2017
export TMPDIR=/data/neurogen/

config_file=config.txt
[ -e "$config_file" ] || config_file=$HOME/neurogen/pipeline/RNAseq/config.txt
echo "Using configuration file at:" $config_file;
source $config_file

list_bw_file=$1
group_lable=$2

format=`head -n1 $list_bw_file`
format=${format##*.} # get extension

N=`cat $list_bw_file | wc -l`
echo "N = $N"

# the nuclear genome size:
# curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes | grep -v "_" | grep -v chrM | awk '{s+=$2}END{print s}'  # 3095677412
GENOME_SIZE=3095677412

if [ "$N" == "1" ]
then
  if [[ $format =~ bigwig|BIGWIG|bw|BW ]];
  then
    awk -v GENOME_SIZE=$GENOME_SIZE '{OFS="\t";S+=$4*($3-$2); }END{TOTAL=GENOME_SIZE; bc=S/TOTAL;print "#basalCoverage="bc, "#"S"/"TOTAL;}' `sed 's/bw/bedGraph/g' $list_bw_file` > mean.$group_lable.bedGraph
    ln -fs `cut -f2 $list_bw_file` mean.$group_lable.bw
  elif [[ $format =~ bam|BAM|cram|CRAM ]];
  then
    genomeCoverageBed -ibam `cut -f2 $list_bw_file` -g $GENOME/Annotation/Genes/ChromInfo.txt -bg | awk -v GENOME_SIZE=$GENOME_SIZE '{OFS="\t";S+=$4*($3-$2); print;} END{TOTAL=GENOME_SIZE; bc=S/TOTAL;print "#basalCoverage="bc, "#"S"/"TOTAL;}'> mean.$group_lable.bedGraph
    ln -fs `cut -f2 $list_bw_file` mean.$group_lable.bw
  else
    echo "the input files must be bam, cram, or bigwig format."
  fi
  exit
fi
    
#-------------------
echo "["`date`"] computing mean of bedGraph (using up to 100 samples randomly picked)"
#-------------------

cut -f2 $list_bw_file | sort --random-sort --random-source=$list_bw_file | head -n100 > $list_bw_file.top100
N=`cat $list_bw_file.top100 | wc -l`

if [[ $format =~ bigwig|BIGWIG|bw|BW ]];
then
  ## solution2: using ucsc-bigWigMerge to add signal values of multiple bigWigs together into a single output bedGraph
  # merge, sort, and add the last line, and also convert sum to mean
  [ -e mean.$group_lable.bedGraph ] || bigWigMerge -inList $list_bw_file.top100 stdout | LC_ALL=C sort -S 1G -T $TMPDIR -k1,1 -k2,2n | awk -v N=$N -v GENOME_SIZE=$GENOME_SIZE '{OFS="\t";$4=$4/N; S+=$4*($3-$2); if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id; TOTAL=GENOME_SIZE; bc=S/TOTAL;print "#basalCoverage="bc, "#"S"/"TOTAL;}' > mean.$group_lable.bedGraph
elif [[ $format =~ bam|BAM|cram|CRAM ]];
then
  ## solution3: using samtools-merge to merge multiple BAM/CRAM files together 
  ## Caution: In the way, the sequencing depth factor is not normalized. This might be a problem if the sequencing depth varies a lot.
  [ -e mean.$group_lable.bedGraph ] || samtools merge -b $list_bw_file.top100 -cp - | genomeCoverageBed -ibam stdin -g $GENOME/Annotation/Genes/ChromInfo.txt -bg -split | LC_ALL=C sort -S 1G -T $TMPDIR -k1,1 -k2,2n | awk -v N=$N -v GENOME_SIZE=$GENOME_SIZE '{OFS="\t";$4=$4/N; S+=$4*($3-$2); print;} END{TOTAL=GENOME_SIZE; bc=S/TOTAL;print "#basalCoverage="bc, "#"S"/"TOTAL;}'> mean.$group_lable.bedGraph
else
  echo "the input files must be bam, cram, or bigwig format."
fi

#-------------------
echo "["`date`"] sorting bedgraph"
#-------------------

LC_ALL=C sort -S 1G -T $TMPDIR -k1,1 -k2,2n mean.$group_lable.bedGraph > mean.$group_lable.bedGraph2
mv mean.$group_lable.bedGraph2 mean.$group_lable.bedGraph

#-------------------
echo "["`date`"] bedgraph --> bigwig"
#-------------------
bedGraphToBigWig mean.$group_lable.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt mean.$group_lable.bw

#-------------------
echo "["`date`"] Done!"
#-------------------
