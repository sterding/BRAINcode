#!/bin/bash
#############################################
# Author: Xianjun Dong
# Email: xdong@rics.bwh.harvard.edu
# Version: 0.1
#############################################
# Bash script convert BAM to bigwig and bedgraph
# requirements:
# 1) install bedtools and Jim Kent's utility and set path in $PATH
# 2) require sam2bed awk script and set its path in $PATH
# 3) Download Annotation and set the path of $ANNOTATION (recommend to use annotation files downloaded from iGenome)

# current usage: bam2bigwig <in.bam|sam|bed> <-split|-nosplit>
# optimal usage: bam2bigwig [-sB] -I <mm9> -i <in.bam|sam> -o <out.bw>
# -B: output bedGraph (instead of bigWig)
# -I: index for annnotated genome (e.g. hg19, mm8 etc.)
# -i: input file
# -o: output file

## Log:
#1. change to use "bamToBed -bed12" and "bedItemOverlapCount -bed12", instead of "bamToBed -split"
#2. 10/10/13: change to use hg19 as default reference 
#3. 10/15/13: allow to convert reads from specific region (only works for BAM format now)
#4. 11/21/14: also generate RPM normalized bigwig
#5. 12/16/15: add normalizationFactor. If ==0, then no normalization; otherwise, use it as the M (in RPM) for normalization
 
export TMPDIR=/data/neurogen/

config_file=config.txt
[ -e "$config_file" ] || config_file=$HOME/neurogen/pipeline/RNAseq/config.txt
echo "Using configuration file at:" $config_file;
source $config_file

inputfile=$1  # bam or cram
split=$2

[ -e "$1" ] || {
	echo "
# Script convert BAM to normalized bigwig and bedgraph
# Usage: bam2bigwig.v2 in.bam <-split|-nosplit>
# -split is to split input file in +/- strand (Default)
# -nosplit is to not split input file by strand
"; exit 0;}

[[ $# -eq 2 ]] || split="-nosplit"

ext=${inputfile##*.}
bname=${inputfile%.*}

#echo "split=" $split "="$2"=";

case $ext in
    bam)
        echo "Input is a $ext file. Converting bam-->bedgraph-->bigwig ...";;
    cram)
        echo "Input is a $ext file. Converting cram--> bam --> bedgraph-->bigwig ...";;
    *)
        echo "Unsupported input format. Exit!"
        exit;;
esac

RPMscale=$(bc <<< "scale=6;1000000/$(samtools view -F 0x100 -c $inputfile)") # reads per million total reads

if [ "$split" == "-split" ]
then
    echo "bam-->bw, by strand..."
    [ "$ext" == "cram" ] && samtools view -b -T $GENOME/Sequence/WholeGenomeFasta/genome.fa $inputfile | bedtools genomecov -ibam stdin -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split -strand + | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.plus.normalized.bedGraph
    [ "$ext" == "bam" ] && bedtools genomecov -ibam $inputfile -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split -strand + | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.plus.normalized.bedGraph
    bedGraphToBigWig $bname.plus.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.plus.normalized.bw
    
    [ "$ext" == "cram" ] && samtools view -b -T $GENOME/Sequence/WholeGenomeFasta/genome.fa $inputfile | bedtools genomecov -ibam stdin -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split -strand - | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.minus.normalized.bedGraph
    [ "$ext" == "bam" ] && bedtools genomecov -ibam $inputfile -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split -strand - | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.minus.normalized.bedGraph
    bedGraphToBigWig $bname.minus.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.minus.normalized.bw
fi

if [ "$split" == "-nosplit" ]
then
    echo "bam-->bw, without spliting by strand..."
    [ "$ext" == "cram" ] && samtools view -b -T $GENOME/Sequence/WholeGenomeFasta/genome.fa $inputfile | bedtools genomecov -ibam stdin -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.normalized.bedGraph
    [ "$ext" == "bam" ] && bedtools genomecov -ibam $inputfile -bg -scale $RPMscale -g $GENOME/Annotation/Genes/ChromInfo.txt -split | LC_COLLATE=C sort -k1,1 -k2,2n > $bname.normalized.bedGraph
    bedGraphToBigWig $bname.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.normalized.bw
fi
