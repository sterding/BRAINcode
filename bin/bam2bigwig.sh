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
 
species=hg19

ANNOTATION=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Annotation/Genes

inputfile=$1
split=$2
normalizationFactor=$3
region=$4 # in the following format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ (region between 1,000,000 and 2,000,000bp including the end points). The coordinate is 1-based.

#[[ -e "$1" ] && 
[ "$2" != "" -a -e "$1" ] || {
	echo "
# Script convert BAM to bigwig and bedgraph
# Usage: bam2bigwig <in.bam|sam|bed> <-split|-nosplit> [chr4|chr4:232,2323-343,5454]
# -split is to split input file in +/- strand (Default)
# -nosplit is to not split input file by strand
"; exit 0;}

ext=${inputfile##*.}
bname=${inputfile%.*}

#echo $split;

case $ext in
    bam|BAM|Bam)
        echo "Input is a BAM file. Converting bam-->sam-->bed ..."
        # use the XS:A for strand
        samtools view $inputfile $region | sam2bed -v bed12=T -v sCol=NH > $bname.bed;;
        # use the FLAG for strand
        #bamToBed -i $inputfile -bed12 > $bname.bed;;
    sam|SAM|Sam)
        echo "Input is a SAM file. Converting sam->bed ..."
        sam2bed -v bed12=T -v sCol=NH $inputfile > $bname.bed;;
    bed|BED|Bed|bed6|bed12)
        echo "Input is a BED file.";;
    *)
        echo "Unsupported input format. Exit!"
        exit;;
esac

# remove reads mapped to contigs other than chr1/2/.../22/M/X/Y
awk '{if($1!~/[_#]/ && $1~/^chr/) {n++; print}}END{print "#total_mapped_reads="n;}' $bname.bed > $bname.bed2
mv $bname.bed2 $bname.bed

# get total mapped reads
total_mapped_reads=`tail -n1 $bname.bed | cut -f2 -d'='`

if [ "$2" == "-split" ]
then
    echo "bed-->bw, by strand..."
    echo "bedToBedGraph..."
    > $bname.bed+; > $bname.bed-;
    awk -v OUTPUT=$bname '{print $0 >> OUTPUT".bed"$6}' $bname.bed
    ncol=`head $bname.bed | awk 'END{print NF}'`
    if [ "$ncol" -eq  12 ]
    then
        echo "in bed12 format ..."
        sort -k1,1 $bname.bed+ | bedItemOverlapCount $species -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.plus.bedGraph
        sort -k1,1 $bname.bed- | bedItemOverlapCount $species -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $bname.minus.bedGraph
    else
        echo "in bed6 format ..."
        sort -k1,1 $bname.bed+ | bedItemOverlapCount $species -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.plus.bedGraph
        sort -k1,1 $bname.bed- | bedItemOverlapCount $species -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $bname.minus.bedGraph
    fi
    echo "bedGraph2bw..."
    bedGraphToBigWig $bname.plus.bedGraph $ANNOTATION/ChromInfo.txt $bname.plus.bw
    bedGraphToBigWig $bname.minus.bedGraph $ANNOTATION/ChromInfo.txt $bname.minus.bw
    
    [ "$normalizationFactor" -gt 0 ] && \
    echo "generated normalized bg..."  && \
    awk -v tmr=$normalizationFactor 'BEGIN{OFS="\t"; print "# total_mapped_reads="tmr;}{$4=$4*1e6/tmr; print}' $bname.plus.bedGraph > $bname.plus.normalized.bedGraph && \
    awk -v tmr=$normalizationFactor 'BEGIN{OFS="\t"; print "# total_mapped_reads="tmr;}{$4=$4*1e6/tmr; print}' $bname.minus.bedGraph > $bname.minus.normalized.bedGraph && \
    bedGraphToBigWig $bname.plus.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.plus.normalized.bw && \
    bedGraphToBigWig $bname.minus.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.minus.normalized.bw

    #echo "removing temp files..."
    #rm $bname.*.bedGraph $bname.bed?
fi

if [ "$2" == "-nosplit" ]
then
    echo "bed-->bw, without spliting by strand..."
    echo "bedGraph2bw..."
    ncol=`head $bname.bed | awk 'END{print NF}'`
    if [ "$ncol" -eq  12 ]
    then
        echo "in bed12 format ..."
        sort -k1,1 $bname.bed | bedItemOverlapCount $species -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.bedGraph
    else
        echo "in bed6 format ..."
        sort -k1,1 $bname.bed | bedItemOverlapCount $species -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.bedGraph
    fi
    echo "bedGraph2bw..."
    bedGraphToBigWig $bname.bedGraph $ANNOTATION/ChromInfo.txt $bname.bw

    [ "$normalizationFactor" -gt 0 ] && \
    echo "generated normalized bg..." && \
    awk -v tmr=$normalizationFactor 'BEGIN{OFS="\t"; print "# total_mapped_reads="tmr;}{$4=$4*1e6/tmr; print}' $bname.bedGraph > $bname.normalized.bedGraph && \
    bedGraphToBigWig $bname.normalized.bedGraph $ANNOTATION/ChromInfo.txt $bname.normalized.bw
    
    #echo "removing temp files..."
    #rm $bname.bedGraph
fi
