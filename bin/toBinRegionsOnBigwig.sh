#!/bin/bash
# ===============================================================
# Script for binning a list of regions into N bins and get bigwig signal on the bins
# Author: Xianjun Dong
# Date: Jun 3, 2014
# Usage:
# toBinRegionsOnBigwig.sh signal.bw input.bed 100
# or
# toBinRegionsOnBigwig.sh http://path.com/signal.bw input.bed 100
# or
# cat input.bed | toBinRegionsOnBigwig.sh signal.bw - 1 max
# Requirement:
# 1. input bed file at least have chr, start, and end fields
# 2. bigWigSummary (from Jim Kent) should be in the $PATH
# ===============================================================

bigwig=$1
bedfile=$2
N=$3
type=$4

[[ $N == "" ]] && {
  echo "Usage: toBinRegionsOnBigwig.sh signal.bw input.bed 100"; 
  exit 0;}
[[ $type == "" ]] && type="mean"

### for N=1
#[[ $type == "sum" ]] && ncol=4
#[[ $type == "mean" ]] && ncol=5
#[[ $type == "min" ]] && ncol=7
#[[ $type == "max" ]] && ncol=8

## NOTE:
## bug in bigwigsummary: https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/pWjcov-xQyQ
## TODO: switch to use bigWigAverageOverBed $bigwig $bedfile stdout -minMax | cut -f1,$ncol

cat $bedfile | while read chr start end name score strand rest
do
    s=`bigWigSummary $bigwig -udcDir=/tmp -type=mean $chr $start $end $N 2>/dev/null | sed 's/n\/a/0/g'`;
    [[ $s == "" ]] && s=`yes 0 | head -n$N | tr '\n' ' '`
    [[ $strand == "-" ]] && s=`echo $s | tr ' ' '\n' | tac | tr '\n' ' '`
    [[ $name == "" ]] && name="$chr_$start_$end";
    echo $name $s;
done