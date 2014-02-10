# ===============================================================
# Pipeline for calculating binned density of any features (in BigWig format)
# Author: Xianjun Dong
# Date: Mon Feb 10
# Usage:
# bigWigAverageOverBed_81bins.sh input.bed12.mRNA.81bins.bed 81 accepted_hits.bw
# or in cluster
# for i in *.bw; do echo bsub -q normal -u sterding.hpcc@gmail.com -N -M 4000 -R rusage[mem=4000] bigWigAverageOverBed_81bins.sh $GENOME/Annotation/Genes/gencode.v19.annotation.bed12.mRNA.81bins.bed 81 $i; done
# Requirement:
# Step1: get bin file in bed format (e.g. awk -v option=mRNA -f bigWigAverageOverBed_generate_81bins.awk input.bed12 > input.bed12.mRNA.81bins.bed)
# Step2: get bigwig file
#
## CHANGE: use cut -f1,5 to get mean0 (instead of mean) from bigWigAverageOverBed- average over bases with non-covered bases counting as zeroes
# ===============================================================

#!/bin/sh

inputfile=$1
option=$2 # number of bins: 80 or 41  (80bins: 100bp bins for both TSS and TTS; 41bins: 100bp bins for TSS-flanking, [2k, TTS] as one big bin;)
filename=$3

export JOBOUTPUT_DIR=~/neurogen/rnaseq_PD/results/aggregationplot
export JOBOUTPUT=$JOBOUTPUT_DIR/output_`basename $inputfile`_of_`basename $filename`
export JOBOUTPUT_HEADER=$JOBOUTPUT_DIR/header_`basename $inputfile`

# Make the directory for the job ID you are running if it does not exist
[ -d $JOBOUTPUT_DIR ] || mkdir -p $JOBOUTPUT_DIR

[ -f $JOBOUTPUT_HEADER ] || cut -f4 -d' ' $inputfile | sed 's/.bin.*//g' | sort -u > $JOBOUTPUT_HEADER

if [[ "$option" == '80' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%80==0) printf $1"\n";}' | sort -k81,81 | cut -f1-80 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi
if [[ "$option" == '41' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%41==0) printf $1"\n";}' | sort -k42,42 | cut -f1-41 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi
if [[ "$option" == '81' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%81==0) printf $1"\n";}' | sort -k82,82 | cut -f1-81 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi