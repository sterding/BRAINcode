###########################################
# script to generate coverage for combined samples
# Usage: $0 HC_SNDA 5 0.5
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
#!/bin/bash

group_lable=$1
N=$2  # at least N number of unique reads (if N is a int) or N rpm values (if N is float)
randN=$3 # randomly selected N samples

[ "$N" = "" ] && N=0.05;

n=`cat /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable | wc -l`
echo "total samples n = $n"

if [[ ! "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with reads cutoff"
    
    unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bw.gt5reads.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}reads.$group_lable.bed;
    sort -k4n,4n covered.${N}reads.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}reads.$group_lable.txt
fi

if [[ "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with RPM cutoff"
    
    if [[ "$randN" = "" ]]; then
      unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}RPM.$group_lable.bed;
      sort -k4,4n covered.${N}RPM.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}RPM.$group_lable.txt
    fi
    if [[ ! "$randN" = "" ]]; then
      TMPFILE=`mktemp random.covered.${N}RPM.$group_lable.XXXXXXXXXX` || exit 1
      unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable | shuf -n $randN` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > $TMPFILE.bed;
      sort -k4,4n $TMPFILE.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > $TMPFILE.txt
    fi  
fi


# for raw reads: at least N unique reads in at least M samples
#[[ ! "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"reads"M"samples:", s;}' > $group_lable.covered.${N}reads${M}samples.bed

# for normalized RPM: at least N RPM in at least M samples
#[[ "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"RPM"M"samples:", s;}' > $group_lable.covered.${N}RPM${M}samples.bed