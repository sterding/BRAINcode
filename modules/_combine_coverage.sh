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
M=$3  # at least M number of samples (i.e. M is a int value between [1, (NF-4)]), or M percent of samples (i.e. M is a float value between (0, 1.0]) 

[ "$N" = "" ] && N=5;
[ "$M" = "" ] && M=0.5;

# for debug
# group_lable="HC_TCPY"; pattern="(HC|ND)_.*_TCPY_[2345]"; N=5; M=0.5

[ "$group_lable" = "HC_TCPY" ]  && pattern="(HC|ND)_.*_TCPY_[^1]";
[ "$group_lable" = "HC_MCPY" ]  && pattern="(HC|ND)_.*_MCPY_[^1]";
[ "$group_lable" = "HC_SNDA" ]  && pattern="(HC|ND)_.*_SNDA_[^1]_[^s]*/";
[ "$group_lable" = "ILB_SNDA" ] && pattern="ILB_.*_SNDA_[^1]";
[ "$group_lable" = "PD_SNDA" ]  && pattern="PD_.*_SNDA_[^1]";
[ "$group_lable" = "HCILB_SNDA" ]  && pattern="(HC|ND|ILB)_.*_SNDA_[^1]_[^s]*/";  # excep batch1 and stranded libs
# control
[ "$group_lable" = "HC_SN" ]  && pattern="HC_.*_SN_[^u]*/";
[ "$group_lable" = "HC_SNDAstranded" ]  && pattern="HC_.*_SNDA_.*stranded";
[ "$group_lable" = "HC_PBMC" ]  && pattern="HC_.*_PBMC_[^u]*/";
[ "$group_lable" = "HC_FB" ]  && pattern="HC_.*_FB_";

echo $group_lable, "$pattern";
ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern"
    
#-------------------

if [[ ! "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with reads cutoff"
    
    #unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bedGraph | grep -E "$pattern"` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}reads.$group_lable.bed;
    unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bw.gt5reads.bedGraph | grep -E "$pattern"` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}reads.$group_lable.bed;
    sort -k4n,4n covered.${N}reads.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}reads.$group_lable.txt
    #intersectBed -a covered.${N}reads.$group_lable.bed -b $GENOME/Annotation/Genes/exons.bed | sort -k4n,4n  | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}reads.$group_lable.txt
fi

#-------------------

if [[ "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with RPM cutoff"
    
    unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern"` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}RPM.$group_lable.bed;
    sort -k4,4n covered.${N}RPM.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}RPM.$group_lable.txt
fi

# for raw reads: at least N unique reads in at least M samples
#[[ ! "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"reads"M"samples:", s;}' > $group_lable.covered.${N}reads${M}samples.bed

# for normalized RPM: at least N RPM in at least M samples
#[[ "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"RPM"M"samples:", s;}' > $group_lable.covered.${N}RPM${M}samples.bed